#include "seq/dirichlet.h"
#include "util/math_fn.h"
#include "util/vector_output.h"

Dirichlet_mixture::Dirichlet_mixture (const vector<Prob>& pseudocounts)
  : log_cpt_prior (1, (Loge) 0),
    alpha (1, vector<vector<Prob> > (1, pseudocounts))
{ }

Dirichlet_mixture::Dirichlet_mixture (int group_size, int components)
  : log_cpt_prior (components, -Prob2Nats (components)),
    alpha (components, vector<vector<Prob> > (1, vector<Prob> (group_size, (Prob) 0)))
{ }

Dirichlet_mixture::Dirichlet_mixture (const vector<int>& group_size, int components)
  : log_cpt_prior (components, -Prob2Nats (components)),
    alpha (components, vector<vector<Prob> > (group_size.size()))
{
  for (int cpt = 0; cpt < components; ++cpt)
    for (int g = 0; g < (int) group_size.size(); ++g)
      alpha[cpt][g] = vector<Prob> (group_size[g], (Prob) 0);
}

Dirichlet_mixture::Dirichlet_mixture (const Multigroup& multigroup, int components)
  : log_cpt_prior (components, -Prob2Nats (components)),
    alpha (components, vector<vector<Prob> > (multigroup.size()))
{
  for (int cpt = 0; cpt < components; ++cpt)
    for (int g = 0; g < (int) multigroup.size(); ++g)
      alpha[cpt][g] = vector<Prob> (multigroup[g].group_size, (Prob) 0);
}

void Dirichlet_mixture::assert_consistent() const
{
  for (int cpt = 0; cpt < components(); ++cpt)
    {
      if ((int) alpha[cpt].size() != groups())
	THROW Standard_exception ("Dirichlet mixture structure internally inconsistent");
      for (int g = 0; g < groups(); ++g)
	if ((int) alpha[cpt][g].size() != vars(g))
	  THROW Standard_exception ("Dirichlet mixture structure internally inconsistent");
    }
}

void Dirichlet_mixture::calc_log_evidence (const PCounts& pcounts,
					   const vector<PGroup>& pgroups,
					   vector<Loge>& log_evidence,
					   Loge& log_total_evidence) const
{
  log_evidence = vector<Loge> (components());
  log_total_evidence = -InfinityScore;

  vector<vector<Prob> > counts_plus_pseudocounts = alpha[0];  // use alpha[0] just to get dimensions right
  for (int cpt = 0; cpt < components(); ++cpt)
    {
      log_evidence[cpt] = log_cpt_prior[cpt];
      for (int group_idx = 0; group_idx < groups(); ++group_idx)
	{
	  for (int var_idx = 0; var_idx < vars(group_idx); ++var_idx)
	    counts_plus_pseudocounts[group_idx][var_idx] = abs (pcounts[pgroups[group_idx]][var_idx]) + alpha[cpt][group_idx][var_idx];
	  NatsPMulAcc (log_evidence[cpt],
		       Math_fn::log_dirichlet_normaliser (counts_plus_pseudocounts[group_idx])
		       - Math_fn::log_dirichlet_normaliser (alpha[cpt][group_idx]));
	}
      NatsPSumAcc (log_total_evidence, log_evidence[cpt]);
    }
  if (CTAGGING(2,DIRICHLET)) CL << "Dirichlet log-evidence = (" << log_evidence << "), total " << log_total_evidence << "\n"; 
}

void Dirichlet_mixture::postcalc_optimise_scores (const PCounts& pcounts,
						  const vector<PGroup>& pgroups,
						  const vector<Loge>& log_evidence,
						  const Loge& log_total_evidence,
						  PScores& pscores) const
{
  vector<Prob> posterior (components());
  for (int cpt = 0; cpt < components(); ++cpt)
    posterior[cpt] = Nats2Prob (log_evidence[cpt] - log_total_evidence);
  
  for (int group_idx = 0; group_idx < groups(); ++group_idx)
    {
      vector<Prob> avg_val (vars(group_idx), (Prob) 0);
      vector<Prob> cpt_val (vars(group_idx));
      for (int cpt = 0; cpt < components(); ++cpt)
	{
	  double norm = 0;
	  for (int var_idx = 0; var_idx < vars(group_idx); ++var_idx)
	    {
	      cpt_val[var_idx] = abs (pcounts[pgroups[group_idx]][var_idx]) + alpha[cpt][group_idx][var_idx];
	      norm += cpt_val[var_idx];
	    }
	  //	  CLOG(6) << "cpt_val#" << cpt << ": (" << cpt_val << ")\n";
	  for (int var_idx = 0; var_idx < vars(group_idx); ++var_idx)
	    avg_val[var_idx] += posterior[cpt] * cpt_val[var_idx] / norm;
	}
      //      CLOG(6) << "avg_val" << ": (" << avg_val << ")\n";
      for (int var_idx = 0; var_idx < vars(group_idx); ++var_idx)
	pscores[pgroups[group_idx]][var_idx] = Prob2Score (avg_val[var_idx]);
    }
}

void Dirichlet_mixture::initialise_default_scores (PScores& pscores, const vector<PGroup>& pgroups) const
{
  for (int g = 0; g < groups(); ++g)
    {
      vector<Prob> avg_val (vars(g), (Prob) 0);
      for (int cpt = 0; cpt < components(); ++cpt)
	{
	  const double norm = accumulate (alpha[cpt][g].begin(), alpha[cpt][g].end(), 0.0);
	  if (norm > 0)
	    for (int var_idx = 0; var_idx < vars(g); ++var_idx)
	      avg_val[var_idx] += Nats2Prob (log_cpt_prior[cpt]) * alpha[cpt][g][var_idx] / norm;
	}
      for (int var_idx = 0; var_idx < vars(g); ++var_idx)
	pscores.group[pgroups[g].group_idx][var_idx] = Prob2Score (avg_val[var_idx]);
    }
}

// The pseudocount for PVar V in a PGroup of size G (with 0<=V<G) is k*pow(krange,V/(G-1))
Laplace_prior::Laplace_prior (const Multigroup& multigroup, double k, double krange)
  : Dirichlet_mixture (multigroup, 1)
{
  for (int g = 0; g < (int) multigroup.size(); ++g)
    for (int v = 0; v < multigroup[g].group_size; ++v)
      alpha[0][g][v] = k * Math_fn::math_pow (krange, v / max (multigroup[g].group_size - 1, 1));
}

Laplace_prior::Laplace_prior (const PGroup& pgroup, double k, double krange)
  : Dirichlet_mixture (pgroup.group_size, 1)
{
  for (int v = 0; v < pgroup.group_size; ++v)
    alpha[0][0][v] = k * Math_fn::math_pow (krange, v / max (pgroup.group_size - 1, 1));
}

Dirichlet_assignment::Dirichlet_assignment (const Dirichlet_mixture& mixture) :
  mixture (mixture),
  pgroups (vector<PGroup> (mixture.groups(), PGroup(-1)))
{
  mixture.assert_consistent();
}

Dirichlet_assignment::Dirichlet_assignment (const Dirichlet_mixture& mixture, const vector<PGroup>& pgroups) : mixture(mixture), pgroups(pgroups)
{
  mixture.assert_consistent();
  if (mixture.groups() != (int) pgroups.size())
    THROW Standard_exception ("Wrong number of groups for Dirichlet mixture assignment");
  for (int g = 0; g < mixture.groups(); ++g)
    if (pgroups[g].group_size != mixture.vars(g))
      THROW Standard_exception ("Group size mismatch in Dirichlet mixture assignment");
}

void Dirichlet_assignment::calc_log_evidence (const PCounts& pcounts)
{
  mixture.calc_log_evidence (pcounts, pgroups, log_evidence, log_total_evidence);
}

void Dirichlet_assignment::optimise_scores (const PCounts& pcounts, PScores& pscores)
{
  if (CTAGGING(2,PCOUNTS))
    pcounts.show (CL);
  calc_log_evidence (pcounts);
  mixture.postcalc_optimise_scores (pcounts, pgroups, log_evidence, log_total_evidence, pscores);
}

void Dirichlet_assignment::initialise_default_scores (PScores& pscores) const
{
  mixture.initialise_default_scores (pscores, pgroups);
}

void Dirichlet_prior::clear()
{
  assignment.clear();
}

void Dirichlet_prior::assign (const vector<PGroup>& pgroups, const Dirichlet_mixture& mixture)
{
  assignment [Multigroup (pgroups)] = Dirichlet_assignment (mixture, pgroups);
}

void Dirichlet_prior::assign (const PGroup& pgroup, const Dirichlet_mixture& mixture)
{
  vector<PGroup> pgroups (1, pgroup);
  assign (pgroups, mixture);
}

void Dirichlet_prior::assign_Laplace (const PGroup& pgroup, const double laplace_pseudocount, const double pseudocount_range_ratio)
{
  assign (pgroup, Laplace_prior (pgroup, laplace_pseudocount, pseudocount_range_ratio));
}

void Dirichlet_prior::assign_Laplace (PScores& pscores_ref, const double laplace_pseudocount, const double pseudocount_range_ratio)
{
  assignment.clear();
  for (int g = 0; g < pscores_ref.groups(); ++g)
    {
      PGroup pg (g, pscores_ref.group_size (g));
      assign_Laplace (pg, laplace_pseudocount, pseudocount_range_ratio);
    }
}

void Dirichlet_prior::assign_Laplace (PScores& pscores_ref, set<int>& mutable_pgroups, const double laplace_pseudocount, const double pseudocount_range_ratio)
{
  assignment.clear();
  for_const_contents (set<int>, mutable_pgroups, g)
    {
      PGroup pg (*g, pscores_ref.group_size (*g));
      assign_Laplace (pg, laplace_pseudocount, pseudocount_range_ratio);
    }
}

void Dirichlet_prior::add (const Dirichlet_prior& prior)
{
  for_const_contents (Assignment_set, prior.assignment, mg_ass)
    assignment[mg_ass->first] = mg_ass->second;
}

void Dirichlet_prior::optimise (const PCounts& pcounts)
{
  for_contents (Assignment_set, assignment, ass) (*ass).second.optimise_scores (pcounts, *pscores);
}

void Dirichlet_prior::initialise() const
{
  for_const_contents (Assignment_set, assignment, ass) (*ass).second.initialise_default_scores (*pscores);
}

Sjolander_prior::Sjolander_prior (const Alphabet& alphabet) : Alphabet_prior(alphabet,9)
{
  // Taken from Sean Eddy's HMMER, which credits Kimmen Sjolander (Blocks9)

				/* mixture coefficients */
  static float defmq[9] = {
    0.178091, 0.056591, 0.0960191, 0.0781233, 0.0834977, 
    0.0904123, 0.114468, 0.0682132, 0.234585 };

				/* mixture Dirichlet components */
  static float defm[9][20] = {
    { 0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 
      0.131916, 0.012391, 0.022599, 0.020358, 0.030727, 
      0.015315, 0.048298, 0.053803, 0.020662, 0.023612,
      0.216147, 0.147226, 0.065438, 0.003758, 0.009621 },
    { 0.021465, 0.010300, 0.011741, 0.010883, 0.385651, 
      0.016416, 0.076196, 0.035329, 0.013921, 0.093517, 
      0.022034, 0.028593, 0.013086, 0.023011, 0.018866, 
      0.029156, 0.018153, 0.036100, 0.071770, 0.419641 },
    { 0.561459, 0.045448, 0.438366, 0.764167, 0.087364,
      0.259114, 0.214940, 0.145928, 0.762204, 0.247320,
      0.118662, 0.441564, 0.174822, 0.530840, 0.465529, 
      0.583402, 0.445586, 0.227050, 0.029510, 0.121090 },
    { 0.070143, 0.011140, 0.019479, 0.094657, 0.013162, 
      0.048038, 0.077000, 0.032939, 0.576639, 0.072293, 
      0.028240, 0.080372, 0.037661, 0.185037, 0.506783, 
      0.073732, 0.071587, 0.042532, 0.011254, 0.028723 },
    { 0.041103, 0.014794, 0.005610, 0.010216, 0.153602, 
      0.007797, 0.007175, 0.299635, 0.010849, 0.999446, 
      0.210189, 0.006127, 0.013021, 0.019798, 0.014509, 
      0.012049, 0.035799, 0.180085, 0.012744, 0.026466 },
    { 0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 
      0.017255, 0.004911, 0.796882, 0.017074, 0.285858, 
      0.075811, 0.014548, 0.015092, 0.011382, 0.012696, 
      0.027535, 0.088333, 0.944340, 0.004373, 0.016741 },
    { 0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 
      0.105877, 0.049776, 0.014963, 0.094276, 0.027761, 
      0.010040, 0.187869, 0.050018, 0.110039, 0.038668, 
      0.119471, 0.065802, 0.025430, 0.003215, 0.018742 },
    { 0.452171, 0.114613, 0.062460, 0.115702, 0.284246,
      0.140204, 0.100358, 0.550230, 0.143995, 0.700649, 
      0.276580, 0.118569, 0.097470, 0.126673, 0.143634, 
      0.278983, 0.358482, 0.661750, 0.061533, 0.199373 },
    { 0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 
      0.016931, 0.003647, 0.002184, 0.005019, 0.005990, 
      0.001473, 0.004158, 0.009055, 0.003630, 0.006583, 
      0.003172, 0.003690, 0.002967, 0.002772, 0.002686 },
  };

  sstring eddy_alphabet ("acdefghiklmnpqrstvwy");
  
  if (alphabet.size() != vars(0)) THROWEXPR ("Sjolander prior expects a 20-residue alphabet");
  for (int cpt = 0; cpt < components(); ++cpt)
    {
      log_cpt_prior[cpt] = Prob2Nats(defmq[cpt]);
      for (int a = 0; a < vars(0); ++a)
	{
	  const int b = alphabet.char2int_deg(eddy_alphabet[a]);
	  if (b < 0) THROWEXPR ("Non-protein alphabet passed to Sjolander prior");
	  alpha[cpt][0][b] = defm[cpt][a];
	}
    }
}

Alphabet_prior::Alphabet_prior (const Alphabet& _alphabet, int components)
  : Dirichlet_mixture (_alphabet.size(), components),
    alphabet (_alphabet)
{ }
