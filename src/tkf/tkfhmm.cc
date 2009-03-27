#include "tkf/tkfhmm.h"
#include "util/math_fn.h"

TKF_joint_pair_HMM_scores::TKF_joint_pair_HMM_scores (const TKF_params& params, double time, bool use_indels, bool use_subst)
  :  Pair_HMM_scores (3, &params.submat_factory.alphabet())
{
  time = max(time,TINY);

  state_type[0] = EmitY;
  state_type[1] = EmitX;
  state_type[2] = EmitXY;

  TKF_branch_transition_scores tb (params, time);
  if (use_indels)
    {
      for (int s = 0; s < 3; s++)
	{
	  start[s] = tb.joint_pair_HMM_transition_score(0,1+s);
	  end[s] = tb.joint_pair_HMM_transition_score(1+s,4);
	  for (int t = 0; t < 3; t++)
	    transition(s,t) = tb.joint_pair_HMM_transition_score(1+s,1+t);
	}
      start_to_end() = tb.joint_pair_HMM_transition_score(0,4);
    }
  else
    {
      for (int s = 0; s < 3; s++)
	{
	  start[s] = end[s] = 0;
	  for (int t = 0; t < 3; t++)
	    transition(s,t) = 0;
	}
      start_to_end() = 0;
    }

  const Alphabet& alphabet = params.submat_factory.alphabet();
  if (use_subst)
    {
      single_emit[0] = single_emit[1] = tb.prior;
      pair_emit[2] = Prob2ScoreArray2d (params.submat_factory.create_joint_substitution_matrix (time));
    }
  else
    {
      single_emit[0] = single_emit[1] = vector<Score> (alphabet.size(), (Score) 0);
      pair_emit[2] = array2d<Score> (alphabet.size(), alphabet.size(), (Score) 0);
    }

  if (CTAGGING(4,TKFHMM TKFHMM_SCORES)) { CL << "Pair HMM scores:\n"; show(CL); }
}


TKF_joint_pair_HMM_dt::TKF_joint_pair_HMM_dt (const TKF_params& p, double t, bool use_indels, bool use_subst)
  :  Pair_HMM_derivatives (3, "time", &p.submat_factory.alphabet())
{
  t = max(t,TINY);

  state_type[0] = EmitY;
  state_type[1] = EmitX;
  state_type[2] = EmitXY;

  if (use_indels)
    {
      start[0] = p.ddes_dt(t);
      start[1] = p.eqm() * (p.dndes_dt(t) * p.nanc(t) + p.ndes(t) * p.dnanc_dt(t));
      start[2] = p.eqm() * (p.dndes_dt(t) * p.anc(t) + p.ndes(t) * p.danc_dt(t));

      transition(0,0) = transition(2,0) = p.ddes_dt(t);
      transition(0,1) = transition(2,1) = p.eqm() * (p.dndes_dt(t) * p.nanc(t) + p.ndes(t) * p.dnanc_dt(t));
      transition(0,2) = transition(2,2) = p.eqm() * (p.dndes_dt(t) * p.anc(t) + p.ndes(t) * p.danc_dt(t));

      transition(1,0) = p.dorp_dt(t);
      transition(1,1) = p.eqm() * (p.dnorp_dt(t) * p.nanc(t) + p.norp(t) * p.dnanc_dt(t));
      transition(1,2) = p.eqm() * (p.dnorp_dt(t) * p.anc(t) + p.norp(t) * p.danc_dt(t));
  
      end[0] = p.neqm() * p.dndes_dt(t);
      end[1] = p.neqm() * p.dnorp_dt(t);
      end[2] = p.neqm() * p.dndes_dt(t);

      start_to_end() = p.neqm() * p.dndes_dt(t);
    }
 
  const Alphabet& alphabet = p.submat_factory.alphabet();
  single_emit[0] = single_emit[1] = vector<double> (alphabet.size(), 0.);
  if (use_subst)
    pair_emit[2] = p.submat_factory.differentiate_joint_substitution_matrix(t);
  else
    pair_emit[2] = array2d<double> (alphabet.size(), alphabet.size(), 0.);

  if (CTAGGING(4,TKFHMM)) { CL << "Pair HMM derivatives:\n"; show(CL); }
}


TKF_counts_function::~TKF_counts_function() { }

Pair_HMM_counts TKF_unaligned_counts_function::operator() (const Pair_HMM_scores& scores)
{
  Pair_HMM_counts counts (scores);
  counts.add_counts_from_unaligned_sequences (scores, xseq, yseq);
  return counts;
}

Pair_HMM_counts TKF_aligned_counts_function::operator() (const Pair_HMM_scores& scores)
{
  Pair_HMM_counts counts (scores);
  counts.add_counts_from_aligned_sequences (scores, xseq, yseq, path);
  return counts;
}
