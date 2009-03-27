#include <iomanip>

#include "kimono/secluster.h"

SE_cluster_model::SE_cluster_model (const Gibbs_alignment_factory& gibbs_factory, const Gaussian_cluster_factory& gauss_factory, int k)
  : gibbs_factory (gibbs_factory),
    gauss_factory (gauss_factory),
    k (k),
    gibbs (k, (Gibbs_alignment*) 0),
    gauss (k, (Gaussian_cluster*) 0),
    model_prior_sc (k+1, (int) Prob2Score (1.0 / (double) (k+1))),
    multi_sc (gauss_factory.params.table.probes(), vector<int> (k+1, (int) 0)),
    gibbs_sc (gauss_factory.params.table.probes(), vector<int> (k+1, (int) 0)),
    gauss_sc (gauss_factory.params.table.probes(), vector<int> (k+1, (int) 0)),
    gibbs_param_sc (k, 0),
    gauss_param_sc (k, 0)
{
  for (int m = 0; m < k; ++m)
    {
      gibbs[m] = gibbs_factory.new_model();
      gauss[m] = gauss_factory.new_model();
    }
}

SE_cluster_model::SE_cluster_model (const SE_cluster_model& model)
  : gibbs_factory (model.gibbs_factory),
    gauss_factory (model.gauss_factory),
    k (model.k),
    gibbs (model.k, (Gibbs_alignment*) 0),
    gauss (model.k, (Gaussian_cluster*) 0),
    model_prior_sc (model.model_prior_sc),
    multi_sc (model.multi_sc),
    gibbs_sc (model.gibbs_sc),
    gauss_sc (model.gauss_sc),
    gibbs_param_sc (model.k, 0),
    gauss_param_sc (model.k, 0)
{
  for (int m = 0; m < k; ++m)
    {
      gibbs[m] = new Gibbs_alignment (*model.gibbs[m]);
      gauss[m] = new Gaussian_cluster (*model.gauss[m]);
    }
}

SE_cluster_model::~SE_cluster_model()
{
  for (int m = 0; m < k; ++m) { delete gibbs[m]; delete gauss[m]; }
}

SE_cluster_model& SE_cluster_model::operator= (const SE_cluster_model& model)
{
  k = model.k;
  gibbs = model.gibbs;
  gauss = model.gauss;
  model_prior_sc = model.model_prior_sc;
  multi_sc = model.multi_sc;
  gibbs_sc = model.gibbs_sc;
  gauss_sc = model.gauss_sc;
  gibbs_param_sc = model.gibbs_param_sc;
  gauss_param_sc = model.gauss_param_sc;
  
  for (int m = 0; m < k; ++m)
  {
    gibbs[m] = new Gibbs_alignment (*model.gibbs[m]);
    gauss[m] = new Gaussian_cluster (*model.gauss[m]);
  }
  
  return *this;
}

void SE_cluster_model::read_model_prior_from_file (const char* model_prior_file)
{
  CLOG(7) << "Reading model prior from file '" << model_prior_file << "'\n";
  ifstream model_prior_stream (model_prior_file);
  if (!model_prior_stream) { sstring e; e << "Model prior file '" << model_prior_file << "' not found"; THROW Standard_exception (e); }
  sstring s;
  s.getline (model_prior_stream);
  for_tmp_contents (vector<sstring>, s.split(), f) model_prior_sc.push_back (Bits2Score (atof ((*f).c_str())));
  k = model_prior_sc.size();
  NormaliseSc (model_prior_sc);
}

double SE_cluster_model::log_gene_model_likelihood (int g, int m) const
{
  return Score2Nats (ScorePMul (gibbs_sc[g][m], gauss_sc[g][m]));
}

double SE_cluster_model::log_gene_likelihood (int g) const
{
  double log_sum_model_likelihoods = -InfinityLoge;
  for (int m = 0; m <= k; ++m)
    NatsPSumAcc (log_sum_model_likelihoods, NatsPMul (Score2Nats (model_prior_sc[m]), log_gene_model_likelihood (g, m)));
  return log_sum_model_likelihoods;
}

double SE_cluster_model::log_likelihood() const
{
  double log_product_gene_likelihoods = 0;
  for (int g = 0; g < genes(); ++g)
    NatsPMulAcc (log_product_gene_likelihoods, log_gene_likelihood(g));
  return log_product_gene_likelihoods;
}

double SE_cluster_model::log_parameter_prior() const
{
  double log_product_parameter_priors = 0;
  for (int m = 0; m < k; ++m)
    NatsPMulAcc (log_product_parameter_priors, Score2Nats (ScorePMul (gibbs_param_sc[m-1], gauss_param_sc[m-1])));
  return log_product_parameter_priors;
}

double SE_cluster_model::log_probability() const
{
  return NatsPMul (log_likelihood(), log_parameter_prior());
}

double SE_cluster_model::log_gene_model_posterior (int g, int m) const
{
  return NatsPMul (NatsPMul (Score2Nats (model_prior_sc[m]), log_gene_model_likelihood (g, m)), -log_gene_likelihood (g));
}

void SE_cluster_model::assign_gene_to_model (int g, int m)
{
  CLOG(5) << "Assigning gene '" << gibbs_factory.dataset.row_name[g] << "' to model #" << m << "\n";
  for (int n = 0; n <= k; ++n)
    multi_sc[g][n] = (m == n) ? 0 : -InfinityScore;
}

void SE_cluster_model::assign_genes_to_random_models()
{
  vector<double> model_prior = Score2ProbVecNorm (model_prior_sc);
  for (int g = 0; g < genes(); ++g)
    {
      int m = 0;
      while (m == 0) m = Rnd::choose (model_prior);
      assign_gene_to_model (g, m);
    }
  update_model_multiplicities();
}

void SE_cluster_model::read_gff_motifs (const char* gff_filename)
{
  ifstream initgff_s (gff_filename);
  sstring g;
  int n = 0;
  vector<int> seen_model (k, (int) 0);
  while (!initgff_s.eof()) {
    g.getline (initgff_s);
    g.chomp();
    ++n;
    vector<sstring> f = g.split("\t");
    if (f.size() == 0) continue;
    if (f.size() < 9) { CLOGERR << "Line " << n << " of '" << gff_filename << "' doesn't have 9 fields; skipping\n"; continue; }
    int row = gibbs_factory.dataset.row_index (f[0]);
    if (row >= 0) {
      int m = atoi (f[2].c_str());
      if (m < 0 || m >= k) { CLOGERR << "Line " << n << " of '" << gff_filename << "': model number out of range; skipping\n"; continue; }
      if (!seen_model[m]) {  // evict randomly allocated genes from this model if this is the first time we've encountered it
	for (int i = 0; i < genes(); ++i)
	  {
	    ScorePSumAcc (multi_sc[i][0], multi_sc[i][m+1]);
	    multi_sc[i][m+1] = -InfinityScore;
	  }
	update_model_multiplicities();
	gibbs[m]->motif_start = 0;
	gibbs[m]->motif_end = gibbs[m]->motif_length - 1;
	gibbs[m]->recalculate_counts();
	seen_model[m] = 1;
      }
      bool rev = f[6] == sstring("-");
      int start;
      if (rev)
	start = gibbs_factory.dataset.profile[row]->size() - atoi (f[4].c_str());
      else
	start = atoi (f[3].c_str()) - 1;
      if (CLOGGING(5))
	CL << "read_gff_motifs: sequence '" << f[0] << "', " << (rev ? "reverse" : "forward") << " strand, offset " << start << ", motif '" << f[2] << "', model #" << m+1 << "\n";
      assign_gene_to_model (row, m+1);
      update_model_multiplicities();
      gibbs[m]->align_row (row, -start, rev);
    }
  }
  optimise_gauss("read_gff_motifs: ");
  if (CLOGGING(3)) display_multiplicities(CL);
  if (CLOGGING(6)) {
    for (int i = 0; i < k; ++i) {
      CL << "Sequence model #" << i+1 << "\n";
      gibbs[i]->display (CL);
    }
  }
}

void SE_cluster_model::calculate_expected_multiplicities (double& rms_delta_ret)
{
  double sum_delta_squared = 0.0;
  int n = 0;
  for (int g = 0; g < genes(); ++g)
    for (int m = 0; m <= k; ++m)
      {
	const double old_multi_sc = multi_sc[g][m];
	multi_sc[g][m] = Nats2Score (log_gene_model_posterior (g, m));
	const double delta = multi_sc[g][m] - old_multi_sc;
	sum_delta_squared += delta * delta;
	++n;
      }
  rms_delta_ret = n > 0 ? sqrt (sum_delta_squared / n) : 0.0;
}

void SE_cluster_model::update_model_multiplicities()
{
  for (int g = 0; g < genes(); ++g)
    for (int m = 1; m <= k; ++m)
      {
	gibbs[m-1]->set_membership_probability (g, multi_sc[g][m]);
	gauss[m-1]->set_membership_probability (g, multi_sc[g][m]);
      }
}

void SE_cluster_model::sample_alignments (double kT, double gibbs_resample_period)
{
  const double effectively_zero = .0000001;

  vector<vector<double> > alignment_sample_prob;
  for (int g = 0; g < genes(); ++g)
    {
      alignment_sample_prob.push_back (Score2ProbVecUnnorm (multi_sc[g]));
      for_contents (vector<double>, alignment_sample_prob.back(), i) *i += 1.0 / max (gibbs_resample_period, effectively_zero);
    }
  
  for (int j = 0; j < k; ++j)
    {
      vector<int> sequence_order (genes());
      for (int g = 0; g < genes(); ++g) sequence_order[g] = g;
      for (int i = 0; i < (int) sequence_order.size(); ++i)
	{
	  swap (sequence_order[i], sequence_order [Rnd::rnd_int (sequence_order.size() - i) + i]);
	  const int g = sequence_order[i];
	  const int m = gibbs_resample_period <= effectively_zero ? j+1 : Rnd::choose (alignment_sample_prob[g]);
	  if (m != 0) {
	    if (CLOGGING(5))
	      CL << "Sampling alignment of gene '" << gibbs_factory.dataset.row_name[g] << "' to model #" << m << "\n";
	    gibbs[m-1]->sample_row (g, kT);
	  }
	}
    }
  if (CLOGGING(4)) { CL << "After sampling alignments:\n"; display_terse(CL); }
}

void SE_cluster_model::sample_gibbs_indent (double kT, const char* log_prefix)
{
  for (int m = 0; m < k; ++m)
    {
      CLOG(4) << log_prefix << "Sampling motif indentation then optimising parameters for sequence model #" << m+1 << "\n";
      gibbs[m]->sample_motif (kT);
      gibbs[m]->flush_left();
      if (CLOGGING(5))
	{
	  CL << "New sequence model #" << m+1 << ":\n";
	  gibbs[m]->display (CL);
	}
    }
}

void SE_cluster_model::optimise_gibbs_constrained()
{
  for (int m = 0; m < k; ++m)
    gibbs[m]->flush_left();
}

void SE_cluster_model::optimise_gauss (const char* log_prefix)
{
  for (int m = 0; m < k; ++m)
    {
      CLOG(4) << log_prefix << "Optimising parameters for expression model #" << m+1 << "\n";
      gauss[m]->optimise_parameters();
      if (CLOGGING(5))
	{
	  CL << "New expression model #" << m+1 << ":\n";
	  gauss[m]->display (CL);
	}
    }
}

void SE_cluster_model::update_gibbs_scores()
{
  for (int g = 0; g < genes(); ++g)
    {
      // null model scores
      //
      gibbs_sc[g][0] = gibbs_factory.null_model.profile_null_score[g];

      // non-null models
      //
      for (int m = 1; m <= k; ++m)
	gibbs_sc[g][m] = gibbs[m-1]->log_likelihood_sc (g);
    }

  // parameter priors
  //
  for (int m = 0; m < k; ++m)
    gibbs_param_sc[m] = gibbs[m]->log_parameter_prior_sc();
}

void SE_cluster_model::update_gauss_scores (double kappa)
{
  for (int g = 0; g < genes(); ++g)
    {
      // null model scores
      //
      gauss_sc[g][0] = gauss_factory.null.profile_null_score (g, kappa);

      // non-null models
      //
      for (int m = 1; m <= k; ++m)
	gauss_sc[g][m] = gauss[m-1]->log_likelihood_sc (g, kappa);
    }

  // parameter priors
  //
  for (int m = 0; m < k; ++m)
    gauss_param_sc[m] = gauss[m]->log_parameter_prior_sc (kappa);
}

void SE_cluster_model::greedy_gauss (const char* log_prefix, double kappa, double min_rms_delta, int max_iter)
{
  bool converged = 0;
  for (int i = 0; i < max_iter; ++i)
    {
      optimise_gibbs_constrained();
      optimise_gauss (log_prefix);
      update_scores (kappa);

      double rms_delta = 0.0;
      calculate_expected_multiplicities (rms_delta);
      update_model_multiplicities();

      if (CLOGGING(4)) { CL << "greedy_gauss iteration " << i << ":\n"; display_terse(CL); display_scores (CL, kappa); display_multiplicities(CL); }

      if (rms_delta < min_rms_delta)
	{
	  CLOG(4) << log_prefix << "Greedy optimisation converged in " << i+1 << " steps (rms delta = " << rms_delta << ")\n";
	  converged = 1;
	  break;
	}
    }
  if (!converged) CLOG(4) << log_prefix << "Greedy optimisation did not converge after " << max_iter << " steps\n";
}

void SE_cluster_model::display_verbose (ostream& o)
{
  for (int m = 0; m < k; ++m)
    {
      o << "-----------------------------------------------------------------------------------------------------------------\n";
      o << "Model #" << m+1 << "\n";
      o << "-----------------------------------------------------------------------------------------------------------------\n";
      gibbs[m]->display (o);
      gauss[m]->display (o);
    }
  o << "-----------------------------------------------------------------------------------------------------------------\n";
  o << "Log-probability = " << Nats2Bits (log_probability()) << " bits\n";
  o << "-----------------------------------------------------------------------------------------------------------------\n";
}

void SE_cluster_model::display_terse (ostream& o)
{
  for (int m = 0; m < k; ++m)
    {
      o << "Model #" << m+1 << ":\n";
      gibbs[m]->display (o);
      gauss[m]->display (o);
    }
}

void SE_cluster_model::display_scores (ostream& o, double kappa)
{
  Gibbs_alignment::display_all_scores (o, gibbs_factory.null_model, gibbs);
  Gaussian_cluster::display_all_scores (o, gauss_factory.null, gauss, kappa);
}


void SE_cluster_model::display_mixture_consensi (ostream& o)
{
  Gibbs_alignment::display_all_mixture_consensi (o, gibbs);
}


void SE_cluster_model::display_alphabet_consensi (ostream& o)
{
  Gibbs_alignment::display_all_alphabet_consensi (o, gibbs);
}


void SE_cluster_model::display_multiplicities (ostream& o)
{
  int old_prec = o.precision(3);
  o.precision(2);
  save_flags (o);
  left_align (o);

  for (int g = 0; g < genes(); ++g)
    {
      sstring n = gibbs_factory.dataset.row_name[g].substr(0,10);
      n << "':";
      o << "Model posteriors for '" << setw(13) << n;
      for (int m = 0; m <= k; ++m)
	{
	  o << " P(";
	  if (m > 0) o << m; else o << "null";
	  o << ")=";
	  o << setw(5) << Score2Prob (multi_sc[g][m]);
	}
      o << "\n";
    }

  restore_flags (o);
  o.precision (old_prec);
}

void SE_cluster_model::save (ostream& o)
{
  for (int g = 0; g < genes(); ++g)
    {
      o << gibbs_factory.dataset.row_name[g];
      for (int m = 0; m <= k; ++m)
	o << ' ' << multi_sc[g][m];
      for (int m = 0; m < k; ++m)
	{
	  o << ' ';
	  int offset = gibbs[m]->row_offset[g];
	  if (offset < 0 || !gibbs[m]->row_aligned[g])
	    o << '*';
	  else
	    {
	      if (gibbs[m]->row_reversed[g])
		o << '-' << offset;
	      else
		o << '+' << offset;
	    }
	}
      o << "\n";
    }
  for (int m = 0; m < k; ++m)
    o << gibbs[m]->motif_start << ' ' << gibbs[m]->motif_end << "\n";
}

void SE_cluster_model::load (istream& i)
{
  for (int g = 0; g < genes(); ++g)
    for (int m = 0; m < k; ++m)
      gibbs[m]->unalign_row (g);
  for (int g = 0; g < genes(); ++g)
    {
      sstring s;
      s.getline (i);
      vector<sstring> f = s.split();
      if ((int) f.size() != 2*k+2 || f[0] != gibbs_factory.dataset.row_name[g])
	THROW Standard_exception ("Bad input format for SE cluster model");
      for (int m = 0; m <= k; ++m)
	multi_sc[g][m] = atoi (f[m+1].c_str());
      for (int m = 0; m < k; ++m)
	{
	  sstring& offset = f[m+k+2];
	  if (offset != "*")
	    gibbs[m]->align_row (g, abs (atoi (offset.c_str())), offset[0] == '-');
	}
    }
  update_model_multiplicities();
  for (int m = 0; m < k; ++m)
    {
      sstring s;
      s.getline (i);
      vector<sstring> f = s.split();
      gibbs[m]->motif_start = atoi (f[0].c_str());
      gibbs[m]->motif_end = atoi (f[1].c_str());
      gibbs[m]->recalculate_counts();
    }
}
