#include <iomanip>
#include "kimono/trues.h"
#include "util/logfile.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  Opts_list opts (argc, argv);
  opts.short_description = "k-means integrated models for oligonucleotide arrays";
  opts.syntax = "[options] <sequence file> <expression file>";
  Rnd::add_opts (opts);
  opts.newline();
  Log_stream::add_opts (opts);

  sstring mu_file;
  sstring sigma_file;
  sstring cov_file;
  double  noise_scale;
  bool    normalise;
  int     motif_len;
  bool    allow_revcomp;
  double  motif_del_prob;
  int     context_order;
  int     k;
  sstring model_prior_file;
  double  p_degen;
  bool    iupac;
  double  iupac_match;
  double  iupac_mismatch;
  int     max_seeds;
  int     max_rounds;
  int     samples_per_round;
  double  gibbs_resample_period;
  double  kT_start;
  double  kT_end;
  double  kappa_start;
  double  kappa_end;
  double  kappa_rise;
  bool    deja_vu_enabled;
  sstring trues_file;
  double  truth_level;
  sstring init_file;
  sstring initgff;
  sstring best_file;

  opts.newline();
  opts.add ("signalmean", mu_file = "",               "file containing signal expectation vector, for each cluster's mean expression profile (default is zero)", 0);
  opts.add ("signalcov",  sigma_file = "",            "file containing signal covariance matrix, for each cluster's mean expression profile (default is identity matrix)", 0);
  opts.add ("noisecov",   cov_file = "",              "file containing noise covariance matrix, for deviations from the mean expression profile within each cluster (default is identity matrix)", 0);
  opts.add ("nscale",     noise_scale = 1,            "scaling factor for noise covariance matrix");
  opts.add ("normalise",  normalise = 1,              "\tnormalise expression profile vectors");
  opts.add ("motiflen",   motif_len = 10,             "motif length");
  opts.add ("revcomp",    allow_revcomp = 1,          "\tallow reverse strand alignments for DNA");
  opts.add ("mdelprob",   motif_del_prob = 1e-60,     "probability of deleting a motif residue");
  opts.add ("k",          k = 10,                     "\tnumber of models (the 'k' in k-means clustering)");
  opts.add ("modelprior", model_prior_file = "",      "file whose first line contains whitespace-delimited prior distribution over models in bits, null model first (overrides -k)", 0);
  opts.add ("context",    context_order = 0,          "order of Markov chain for null model");
  opts.add ("iupac",      iupac = 1,                  "use IUPAC degeneracy classes as Dirichlet mixture prior components");
  opts.add ("iupacmat",   iupac_match = 100,          "pseudocounts for matches to IUPAC classes");
  opts.add ("iupacmis",   iupac_mismatch = 1,         "pseudocounts for mismatches to IUPAC classes");
  opts.add ("pdegen",     p_degen = 1e-10,            "combined probability of degenerate IUPAC classes");
  opts.add ("seeds",      max_seeds = 1,              "\tnumber of times clustering will be (re-)started from an initial random seed");
  opts.add ("rounds",     max_rounds = 100,           "\tnumber of rounds of sampling algorithm");
  opts.add ("smplround",  samples_per_round = 1,      "number of alignment sampling steps per round of sampling algorithm");
  opts.add ("resample",   gibbs_resample_period = 10, "average number of rounds between sampling unfavourable alignments");
  opts.add ("ktstart",    kT_start = 1,               "initial value of kT for simulated annealing (use ktstart=ktend=1 for regular Gibbs sampling)");
  opts.add ("ktend",      kT_end = 1,                 "final value of kT for simulated annealing");
  opts.add ("kappastart", kappa_start = 1,            "initial value of kappa (the expression variance scaling parameter)");
  opts.add ("kappaend",   kappa_end = 1,              "final value of kappa");
  opts.add ("kapparise",  kappa_rise = .1,            "fractional number of sampling rounds before kappaend is reached");
  opts.add ("dejavu",     deja_vu_enabled = 1,        "\tstop if a previously-seen optimum is encountered (only makes sense for multiple -seeds)");
  opts.add ("trues",      trues_file = "",            "file containing true cluster assignments for test purposes (one gene per line, format '<name> <cluster> <motif start> <motif end>')\n\t\t\t (NB residue numbering starts at 1, not 0)", 0);
  opts.add ("truthlevel", truth_level = .9,           "for a model to be 'true', all the assignments must have posterior probabilities greater than (truthlevel)\n\t\t\t and, additionally, a proportion (truthlevel) of these must be correctly aligned");
  opts.add ("initfile",   init_file = "",             "file to load initial clusters from", 0);
  opts.add ("initgff",    initgff = "",               "file to load initial motif positions from", 0);
  opts.add ("bestfile",   best_file = "",             "file to write best intermediate results into", 0);
  
  opts.newline();
  opts.add ("exprformathelp", &Expr_table::format_help,              "\tdescribe expression file format");
  opts.add ("gpformathelp",   &Gaussian_cluster_params::format_help, "\t\tdescribe Gaussian prior file format");

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 2) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // get parameters

      sstring seq_file = opts.args[0];
      sstring expr_file = opts.args[1];

      // read in the sequences; detect alphabet; digitise; revcomp if allowed

      FASTA_sequence_database seq_db (seq_file.c_str());
      const Alphabet& alphabet = seq_db.detect_alphabet();
      seq_db.seqs2dsqs (alphabet);

      Sequence_database revcomp_db;
      if (alphabet.has_complement() && allow_revcomp) revcomp_db = seq_db.revcomp (alphabet);

      // set up database indices

      Sequence_database_index db_index (seq_db);
      Sequence_database_index revcomp_index (revcomp_db);

      // read in the expression vectors

      Scaled_expr_table expr_tab (expr_file.c_str());
      if (normalise) expr_tab.normalise();

      // set up Dirichlet mixture prior for motifs

      int degen_classes = iupac ? (1 << alphabet.size()) - 1 : 1;

      vector<double> log_dirichlet_component_weight (degen_classes, -InfinityLoge);
      vector<vector<double> > dirichlet_mixture (degen_classes, vector<double> (alphabet.size()));
      vector<sstring> dirichlet_component_name (degen_classes, sstring());

      if (iupac) {

	double log_p_degen_indiv = log (p_degen / ((double) (degen_classes - alphabet.size())));
	double log_p_nondegen_indiv = log ((1 - p_degen) / ((double) alphabet.size()));
	for (int c = 0; c < degen_classes; ++c) {
	  int members_of_class = 0;
	  for (int b = 0; b < alphabet.size(); ++b) {
	    const bool in_class = (c+1) & (1<<b);
	    if (in_class) { ++members_of_class; dirichlet_component_name[c] << alphabet.int2char(b); }
	  }
	  for (int b = 0; b < alphabet.size(); ++b) {
	    const bool in_class = (c+1) & (1<<b);
	    dirichlet_mixture[c][b] = iupac_mismatch + (in_class ? iupac_match / members_of_class : 0);
	  }
	  log_dirichlet_component_weight[c] = members_of_class==1 ? log_p_nondegen_indiv : log_p_degen_indiv;
	}

      } else {
	
	log_dirichlet_component_weight[0] = 0;
	dirichlet_mixture[0] = vector<double> (alphabet.size(), iupac_mismatch);
	dirichlet_component_name[0] = "*";

      }
      
      if (CLOGGING(5))
	for (int c = 0; c < degen_classes; ++c) {
	  CL << "Mixture component #" << c << " " << dirichlet_component_name[c] << ": ";
	  CL << "log_dirichlet_component_weight = " << log_dirichlet_component_weight[c] << " ";
	  CL << "dirichlet_mixture = (" << dirichlet_mixture[c] << ")\n";
	}
      
      // set up the Gibbs_alignment_factory
      
      Gibbs_dataset     gibbs_dataset (db_index, revcomp_index, expr_tab.probe_name, alphabet);
      Gibbs_null_model  gibbs_null (gibbs_dataset, alphabet, 1, motif_del_prob, context_order);

      Gibbs_alignment_factory  gibbs_factory (gibbs_dataset, alphabet, gibbs_null, log_dirichlet_component_weight, dirichlet_mixture, dirichlet_component_name, motif_len);
      
      // set up the Gaussian_cluster_factory

      Gaussian_cluster_params gauss_params (expr_tab);
      gauss_params.read (mu_file, sigma_file, cov_file, noise_scale);

      Gaussian_cluster_null_model gauss_null (expr_tab, gauss_params);
      Gaussian_cluster_factory    gauss_factory (gauss_params, gauss_null);

      // need some variables to store best models to date

      SE_cluster_model best_se (gibbs_factory, gauss_factory, k);
      if (model_prior_file != "")
	best_se.read_model_prior_from_file (model_prior_file.c_str());
      
      double best_loglike = -InfinityLoge;
      int best_seed = -1;

      // read in trues file (if there is one)

      True_cluster_set trues (trues_file, gibbs_factory, gauss_factory, best_se.model_prior_sc);
      
      CLOG(7) << "Ready to start.\n";
      
      // start the outer seeding loop
      
      for (int seed = 0; seed < max_seeds; ++seed)
	{
	  // initialisation
	  
	  double best_seed_loglike = -InfinityLoge;

	  SE_cluster_model se (gibbs_factory, gauss_factory, k);
	  se.model_prior_sc = best_se.model_prior_sc;         // ugh
	  
	  CLOG(7) << "Seed #" << seed+1 << ": ";
	  
	  if (init_file != "") {
	    
	    CL << "reading initial clusters from file '" << init_file << "'\n";
	    
	    ifstream init (init_file.c_str());
	    se.load (init);
	    
	  } else {

	    CL << "assigning each gene to a random cluster\n";
	    
	    se.assign_genes_to_random_models();
	    if (initgff != "") se.read_gff_motifs (initgff.c_str());
	    else se.greedy_gauss ("", kappa_start);
	  }

	  // Sampling loop is as follows:
	  // (....THIS DESCRIPTION NOW OUT OF DATE....)
	  //  +->  For each model
	  //  | |->  For each gene (taken in a random order)
	  //  | | |->  Pick a model at random according to the model membership probabilities plus (1 / gibbs_resample_period) `pseudocounts'
	  //  | | |->  Optimise sequence model with this alignment removed (implicit step)
	  //  | | \->  Sample alignment of sequence to model
	  //  | |->  Optimise sequence model for all alignments
	  //  | |->  Optimise expression model
	  //  | \->  For each gene
	  //  |   \->  Calculate sequence & expression likelihoods
	  //  |->  For each gene
	  //  | \->  Estimate new model membership probabilities
	  //  |->  Calculate expected total log-likelihood
	  //  \->  If best score so far, save all parameters
	  //
	  //  In summary: (i) sample alignments; (ii) optimise models; (iii) reestimate membership probabilities; (iv) remember if best.

	  // (i) & (ii): sample alignments, optimise models

	  for (int round = 0; round < max_rounds; ++round)         // main sampling loop
	    {
	      double round_frac = (double) round / (double) max_rounds;
	      double kT = kT_start + round_frac * (kT_end - kT_start);
	      double kappa = kappa_start + max (round_frac / kappa_rise, (double) 1) * (kappa_end - kappa_start);
	      
	      sstring round_str;
	      round_str << "Round " << seed+1 << "/" << round+1 << ": ";
	      CLOG(6) << round_str << "Sampling sequence alignments  (kT = " << kT << ", kappa = " << kappa << ")\n";

	      for (int sample = 0; sample < samples_per_round; ++sample)
		se.sample_alignments (kT, gibbs_resample_period);

	      se.sample_gibbs_indent (kT, round_str.c_str());

	      if (CLOGGING(4)) se.display_scores (CL, kappa);
	      
	      // (iii): reestimate membership probabilities; calculate expected total log-likelihood
	      
	      CLOG(6) << round_str << "Estimating model membership probabilities P(m|g) and overall likelihood P({g},{theta})\n";
	      
	      se.greedy_gauss (round_str.c_str(), kappa);
	      if (CLOGGING(6)) { if (iupac) se.display_mixture_consensi (CL); else se.display_alphabet_consensi (CL); }
	      if (CLOGGING(3)) se.display_scores (CL, kappa);
	      if (CLOGGING(4)) se.display_multiplicities (CL);

	      double loglike = se.log_probability();
	      
	      CLOG(7) << round_str << "Log-likelihood = " << Nats2Bits (loglike) << " bits, seed best = " << Nats2Bits (best_seed_loglike) << ", overall best = " << Nats2Bits (best_loglike) << " bits\n";
	      
	      // SHOULD resample expression profile normalisation parameters here
	      
	      // (iv): store models if best so far

	      if (loglike > best_seed_loglike)
		{
		  best_seed_loglike = loglike;
		  CLOG(7) << round_str << "Best score for this seed; log-likelihood = " << Nats2Bits (loglike) << " bits\n";
		}

	      if (loglike > best_loglike)
		{
		  best_se = se;
		  best_loglike = loglike;
		  best_seed = seed;
		  CLOG(7) << round_str << "New hi-score! Log-likelihood = " << Nats2Bits (loglike) << " bits\n";
		  if (best_file != "")
		    {
		      ofstream best_out (best_file.c_str());
		      best_se.save (best_out);

		      sstring disp_file = best_file;
		      disp_file << ".display";
		      ofstream disp_out (disp_file.c_str());
		      best_se.display_verbose (disp_out);
		    }
		}

	      // compare with true clusters, if we know them
	      
	      if (trues.true_enough (se.gibbs, truth_level)) { CLOG(7) << "Matched trues set specified in '" << trues_file << "'; stopping\n"; break; }
	      
	      // end of sampling loop
	    }

	  // clean up
	  //
	  if (deja_vu_enabled && abs(best_seed_loglike-best_loglike) < .0001 && seed != best_seed && best_loglike > -InfinityLoge)
	    {
	      CLOG(7) << "Deja vu! Got to this score from seed #" << best_seed+1 << ". Time to retire...\n";
	      break;
	    }
	  
	  // end of seeding loop
	}

      if (best_loglike > -InfinityLoge)
	best_se.display_verbose (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
