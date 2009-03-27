#include "kimono/gibbs.h"
#include "kimono/gp.h"
#include "util/logfile.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  Opts_list opts (argc, argv, CLOGERR);
  opts.short_description = "integrated sequence / expression modelling software";
  opts.syntax = "[options] <sequence file> <expression file>";
  Rnd::add_opts (opts);
  opts.newline();
  Log_stream::add_opts (opts);

  sstring mu_file;
  sstring sigma_file;
  sstring cov_file;
  bool    normalise;
  int     motif_len;
  bool    allow_revcomp;
  sstring k_prior_file;
  int     max_models;
  int     max_seeds;
  int     max_samples;
  bool    deja_vu_enabled;

  opts.newline();
  opts.add ("signalmean", mu_file = "",        "file containing expectation vector for each cluster's mean expression profile (default is zero)", 0);
  opts.add ("signalcov",  sigma_file = "",     "file containing covariance matrix for each cluster's mean expression profile (default is identity matrix)", 0);
  opts.add ("noisecov",   cov_file = "",       "file containing covariance matrix for deviations from the mean expression profile within each cluster (default is identity matrix)", 0);
  opts.add ("normalise",  normalise = 1,       "\tnormalise expression profile vectors");
  opts.add ("motiflen",   motif_len = 10,      "\texpected prior length for ungapped motif");
  opts.add ("revcomp",    allow_revcomp = 1,   "\tallow reverse strand alignments for DNA");
  opts.add ("modelprior", k_prior_file = "",   "file containing prior distribution for number of models (one line of whitespace-separated probabilities; overrides -maxmodels)");
  opts.add ("models",     max_models = 10,     "number of models (the 'k' in k-means clustering)");
  opts.add ("seeds",      max_seeds = 1,       "\tnumber of times clustering will be (re-)started from an initial random seed");
  opts.add ("samples",    max_samples = 100,   "\tnumber of rounds of sampling algorithm");
  opts.add ("dejavu",     deja_vu_enabled = 1, "\t\tstop if a previously-seen optimum is encountered (only makes sense for multiple -seeds)");
  
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
      seq_db.clear_seqs();

      Sequence_database revcomp_db;
      if (alphabet.has_complement() && allow_revcomp) revcomp_db = seq_db.revcomp (alphabet);

      // set up database indices

      Sequence_database_index db_index (seq_db);
      Sequence_database_index revcomp_index (revcomp_db);

      // read in the expression vectors

      Stanford_expression_profile_table expr_tab (expr_file.c_str());
      if (normalise) expr_tab.normalise();

      // set up the Gibbs_alignment_factory
      
      Gibbs_dataset     gibbs_dataset (db_index, revcomp_index, expr_tab.sequence_name, alphabet);
      vector<double>    pseudocounts (alphabet.size(), (double) 1);
      Gibbs_null_model  gibbs_null (gibbs_dataset, alphabet, pseudocounts);
      Gibbs_alignment_factory  gibbs_factory (gibbs_dataset, alphabet, gibbs_null, pseudocounts, motif_len);
      
      // set up the Gaussian_cluster_factory

      Gaussian_cluster_params gauss_params (expr_tab);
      gauss_params.read (mu_file, sigma_file, cov_file);

      Gaussian_cluster_null_model gauss_null (expr_tab);
      Gaussian_cluster_factory    gauss_factory (gauss_params, gauss_null);

      // set up the prior distribution for the number of models (call this k, as in k-means clustering)
      // NB the size of the k_prior vector is (max_models + 1) to account for the probability that there are no models (except the null model)
      
      vector<double> k_prior;
      if (k_prior_file = "")
	{
	  CLOG(7) << "Using default prior for number of models (flat)\n";
	  k_prior = vector<double> (max_models + 1, 1.0 / (double) max_models);
	  k_prior[0] = 0;
	}
      else
	{
	  CLOG(7) << "Reading number-of-models prior from file '" << k_prior_file << "'\n";
	  ifstream k_prior_stream (k_prior_file.c_str());
	  sstring s;
	  s.getline (k_prior_stream);
	  for_tmp_contents (vector<sstring>, s.split(), f) k_prior.push_back (atof ((*f).c_str()));
	  max_models = k_prior.size();

	  double k_sum = accumulate (k_prior.begin(), k_prior.end(), 0.0);
	  if (k_sum > 1.001) { sstring e; e << "Vector in file '" << k_prior_file << "' is not a probability distribution"; THROW Standard_exception (e); }
	  k_prior.insert (k_prior.begin(), max (1.0 - k_sum, 0.0));
	}
      vector<int> log_k_prior = Prob2Score_vector (k_prior);
      
      // need some variables to store best models to date
      
      vector<Gibbs_alignment*>  best_gibbs (max_models, (Gibbs_alignment*) 0);
      vector<Gaussian_cluster*> best_gauss (max_models, (Gaussian_cluster*) 0);

      int  best_score = -InfinityScore;
      int  best_seed = -1;

      CLOG(7) << "Ready to start.\n";
      
      // start the outer seeding loop
      
      for (int seed = 0; seed < max_seeds; ++seed)
	{
	  CLOG(7) << "Seed #" << seed+1 << ": assigning each gene to a random cluster\n";

	  int best_seed_score = -InfinityScore;
	  
	  // initialisation - create (max_models) models
	  
	  vector<Gibbs_alignment*>  gibbs (max_models);
	  vector<Gaussian_cluster*> gauss (max_models);
	  
	  for (int m = 0; m < max_models; ++m)
	    {
	      gibbs[m] = gibbs_factory.new_model();
	      gauss[m] = gauss_factory.new_model();
	    }
	  
	  // assign each gene randomly to a model, using P(m) = sum_{k=m}^{k=k_max} Prior(k) / (k+1)
	  
	  vector<double> initial_model_probability (max_models + 1, (double) 0);
	  for (int m = 0; m <= max_models; ++m)
	    for (int k = m; k <= max_models; ++k)
	      initial_model_probability[m] += k_prior[k] / (k + 1);
	  
	  for (int g = 0; g < expr_tab.sequences(); ++g)
	    {
	      int initial_model = Rnd::choose (initial_model_probability);
	      for (int m = 0; m < max_models; ++m)
		if (initial_model == m+1)
		  {
		    gibbs[m]->set_membership_probability (g, 0);
		    gauss[m]->set_membership_probability (g, 0);
		  }
		else
		  {
		    gibbs[m]->set_membership_probability (g, -InfinityScore);
		    gauss[m]->set_membership_probability (g, -InfinityScore);
		  }
	    }
	  
	  // Sampling loop is as follows:
	  //  +->  For each model
	  //  | |->  For each gene (taken in a random order)
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

	  for (int sample = 0; sample < max_samples; ++sample)         // main sampling loop
	    {
	      sstring round;
	      round << "Round " << seed+1 << "/" << sample+1 << ": ";
	      CLOG(6) << round << "Sampling sequence alignments\n";
	      
	      for (int m = 0; m < max_models; ++m)
		{
		  CLOG(4) << round << "Sampling alignments to sequence model #" << m+1 << "\n";
		  vector<int> sequence_order (expr_tab.sequences());
		  for (int g = 0; g < expr_tab.sequences(); ++g) sequence_order[g] = g;
		  for (int i = 0; i < sequence_order.size(); ++i)
		    {
		      swap (sequence_order[i], sequence_order [Rnd::rnd_int (sequence_order.size() - i) + i]);
		      const int g = sequence_order[i];
		      gibbs[m]->sample_row (g, 0, 1.0);
		    }

		  CLOG(4) << round << "Optimising parameters for sequence model #" << m+1 << "\n";
		  gibbs[m]->find_MAP_length_motif();
		  if (CLOGGING(5))
		    {
		      CL << "New sequence model #" << m+1 << ":\n";
		      gibbs[m]->display (CL);
		      if (CLOGGING(2)) gibbs[m]->display_scores (CL);
		    }

		  CLOG(4) << round << "Optimising parameters for expression model #" << m+1 << "\n";
		  gauss[m]->optimise_parameters();
		  if (CLOGGING(5))
		    {
		      CL << "New expression model #" << m+1 << ":\n";
		      gauss[m]->display (CL);
		      if (CLOGGING(2)) gauss[m]->display_scores (CL);
		    }
		}
	      
	      // (iii): reestimate membership probabilities; calculate expected total log-likelihood
	      //
	      // relevant formulae are:  (NB {g} = all genes, {theta} = all model params, g = individual gene, theta_m = params for model m)
	      //
	      // P(m|g) = sum_{k>=m} P(m|g,k) . P(k)                                                                           (implicitly all conditioned on {theta}; P(k) = prior for k)
	      //        = sum_{k>=m} P(g|m) . P(m|k) / ( sum_{m'<=k} P(g|m') . P(m'|k) )
	      //        = P(g|m) . sum_{k>=m} P(k) / sum_{m'<=k} P(g|m')                                                       (assuming P(m|k) does not depend on m, i.e. P(m|k) = 1/(k+1))
	      //
	      // P({g},{theta}) = P({theta}) . P({g} | {theta})
	      //                = (\prod_{m} P(theta_m)) . \sum_k ( P(k) . \prod_g \sum_{m<=k} ( P(m|k) . P(g|m,theta_m) ) )
	      //                = (\prod_{m} P(theta_m)) . \sum_k ( (P(k)/k) . \prod_g \sum_{m<=k} P(g|m,theta_m) )            (if P(m|k) = 1/k)
	      //
	      // NB model #0 is the null model.
	      //
	      CLOG(6) << round << "Estimating model membership probabilities P(m|g) and overall likelihood P({g},{theta})\n";

	      vector<int> log_g_prod (max_models + 1, 0);                       // log_g_prod[k] = log \prod_g \sum_{m<=k} P(g|m,theta_m)
	      for (int g = 0; g < expr_tab.sequences(); ++g)
		{
		  vector<int> log_joint_likelihood (max_models + 1);            // log_joint_likelihood[m] = log P(g|m) = log ( P(seq_g|seqmodel_m) . P(expr_g|exprmodel_m) )
		  vector<int> log_sum_joint_likelihood (max_models + 1);        // log_sum_joint_likelihood[k] = log sum_{m'=0}^k P(g|m')
		  log_joint_likelihood[0] = log_sum_joint_likelihood[0] = ScorePMul (gibbs_null.profile_null_score[g], gauss_null.profile_null_score[g]);
		  for (int m = 0; m < max_models; ++m)
		    {
		      log_joint_likelihood[m+1]     = ScorePMul (gibbs[m]->log_likelihood_sc (g), gauss[m]->log_likelihood_sc (g));
		      log_sum_joint_likelihood[m+1] = ScorePSum (log_sum_joint_likelihood[m], log_joint_likelihood[m+1]);
		    }
		  int log_k_sum = -InfinityScore;                             // log_k_sum = log sum_{k=m}^{k_max} (P(k) / sum_{m'=0}^k P(g|m'))
		  sstring log_string;
		  if (CLOGGING(4))
		    {
		      log_string.flags ((log_string.flags() & ~ios::right & ~ios::fixed) | ios::left | ios::scientific);
		      sstring n = gibbs_dataset.row_name[g].substr(0,10);
		      n << "':";
		      log_string << "Model posteriors for '" << setw(10) << n;
		    }
		  for (int m = max_models; m >= 0; --m)
		    {
		      ScorePSumAcc (log_k_sum, ScorePMul (log_k_prior[m], -log_sum_joint_likelihood[m]));
		      int log_new_membership_prob = ScorePMul (log_joint_likelihood[m], log_k_sum);     // log_new_membership_prob = log P(m|g)
		      if (m > 0)
			{
			  gibbs[m-1]->set_membership_probability (g, log_new_membership_prob);
			  gauss[m-1]->set_membership_probability (g, log_new_membership_prob);
			}
		      if (CLOGGING(4))
			{
			  log_string << " P(";
			  if (m > 0) log_string << m; else log_string << "null";
			  log_string << ")=";
			  log_string.precision(3);
			  log_string.width(5);
			  log_string << Score2Prob (log_new_membership_prob);
			}
		    }
		  for (int k = 0; k <= max_models; ++k) ScorePMulAcc (log_g_prod[k], log_sum_joint_likelihood[k]);
		  if (CLOGGING(4)) CL << log_string << "\n";
		}
	      
	      multiset<int> k_series;                   // score_pr_sum(k_series) = log \sum_k ( (P(k)/k) . \prod_g \sum_{m<=k} P(g|m,theta_m) )
	      for (int k = 0; k <= max_models; ++k)
		k_series.insert (log_k_prior[k] - Prob2Score (k+1) + log_g_prod[k]);

	      int score = ScorePSumSet (k_series);        // score = log P({g},{theta})
	      for (int m = 0; m < max_models; ++m)
		ScorePMulAcc (score, ScorePMul (gibbs[m]->log_parameter_prior_sc(), gauss[m]->log_parameter_prior_sc()));
	      
	      CLOG(7) << round << "Log-likelihood = " << Score2Bits (score) << " bits, best so far = " << Score2Bits (best_score) << " bits\n";
	      
	      // (iv): store models if best so far

	      if (score > best_score)
		{
		  for (int m = 0; m < max_models; ++m)
		    {
		      if (best_score > -InfinityScore) { delete best_gibbs[m]; delete best_gauss[m]; }
		      best_gibbs[m] = new Gibbs_alignment (*gibbs[m]);
		      best_gauss[m] = new Gaussian_cluster (*gauss[m]);
		    }
		  best_score = score;
		  best_seed = seed;
		  CLOG(7) << round << "New hi-score! Log-likelihood = " << Score2Bits (score) << " bits\n";
		}
	      
	      if (score > best_seed_score) best_seed_score = score;

	      // end of sampling loop
	    }

	  // clean up
	  //
	  for (int m = 0; m < max_models; ++m) { delete gibbs[m]; delete gauss[m]; }

	  if (deja_vu_enabled && best_seed_score == best_score && seed != best_seed && best_score > -InfinityScore)
	    {
	      CLOG(7) << "Deja vu! Got to this score from seed #" << best_seed+1 << ". Time to retire...\n";
	      break;
	    }
	  
	  // end of seeding loop
	}

      if (best_score > -InfinityScore)
	{
	  for (int m = 0; m < max_models; ++m)
	    {
	      cout << "-----------------------------------------------------------------------------------------------------------------\n";
	      cout << "Model #" << m+1 << "\n";
	      cout << "-----------------------------------------------------------------------------------------------------------------\n";
	      best_gibbs[m]->display (cout);
	      best_gauss[m]->display (cout);
	      delete best_gibbs[m];
	      delete best_gauss[m];
	    }
	  cout << "-----------------------------------------------------------------------------------------------------------------\n";
	  cout << "Expected log-likelihood = " << Score2Bits (best_score) << " bits\n";
	  cout << "-----------------------------------------------------------------------------------------------------------------\n";
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
