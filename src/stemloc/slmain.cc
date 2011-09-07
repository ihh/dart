#include "stemloc/slmain.h"
#include "scfg/covmodel.h"
#include "scfg/paircfgem.h"
#include "scfg/cfgdotplot.h"
#include "amap/amap_adapter.h"
#include "stemloc/slkeywords.h"

void Stemloc::init_opts()
{
  INIT_CONSTRUCTED_OPTS_LIST (opts, -1, "[options] <sequence file(s)>",
			      "comparative RNA structure finder using pairwise SCFGs");

  opts.newline();
  opts.print_title ("General options");

  sstring ps_str;
  ps_str << "use this parameterisation instead of trying to guess (valid values: " << map_keys(gramset.phmm_cfg_map) << ")";

  opts.newline();
  opts.add ("g -global", global = true, "force global alignment", true,
	    "l -local", "force local alignment");
  opts.add ("ps -paramset", param_set, ps_str.c_str(), false);
  opts.add ("to -trainonly", train_only = false, "don't do search, just do any specified training");

  opts.newline();
  opts.print_title ("Training options");

  opts.newline();
  opts.add ("tf -trainfold", trainfold = "", "train folding parameters from file (Stockholm format)", false);
  opts.add ("ta -trainalign", trainalign = "", "train alignment parameters from file (Stockholm format)", false);
  opts.add ("tc -traincmp", traincmp = "", "train sequence comparison parameters from file (Stockholm format)", false);

  opts.newline();
  opts.add ("lg -loadgramset", loadgramset = "", "load grammars and parameters from file, pre-training (stemloc grammar format)", false);
  opts.add ("sg -savegramset", savegramset = "", "save grammars and parameters to file, post-training (stemloc grammar format)", false);

  opts.newline();
  opts.print_title ("Expert training options");

  opts.add ("mr -maxrounds", em_max_iter = -1, "max number of rounds (iterations) of EM", false);
  opts.add ("mi -mininc", em_min_inc = .001, "minimum fractional increase in log-likelihood per round of EM");
  opts.add ("f -forgive", em_forgive = 0, "number of non-increasing rounds of EM to forgive");

  opts.newline();
  opts.print_title ("Fold envelope generation options");

  opts.newline();
  opts.add ("nf -nfold", pairwise_nfold = -1, "number of folds to sample; -1 to unlimit");
  opts.add ("rf -rndfold", random_fold = false, "sample <nf> folds stochastically, rather than using <nf> best hits");
  opts.add ("cf -cachefold", cache_fold = true, "cache fold envelopes: memory-hungry but fast");
  opts.add ("len -maxstemlen", max_subseq_len = -1, "maximum length of envelope subsequences; -1 to unlimit");
  opts.add ("loop -minlooplen", min_loop_len = 3, "minimum length of loop regions");

  opts.newline();
  opts.print_title ("Alignment envelope generation options");

  opts.newline();
  opts.add ("na -nalign", pairwise_nalign = 100, "number of alignments to sample; -1 to unlimit");
  opts.add ("mb -megabytes", max_megabytes = 0, "truncate <na> so Pair-SCFG uses at most this many megabytes", false);
  opts.add ("ra -rndalign", random_align = false, "sample <na> alignments stochastically, rather than using <na> best hits");
  opts.add ("band -bandalign", align_band = -1, "band size for pairwise alignment; -1 to unlimit");

  opts.newline();
  opts.print_title ("Pairwise alignment options");

  opts.newline();
  opts.add ("ms -minscore", min_bitscore = Score2Bits(-InfinityScore), "minimum bit-score for reporting a match", false);
  opts.add ("mh -maxhits", max_hits = 1, "maximum number of matches to return per sequence-pair");
  opts.add ("xbp -extrabasepairs", allow_extra_basepairs = false, "allow extra basepairs in pre-specified structures");

  opts.newline();
  opts.print_title ("Multiple alignment options: progressive");

  opts.newline();
  opts.add ("m -multiple", multi = true, "build progressive multiple alignment", true,
	    "p -pairwise", "don't build multiple alignment: stop after pairwise stage");

  // disable profile mode for now (pending resolution of RT ticket #143, "fix or remove profile mode") -- IH, 12/18/2007
  progressive = true;
  /*
  opts.add ("prog -progressive", progressive = true, "align by progressive pairwise comparisons, like CLUSTAL", true,
	    "prof -profile", "align everything to EM-optimized consensus profile [EXPERIMENTAL]");
  */

  opts.add ("lib -liberal", liberal = true, "for progressive alignment:  use liberal strategy", true,
	    "con -conservative", "for progressive alignment:  use conservative strategy");
  opts.add ("si -show-intermediates", show_intermediates = true, "show all intermediate alignments");
  opts.add ("ma -multi-nalign", multi_nalign = -2, "number of alignments in envelope; -1 to unlimit, -2 to use '-na' value");
  opts.add ("mf -multi-nfold", multi_nfold = -2, "number of folds in envelope; -1 to unlimit, -2 to use '-nf' value");
  opts.add ("mxs -minextendscore", min_extend_bitscore = Score2Bits(-InfinityScore), "minimum bit-score for extending hits", false);
  opts.add ("mxl -minextendlen", min_overlap_len = 10, "minimum overlap for extending hits");

  opts.newline();
  opts.print_title ("Multiple alignment options: sequence annealing");
  opts.print ("(Use --noshow-intermediates flag to disable display of intermediate alignments.)\n");
  opts.add ("ama -ama", use_ama = false, "build multiple alignment by 'sequence annealing'", false);
  opts.add ("gf -gapfactor", gap_factor = 1, "AMA [Alignment Metric Accuracy] gap factor", true);
  opts.add ("ewt -edgeweightthreshold", edge_weight_threshold = 0, "edge weight threshold", true);
  opts.add ("cr -consistencyreps", num_consistency_reps = 0, "number of PROBCONS-style consistency transformations", true);
  opts.add ("gui -gui", output_for_gui = false, "output AMA intermediate alignments for interactive Java GUI", false);


  opts.newline();
  opts.print_title ("Posterior probability matrix options");

  opts.newline();
  opts.add ("paf -postalignfile", palign_prefix = "", "filename prefix for alignment post.prob. matrices", false);
  opts.add ("pff -postfoldfile", pfold_prefix = "", "filename prefix for fold post.prob. matrices", false);

  opts.newline();
  opts.print_title ("Common shorthand options");

  opts.newline();
  opts.add ("verbose", "-log 8 -log ALLOC -log CFGDP", "show progress of dynamic programming recursions");
  opts.add ("utr", "--multiple --local --liberal --maxstemlen 100", "analysis of UTRs & introns");
  opts.add ("gene", "--multiple --global", "analysis of RNA genes");
  opts.add ("fast", "--nfold 100 --norndfold", "fast analysis: considers best 100 RNA structures");
  opts.add ("slow", "--nfold 1000 --norndfold", "slow analysis: considers best 1000 RNA structures");
  opts.add ("tedious", "--nfold -1", "tediously slow analysis: considers all RNA structures");
  opts.add ("pleasekillme", "--nfold -1 --nalign -1", "deathly slow analysis: considers all structures & all alignments");

  opts.newline();
  opts.newline();
  opts.print ("Input sequence file can be in either Stockholm or FASTA format.\n");
  opts.print ("If in Stockholm format, the corresponding multiple alignment is taken as fixed.");
  opts.newline();

  opts.newline();
  opts.print ("References\n");
  opts.print ("----------\n");
  opts.print ("Holmes I, 2005.  Accelerated probabilistic inference of RNA structure evolution.\n");
  opts.print ("BMC Bioinformatics 6:73.   http://www.biomedcentral.com/1471-2105/6/73\n");

}

void Stemloc::parse_opts()
{
  opts.parse_or_die();
  // adjust some options
  if (multi_nfold == -2) multi_nfold = pairwise_nfold;
  if (multi_nalign == -2) multi_nalign = pairwise_nalign;
  // do some extra checking on certain options
  if ((trainalign.size() || traincmp.size()) && param_set.empty())
    THROWEXPR ("You must specify a parameter set when training comparative models.");
}

void Stemloc::input_data()
{
  // get sequence filenames
  const vector<sstring>& seq_filename = opts.args;
  
  // initialise null model
  null_emit_prob = vector<Prob> (alphabet.size(), 1. / (double) alphabet.size());  // default null model

  // if not just training, read in sequences & estimate null model
  if (train_only)
    {
      if (seq_filename.size())
	CLOGERR << "Warning: the following sequence files will be ignored, since the --trainonly flag is set: " << seq_filename << "\n";
    }
  else
    {
      // check for at least one sequence file
      if (seq_filename.size() == 0)
	THROWEXPR (opts.short_help() << "You need to specify at least one sequence file");
      // read everything in as a Stockholm alignment
      for_const_contents (vector<sstring>, seq_filename, sf)
	{
	  ifstream in (sf->c_str());
	  if (!in) THROWEXPR ("Sequence file not found: '" << *sf << "'");
	  align_db.read_Stockholm_or_FASTA (in, seq_db);
	}
      // check names are unique
      seq_db.index.assert_names_unique (true);
      // conservatively propagate consensus folds
      align_db.propagate_consensus_folds (false);
      // update the sequence database index
      seq_db.update_index();
      seq_db.update_alphabet (&alphabet);
      seq_db.seqs_update (alphabet, (Profile_flags_enum::Profile_flags) (Profile_flags_enum::DSQ | Profile_flags_enum::SCORE));
      
      // get the null model
      null_emit_prob = seq_db.get_null_model (alphabet.size());
    }
}

void Stemloc::load_gramset()
{
  if (loadgramset.size())
    {
      SExpr_file sexpr_file (loadgramset.c_str());
      PCFG_builder::init_gramset (sexpr_file.sexpr.find_or_die (CFG_GRAMSET), gramset);
    }

  const sstring local_ss (STEMLOC_LOCAL), global_ss (STEMLOC_GLOBAL);
  for_contents (PHMM_CFG_map, gramset.phmm_cfg_map, phmm_cfg)
    {
      PScores& psc = phmm_cfg->second.hmm_pscores;
      if (psc.contains_suffix (local_ss))
	psc[local_ss] = global ? -InfinityScore : 0;
      if (psc.contains_suffix (global_ss))
	psc[global_ss] = global ? 0 : -InfinityScore;
    }
}

void Stemloc::save_gramset()
{
  if (savegramset.size())
    {
      ofstream gramset_stream (savegramset.c_str());
      Stemloc_gramset::prefix_comments (gramset_stream);
      PCFG_builder::gramset2stream (gramset_stream, gramset);
    }
}

void Stemloc::init_fold_params()
{
  // train params, if asked
  if (trainfold.size())
    {
      CTAG(6,STEMLOC) << "Training pre-folding parameters from alignments in file '" << trainfold << "'\n";
      // read in alignments
      Sequence_database stock_seq_db;
      Stockholm_database stock_db;
      ifstream stock_db_file (trainfold.c_str());
      if (!stock_db_file)
	THROWEXPR ("Secondary structure database '" << trainfold << "' not found");
      stock_db.read (stock_db_file, stock_seq_db);

      // digitize sequences
      stock_seq_db.seqs2dsqs (alphabet);  // for Pair_CFG DP

      // fit null model to training sequences before starting EM
      // (this prevents the first round of training being dominated by the null model)
      gramset.single_scfg_pscores.set_null_model (stock_seq_db, gramset.single_scfg.null_emit, gramset.single_scfg.null_extend);

      // do EM
      Pair_CFG_trainer trainer (stock_db, gramset.single_scfg, gramset.single_scfg_prior, gramset.single_scfg_pscores, false);
      trainer.train_single (em_min_inc, em_max_iter, em_forgive);
    }

  // initialise CFG scores
  single_scfg = gramset.single_scfg.eval_cfg_scores (gramset.single_scfg_pscores);
}

void Stemloc::init_pair_fold_env (Fold_envelope& env, const Named_profile& np, int nfold)
{
  // print log message
  CTAG(6,STEMLOC) << "Finding fold envelope for sequence '" << np.name << "'\n";

  // keep a few flags saying whether we've got an envelope yet
  bool got_foldenv = false;
  bool got_from_cache = false;

  // check in pairwise cache
  if (cache_fold)
    if (foldenv_db.find (np.name) != foldenv_db.end())
      {
	CTAG(3,STEMLOC) << "Retrieving fold envelope '" << np.name << "' from cache\n";
	env = foldenv_db[np.name];  // this assignment is costly; env should probably be a pointer-to-pointer instead
	got_foldenv = got_from_cache = true;
      }

  // check if a structure was specified for this sequence; if so, don't try to pre-fold it
  if (!got_foldenv)
    got_foldenv = attempt_init_from_fold_string (env, np);

  // if no structure was specified, do the pre-folding thing
  if (!got_foldenv)
    {
      // initialise full, banded fold envelope
      env.initialise_3prime_banded (np.size(), max_subseq_len, min_loop_len);
      // show dotplot
      if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_FOLDENV STEMLOC_BANDED_FOLDENV))
	{
	  CL << "Banded folding pre-envelope for '" << np.name << "':\n";
	  env.render_dotplot (CL, np.seq, true);
	}

      // sample structures
      sample_folds (env, np, nfold);
    }
  // print no. of subseqs to log
  if (CTAGGING(6,STEMLOC_FOLDENV STEMLOC_ENVSIZE))
    CL << "Fold envelope for " << np.name << " has " << env.subseqs() << " subseqs\n";

  // store in cache
  if (cache_fold && !got_from_cache)
    {
      CTAG(2,STEMLOC) << "Storing fold envelope '" << np.name << "' in cache\n";
      foldenv_db[np.name] = env;  // another costly assignment; there are better ways to do this
    }
}

void Stemloc::sample_folds (Fold_envelope& env, const Named_profile& np, int nfold)
{
  // create dummy y-sequence & envelope
  const Fold_envelope dummy_yenv;
  const Named_profile dummy_npy = Named_profile();

  // if unlimited structures, stick with banded envelope; otherwise, sample structures
  if (nfold < 0)
    CTAG(6,STEMLOC_PREFOLD) << "Using banded fold envelope for sequence '" << np.name << "'\n";
  else
    {
      // find the sampled fold envelope
      CTAG(6,STEMLOC_PREFOLD) << "Pre-folding sequence '" << np.name << "'\n";
      Sampled_fold_envelope sample_env;
      Subseq_coords_count subseq_count;
      // stochastic or deterministic?
      if (random_fold)
	{
	  // create Inside DP matrix
	  const Pair_inside_matrix inside (np, dummy_npy, env, dummy_yenv, single_scfg, false);
	  if (CTAGGING(1,STEMLOC_DP))
	    inside.show (CL);
	  // sample the envelope
	  subseq_count = sample_env.initialise_sampled (inside, nfold, min_loop_len);
	}
      else
	{
	  // create CYK-KYC DP matrix
	  const Pair_CYK_KYC_matrix cyk_kyc (np, dummy_npy, env, dummy_yenv, single_scfg, false);
	  if (CTAGGING(1,STEMLOC_DP))
	    cyk_kyc.show (CL);
	  // sample the envelope
	  subseq_count = sample_env.initialise_best (cyk_kyc, nfold, min_loop_len);
	}
      // assign the envelope
      env.swap (sample_env);
      // show dotplot
      if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_FOLDENV STEMLOC_SAMPLED_FOLDENV))
	{
	  CL << "Sampled fold envelope for '" << np.name << "':\n";
	  Fold_envelope::render_dotplot_from_counts (CL, subseq_count, np.seq, nfold, true);
	}
    }
}

bool Stemloc::attempt_init_from_fold_string (Fold_envelope& env, const Named_profile& np)
{
  const Stockholm& stock = *align_db.align_index[align_db.align_row[np.name].nAlign];
  const sstring fold_string = stock.get_fold_string (np.name);
  const bool got_ss = fold_string.size() > 0;
  if (got_ss)
    {
      // initialise envelope from fold string
      env.initialise_from_fold_string (fold_string);
      if (allow_extra_basepairs)
	env.connect_all (min_loop_len);
      // show dotplot
      if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PREFOLD))
	{
	  CL << "Pre-specified fold envelope for sequence '" << np.name << "':\n";
	  env.render_dotplot (CL, np.seq, true);
	}
    }
  return got_ss;
}

void Stemloc::init_multi_fold_env (Fold_envelope& env, const Named_profile& np, const Local_fold_string& local_fold)
{
  // print log message
  CTAG(6,STEMLOC) << "Finding fold envelope for local fold '" << local_fold.fold << "' starting at " << local_fold.start << " in sequence '" << np.name << "' (length " << local_fold.seqlen << ")\n";
  
  // check in multiple alignment cache
  bool got_from_cache = false;
  if (cache_fold)
    if (multi_foldenv_db.find (np.name) != multi_foldenv_db.end())
      {
	CTAG(2,STEMLOC) << "Retrieving fold envelope '" << np.name << "' from cache\n";
	env = multi_foldenv_db[np.name];  // this assignment is costly; env should probably be a pointer-to-pointer instead
	got_from_cache = true;
      }
  
  // create envelope
  // if a fold string was specified in the input data, this overrides the local fold string.
  // FIXME: this is not quite right, because '-xbp' option might have added extra basepairs;
  // what we should really do is merge the local fold string with the global one when we create
  // the local fold string in the first place, but i'm feeling too lazy to fix this today (NYE'03).
  if (!got_from_cache)
    if (!attempt_init_from_fold_string (env, np))
      {
	// initialise local envelope
	if (global || multi_nfold >= 0)
	  env.initialise_local_fold_string_constrained (local_fold, max_subseq_len, min_loop_len);
	else
	  env.initialise_local_fold_string_neighbourhood (local_fold, max_subseq_len, min_loop_len, min (min_overlap_len, (int) local_fold.fold.size()), false);
	// show annotated dotplot
	if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_FOLDENV STEMLOC_CONSTRAINED_FOLDENV))
	  {
	    CL << "Constrained pre-envelope for '" << np.name << "':\n";
	    env.render_annotated_dotplot (CL, np.seq, true);
	  }
	
	// if there is room for more basepairing (i.e. local fold string is not global), then sample structures
	if (!local_fold.is_global())
	  sample_folds (env, np, multi_nfold);

	// make it local
	if (!global && multi_nfold >= 0)
	  env.make_local_overlap (local_fold.start, local_fold.fold.size(), min_overlap_len);
      }

  // show dotplot
  if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_FOLDENV STEMLOC_MULTI_FOLDENV))
    {
      CL << "Local fold envelope for '" << np.name << "':\n";
      env.render_annotated_dotplot (CL, np.seq, true);
    }

  // print no. of subseqs to log
  if (CTAGGING(6,STEMLOC_FOLDENV STEMLOC_ENVSIZE))
    CL << "Fold envelope for " << np.name << " has " << env.subseqs() << " subseqs\n";

  // store in cache
  if (cache_fold && !got_from_cache)
    {
      CTAG(2,STEMLOC) << "Storing fold envelope '" << np.name << "' in cache\n";
      multi_foldenv_db[np.name] = env;  // another costly assignment; there are better ways to do this
    }
}

void Stemloc::init_align_params()
{
  // initialise null emit scores
  for_contents (PHMM_CFG_map, gramset.phmm_cfg_map, name_phmm_cfg)
    name_phmm_cfg->second.hmm_pscores[name_phmm_cfg->second.hmm_null_emit] = Prob2FScoreVec (null_emit_prob);

  // train params, if asked
  if (trainalign.size())
    {
      CTAG(6,STEMLOC) << "Training pre-alignment parameters from alignments in file '" << trainalign << "'\n";
      // read in alignments
      Sequence_database stock_seq_db;
      Stockholm_database stock_db;
      ifstream stock_db_file (trainalign.c_str());
      if (!stock_db_file)
	THROWEXPR ("Pairwise alignment database '" << trainalign << "' not found");
      stock_db.read (stock_db_file, stock_seq_db);

      // digitize sequences
      stock_seq_db.seqs2dsqs (alphabet);  // for counting symbol frequencies to get null model (see below)
      stock_seq_db.seqs2scores (alphabet);  // for Pair_HMM DP

      // fit null model to training sequences before starting EM
      // (this prevents the first round of training being dominated by the null model)
      PHMM_CFG& phmm_cfg = gramset.phmm_cfg_map[param_set];
      Boolean_group dummy_null_extend = phmm_cfg.hmm_pscores.new_boolean_group ("dummyNullExtend");   // create a dummy null-extend group... this is truly sad
      phmm_cfg.hmm_pscores.set_null_model (stock_seq_db, phmm_cfg.hmm_null_emit, dummy_null_extend);
      phmm_cfg.hmm_pscores.delete_group (dummy_null_extend.group_idx);   // THE SHAME!!! it burns

      // do EM
      Pair_HMM_trainer trainer (stock_db, phmm_cfg.hmm, phmm_cfg.hmm_prior, phmm_cfg.hmm_pscores);
      trainer.train (em_min_inc, em_max_iter, em_forgive);
    }
  
  // create HMM scores, get valid (non-padding) states
  for_contents (PHMM_CFG_map, gramset.phmm_cfg_map, name_phmm_cfg)
    {
      pair_hmm[name_phmm_cfg->first] = name_phmm_cfg->second.hmm.eval_hmm_scores (name_phmm_cfg->second.hmm_pscores);

      set<int> non_pad_states;
      for (int s = 0; s < pair_hmm[name_phmm_cfg->first].states(); ++s)
	if (name_phmm_cfg->second.hmm_pad_states.find(s) == name_phmm_cfg->second.hmm_pad_states.end())
	  non_pad_states.insert (s);
      pair_hmm_states[name_phmm_cfg->first] = non_pad_states;
    }
}

void Stemloc::init_cmp_params()
{
  // train params, if asked
  if (traincmp.size())
    {
      CTAG(6,STEMLOC) << "Training sequence comparison parameters from alignments in file '" << traincmp << "'\n";
      // read in alignments
      Sequence_database stock_seq_db;
      Stockholm_database stock_db;
      ifstream stock_db_file (traincmp.c_str());
      if (!stock_db_file)
	THROWEXPR ("Alignment database '" << traincmp << "' not found");
      stock_db.read (stock_db_file, stock_seq_db);

      // digitize sequences
      stock_seq_db.seqs2dsqs (alphabet);  // for Pair_CFG DP

      // fit null model to training sequences before starting EM
      // (this prevents the first round of training being dominated by the null model)
      PHMM_CFG& phmm_cfg = gramset.phmm_cfg_map[param_set];
      phmm_cfg.cfg_pscores.set_null_model (stock_seq_db, phmm_cfg.cfg.null_emit, phmm_cfg.cfg.null_extend);

      // do EM
      // we enforce global training, even though alignment may be local
      Pair_CFG_trainer trainer (stock_db, phmm_cfg.cfg, phmm_cfg.cfg_prior, phmm_cfg.cfg_pscores, false);
      trainer.train_pairwise (em_min_inc, em_max_iter, em_forgive);
    }

  // create the Pair_CFG_scores objects (the actual numerical scores for the SCFG)
  for_contents (PHMM_CFG_map, gramset.phmm_cfg_map, name_phmm_cfg)
    pair_cfg[name_phmm_cfg->first] = name_phmm_cfg->second.cfg.eval_cfg_scores (name_phmm_cfg->second.cfg_pscores);
}

void Stemloc::build_pairwise_alignments()
{
  // skip if train_only flag is set
  if (!train_only)
    {
      // do alignments
      CTAG(8,STEMLOC) << "Building pairwise alignments\n";
      // loop through sequence database, do CYK algorithm on all sequence pairs
      for (int i = 0; i < seq_db.size(); ++i)
	for (int j = i + 1; j < seq_db.size(); ++j)
	  {
	    // get sequences
	    const Named_profile& npx = *seq_db.index.profile[i];
	    const Named_profile& npy = *seq_db.index.profile[j];

	    // retrieve fold envelopes
	    Fold_envelope envx, envy;
	    init_pair_fold_env (envx, npx, pairwise_nfold);
	    init_pair_fold_env (envy, npy, pairwise_nfold);
	    
	    // create the sampled alignment envelope
	    Pair_envelope pair_env;
	    sstring best_param_set;
	    init_align_env (pair_env, npx, npy, envx, envy, best_param_set, pairwise_nalign);

	    // if doing dotplots, create Inside-Outside matrix & print to file
	    if (palign_prefix.size() || pfold_prefix.size())
	      {
		// print log message
		CTAG(6,STEMLOC) << "Summing over alignments and folds of '" << npx.name << "' and '" << npy.name << "' to get posterior probs.\n";

		const Pair_inside_outside_matrix in_out (npx, npy, envx, envy, pair_cfg[best_param_set], pair_env, !global);
		if (palign_prefix.size())
		  {
		    PairCFG_alignment_dotplot align_dotplot (in_out);
		    align_dotplot.write_dotplot (palign_prefix, npx.name, npy.name);
		  }
 		if (pfold_prefix.size())
		  {
		    const PairCFG_fold_dotplot xfold_dotplot (in_out, 0);
		    xfold_dotplot.write_dotplot (pfold_prefix, npx.name);

		    const PairCFG_fold_dotplot yfold_dotplot (in_out, 1);
		    yfold_dotplot.write_dotplot (pfold_prefix, npy.name);
		  }
	      }

	    // print log message
	    CTAG(6,STEMLOC) << "Aligning & folding '" << npx.name << "' and '" << npy.name << "'\n";
		
	    // do CYK with Waterman-Eggert style blocking
	    int hit = 0;
	    while (1)
	      {
		if (++hit > max_hits)
		  {
		    CTAG(5,STEMLOC_CMP) << "Reported " << max_hits << " alignments for this sequence-pair; stopping\n";
		    break;
		  }
		CTAG(5,STEMLOC_CMP) << "Finding alignment #" << hit << " for '" << npx.name << "' and '" << npy.name << "'\n";
		  
		// if this isn't the first hit, print out the alignment envelope to show Waterman-Eggert in progress
		if (hit > 1 && CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PAIRENV))
		  {
		    CL << "Alignment envelope for '" << npx.name << "' vs '" << npy.name << "', alignment #" << hit << ":\n";
		    pair_env.render_dotplot (CL, npx.seq, npy.seq, 1, true);
		  }
		    
		// build the DP matrix
		const Pair_CYK_matrix cyk (npx, npy, envx, envy, pair_cfg[best_param_set], pair_env, !global);
		if (CTAGGING(1,STEMLOC_DP))
		  CL << cyk.cfg_dump();
		    
		// test for -inf score
		if (cyk.final_score == -InfinityScore && hit == 1)
		  {
		    CLOGERR << "Warning: alignment score for '" << npx.name << "' and '" << npy.name << "' is -infinity; skipping.\n";
		    CLOGERR << " (this usually means there is no valid alignment of the sequences, which often means that the precomputed structures or alignments are too restrictive;\n";
		    CLOGERR << "  try setting -nf or -na higher?)\n";
		    break;
		  }
		
		// get alignment, store
		const Pair_CFG_alignment cyk_alignment = cyk.alignment();
		pair_align_db.push_back (cyk_alignment);

		// check if hit is worth reporting
		if (Score2Bits(cyk.final_score) < min_bitscore)
		  {
		    CTAG(5,STEMLOC_CMP) << "Score of " << Score2Bits(cyk.final_score) << " bits is below threshold of " << min_bitscore  << " bits; stopping\n";
		    break;
		  }
		    
		// output
		if (show_intermediates || !multi)
		  cyk_alignment.show (cout);
		  
		// do Waterman-Eggert style blocking-out
		cyk_alignment.add_pairwise_path (pair_env, false);
	      }
	  }
    }
}

void Stemloc::build_multiple_alignments_ama()
{

  CTAG(8,STEMLOC) << "Collecting posterior probabilities for AMA.\n";
      
  // to hold post prob matrices for AMA
  Dotplot_map dotplots;

  // loop through sequence database, get post prob alignment matrices for each pair
  for (int i = 0; i < seq_db.size(); ++i)
    for (int j = i + 1; j < seq_db.size(); ++j)
      {
	// get sequences
	const Named_profile& npx = *seq_db.index.profile[i];
	const Named_profile& npy = *seq_db.index.profile[j];

	// retrieve fold envelopes
	Fold_envelope envx, envy;
	init_pair_fold_env (envx, npx, pairwise_nfold);
	init_pair_fold_env (envy, npy, pairwise_nfold);
	    
	// create the sampled alignment envelope
	Pair_envelope pair_env;
	sstring best_param_set;
	init_align_env (pair_env, npx, npy, envx, envy, best_param_set, pairwise_nalign);

	// print log message
	CTAG(6,STEMLOC) << "Summing over alignments and folds of '" << npx.name << "' and '" << npy.name << "' to get posterior probs.\n";

	// create Inside-Outside matrix and store alignment dotplot
	const Pair_inside_outside_matrix in_out (npx, npy, envx, envy, pair_cfg[best_param_set], pair_env, !global);
	PairCFG_alignment_dotplot align_dotplot (in_out);
	dotplots[i][j] = align_dotplot;

	// write alignment post prob matrices if requested
	if (palign_prefix.size())
	  align_dotplot.write_dotplot (palign_prefix, npx.name, npy.name);
	// write fold post prob matrices if requested
	if (pfold_prefix.size())
	  {
	    const PairCFG_fold_dotplot xfold_dotplot (in_out, 0);
	    xfold_dotplot.write_dotplot (pfold_prefix, npx.name);
	    
	    const PairCFG_fold_dotplot yfold_dotplot (in_out, 1);
	    yfold_dotplot.write_dotplot (pfold_prefix, npy.name);
	  }
      }

  AMAP_adapter amap_adapter (seq_db, dotplots);
  Stockade stockade = amap_adapter.get_alignment (gap_factor, edge_weight_threshold, num_consistency_reps, show_intermediates, output_for_gui);

  if (global)
    stockade.align.write_Stockholm (cout);
  else
    stockade.align.write_Stockholm_NSE (cout);

}

void Stemloc::init_align_env (Pair_envelope& pair_env, const Named_profile& npx, const Named_profile& npy, const Fold_envelope& xenv, const Fold_envelope& yenv, sstring& best_param_set, int nalign)
{
  // print log message
  CTAG(6,STEMLOC) << "Finding alignment envelope for sequences '" << npx.name << "' and '" << npy.name << "'\n";

  // figure out how much memory is available
  unsigned long max_bytes = max_megabytes << 20;
  const bool need_inside_outside =
    use_ama
    ? true
    : (palign_prefix.size() || pfold_prefix.size());
  if (need_inside_outside)
    max_bytes >>= 1;

  // need to choose a parameterisation
  set<sstring> param_sets;
  if (param_set.size())
    param_sets.insert (param_set);
  else
    {
      vector<sstring> ps = map_keys (gramset.phmm_cfg_map);
      param_sets = set<sstring> (ps.begin(), ps.end());
    }

  // keep a flag saying whether we've got the alignment yet
  bool got_align = false;

  // if an alignment between the two sequences was prespecified, then don't bother pre-aligning
  if (!got_align)
    {
      const Stockholm_database::Align_row& arx = align_db.align_row[npx.name];
      const Stockholm_database::Align_row& ary = align_db.align_row[npy.name];
      if (arx.nAlign == ary.nAlign)
	{
	  const Pairwise_path path (align_db.align_index[arx.nAlign]->path, arx.nRow, ary.nRow, true);
	  pair_env = Pair_envelope (npx.size(), npy.size(), 0);
	  pair_env.add_pairwise_path (path);
	  got_align = true;
	  if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PREALIGN))
	    {
	      CL << "Using pre-specified alignment for '" << npx.name << "' vs '" << npy.name << "':\n";
	      pair_env.render_dotplot (CL, npx.seq, npy.seq, 1, true);
	    }
	}
    }

  // if no alignment envelope has been found yet, do the pre-alignment thing
  if (!got_align)
    {
      // initialise full, banded alignment envelope
      Pair_envelope banded_pair_env;
      banded_pair_env.initialise_banded (npx.size(), npy.size(), align_band);
      if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PAIRENV STEMLOC_BANDED_PAIRENV))
	{
	  CL << "Banded alignment pre-envelope for '" << npx.name << "' vs '" << npy.name << "':\n";
	  banded_pair_env.render_dotplot (CL, npx.seq, npy.seq, 1, true);
	}

      // if unlimited structures, copy banded envelope to sampled envelope; otherwise, sample structures
      if (nalign < 0)
	{
	  CTAG(6,STEMLOC_PREALIGN) << "Using banded alignment envelope for '" << npx.name << "' and '" << npy.name << "'\n";
	  pair_env = banded_pair_env;
	}
      else   // nalign >= 0
	{
	  // find the sampled alignment envelope
	  CTAG(6,STEMLOC_PREALIGN) << "Pre-aligning '" << npx.name << "' and '" << npy.name << "'\n";
	  // stochastic or deterministic?
	  if (random_align)
	    {
	      // choose a parameterisation
	      const Pair_forward_DP_matrix* best_forward = 0;
	      Score best_forward_score = -InfinityScore;
	      for_const_contents (set<sstring>, param_sets, n)
		{
		  CTAG(5,STEMLOC_PREALIGN) << "Trying parameterisation " << *n << "\n";
		  // create Forward DP matrix
		  const Pair_forward_DP_matrix* forward = new Pair_forward_DP_matrix (pair_hmm[*n], npx.prof_sc, npy.prof_sc, banded_pair_env.allow_cut);
		  // log score & possibly matrix
		  CTAG(6,STEMLOC_PREALIGN) << "Parameterisation " << *n << " scored " << Score2Bits(forward->final_score) << " bits\n";
		  if (CTAGGING(1,STEMLOC_DP))
		    forward->show (CL);
		  // check if this is the best score
		  if (forward->final_score > best_forward_score)
		    {
		      best_forward_score = forward->final_score;
		      if (best_forward) delete best_forward;
		      best_forward = forward;
		      best_param_set = *n;
		    }
		  else
		    delete forward;
		}
	      // check we got something
	      if (best_forward == 0)
		THROWEXPR ("Failed to find an optimal parameterisation");
	      // sample the envelope
	      Sampled_pair_envelope sampled_pair_env (npx, npy, xenv, yenv);
	      sampled_pair_env.cells_available = max_bytes / (pair_cfg[best_param_set].states() * sizeof(Score));
	      sampled_pair_env.initialise_sampled (*best_forward, nalign, pair_hmm_states[best_param_set]);
	      pair_env = sampled_pair_env;
	      // delete the best Forward matrix
	      delete best_forward;
	    }
	  else  // !random_align
	    {
	      // create HMM-emulating fold envelopes
	      Fold_envelope hmm_xenv, hmm_yenv;
	      hmm_xenv.initialise_3prime (npx.size());
	      hmm_yenv.initialise_3prime (npy.size());
	      // create a vector of Pair HMM-emulating Pair_CFG_scores
	      map<sstring,Pair_CFG_scores> hmm_cfg;
	      for_const_contents (set<sstring>, param_sets, n)
		hmm_cfg[*n] = Pair_CFG_scores (pair_hmm[*n]);
	      // choose a parameterisation
	      Pair_CYK_KYC_matrix* best_cyk_kyc = 0;
	      Score best_cyk_score = -InfinityScore;
	      for_const_contents (set<sstring>, param_sets, n)
		{
		  CTAG(5,STEMLOC_PREALIGN) << "Trying parameterisation " << *n << "\n";
		  // create CYK-KYC DP matrix, only filling the CYK part at this stage
		  Pair_CYK_KYC_matrix* cyk_kyc = new Pair_CYK_KYC_matrix (npx, npy, hmm_xenv, hmm_yenv, hmm_cfg[*n], banded_pair_env, false, false);
		  cyk_kyc->cyk.fill();
		  // log score & possibly matrix
		  CTAG(6,STEMLOC_PREALIGN) << "Parameterisation " << *n << " scored " << Score2Bits(cyk_kyc->cyk.final_score) << " bits\n";
		  if (CTAGGING(1,STEMLOC_DP))
		    cyk_kyc->cyk.show (CL);
		  // check if this is the best score
		  if (cyk_kyc->cyk.final_score > best_cyk_score)
		    {
		      best_cyk_score = cyk_kyc->cyk.final_score;
		      if (best_cyk_kyc) delete best_cyk_kyc;
		      best_cyk_kyc = cyk_kyc;
		      best_param_set = *n;
		    }
		  else
		    delete cyk_kyc;
		}
	      // check we got something
	      if (best_cyk_kyc == 0)
		THROWEXPR ("Failed to find an optimal parameterisation");
	      // fill the KYC bit and log
	      best_cyk_kyc->kyc.fill();
	      if (CTAGGING(1,STEMLOC_DP))
		best_cyk_kyc->show (CL);
	      // sample the envelope
	      Sampled_pair_envelope sampled_pair_env (npx, npy, xenv, yenv);
	      sampled_pair_env.cells_available = max_bytes / (pair_cfg[best_param_set].states() * sizeof(Score));
	      sampled_pair_env.initialise_best (*best_cyk_kyc, pair_hmm[best_param_set], nalign, pair_hmm_states[best_param_set]);
	      pair_env = sampled_pair_env;
	      // delete the best CYK-KYC matrix
	      delete best_cyk_kyc;
	    }
	  // show dotplot
	  if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PAIRENV STEMLOC_SAMPLED_PAIRENV))
	    {
	      CL << "Sampled alignment envelope for '" << npx.name << "' vs '" << npy.name << "':\n";
	      pair_env.render_dotplot (CL, npx.seq, npy.seq, nalign, true);
	    }
	}
    }
  if (CTAGGING(6,STEMLOC_PAIRENV STEMLOC_ENVSIZE))
    CL << "Alignment envelope for '" << npx.name << "' vs '" << npy.name << "' has " << pair_env.cells() << " cells\n";

  // if we haven't got a best parameter set, then we need to choose one now.
  if (best_param_set.empty())
    {
      Score best_forward_score = -InfinityScore;
      for_const_contents (set<sstring>, param_sets, n)
	{
	  CTAG(6,STEMLOC_PREALIGN) << "Trying parameterisation " << *n << "\n";
	  // create Forward DP matrix
	  const Pair_forward_DP_matrix forward (pair_hmm[*n], npx.prof_sc, npy.prof_sc, pair_env.allow_cut);
	  // log score & possibly matrix
	  CTAG(6,STEMLOC_PREALIGN) << "Parameterisation " << *n << " scored " << Score2Bits(forward.final_score) << " bits\n";
	  if (CTAGGING(1,STEMLOC_DP))
	    forward.show (CL);
	  // check if this is the best score
	  if (forward.final_score > best_forward_score || best_param_set.empty())
	    {
	      best_forward_score = forward.final_score;
	      best_param_set = *n;
	    }
	}
    }
  CTAG(6,STEMLOC_PREALIGN) << "Using parameterisation " << best_param_set << "\n";
}

void Stemloc::remove_seeds (const Alignment& align)
{
  typedef list<const Pair_CFG_alignment*>::iterator Seed_iterator;

  set<sstring> row_name (align.row_name.begin(), align.row_name.end());

  vector<Seed_iterator> to_erase;
  for_contents (list<const Pair_CFG_alignment*>, seed_db, seed_iter)
    {
      const sstring& xname = (*seed_iter)->npx.name;
      const sstring& yname = (*seed_iter)->npy.name;
      if (row_name.find (xname) != row_name.end()
	  || row_name.find (yname) != row_name.end())
	to_erase.push_back (seed_iter);
      // unaligned_seqname should already be cleared, but just in case
      if (unaligned_seqname.find (xname) != unaligned_seqname.end())
	unaligned_seqname.erase (xname);
      if (unaligned_seqname.find (yname) != unaligned_seqname.end())
	unaligned_seqname.erase (yname);
    }

  for_contents (vector<Seed_iterator>, to_erase, seed_iter_ptr)
    seed_db.erase (*seed_iter_ptr);
}

void Stemloc::build_multiple_alignments()
{
  if (multi && !train_only && seq_db.size() > 2)
    {
      // flush fold envelope database if 'nfold' parameter is different for multiple alignment than pairwise alignment
      if (multi_nfold != pairwise_nfold)
	foldenv_db.clear();

      // populate list of pointers to pairwise alignments for unaligned sequences
      seed_db.clear();
      for_const_contents (list<Pair_CFG_alignment>, pair_align_db, pa)
	seed_db.push_back (&*pa);

      // populate list of unaligned sequence names
      unaligned_seqname.insert (seq_db.index.name.begin(), seq_db.index.name.end());

      // multiple alignment generation loop
      for (int n_multi = 0; true; ++n_multi)
	{
	  CTAG(7,STEMLOC STEMLOC_MULTI) << "Building multiple alignment (looking for best pairwise seed...)\n";
	  // find highest-scoring pairwise seed
	  Score best_sc = -InfinityScore;
	  const Pair_CFG_alignment* best_pair_align = 0;
	  const Score min_ext_sc = Bits2Score (min_extend_bitscore);
	  for_const_contents (list<const Pair_CFG_alignment*>, seed_db, pa)
	    if ((*pa)->score > min_ext_sc && (best_pair_align == 0 || (*pa)->score > best_sc))
	      {
		best_pair_align = *pa;
		best_sc = (*pa)->score;
	      }
	  if (best_pair_align == 0)  // if no good pairwise seeds, then bail
	    {
	      if (n_multi == 0)
		CLOGERR << "Warning: couldn't find high-scoring pairwise seed alignment with score over ";
	      else
		CTAG(7,STEMLOC STEMLOC_MULTI) << "No more pairwise alignments with scores over ";
	      if (min_ext_sc > -InfinityScore)
		CL << min_extend_bitscore;
	      else
		CL << "-oo";
	      CL << " bits\n";
	      break;  // exit multiple alignment generation loop
	    }
	  else
	    {
	      // print log message saying what seed we're using
	      CTAG(7,STEMLOC STEMLOC_MULTI) << "Best pairwise seed alignment is between '" << best_pair_align->npx.name << "' and '" << best_pair_align->npy.name << "'\n";

	      // remove seed sequence names from unaligned_seqname
	      unaligned_seqname.erase (best_pair_align->npx.name);
	      unaligned_seqname.erase (best_pair_align->npy.name);

	      // delegate to specified build algorithm & display appropriately
	      if (progressive)
		if (liberal)
		  {
		    // progressive, liberal
		    const Stockholm multi_align = build_liberal_multiple_alignment (*best_pair_align);
		    if (global)
		      multi_align.write_Stockholm (cout);
		    else
		      multi_align.write_Stockholm_NSE (cout);
		    remove_seeds (multi_align);  // trim the database of pairwise seed alignments
		  }
		else
		  {
		    // progressive, conservative
		    const Stockade multi_stockade = build_conservative_multiple_alignment (*best_pair_align);
		    multi_stockade.align.write_Stockholm (cout);
		    remove_seeds (multi_stockade.align);  // trim the database of pairwise seed alignments
		  }
	      else
		{
		  // consensus motif
		  const Stockholm multi_align = build_motif_alignment (*best_pair_align);
		  multi_align.write_Stockholm_NSE (cout);
		  remove_seeds (multi_align);  // trim the database of pairwise seed alignments
		}
	    }
	}
    }
}

Stockade Stemloc::build_conservative_multiple_alignment (const Pair_CFG_alignment& best_pair_align)
{
  // print log message
  CTAG(8,STEMLOC) << "Doing conservative progressive multiple alignment\n";
  // build Stockade
  // the Stockade is just used to hold the Named_profile subseqs for local alignment
  // the path and the fold strings are continually updated by the code
  Stockade seed = best_pair_align.stockade();
  // build the fold string & fold envelope vectors for rows of seed
  vector<sstring> seed_fold (2);
  vector<Fold_envelope> seed_foldenv (2);
  for (int i = 0; i < 2; ++i)
    {
      seed_fold[i] = seed.align.get_fold_string (seed.align.row_name[i]);
      seed_foldenv[i].initialise_from_fold_string (seed_fold[i]);
    }
  seed.align.make_consensus_fold();
  
  // initialise Decomposition
  Decomposition decomp;
  decomp [Row_pair (0,1)] = seed.align.path;
  // initialise multiple alignment score
  // (this isn't quite the log-odds score, since pairwise scores aren't conditionally normalised)
  Score multi_align_sc = best_pair_align.score;

  // min score for extending a hit
  const Score min_ext_sc = Bits2Score (min_extend_bitscore);
  
  // for each unaligned sequence, this map will contain best pairwise match to seed
  Best_match_table best_match_by_seq;
  int last_row = -1;  // last row scanned for best_match
  
  // greedily grow the alignment
  while (!unaligned_seqname.empty())   // NB assumes EVERY sequence will go into the multiple alignment
    {
      // update best_match_by_seq
      while (last_row < ((int) seed_foldenv.size()) - 1)
	{
	  // update last_row
	  ++last_row;
	  // get sequence for last row
	  const Named_profile& npx = seed.np[last_row];
	  // get fold envelope for last row
	  const Fold_envelope& xenv = seed_foldenv[last_row];
	  // loop through unaligned seqs
	  for_contents (set<sstring>, unaligned_seqname, yname)
	    {
	      // get sequence
	      const Named_profile& npy (*seq_db.index.name2profile (*yname));
	      // print log message
	      CTAG(4,STEMLOC) << "Comparing aligned sequence '" << npx.name << "' and unaligned sequence '" << npy.name << "'\n";
	      // get fold envelope
	      Fold_envelope yenv;
	      init_pair_fold_env (yenv, npy, multi_nfold);
	      // get alignment envelope
	      Pair_envelope pair_env;
	      sstring best_param_set;
	      init_align_env (pair_env, npx, npy, xenv, yenv, best_param_set, multi_nalign);
	      // build the DP matrix
	      const Pair_CYK_matrix cyk (npx, npy, xenv, yenv, pair_cfg[best_param_set], pair_env, !global);
	      if (CTAGGING(1,STEMLOC_DP))
		CL << cyk.cfg_dump();
	      // get CYK score
	      const Score cyk_sc = cyk.final_score;
	      // print log message
	      CTAG(5,STEMLOC STEMLOC_MULTI) << "Score for '" << npx.name << "' vs '" << npy.name << "' is " << Score2Bits(cyk_sc) << " bits\n";
	      // check for low score; if so, don't attempt traceback
	      if (cyk_sc <= -InfinityScore || cyk_sc < min_ext_sc)
		{
		  CTAG(4,STEMLOC STEMLOC_MULTI) << "Score for '" << npx.name << "' vs '" << npy.name << "' is too low; skipping.\n";
		  continue;
		}
	      // get CYK alignment
	      const Pair_CFG_alignment cyk_alignment = cyk.alignment();
	      // show alignment
	      if (show_intermediates)
		cyk_alignment.show (cout);
	      // store
	      if (best_match_by_seq.find (*yname) == best_match_by_seq.end()
		  ? true
		  : best_match_by_seq[*yname].score < cyk_sc)
		best_match_by_seq[*yname] = Multi_match (cyk_alignment, last_row);
	    }
	}
      // if there are no matches (i.e. all scores were -infinity and/or below threshold), then bail out now
      if (best_match_by_seq.size() == 0)
	{
	  CTAG(5,STEMLOC) << "No more matches; halting multiple alignment\n";
	  break;
	}
      // find best match score
      Score best_match_sc = -InfinityScore;
      Best_match_table::const_iterator best_match = best_match_by_seq.end();
      for_const_contents (Best_match_table, best_match_by_seq, match)
	if (match->second.score > best_match_sc)
	  {
	    best_match_sc = match->second.score;
	    best_match = match;
	  }
      // recover best match
      const sstring& match_name (best_match->first);
      CTAG(4,STEMLOC STEMLOC_MULTI) << "Adding sequence '" << match_name << "' to multiple alignment with score " << Score2Bits(best_match_sc) << " bits\n";
      const Multi_match& match (best_match->second);
      const int seed_row = match.seed_row;
      const int new_row = seed_fold.size();
      const Named_profile& seed_np = seed.np[seed_row];
      // since this is a local hit, we have to extend the path so it's global w.r.t. the seed row
      Pairwise_path path;
      for (int i = 0; i < match.seed_lpad; ++i)
	path.append_column (1, 0);
      for (int i = 0; i < match.path.columns(); ++i)
	path.append_column (match.path(0,i), match.path(1,i));
      for (int i = path.count_steps_in_row(0); i < seed_np.size(); ++i)
	path.append_column (1, 0);
      decomp [Row_pair (seed_row, new_row)] = path;
      // add best match to seed
      seed.add_np (match.subseq);
      seed_fold.push_back (match.fold);
      seed_foldenv.push_back (Fold_envelope());
      seed_foldenv.back().initialise_from_fold_string (match.fold);
      ScorePMulAcc (multi_align_sc, best_match_sc);
      // compose multiple alignment
      seed.align.path.compose (decomp, true);
      // update fold strings
      for (int i = 0; i <= new_row; ++i)
	seed.align.set_fold_string (seed.align.row_name[i], seed_fold[i]);
      // set consensus fold
      seed.align.make_consensus_fold();
      // show an intermediate alignment, if this is not the last unaligned sequence
      const int seqs_remaining = unaligned_seqname.size() - 1;
      if (show_intermediates && seqs_remaining > 0)
	seed.align.write_Stockholm (cout);
      // print log message(s)
      CTAG(6,STEMLOC) << "Aligned '" << seed.align.row_name.back() << "' to '" << seed.align.row_name[match.seed_row] << "' (" << new_row+1 << " sequences aligned, " << seqs_remaining << " unaligned)\n";
      if (CTAGGING(2,STEMLOC STEMLOC_MULTI))
	seed.align.write_Stockholm (CL);
      // remove best match from unaligned_seqname and best_match_by_seq
      unaligned_seqname.erase (match_name);
      best_match_by_seq.erase (match_name);  // MUST do this last of all, as it invalidates best_match, match & match_name
    }

  // set score
  sstring sc_string;
  sc_string << Score2Bits (multi_align_sc);
  seed.align.add_gf_annot (sstring (Stockholm_bit_score_tag), sc_string);

  // return
  return seed;
}

Stockholm Stemloc::build_liberal_multiple_alignment (const Pair_CFG_alignment& seed)
{
  // print log message
  CTAG(8,STEMLOC) << "Doing liberal progressive multiple alignment\n";

  // get seed path
  const Pairwise_path seed_path = seed.pairwise_path (true);

  // initialise Decomposition
  Decomposition decomp;
  decomp [Row_pair (0,1)] = seed_path;

  // initialise multi_foldenv_db and multi_local_fold
  Local_fold_string_database multi_local_fold;
  const Local_fold_string& seed_xfold (multi_local_fold[seed.npx.name] = seed.local_xfold());
  const Local_fold_string& seed_yfold (multi_local_fold[seed.npy.name] = seed.local_yfold());
  multi_foldenv_db.clear();

  // initialise Stockholm
  Stockholm multi_align (2);
  multi_align.set_np (0, seed.npx);
  multi_align.set_np (1, seed.npy);
  multi_align.path = seed_path;
  multi_align.set_fold_string (seed.npx.name, seed_xfold.global_fold());
  multi_align.set_fold_string (seed.npy.name, seed_yfold.global_fold());
  multi_align.make_consensus_fold();

  // initialise multiple alignment score
  // (this isn't quite the log-odds score, since pairwise scores aren't conditionally normalised)
  Score multi_align_sc = seed.score;

  // min score for extending a hit
  const Score min_ext_sc = Bits2Score (min_extend_bitscore);
  
  // for each unaligned sequence, this map will contain best pairwise match to seed
  Extension_table extension_by_seq;  // indexing: extension_by_seq[seed_name][new_name]
  set<int> unscanned_rows;  // indices of rows that need to be (re)scanned for extension_by_seq
  for (int r = 0; r < multi_align.rows(); ++r)
    unscanned_rows.insert (r);
  
  // greedily grow the alignment
  while (!unaligned_seqname.empty())   // NB assumes EVERY sequence will go into the multiple alignment
    {
      // update extension_by_seq
      for_const_contents (set<int>, unscanned_rows, scan_row)
	{
	  // get sequence for newly scanned row
	  const Named_profile& npx = *multi_align.np[*scan_row];
	  // get fold envelope for newly scanned row
	  Fold_envelope xenv;
	  init_multi_fold_env (xenv, npx, multi_local_fold[npx.name]);
	  if (CTAGGING(1,STEMLOC_DEBUG))
	    {
	      CL << "init_multi_fold_env fold envelope for " << npx.name << ", local fold string [" << multi_local_fold[npx.name].start << ']' << multi_local_fold[npx.name].fold << '\n';
	      xenv.dump (CL);
	    }
	  // loop through unaligned seqs
	  for_contents (set<sstring>, unaligned_seqname, yname)
	    {
	      // get sequence
	      const Named_profile& npy (*seq_db.index.name2profile (*yname));
	      // print log message
	      CTAG(4,STEMLOC) << "Comparing aligned sequence '" << npx.name << "' and unaligned sequence '" << npy.name << "'\n";
	      // get fold envelope
	      Fold_envelope yenv;
	      init_pair_fold_env (yenv, npy, multi_nfold);
	      // make fold envelope local.
	      // FIXME: should do this more consistently throughout this module, rather than using Pair_CFG_DP_matrix::local flag
	      if (!global)
		yenv.make_local();
	      // get alignment envelope
	      Pair_envelope pair_env;
	      sstring best_param_set;
	      init_align_env (pair_env, npx, npy, xenv, yenv, best_param_set, multi_nalign);
	      // build the DP matrix
	      const Pair_CYK_matrix cyk (npx, npy, xenv, yenv, pair_cfg[best_param_set], pair_env, false);
	      if (CTAGGING(1,STEMLOC_DP))
		CL << cyk.cfg_dump();
	      // get CYK score
	      const Score cyk_sc = cyk.final_score;
	      // print log message
	      CTAG(5,STEMLOC STEMLOC_MULTI) << "Score for '" << npx.name << "' vs '" << npy.name << "' is " << Score2Bits(cyk_sc) << " bits\n";
	      // check for low score; if so, don't attempt traceback
	      if (cyk_sc <= -InfinityScore || cyk_sc < min_ext_sc)
		{
		  CTAG(4,STEMLOC STEMLOC_MULTI) << "Score for '" << npx.name << "' vs '" << npy.name << "' is too low; skipping.\n";
		  continue;
		}
	      // get CYK alignment
	      const Pair_CFG_alignment cyk_alignment = cyk.alignment();
	      // show alignment
	      if (show_intermediates)
		cyk_alignment.show (cout);
	      // store
	      extension_by_seq[*yname][npx.name] = Extension (cyk_alignment, *scan_row);
	    }
	}
      // clear unscanned_rows
      unscanned_rows.clear();
      // if there are no extensions (i.e. all scores were -infinity and/or below threshold), then bail out now
      if (extension_by_seq.size() == 0)
	{
	  CTAG(5,STEMLOC) << "No more matches; halting multiple alignment\n";
	  break;
	}
      // find best extension score
      Score best_extension_sc = -InfinityScore;
      Extension_table_cell::const_iterator best_extension = extension_by_seq.begin()->second.end();
      sstring new_name;
      for_const_contents (Extension_table, extension_by_seq, extension_tbl_cell)
	for_const_contents (Extension_table_cell, extension_tbl_cell->second, extension)
	if (extension->second.align.score > best_extension_sc)
	  {
	    best_extension_sc = extension->second.align.score;
	    best_extension = extension;
	    new_name = extension_tbl_cell->first;
	  }
      // recover best extension
      const Extension& ext (best_extension->second);
      const int seed_row = ext.seed_row;
      const Named_profile& new_np = *seq_db.index.name2profile (new_name);
      const Named_profile& seed_np = *multi_align.np[seed_row];
      const sstring seed_name = seed_np.name;  // make sure this isn't a reference
      // print log message(s)
      CTAG(4,STEMLOC STEMLOC_MULTI) << "Adding sequence '" << new_name << "' to multiple alignment with score " << Score2Bits(best_extension_sc) << " bits\n";
      if (CTAGGING(1,STEMLOC_DEBUG))
	{
	  CL << "Local fold string alignment (starting at position " << ext.align.xstart << " in " << seed_name << " and " << ext.align.ystart << " in " << new_name << ")\n";
	  CL << ext.align.xfold << '\n';
	  CL << ext.align.yfold << '\n';
	  CL << "multi_local_fold[" << seed_name << "] = [" << multi_local_fold[seed_name].start << ']' << multi_local_fold[seed_name].fold << '\n';
	  CL << "ext.align.local_xfold() = [" << ext.align.local_xfold().start << ']' << ext.align.local_xfold().fold << '\n';
	}
      // add best match to multi_align
      const int new_row = multi_align.add_row();
      multi_align.set_np (new_row, new_np);
      ScorePMulAcc (multi_align_sc, best_extension_sc);
      // compose multiple alignment
      decomp [Row_pair (seed_row, new_row)] = ext.align.pairwise_path (true);
      multi_align.path.compose (decomp, true);
      // update multi_local_fold
      multi_local_fold[seed_name].add (ext.align.local_xfold());
      multi_local_fold[new_name] = ext.align.local_yfold();
      // update multi_align fold strings
      for_const_contents (Local_fold_string_database, multi_local_fold, name_mlf)
	multi_align.set_fold_string (name_mlf->first, name_mlf->second.global_fold());
      // set consensus fold
      multi_align.make_consensus_fold();
      // show an intermediate alignment, if this is not the last unaligned sequence
      const int seqs_remaining = unaligned_seqname.size() - 1;
      if (show_intermediates && seqs_remaining > 0)
	{
	  if (global)
	    multi_align.write_Stockholm (cout);
	  else
	    multi_align.write_Stockholm_NSE (cout);
	}
      // print log message(s)
      CTAG(6,STEMLOC) << "Aligned '" << new_name << "' to '" << seed_name << "' (" << new_row+1 << " sequences aligned, " << seqs_remaining << " unaligned)\n";
      if (CTAGGING(2,STEMLOC STEMLOC_MULTI))
	{
	  if (global)
	    multi_align.write_Stockholm (CL);
	  else
	    multi_align.write_Stockholm_NSE (CL);
	}
      // add seed row & new row to the list of rows to scan for the next seed
      unscanned_rows.insert (seed_row);  // must rescan the seed row, since its fold envelope may have changed
      unscanned_rows.insert (new_row);
      // flush seed_row from multi_foldenv_db fold envelope cache, since foldenv may have changed
      if (cache_fold && multi_foldenv_db.find (seed_name) != multi_foldenv_db.end())
	multi_foldenv_db.erase (seed_name);
      // remove best extension from unaligned_seqname and extension_by_seq
      unaligned_seqname.erase (new_name);
      extension_by_seq.erase (new_name);  // MUST do this last of all, as it invalidates extension, ext & new_name
      // finally, invalidate extensions that relied on the old seed row, since its foldenv may have changed
      for_contents (Extension_table, extension_by_seq, ext_tbl_cell)
	if (ext_tbl_cell->second.find (seed_name) != ext_tbl_cell->second.end())
	  ext_tbl_cell->second.erase (seed_name);
    }

  // set score
  sstring sc_string;
  sc_string << Score2Bits (multi_align_sc);
  multi_align.add_gf_annot (sstring (Stockholm_bit_score_tag), sc_string);

  // return
  return multi_align;
}

Stockholm Stemloc::build_motif_alignment (const Pair_CFG_alignment& seed)
{
  CLOGERR << "Warning: motif alignment untested; email Ian Holmes <ihh@berkeley.edu> for beta testing info\n";

  // print log message
  CTAG(8,STEMLOC PROFALIGN) << "Doing motif-based consensus multiple alignment\n";

  // build a Covariance_model from the seed alignment
  CTAG(4,STEMLOC PROFALIGN) << "Building covariance model\n";
  Covariance_model cov;
  cov.build_from_parse_tree (seed.parse);
  if (CTAGGING(3,PROFALIGN))
    cov.show (CL);

  // create a Super_pair model, initialize w/defaults
  PScores sp_pscores;
  Super_pair sp (&sp_pscores);
  Telegraph_PScores_adaptor sp_tgio (sp_pscores);
  Stemloc_defaults::init_superpair_params (sp_tgio);

  // initialise Dirichlet_prior & use it to set initial PScores
  const double pseudocount_multipler = 1.;
  cov.init_prior (sp.single_dinuc, sp.single_nuc, sp_pscores, pseudocount_multipler);
  cov.prior.initialise();  // initialise the PScores

  // build a Pair_CFG_trainer for the seed, and run EM
  CTAG(5,STEMLOC PROFALIGN) << "Training covariance model on seed alignment\n";
  Stockade seed_stockade = seed.stockade (false);
  Stockholm_database seed_db;
  seed_db.add (seed_stockade.align);
  Pair_CFG_trainer seed_trainer (seed_db, cov.single_model, cov.prior, cov.pscores, true);
  seed_trainer.train_single (em_min_inc, em_max_iter, em_forgive);

  // build a Pair_CFG_trainer for the entire input set, and run EM
  CTAG(5,STEMLOC PROFALIGN) << "Training covariance model on all input sequences\n";
  Pair_CFG_trainer trainer (align_db, cov.single_model, cov.prior, cov.pscores, true);
  trainer.train_single (em_min_inc, em_max_iter, em_forgive);

  // initialize the Decomposition & the Stockholm alignment
  // last row is consensus sequence; we'll lop this off at the end
  Decomposition decomp;
  Stockholm multi;

  // initialise DP stuff
  CTAG(5,STEMLOC PROFALIGN) << "Aligning sequences to covariance model\n";
  const int cons_row = seq_db.index.size();
  const Pair_CFG_scores& profile_cfg = cov.single_model.eval_cfg_scores (cov.pscores);
  const Named_profile dummy_npy;
  const Fold_envelope dummy_envy (0);
  Fold_envelope envx;

  // loop over sequences
  for (int r = 0; r < seq_db.index.size(); ++r)
    {
      // do CYK alignment, get local path
      const Named_profile& npx (*seq_db.index.profile[r]);
      init_pair_fold_env (envx, npx, multi_nfold);
      const Pair_CYK_matrix cyk (npx, dummy_npy, envx, dummy_envy, profile_cfg, !global);
      const Pair_CFG_local_path local_path = cyk.traceback_with_coords();

      // align to consensus by calling Covariance_model::convert_local_path_to_consensus_alignment
      const Pairwise_path cons_ppath = cov.convert_local_path_to_consensus_alignment (local_path, npx);

      // add to alignment
      multi.add_row();
      multi.set_np (r, npx);

      // add to Decomposition
      decomp [Row_pair (cons_row, r)] = cons_ppath;
    }

  // create alignment from Decomposition
  multi.path.compose (decomp, true);

  // lop off consensus row, to avoid having to guarantee persistence of consensus Named_profile
  multi.path.erase_rows (cons_row - 1, 1);

  // return the alignment
  return multi;
}

Stemloc::Multi_match::Multi_match() { }
Stemloc::Multi_match::Multi_match (const Pair_CFG_alignment& align, int seed_row)
  : seed_row (seed_row), seed_lpad (align.xstart), score (align.score)
{
  Stockade stockade = align.stockade();
  subseq = stockade.np[1];
  path = stockade.align.path;
  fold = stockade.align.get_fold_string (subseq.name);
}

Stemloc::Extension::Extension() { }
Stemloc::Extension::Extension (const RNA_pair_path& align, int seed_row)
  : seed_row (seed_row), align (align)
{ }

Stemloc::Stemloc (int argc, char** argv) :
  opts (argc, argv),
  alphabet (CFG_alphabet),
  gramset (Stemloc_gramset())
{ }

int Stemloc::run()
{
  try
    {
      init_opts();
      parse_opts();
      input_data();
      load_gramset();
      init_fold_params();
      init_align_params();
      init_cmp_params();
      save_gramset();
      if (use_ama)
	{
	  build_multiple_alignments_ama();
	}
      else
	{
	  build_pairwise_alignments();
	  build_multiple_alignments();
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  return 0;
}
