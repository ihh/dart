#include "ecfg/ecfgmain.h"
#include "ecfg/ecfgsexpr.h"
#include "ecfg/pfold.h"
#include "ecfg/protgrammar.h"
#include "util/vector_output.h"
#include "irrev/irrev_em_matrix.h"

// maximum length for lines in Stockholm input alignments
#define MAX_LINE_LEN 12345678

// branch length resolution & max value
#define BRANCH_LENGTH_RES TINY
#define BRANCH_LENGTH_MAX 10.

// Stockholm tags for posterior probabilities
#define CONFIDENCE_TAG_PREFIX  "PC"
#define LOGPOSTPROB_TAG_PREFIX "PP"

// GFF tags for posterior probabilities
#define CYK_STATE_LABEL "CYK"  /* indicates that GFF postprob line is part of CYK trace */

// tags for training meta-information in grammar
#define TRAINING_INFO  "training-info"
#define TRAINING_TIME  "unix-time"
#define TRAINING_BITS  "final-bits"
#define ALIGN_FILENAME "alignment-filename"
#define ALIGN_STDIN    "<stdin>"

// ECFG_main
ECFG_main::ECFG_main (int argc, char** argv)
  : opts (argc, argv),
    max_subseq_len (0),
    alph (0),
    tree_estimation_grammar (0),
    tree_estimation_chain (0),
    align_db (seq_db)
{ }

ECFG_main::ECFG_main()
  : opts (0, (char**) 0),
    max_subseq_len (0),
    alph (0),
    tree_estimation_grammar (0),
    tree_estimation_chain (0),
    align_db (seq_db)
{
  init_opts("");
}

ECFG_main::~ECFG_main()
{
  delete_loaded_grammars();
  delete_trainers();
}

void ECFG_main::add_grammar (const char* name, ECFG_scores* ecfg)
{
  ecfg_map[sstring (name)] = ecfg;
  ecfg->name = name;
  grammar_list.push_back (name);
}

void ECFG_main::add_standard_grammars (const char* default_grammar)
{
  add_grammar ("rev", new Null_ECFG (DNA_alphabet, true));
  add_grammar ("irrev", new Null_ECFG (DNA_alphabet, false));
  add_grammar ("dinuc", new Nearest_neighbor_ECFG());
  add_grammar ("codon", new Codon_ECFG (true));
  add_grammar ("aa", new Protein_grammar(1,1));
  add_grammar ("aa2", new Protein_grammar(1,2));
  add_grammar ("aa3", new Protein_grammar(1,3));
  add_grammar ("aa4", new Protein_grammar(1,4));
  add_grammar ("pfold", new PFOLD_ECFG());

  default_grammars = default_grammar;
}

void ECFG_main::delete_grammars()
{
  for_contents (ECFG_map, ecfg_map, sg)
    delete (*sg).second;
}

void ECFG_main::init_opts (const char* desc)
{
  INIT_CONSTRUCTED_OPTS_LIST (opts, -1, "[options] [<alignment database(s) in Stockholm format>]", desc);

  opts.print_title ("Model selection");

  opts.add ("g -grammar", grammars_filename, "filename of grammar to use for training & annotation", false);

  if (!(preset = default_grammars).size() && grammar_list.size() > 0)
    preset = grammar_list.front();
  if (preset.size())
    {
      sstring preset_help;
      preset_help << "use one of the presets: " << sstring::join (grammar_list, ",");
      if (grammar_list.size() > 1)
	opts.add ("p -preset", preset, preset_help.c_str());
    }

  opts.add ("x -expand", dump_expanded = "", "dump macro-expanded grammar to file (prior to any training)", false);
  opts.add ("e -tree-grammar", tree_grammar_filename = "", "grammar to use for tree estimation (defaults to the training/annotation grammar)", false);
  opts.add ("l -length", max_subseq_len = -1, "limit maximum length of infix subseqs (context-free grammars only)", false);

  sstring gap_chars_help;
  gap_chars_help << "set character(s) that are used to denote gaps in alignment (default is \"" << DEFAULT_GAP_CHARS << "\")";
  opts.add ("gc -gapchar", gap_chars = "", gap_chars_help.c_str(), false);

  opts.newline();
  opts.print_title ("Tree estimation algorithms");

  opts.add ("nj -neighbor-joining", do_neighbor_joining = false, "do neighbor-joining to estimate trees for alignments which have none already annotated");
  opts.add ("obl -optimize-branch-lengths", do_branch_length_EM = false, "optimize branch lengths in trees, using EM");
  opts.add ("ps -point-sub", avoid_ECFG_for_branch_length_EM = false, "modifies --optimize-branch-lengths option to use point-substitution model, rather than entire phylogrammar");
  opts.add ("at -attach", attach_rows = false, "attempt to place unattached alignment rows on the tree by adding extra branches");

  opts.newline();
  opts.print_title ("Precision of tree estimation");

  opts.add ("bmin -branch-min", min_branch_len = .0001, "minimum branch length in phylogenies");
  opts.add ("bres -branch-resolution", tres = .0001, "resolution of branch lengths in phylogenies");

  opts.newline();
  opts.print_title ("Parameter estimation algorithms");

  opts.add ("t -train", train, "use EM algorithm to \"train\" model (i.e. fit it to alignment data), saving trained grammar to file", false);
  opts.add ("tlog -training-log", training_log_filename = "", "use EM to train model, dumping every intermediate grammar to this file", false);

  opts.newline();
  opts.print_title ("Pseudocounts for parameter estimation");
  opts.print ("(Pseudocount options are obsolescent; they are superceded by pseudocount tags within the grammar. See biowiki.org/XrateFormat)\n\n");

  opts.add ("pi -pseudinitial", pseud_init = 1e-9, "pseudocount for initial state occupancies and probability parameters");
  opts.add ("pm -pseudmutate", pseud_mutate = 0., "pseudocount for mutation rates and rate parameters");
  opts.add ("pt -pseudtime", pseud_wait = 1e-4, "pseudo-wait time for mutation rates and rate parameters");

  opts.newline();
  opts.print_title ("Convergence criteria for tree/parameter estimation");

  opts.add ("mr -maxrounds", em_max_iter = -1, "max number of \"rounds\" (iterations) of EM", false);
  opts.add ("mi -mininc", em_min_inc = .001, "minimum fractional increase in log-likelihood per round of EM");
  opts.add ("f -forgive", em_forgive = 0, "number of consecutive non-increasing rounds of EM to \"forgive\" before stopping");

  opts.newline();
  opts.print_title ("Annotation algorithms");

  opts.add ("a -annotate", annotate = true, "generate #=GC, GFF and/or WIG annotations, running CYK/Inside/Outside algorithms as appropriate");
  opts.add ("ms -maxscore", report_maxscore = false, "report CYK log-likelihood, corresponding to maximum-likelihood parse tree");
  opts.add ("s -score", report_sumscore = false, "report Inside log-likelihood, corresponding to a sum over all parse trees");
  opts.add ("c -confidence", report_confidence = false, "report Inside-Outside posterior log-probabilities of nodes in CYK parse tree");
  opts.add ("pp -postprob", report_postprob = false, "report Inside-Outside posterior log-probabilities of all possible parse tree nodes");
  opts.add ("hc -hidden-classes", report_hidden_classes = false, "impute ML hidden classes at each site (for substitution models with hidden classes)");

  opts.newline();
  opts.print_title ("Annotation output");
  opts.add ("gff", gff_filename, "save GFF annotations to file, rather than interleaving into Stockholm output", false);
  opts.add ("wig -wiggle", wiggle_filename, "save Wiggle annotations to file, rather than interleaving into Stockholm output", false);

  opts.newline();
  opts.print_title ("Acceleration of annotation DP algorithms (experimental)");

  opts.add ("fp -fast-prune", use_fast_prune = false, "attempt pruning algorithm in probability-space, rather than log-space (caveat: prone to underflow)");
#ifdef BEAGLE_INCLUDED
  opts.add ("bgl -beagle", use_beagle = false, "use Beagle GPU library to do pruning");
#else /* BEAGLE_INCLUDED */
  opts.add ("bgl -beagle", use_beagle = false);   // if Beagle not compiled, then allow this option, but don't advertise it (using it will cause a warning)
#endif /* BEAGLE_INCLUDED */

  opts.newline();
  opts.print_title ("Ancestral reconstruction");

  opts.add ("ar -ancrec-cyk-map", ancrec_CYK_MAP = false, "reconstruct Maximum-A-Posteriori ancestral sequences, using CYK parse tree");
  opts.add ("arpp -ancrec-postprob", ancrec_postprob = false, "report posterior probabilties of alternate reconstructions on CYK parse tree");
  opts.add ("marp -min-ancrec-prob", min_ancrec_postprob = .01, "minimum probability to report for --ancrec-postprob option");

}

void ECFG_main::annotate_loglike (Stockholm& stock, const char* tag, const sstring& ecfg_name, Loge loglike) const
{
  sstring score_tag, score_string;
  score_tag << tag << "_" << ecfg_name;
  score_string << Nats2Bits (loglike) << " bits";
  stock.add_gf_annot (score_tag, score_string);
}

void ECFG_main::parse_opts()
{
  // parse the command line
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }
}

void ECFG_main::read_alignments (const Stockholm& stock)
{
  Stockholm_database dummy_stock_db;
  dummy_stock_db.add (stock);
  read_alignments (&dummy_stock_db);
}

void ECFG_main::read_alignments (const Stockholm_database* stock)
{
  // set up the gap characters
  if (gap_chars.size())
    Alignment::set_gap_chars (gap_chars);

  // initialise Stockholm_database
  if (stock != 0)
    {
      stock_db.add (*stock);
      training_alignment_filename.push_back ("-");  // default "filename" for alignments supplied through the API
    }
  else if (!opts.args.size())
    {
      // no alignment filenames specified; read from stdin
      CLOGERR << "[waiting for alignments on standard input]\n";
      stock_db.read (cin, seq_db, MAX_LINE_LEN);
      training_alignment_filename.push_back (sstring (ALIGN_STDIN));  // default "filename" for alignments read from stdin
    }
  else
    for_const_contents (vector<sstring>, opts.args, align_db_filename)
      {
	ifstream align_db_in ((const char*) align_db_filename->c_str());
	if (!align_db_in) THROWEXPR ("Couldn't open alignment file '" << *align_db_filename << "'");
	stock_db.read (align_db_in, seq_db, MAX_LINE_LEN);
	training_alignment_filename.push_back (*align_db_filename);
      }

  // initialise Tree_alignment_database
  align_db.initialise_from_Stockholm_database (stock_db, false);
}

void ECFG_main::estimate_trees (SExpr* grammar_alphabet_sexpr, Sequence_database* seq_db_ptr) {
  // HACK REPORT
  // Score_profile's are needed by the Tree_alignment tree estimation routines
  // (neighbor joining & branch-length EM).
  // thus, if any trees are missing, or if we are optimising branch lengths,
  // we first convert sequences to Score_profile's.
  // this is hacky, but appeases the conflicting needs of keeping memory low (for genomic seqs)
  // and being able to supply alignments without trees (for protein seqs and small nucleotide alignments).
  //
  // ADDENDUM: July 10, 2006; IH. Situation getting hackier.
  // Since the rate matrix for tree estimation is now in a separate grammar file,
  // it's conceivable that the Score_profile's could end up getting converted to the "wrong" alphabet here.
  // That is, it'll be the right Alphabet for tree estimation,
  // but wrong for any later operations involving the main grammar.
  // Workaround: call seq_db_ptr->clear_scores() after tree algorithms.

  // use our Sequence_database unless caller overrides this
  if (seq_db_ptr == 0)
    seq_db_ptr = &seq_db;

  // read tree-estimation grammar & alphabet from file or SExpr
  bool mark_for_deletion = true;
  if (grammar_alphabet_sexpr)
    ECFG_builder::load_xgram_alphabet_and_grammars (*grammar_alphabet_sexpr, tree_estimation_grammar_alphabet, tree_estimation_grammars);
  else if (tree_grammar_filename.size())
    ECFG_builder::load_xgram_alphabet_and_grammars (tree_grammar_filename, tree_estimation_grammar_alphabet, tree_estimation_grammars);
  else if (grammars_filename.size())
    ECFG_builder::load_xgram_alphabet_and_grammars (grammars_filename, tree_estimation_grammar_alphabet, tree_estimation_grammars);
  else
    {
      // initialise tree grammars from preset (and do NOT mark them for later deletion)
      ECFG_map::iterator ecfg_iter = ecfg_map.find (preset);
      if (ecfg_iter == ecfg_map.end())
	THROWEXPR ("Preset grammar '" << preset << "' not known (did you mean to load a grammar from a file?)");
      tree_estimation_grammars.push_back ((*ecfg_iter).second);
      ((Alphabet&) tree_estimation_grammar_alphabet) = tree_estimation_grammars[0]->alphabet;
      mark_for_deletion = false;
    }

  // mark grammars for later deletion
  if (mark_for_deletion)
    grammars_to_delete.insert (grammars_to_delete.begin(), tree_estimation_grammars.begin(), tree_estimation_grammars.end());

  // set the matrix used for tree estimation to be the first single-character matrix in any grammar in the tree grammar file
  for (int g = 0; tree_estimation_chain == 0 && g < (int) tree_estimation_grammars.size(); ++g)
    {
      tree_estimation_grammar = tree_estimation_grammars[g];
      tree_estimation_chain = tree_estimation_grammar->first_single_pseudoterminal_chain();
    }

  // turn on tree estimation algorithms if it looks like the user wanted that
  if (missing_trees() && !do_neighbor_joining)
    {
      CLOGERR << "Warning: you did not request neighbor-joining, and input alignments did not include trees.\n"
	      << "Activating neighbor-joining and branch-length optimization by default.\n";
      do_neighbor_joining = do_branch_length_EM = true;
    }

  if (tree_grammar_filename.size() && !do_neighbor_joining && !do_branch_length_EM)
    {
      CLOGERR << "Warning: you specified a tree estimation grammar, but no tree estimation algorithms.\n"
	      << "Activating neighbor-joining and branch-length optimization by default.\n";
      do_neighbor_joining = do_branch_length_EM = true;
    }

  // figure out what we need to do
  const bool need_tree_estimation_grammar = do_neighbor_joining || do_branch_length_EM;
  const bool need_tree_estimation_chain = do_neighbor_joining || (do_branch_length_EM && avoid_ECFG_for_branch_length_EM);
  const bool convert_seq_db = do_neighbor_joining || do_branch_length_EM || attach_rows;

  updated_trees = do_neighbor_joining || do_branch_length_EM || attach_rows;   // this flag can also be set later if any branch lengths are rounded up

  sstring error_preamble;
  if (missing_trees())
    error_preamble << "Input alignments do not include Newick trees.\n";

  if (need_tree_estimation_grammar && tree_estimation_grammar == 0)
    THROWEXPR (error_preamble << "In order to estimate trees by neighbor-joining, optimize branch lengths,\n"
	       "or place unattached rows on the tree, a tree-estimation grammar\n"
	       "must be specified.\n");

  if (need_tree_estimation_chain && tree_estimation_chain == 0)
    THROWEXPR (error_preamble << "The tree-estimation grammar did not contain a single-terminal chain\n"
	       "(i.e. a point substitution matrix), which is needed for some operations\n"
	       "(neighbor-joining, or branch-length optimization with the --point-sub option).\n");

  if (convert_seq_db)
    {
      // convert sequences to Score_profile's for tree estimation
      tree_estimation_hidden_alphabet.init_hidden (tree_estimation_grammar_alphabet, tree_estimation_chain->class_labels);
      seq_db_ptr->seqs2scores (tree_estimation_hidden_alphabet);
    }

  if (do_neighbor_joining && missing_trees())
    {
      // estimate missing trees by neighbor-joining
      CTAG(6,XRATE) << "Estimating missing trees by neighbor-joining\n";
      Subst_dist_func_factory dist_func_factory (*tree_estimation_chain->matrix);
      align_db.estimate_missing_trees_by_nj (dist_func_factory);
      copy_trees_to_stock_db();
    }

  if (do_branch_length_EM)
    {
      // optimise branch lengths by EM
      CTAG(6,XRATE) << "Optimizing tree branch lengths by EM\n";
      if (avoid_ECFG_for_branch_length_EM)
	{
	  Irrev_EM_matrix nj_hsm (1, 1);  // create temporary EM_matrix
	  nj_hsm.assign (*tree_estimation_chain->matrix);
	  align_db.optimise_branch_lengths_by_EM (nj_hsm, 0., em_max_iter, em_forgive, em_min_inc, BRANCH_LENGTH_RES, BRANCH_LENGTH_MAX, min_branch_len);
	}
      else
	align_db.optimise_branch_lengths_by_ECFG_EM (*tree_estimation_grammar, 0., em_max_iter, em_forgive, em_min_inc, BRANCH_LENGTH_RES, BRANCH_LENGTH_MAX, min_branch_len);
      copy_trees_to_stock_db();
    }

  if (attach_rows && unattached_rows())
    {
      // attach all alignment rows to tree
      CTAG(6,XRATE) << "Placing unattached alignment rows on trees\n";
      align_db.attach_rows (*tree_estimation_grammar, 0., BRANCH_LENGTH_RES, BRANCH_LENGTH_MAX, min_branch_len);
      copy_trees_to_stock_db();
    }

  if (convert_seq_db)
    {
      // clear the (potentially inconsistent with other grammars) Score_profile's
      seq_db_ptr->clear_scores();
    }

  // ensure all tree branches meet minimum length requirement
  for (int n_align = 0; n_align < align_db.size(); n_align++)
    {
      bool adjusted_lengths = false;
      PHYLIP_tree& tree = align_db.tree_align[n_align]->tree;
      for_rooted_branches_post (tree, b)
	if ((*b).length < min_branch_len)
	  {
	    tree.branch_length (*b) = min_branch_len;
	    adjusted_lengths = updated_trees = true;
	  }

      if (adjusted_lengths)
	copy_trees_to_stock_db();
    }
}

void ECFG_main::copy_trees_to_stock_db()
{
  // copy all trees back into the Stockholm alignments (ugh)
  for (int n_align = 0; n_align < align_db.size(); n_align++)
    {
      Tree_alignment& tree_align = *align_db.tree_align[n_align];
      Stockholm& stock = *stock_db.align_index[n_align];
      tree_align.copy_tree_to_Stockholm (stock);
    }
}

void ECFG_main::read_grammars (SExpr* grammar_alphabet_sexpr)
{
  // get list of grammars; either by loading from file, or by using a preset
  if (grammar_alphabet_sexpr) {
    ECFG_builder::load_xgram_alphabet_and_grammars (*grammar_alphabet_sexpr, user_alphabet, grammar, &align_db, max_subseq_len, tres);
    alph = &user_alphabet;

    // mark grammars for later deletion
    grammars_to_delete.insert (grammars_to_delete.begin(), grammar.begin(), grammar.end());

  } else if (grammars_filename.size())
    {
      ECFG_builder::load_xgram_alphabet_and_grammars (grammars_filename, user_alphabet, grammar, &align_db, max_subseq_len, tres);
      alph = &user_alphabet;

      // mark grammars for later deletion
      grammars_to_delete.insert (grammars_to_delete.begin(), grammar.begin(), grammar.end());
    }
  else
    {
      // initialise grammars from preset (and do NOT mark them for later deletion)
      ECFG_map::iterator ecfg_iter = ecfg_map.find (preset);
      if (ecfg_iter == ecfg_map.end())
	THROWEXPR ("Preset grammar '" << preset << "' not known (did you mean to load a grammar from a file?)");
      grammar.push_back ((*ecfg_iter).second);

      // update timepoint_res for loaded preset grammar (hacky)
      for_const_contents (vector<ECFG_chain>, (*ecfg_iter).second->matrix_set.chain, ec)
	ec->matrix->timepoint_res = tres;

      // initialise alphabet from grammar (presets use global singleton alphabet objects)
      alph = &grammar[0]->alphabet;
    }

  // save macro-expanded grammar (prior to training)
  if (dump_expanded.size())
    {
      // save
      ofstream dump_file (dump_expanded.c_str());
      for (int n = 0; n < (int) grammar.size(); ++n)
	ECFG_builder::ecfg2stream (dump_file, *alph, *grammar[n], (ECFG_counts*) 0);
      ECFG_builder::alphabet2stream (dump_file, *alph);
    }
}

void ECFG_main::convert_sequences()
{
  // create Aligned_score_profile's using the alphabet, for training & annotation
  vector<Aligned_score_profile> av (align_db.size());
  asp_vec.swap(av);
  for (int n = 0; n < align_db.size(); ++n)
    asp_vec[n].init (align_db.tree_align[n]->align, stock_db.align_index[n]->np, *alph);
}

void ECFG_main::train_grammars()
{
  // do training
  if (missing_trees())
    THROWEXPR("All alignments must be annotated with trees before training can begin");

  // train
  trainer = vector<ECFG_trainer*> (grammar.size(), (ECFG_trainer*) 0);
  if (stock_db.size()) {

    // init training log
    ofstream* training_log = 0;
    if (training_log_filename.size())
      {
	training_log = new ofstream(training_log_filename.c_str());
	ECFG_builder::alphabet2stream (*training_log, *alph);
      }

    // loop over grammars
    for (int n_grammar = 0; n_grammar < (int) grammar.size(); ++n_grammar)
      {
	ECFG_scores& ecfg = *grammar[n_grammar];
	ECFG_trainer*& tr (trainer[n_grammar]);

	tr = new ECFG_trainer (ecfg, stock_db.align_index, asp_vec, ecfg.is_left_regular() || ecfg.is_right_regular() ? 0 : max_subseq_len, training_log);

	tr->pseud_init = pseud_init;
	tr->pseud_mutate = pseud_mutate;
	tr->pseud_wait = pseud_wait;

	tr->em_min_inc = em_min_inc;
	tr->em_max_iter = em_max_iter;

	tr->iterate_quick_EM (em_forgive);
      }

    if (training_log)
      {
	training_log->close();
	delete training_log;
	training_log = 0;
      }
  }

  // save
  if (train.size())
    {
      ofstream train_file (train.c_str());
      for (int n = 0; n < (int) grammar.size(); ++n)
	{
	  // add training meta-info
	  sstring training_meta;
	  training_meta << TRAINING_INFO
			<< " (" << TRAINING_TIME << ' ' << ((unsigned long) time ((time_t*) 0))
			<< ") (" << TRAINING_BITS << ' ' << Nats2Bits (trainer[n]->best_loglike)
			<< ") (" << ALIGN_FILENAME << ' ' << training_alignment_filename
			<< ')';
	  SExpr training_meta_sexpr (training_meta.c_str());
	  grammar[n]->transient_meta.push_back (training_meta_sexpr);
	  // save
	  ECFG_builder::ecfg2stream (train_file, *alph, *grammar[n], stock_db.size() ? &trainer[n]->counts : (ECFG_counts*) 0);
	}
      ECFG_builder::alphabet2stream (train_file, *alph);
    }
}

void ECFG_main::delete_trainers()
{
  for_const_contents (vector<ECFG_trainer*>, trainer, t)
    if (*t)
      delete *t;

  trainer.clear();
}

void ECFG_main::delete_loaded_grammars()
{
  for_contents (vector<ECFG_scores*>, grammars_to_delete, g)
    delete *g;
}

void ECFG_main::annotate_alignments (ostream* align_stream)
{
  if (missing_trees())
    THROWEXPR("All alignments must be annotated with trees before annotation can begin");

  // open GFF file stream
  ostream* gff_stream = 0;
  bool created_gff_stream = false;
  if (gff_filename.size())
    {
      bool has_GFF = false;
      for_const_contents (vector<ECFG_scores*>, grammar, g)
	if ((*g)->has_GFF())
	  has_GFF = true;
      if (!has_GFF)
	THROWEXPR ("GFF filename specified, but no GFF annotations in grammar");
      gff_stream = new ofstream (gff_filename.c_str());
      created_gff_stream = true;
    }

  // annotation loop over alignments
  for (int n_align = 0; n_align < align_db.size(); n_align++)
    {
      // get Tree_alignment
      const Tree_alignment& tree_align = *align_db.tree_align[n_align];
      const int seqlen = tree_align.align.columns();

      // get the original Stockholm alignment
      Stockholm* stock = stock_db.align_index[n_align];

      // get the alignment ID, if it has one
      sstring align_id = stock->get_name();
      if (!align_id.size())
	align_id << "Alignment" << n_align+1;

      // print log message
      CTAG(7,XRATE) << "Processing alignment " << align_id
		    << " (" << n_align+1 << " of " << align_db.size()
		    << "): " << seqlen << " columns\n";

      // get the Aligned_score_profile
      const Aligned_score_profile& asp = asp_vec[n_align];

      // clear the #=GF annotation for the Stockholm alignment,
      // then add the tree as a "#=GF NH" line
      const sstring nh_tag (Stockholm_New_Hampshire_tag);
      stock->clear_gf_annot (nh_tag);
      ostringstream tree_stream;
      tree_align.tree.write (tree_stream, 0);
      sstring tree_string (tree_stream.str());
      tree_string.chomp();
      stock->add_gf_annot (nh_tag, tree_string);

      // loop over grammars
      for (int n_grammar = 0; n_grammar < (int) grammar.size(); ++n_grammar)
	{
	  ECFG_scores& ecfg = *grammar[n_grammar];
	  const sstring& ecfg_name = ecfg.name;

	  // create fold envelope
	  ECFG_auto_envelope env (*stock, ecfg, max_subseq_len);

	  // create GFF container
	  GFF_list gff_list;

	  // create null pointer to persistent inside matrix
	  ECFG_inside_matrix* inside_mx = 0;
	  bool delete_inside_mx = false;

	  // create null pointer to persistent inside-outside matrix
	  ECFG_inside_outside_matrix* inout_mx = 0;

	  // create empty CYK traceback & emit_loglike array
	  ECFG_cell_score_map cyk_trace;
	  ECFG_EM_matrix::Emit_loglike_matrix cyk_emit_loglike;
	  bool cyk_emit_loglike_initialized = false;

	  // decide what annotations we want
	  const bool want_ancestral_reconstruction = ancrec_CYK_MAP || ancrec_postprob;
	  const bool want_hidden_classes = report_hidden_classes && ecfg.has_hidden_classes();
	  const bool want_GC = annotate && ecfg.has_GC();
	  const bool want_GFF = annotate && ecfg.has_GFF();
	  const bool want_wiggle = (annotate || wiggle_filename.size() > 0) && ecfg.has_wiggle();

	  // decide whether we need CYK, Inside, Outside, peeling
	  const bool want_CYK = report_maxscore || want_GC || want_GFF;
	  const bool want_outside = report_postprob || report_confidence || want_wiggle || want_hidden_classes || want_ancestral_reconstruction;
	  const bool want_inside = report_sumscore || want_outside;
	  const bool want_fill_down = want_hidden_classes || want_ancestral_reconstruction;
	  const bool want_fast_prune = use_fast_prune && !want_fill_down;

	  // log what we're doing
	  CTAG(6,XRATE) << "Desired annotations:"
			<< (want_ancestral_reconstruction ? " Ancestral_reconstruction" : "")
			<< (want_hidden_classes ? " Hidden_classes" : "")
			<< (want_GC ? " Stockholm(#=GC)" : "")
			<< (want_GFF ? " GFF" : "")
			<< (want_wiggle ? " Wiggle" : "")
			<< (report_maxscore ? " CYK_score" : "")
			<< (report_sumscore ? " Inside_score" : "")
			<< (report_confidence ? " CYK_postprobs" : "")
			<< (report_postprob ? " All_postprobs" : "")
			<< '\n';

	  CTAG(6,XRATE) << "Selected algorithms:"
			<< (want_CYK ? " CYK" : "")
			<< (want_inside ? " Inside" : "")
			<< (want_outside ? " Outside" : "")
			<< (want_fill_down ? " Peeling" : "")
			<< '\n';

	  // for backward compatibility, issue a warning if inside-outside data unavailable to GFF
	  if (want_GFF && !want_outside)
	    CLOGERR << "Warning: --score, --postprob or --confidence options not specified.\n"
		    << "GFF annotations will not contain Inside likelihoods or Inside-Outside posterior probabilities.\n";

	  // if needed, do CYK algorithm
	  bool zero_likelihood = false;  // this flag becomes set if final likelihood is 0, i.e. no traceback path
	  if (want_CYK)
	    {
	      CTAG(6,XRATE) << "Annotating using grammar '" << ecfg_name << "'\n";

	      // create CYK matrix; get score & traceback; save emit_loglike; & delete matrix
	      ECFG_CYK_matrix* cyk_mx = new ECFG_CYK_matrix (ecfg, *stock, asp, env, want_fast_prune);
	      if (use_beagle)
		cyk_mx->compute_phylo_likelihoods_with_beagle();
	      cyk_mx->fill();
	      const Loge cyk_loglike = cyk_mx->final_loglike;
	      zero_likelihood = (cyk_loglike <= -InfinityLoge);
	      if (zero_likelihood)
		CLOGERR << "Alignment likelihood P(" << align_id << "|" << ecfg_name << ") is zero; skipping annotation step\n";
	      else
		cyk_trace = cyk_mx->traceback();

	      cyk_emit_loglike.swap (cyk_mx->emit_loglike);
	      cyk_emit_loglike_initialized = true;

	      delete cyk_mx;
	      cyk_mx = 0;

	      // if needed, do inside-outside as well as CYK
	      if (want_outside && !zero_likelihood)
		{
		  inout_mx = new ECFG_inside_outside_matrix (ecfg, *stock, asp, env, (ECFG_counts*) 0);
		  inout_mx->inside.use_precomputed (cyk_emit_loglike);  // re-use emit log-likelihoods calculated by CYK
		  cyk_emit_loglike_initialized = false;  // safeguard against trying to use_precomputed twice
		  inout_mx->outside.want_substitution_counts = want_fill_down;  // only call fill_down if necessary
		  inout_mx->fill();

		  zero_likelihood = zero_likelihood || (inout_mx->inside.final_loglike <= -InfinityLoge);

		  if (zero_likelihood)
		    CLOGERR << "Alignment likelihood P(" << align_id << "|" << ecfg_name << ") is zero; skipping posterior probability annotation step\n";
		  else
		    {
		      if (report_postprob)
			{
			  sstring pp_tag;
			  pp_tag << LOGPOSTPROB_TAG_PREFIX << '_' << ecfg_name;
			  sstring trace_tag (CYK_STATE_LABEL);
			  inout_mx->annotate_all_post_state_ll (gff_list, align_id, cyk_trace, trace_tag);
			}

		      if (report_confidence)
			{
			  sstring pp_tag;
			  pp_tag << CONFIDENCE_TAG_PREFIX << '_' << ecfg_name;
			  inout_mx->annotate (*stock, gff_list, align_id, cyk_trace, pp_tag);
			}

		      if (want_hidden_classes)
			inout_mx->annotate_hidden_classes (*stock, cyk_trace);

		      if (want_ancestral_reconstruction)
			inout_mx->inside.reconstruct_MAP (*stock, cyk_trace, CYK_MAP_reconstruction_tag, ancrec_CYK_MAP, ancrec_postprob);
		    }
		}

	      // since we have computed the CYK score, add it as annotation, regardless of whether it was requested
	      annotate_loglike (*stock, ECFG_max_score_tag, ecfg_name, cyk_loglike);
	    }

	  // if needed, do Inside algorithm (or retrieve previously-computed Inside matrix)
	  if (want_inside)
	    {
	      CTAG(6,XRATE) << "Running Inside algorithm using grammar '" << ecfg_name << "'\n";

	      // create Inside matrix, or get score from previously created matrix
	      if (inout_mx)
		inside_mx = &inout_mx->inside;
	      else
		{
		  inside_mx = new ECFG_inside_matrix (ecfg, *stock, asp, env, true);
		  delete_inside_mx = true;

		  if (cyk_emit_loglike_initialized)  // can we re-use emit log-likelihoods calculated by CYK?
		    {
		      inside_mx->use_precomputed (cyk_emit_loglike);
		      cyk_emit_loglike_initialized = false;
		    }
		  else if (use_beagle)
		    inside_mx->compute_phylo_likelihoods_with_beagle();

		  inside_mx->fill();
		}

	      // TODO:
	      // add stochastic traceback method to ECFG_inside_matrix
	      // create vector<ECFG_cell_score_map> for holding N inside tracebacks
	      // get the tracebacks
	      // do the MAP ancestral reconstructions (rename ancrec_CYK_MAP to ancrec_MAP)

	      // since we have computed the Inside score, add it as annotation, regardless of whether it was requested
	      annotate_loglike (*stock, ECFG_sum_score_tag, ecfg_name, inside_mx->final_loglike);
	    }

	  // annotate CYK traceback *after* DP is finished, to make use of inside-outside probs if available
	  if (want_GC)
	    ecfg.annotate (*stock, cyk_trace);

	  if (want_GFF)
	    ecfg.make_GFF (gff_list, cyk_trace, align_id.c_str(), inout_mx, inside_mx, false);

	  if (want_wiggle && inout_mx != 0)
	    {
	      if (wiggle_filename.size())
		{
		  ofstream wiggle_file (wiggle_filename.c_str());
		  ecfg.make_wiggle (wiggle_file, env, *inout_mx, align_id.c_str());
		}
	      else
		{
		  sstring wiggle_string;
		  ecfg.make_wiggle (wiggle_string, env, *inout_mx, align_id.c_str());
		  stock->add_multiline_gf_annot (sstring(Stockholm_WIG_tag), wiggle_string);
		}
	    }

	  // flush GFF
	  if (!gff_list.empty())
	    {
	      if (gff_stream)
		*gff_stream << gff_list;
	      else
		stock->add_gff (gff_list);
	    }

	  // delete matrices
	  if (inout_mx) {
	    delete inout_mx;
	    inout_mx = 0;
	  }

	  if (delete_inside_mx) {
	    delete inside_mx;
	    inside_mx = 0;
	  }
	}

      // print out the annotated Stockholm alignment
      if (align_stream)
	stock->write_Stockholm (*align_stream);

    } // end of annotation loop over alignments
  
  // close the GFF stream
  if (created_gff_stream) {
    delete gff_stream;
    gff_stream = 0;
  }
}

bool ECFG_main::missing_trees() const
{
  for_const_contents (list<Tree_alignment>, align_db.tree_align_list, ta)
    if (ta->tree.nodes() == 0)
      return true;
  return false;
}

bool ECFG_main::unattached_rows() const
{
  return align_db.has_unattached_rows();
}

// run methods
void ECFG_main::run_xrate (ostream& alignment_output_stream)
{
  parse_opts();
  read_alignments();
  estimate_trees();
  read_grammars();
  convert_sequences();

  if (train.size() || training_log_filename.size())
    {
      train_grammars();
      delete_trainers();
    }

  if (annotate || report_sumscore || updated_trees)
    annotate_alignments (&alignment_output_stream);

  // void return value
  // annotated alignments written to alignment_output_stream
  // trained grammars & GFF annotations written to files
}

Stockholm& ECFG_main::run_tree_estimation (Stockholm& stock, Sequence_database& stock_seq_db, SExpr& grammar_alphabet_sexpr)
{
  read_alignments (stock);
  estimate_trees (&grammar_alphabet_sexpr, &stock_seq_db);

  // return
  return *stock_db.align.begin();
}

Stockholm& ECFG_main::run_alignment_annotation (Stockholm& stock, SExpr& grammar_alphabet_sexpr)
{
  read_alignments (stock);
  read_grammars (&grammar_alphabet_sexpr);
  convert_sequences();
  annotate_alignments();

  // return
  return *stock_db.align.begin();
}

void ECFG_main::run_grammar_training (Stockholm_database& stock, SExpr& grammar_alphabet_sexpr, ECFG_scores** grammar_ret, ECFG_counts** counts_ret)
{
  read_alignments (&stock);
  read_grammars (&grammar_alphabet_sexpr);
  convert_sequences();
  train_grammars();

  // set return values
  if (grammar.size() > 0)
    *grammar_ret = *grammar.begin();
  else
    *grammar_ret = 0;

  if (trainer.size() > 0 && *trainer.begin() != 0)
    *counts_ret = &(*trainer.begin())->counts;
  else
    *counts_ret = 0;
}

ECFG_scores* ECFG_main::run_macro_expansion (SExpr& grammar_alphabet_sexpr) {
  read_grammars (&grammar_alphabet_sexpr);
  return grammar.size() > 0 ? grammar[0] : (ECFG_scores*) 0;
}

ECFG_scores* ECFG_main::run_macro_expansion (Stockholm& stock, SExpr& grammar_alphabet_sexpr) {
  read_alignments (stock);
  read_grammars (&grammar_alphabet_sexpr);
  return grammar.size() > 0 ? grammar[0] : (ECFG_scores*) 0;
}
