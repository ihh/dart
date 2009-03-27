#include "seq/stockholm.h"
#include "tree/tree_alignment.h"
#include "scfg/foldenv.h"
#include "scfg/postenv.h"
#include "scfg/paircfgdp.h"
#include "indiegram/tripletscfgdp.h"
#include "indiegram/igmain.h"

Indiegram::Indiegram (int argc, char** argv)
  : opts (argc, argv),
    alphabet (SCFG_alphabet),
    qs_pscore(), // parameters for pre-folding
    qstem (qs_pscore),
    qs_prior (qstem.default_prior())
{ }

void Indiegram::init_opts()
{
  INIT_CONSTRUCTED_OPTS_LIST (opts, -1, "g[options] <sequence file>",
			      "(Global) RNA structural alignment of three sequences");

  opts.newline();
  opts.print_title ("Fold envelope generation options");
  opts.add ("nf -nfold", nfold = 100, "number of folds to sample; -1 to unlimit");
  opts.add ("rf -rndfold", random_fold = false, "sample folds stochastically, rather than using <nf> best hits");
  opts.add ("len -maxstemlen", max_subseq_len = -1, "maximum length of envelope subsequences; -1 to unlimit");
  opts.add ("loop -minlooplen", min_loop_len = 0, "minimum length of loop regions");

  opts.newline();
  opts.print_title ("Posterior probability fold dotplots");
  opts.add ("pff -postfoldfile", postfold_prefix = "", "filename prefix for fold posterior probability matrices", FALSE);

  opts.newline();
  opts.print_title ("Grammar choices");
  opts.add ("test -testgrammar", use_test_scfg = 0, "use Triplet_SCFG from testtripletscfgdp.cc rather than TKFST");

}

void Indiegram::parse_opts()
{
  opts.parse_or_die();
}

// see ECFG_main::run() for similar code.
void Indiegram::input_data()
{

  // initialize Stockholm_database and store sequences in seq_db
  const vector<sstring>& alignment_filenames = opts.args;
  if (!alignment_filenames.size())
    {
      // no alignment filenames give; read from stdin
      CLOGERR << "[waiting for alignments on standard input]\n";
      stock_db.read (cin, seq_db);
    }
  else
    {
      // read in the sequences as stockholm alignment
      for_const_contents (vector<sstring>, alignment_filenames, sf)
	{
	  ifstream in (sf->c_str());
	  if (!in) THROWEXPR ("Sequence file not found: '" << *sf << "'");
	  stock_db.read (in, seq_db); // read in the sequences and store in seq_db
	}
    }
  // check that names are unique
  seq_db.index.assert_names_unique (TRUE);
  // update sequence database index
  seq_db.update_index();
  seq_db.update_alphabet (&alphabet);
  seq_db.seqs_update (alphabet, (Profile_flags_enum::Profile_flags) (Profile_flags_enum::DSQ | Profile_flags_enum::SCORE));

  // get the null model
  null_emit_prob = seq_db.get_null_model (alphabet.size());

}

// this is done a la stemloc, and in fact explicitly relies on such
void Indiegram::init_prefold_scfg()
{
  // create Telegraph I/O adaptor for parameters
  Telegraph_PScores_adaptor qs_tgio (qs_prior);

  // initialize parameters
  Indiegram_defaults::init_superstem_params (qs_tgio);

  // initialize pair SCFG scores
  qs_cfg = qstem.eval_cfg_scores (qs_pscore);

}

Triplet_SCFG Indiegram::init_triplet_scfg (double branch_length_x, double branch_length_y, double branch_length_z)
{
  Triplet_SCFG triplet_scfg;

  if (use_test_scfg) {
    triplet_scfg = Indiegram_defaults::init_test_scfg();
  }
  else {
    CTAG (6, INDIEGRAM) << "Using TKFST triplet grammar with default parameters.\n";
    triplet_scfg = Indiegram_defaults::init_tkfst_scfg (branch_length_x, branch_length_y, branch_length_z);
  }

  return triplet_scfg;
}


void Indiegram::init_foldenv (Fold_envelope& foldenv, const Named_profile& np, int nfold, const sstring& foldstring /* = "" */)
{
  // log message
  CTAG(6, INDIEGRAM) << "Finding fold envelope for sequence '" << np.name << "'\n";
  
  // initialize envelope from foldstring if passed
  if (foldstring != "")
    {
      // create dummy subseq to test foldstring for validity
      Subseq subseq (0, np.size());
      if (!subseq.test_fold_string (foldstring)) THROWEXPR ("Invalid fold string '" << foldstring << "' for sequence np.name\n");
      foldenv.initialise_from_fold_string (foldstring);

      // show dotplot for foldenv
      if (CTAGGING(4,DOTPLOT INDIEGRAM_DOTPLOT INDIEGRAM_FOLDENV))
	{
	  CL << "Fold envelope (initialized from foldstring '" << foldstring << "' for '" << np.name << "'):\n";
	  foldenv.render_dotplot (CL, np.seq, TRUE);
	}
    }
  // else create sampled foldenv
  else
    {
      // first initialize full (well, banded by min_loop_len) fold envelope
      // NB: initialise_full (np.size()) is equivalent to initialise_3prime_banded (np.size(), -1, 0).
      foldenv.initialise_3prime_banded (np.size(), max_subseq_len, min_loop_len);

      // show dotplot for banded foldenv (prior to sampling)
      if (CTAGGING(4,DOTPLOT INDIEGRAM_DOTPLOT INDIEGRAM_FOLDENV))
	{
	  CL << "Banded folding pre-envelope for '" << np.name << "':\n";
	  foldenv.render_dotplot (CL, np.seq, TRUE);
	}

      // sample structures according to parameter nfoldd
      sample_folds (foldenv, np, nfold);

      // show dotplot for sampled foldenv
      if (CTAGGING(4,DOTPLOT INDIEGRAM_DOTPLOT INDIEGRAM_FOLDENV))
	{
	  CL << "Sampled banded fold envelope for '" << np.name << "':\n";
	  foldenv.render_dotplot (CL, np.seq, TRUE);
	}
    }

  // dump foldenv if requested
  if (CTAGGING(2,INDIEGRAM_FOLDENV))
    foldenv.dump (CL);
  
  // log number of subseqs
  if (CTAGGING(6, INDIEGRAM_FOLDENV INDIEGRAM_ENVSIZE))
    CL << "Fold envelope for " << np.name << " has " << foldenv.subseqs() << " subseqs\n";

}


// Use scfg/paircfgdp.h and scfg/postenv.h routines to create a sampled
// fold envelope.  I could write code to accomplish this with the triplet dp that 
// I've written, but it's much simpler to just use the pair code.
// Taken from stemloc.
// NB: We assume that the calling function has initialized foldenv 
// to the full fold envelope.
void Indiegram::sample_folds (Fold_envelope& foldenv, const Named_profile& np, int nfold)
{
  // create dummy Y seq and foldenv
  const Fold_envelope dummy_foldenv_y;
  const Named_profile dummy_np_y;

  if (nfold < 0)  // if full fold envelope, don't sample
    CTAG(6,INDIEGRAM_PREFOLD << "Using banded fold envelope for sequence '" << np.name << "'\n");
  else  // if not, sample folds
    {
      CTAG(6,INDIEGRAM_PREFOLD) << "Pre-foldings sequence '" << np.name << "'\n";
      Sampled_fold_envelope sampled_foldenv;
      Subseq_coords_count subseq_count;

      // stochastic or deterministic?
      if (random_fold) // stochastic
	{
	  // create Inside DP matrix
	  const Pair_inside_matrix inside (np, dummy_np_y, foldenv, dummy_foldenv_y, qs_cfg, false);
	  if (CTAGGING(1,INDIEGRAM_DP))
	    inside.show (CL);
	  // sample the envelope
	  subseq_count = sampled_foldenv.initialise_sampled (inside, nfold, min_loop_len);
	}
      else // deterministic (N-best)
	{
	  // create CYK-KYC DP matrix
	  const Pair_CYK_KYC_matrix cyk_kyc (np, dummy_np_y, foldenv, dummy_foldenv_y, qs_cfg, false);
	  if (CTAGGING(1,INDIEGRAM_DP))
	    cyk_kyc.show (CL);
	  // sample the envelope
	  subseq_count = sampled_foldenv.initialise_best (cyk_kyc, nfold, min_loop_len);
	}

      // assign the envelope
      foldenv.swap (sampled_foldenv);
      // show dotplot
      if (CTAGGING(4,DOTPLOT INDIEGRAM_DOTPLOT INDIEGRAM_FOLDENV INDIEGRAM_SAMPLED_FOLDENV))
	{
	  CL << "Sampled fold envelope for '" << np.name << "':\n";
	  Fold_envelope::render_dotplot_from_counts (CL, subseq_count, np.seq, nfold, true);
	}
    }

}


void Indiegram::build_triplet_alignments()
{

  CTAG(8,INDIEGRAM) << "\nBuilding triplet structural alignments.\n";

  // get Tree_alignment_database
  Tree_alignment_database align_db (seq_db);
  align_db.initialise_from_Stockholm_database (stock_db, true);

  // loop over alignments, performing CYK on each
  for (int n_align = 0; n_align < align_db.size(); n_align++) {

    // get Tree_alignment
    const Tree_alignment& tree_align = *align_db.tree_align[n_align];
    const int len = tree_align.align.columns();

    // get corresponding Stockholm alignment
    Stockholm& stock = *stock_db.align_index[n_align];

    // get alignment ID
    sstring align_id = stock.get_name();
    if (!align_id.size())
      align_id << "Alignment" << n_align+1;

    // log message
    CTAG (7,INDIEGRAM) << "\nProcessing alignment '" << align_id
		       << "' (" << n_align+1 << " of " << align_db.size()
		       << "): " << len << " columns\n";
    
    // sanity checks: must be 3 leaves in the tree
    // and can only 3 seqs in alignment
    if (tree_align.tree.leaves() != 3) THROWEXPR ("There must be 3 leaves in the tree, but input alignment " << align_id << " has " << tree_align.tree.leaves() << " leaves.\n");
    if (tree_align.align.rows() != 3) THROWEXPR ("Alignment " << align_id << " doesn't have 3 sequences.\n");

    // hacky way to get Named_profile's:
    //  (the problem: Triplet_DP code wants Named_profile's, but 
    //  these aren't stored in a Tree_alignment, so we need to grab
    //  them from the corresponding Stockholm object)
    // Require that sequences have unique names to make sure that
    // the below mapping between rows in tree_align and stock doesn't break
    // (no guarantee that indexing will be the same).
    // Note that we could alternately do this by taking the Named_profile's
    // from the FASTA_sequence_database (although in this case we would have to 
    // require that all sequences in all alignments have unique names).
    tree_align.align.assert_names_unique();

    Named_profile np_x, np_y, np_z;
    double depth_x, depth_y, depth_z;
    // X (row 0)
    sstring seqname = tree_align.align.row_name[0];
    np_x = *stock.np[stock.row_index[seqname]];
    int node = tree_align.row2node[0];
    depth_x = tree_align.tree.branch_depth (node);
    // Y (row 1)
    seqname = tree_align.align.row_name[1];
    np_y = *stock.np[stock.row_index[seqname]];
    node = tree_align.row2node[1];
    depth_y = tree_align.tree.branch_depth (node);
    // Z (row 2)
    seqname = tree_align.align.row_name[2];
    np_z = *stock.np[stock.row_index[seqname]];
    node = tree_align.row2node[2];
    depth_z = tree_align.tree.branch_depth (node);

    // initialize the SCFG
    Triplet_SCFG triplet_scfg (init_triplet_scfg (depth_x, depth_y, depth_z));

    // store state_type information
    (*this).state_type.assign (triplet_scfg.state_type.begin(), triplet_scfg.state_type.end());

    // get fold envelopes
    Fold_envelope foldenv_x, foldenv_y, foldenv_z;
    // if SS specified in input Stockholm file, pull it out to initialize the Fold_envelope
    //  (holds empty string if not specified)
    sstring xfold = stock.get_gr_annot (np_x.name, Stockholm_secondary_structure_tag);
    sstring yfold = stock.get_gr_annot (np_y.name, Stockholm_secondary_structure_tag);
    sstring zfold = stock.get_gr_annot (np_z.name, Stockholm_secondary_structure_tag);
    // if fold strings are present, ensure that they're the proper length
    if (xfold.length() && (xfold.length() != np_x.seq.length()))
      THROWEXPR ("Foldstring for sequence '" << np_x.name << "' isn't flush with the sequence.");
    if (yfold.length() && (yfold.length() != np_y.seq.length()))
      THROWEXPR ("Foldstring for sequence '" << np_y.name << "' isn't flush with the sequence.");
    if (zfold.length() && (zfold.length() != np_z.seq.length()))
      THROWEXPR ("Foldstring for sequence '" << np_z.name << "' isn't flush with the sequence.");
    init_foldenv (foldenv_x, np_x, nfold, xfold);
    init_foldenv (foldenv_y, np_y, nfold, yfold);
    init_foldenv (foldenv_z, np_z, nfold, zfold);

    // log message
    CTAG(6,INDIEGRAM) << "Performing structural alignment on '" << np_x.name << "', '" << np_y.name << "', '" << np_z.name << "'.\n";	  

    // do structural alignment via CYK
    const Triplet_CYK_matrix cyk (triplet_scfg, np_x, np_y, np_z, foldenv_x, foldenv_y, foldenv_z, true);
    if (CTAGGING(1,INDIEGRAM_DP)) // dump DP matrix at lowest log level
      {
	CL << "(CYK):\n";
	CL << cyk.scfg_dump();
      }

    // test for -InfinityScore
    if (cyk.final_sc == -InfinityScore)
      {
	CLOGERR << "Warning: structural alignment score for '" << np_x.name << "', '" << np_y.name << "', '" << np_z.name << "' is -infinity; skipping.\n";
	CLOGERR << " (this usually means there is no valid alignment of the sequences, which often means that the precomputed structures are too restrictive;\n";
	CLOGERR << "  try setting -nf higher?)\n";
      }

    // get alignment and store it, along with the score,
    // as well as associated parse_tree
    //  NB: the score must be stored manually because 
    //  Triplet_SCFG_parse_tree::alignment()
    //  doesn't know about the scfg and so can't calculate the score.
    Triplet_SCFG_alignment cyk_alignment = cyk.alignment();
    cyk_alignment.score = cyk.final_sc;
    Triplet_alignment_with_tracebacks tracebacks (cyk_alignment);
    tracebacks.push_back (cyk.parse_tree());

    triplet_tracebacks_db.push_back (tracebacks);

    // perform posterior probability calculations for fold dotplots
    if (postfold_prefix.size())
      {
	const Triplet_inside_outside_matrix in_out (triplet_scfg, np_x, np_y, np_z, foldenv_x, foldenv_y, foldenv_z, true);
	const Triplet_fold_dotplot xfold_dotplot (in_out, 0);
	xfold_dotplot.write_dotplot (postfold_prefix, np_x.name);
	const Triplet_fold_dotplot yfold_dotplot (in_out, 1);
	yfold_dotplot.write_dotplot (postfold_prefix, np_y.name);
	const Triplet_fold_dotplot zfold_dotplot (in_out, 2);
	zfold_dotplot.write_dotplot (postfold_prefix, np_z.name);

	if (CTAGGING(1,INDIEGRAM_DP)) // dump DP matrix at lowest log level
	  {
	    CL << "(Inside):\n";
	    CL << in_out.inside.scfg_dump();
	    CL << "(Outside):\n";
	    CL << in_out.outside.scfg_dump();
	  }
      }


  }

  CTAG(8,INDIEGRAM) << "Completed triplet structural alignments.\n";
}


void Indiegram::show_triplet_tracebacks (ostream& o)
{

  for_const_contents (list<Triplet_alignment_with_tracebacks>, triplet_tracebacks_db, tb)
    {
      (*tb).alignment.show (o);
      for_const_contents (vector<Triplet_SCFG_parse_tree>, (*tb), t)
	(*t).show (o, &(*this).state_type);
      o << "\n";
    }

}

int Indiegram::run()
{
  try
    {
      init_opts();
      parse_opts();
      input_data();
      init_prefold_scfg();
      build_triplet_alignments();
      show_triplet_tracebacks (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit (1);
    }
  return 0;
}

