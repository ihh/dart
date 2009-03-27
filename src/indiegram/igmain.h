#ifndef IGMAIN_INCLUDED
#define IGMAIN_INCLUDED

#include "util/logfile.h" // NB: This includes util/optslist.h
#include "seq/stockholm.h"
#include "tree/tree_alignment.h"
#include "indiegram/igdefaults.h"
#include "indiegram/tripletscfg.h"
#include "indiegram/postprobs.h"
#include "scfg/foldenv.h"
#include "scfg/paircfgdp.h"
#include "stemloc/superstem.h"

/// Encapsulates a running Indiegram program.
struct Indiegram
{
  
  /// command-line options
  int nfold;
  bool random_fold;
  int max_subseq_len;
  int min_loop_len;

  bool use_test_scfg;

  sstring postfold_prefix;

  /// options list
  Opts_list opts;
  
  const Alphabet& alphabet;

  /// Basically a helper to read input Stockholm sequence files.
  Stockholm_database stock_db;

  /// Sequence data.
  /*
   * Used for internal storage of sequence data.
   * Populated by Stockholm_database::read().
   * Note that while stock_db holds each alignment separately,
   * seq_db globs all alignments together.
   */
  FASTA_sequence_database seq_db;

  /// Null model for emissions
  vector<Prob> null_emit_prob;

  /// triplet alignment databases
  /*
   * Used to store the results of structural alignment.
   */
  list<Triplet_alignment_with_tracebacks> triplet_tracebacks_db;


  // Grammars //

  /// Super_stem model.
  /*
   * Single-sequence SCFG for pre-folding to predict the fold envelope.
   * Taken straight from stemloc.
   */
  PScores         qs_pscore;
  Super_stem      qstem;
  Dirichlet_prior qs_prior;
  Pair_CFG_scores qs_cfg;

  /// Set up command-line options.
  void init_opts();

  /// Read command-line options.
  void parse_opts();

  /// Read in input data.
  void input_data();

  /// Initialize parameters for pre-folding the sequences.
  /*
   * Parameters for a single-sequence SCFG for constructing
   * the fold envelopes (see the Super_stem model).
   * A la Stemloc::init_fold_params().
   */
  void init_prefold_scfg();

  /// Initialize parameters for triplet structural alignment.
  /*
   * The Triplet_SCFG will depend on the branch lengths t, u, v.
   */
  Triplet_SCFG init_triplet_scfg (double branch_length_x, double branch_length_y, double branch_length_z);

  /// Build fold envelope.
  void init_foldenv (Fold_envelope& env, const Named_profile& np, int nfold, const sstring& foldstring = "");

  /// Sample structures for pre-folding.
  /*
   * For building the fold envelopes.
   * Currently only deterministic N-best sampling is implemented.
   */
  void sample_folds (Fold_envelope& foldenv, const Named_profile& np, int nfold);

  /// Build a structural alignment of the input sequences.
  /*
   * Uses constrained DP to a triplet SCFG.
   */
  void build_triplet_alignments();

  /// Show the alignments and parse trees.
  void show_triplet_tracebacks (ostream& o);

  /// run indiegram
  int run();

  /// constructor
  Indiegram (int argc, char** argv);

private:

  /// state types
  vector<State_type> state_type;

};

#endif /* IGMAIN_INCLUDED */
