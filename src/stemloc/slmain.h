#ifndef SLMAIN_INCLUDED
#define SLMAIN_INCLUDED

#include <math.h>

#include "stemloc/superstem.h"
#include "stemloc/quickalign.h"
#include "stemloc/superpair.h"
#include "util/logfile.h"
#include "util/vector_output.h"
#include "seq/stockholm.h"
#include "telegraph/tgvar.h"
#include "scfg/paircfgem.h"
#include "hmm/pairhmmem.h"
#include "scfg/postenv.h"
#include "stemloc/sldefaults.h"

// This class encapsulates a running stemloc program
struct Stemloc
{
  // typedefs
  typedef Alignment_path::Row_pair      Row_pair;
  typedef Alignment_path::Decomposition Decomposition;
  typedef map<sstring, set<sstring> >   Pairwise_schedule;

  // helper struct for conservative multiple alignment (yuck)
  struct Multi_match
  {
    // data
    int seed_row;  // the row of the parent path in the seed alignment
    int seed_lpad;  // number of padding columns at 5' end of seed row
    Score score;  // score
    Named_profile subseq;  // local subseq
    Pairwise_path path;  // local path
    sstring fold;  // local fold string
    // constructors
    Multi_match();
    Multi_match (const Pair_CFG_alignment& align, int seed_row);
  };
  typedef map<sstring,Multi_match> Best_match_table;

  // better helper for liberal multiple alignment
  struct Extension
  {
    // data
    int seed_row;  // the row of the parent path in the seed alignment
    RNA_pair_path align;
    // constructors
    Extension();
    Extension (const RNA_pair_path& align, int seed_row);
  };
  typedef map<sstring,Extension> Extension_table_cell;
  typedef map<sstring,Extension_table_cell> Extension_table;

  // option list
  Opts_list opts;

  // command-line options
  sstring param_set;
  bool train_only;
  double em_min_inc;
  int em_max_iter;
  int em_forgive;

  int pairwise_nfold;
  bool random_fold;
  bool cache_fold;
  int max_subseq_len;
  int min_loop_len;

  sstring trainfold;
  sstring trainalign;
  sstring traincmp;
  sstring loadgramset;
  sstring savegramset;

  int pairwise_nalign;
  int max_megabytes;
  bool random_align;
  int align_band;
  bool allow_extra_basepairs;

  bool global;
  double min_bitscore;
  int max_hits;

  bool multi;
  bool show_intermediates;
  bool progressive;
  bool liberal;
  double min_extend_bitscore;
  int multi_nfold;
  int multi_nalign;
  int min_overlap_len;

  /// AMA options
  bool use_ama;  /// use AMA to create a multiple alignment
  double gap_factor;
  double edge_weight_threshold;
  int num_consistency_reps;
  bool output_for_gui;

  sstring palign_prefix; /// prefix for files holding posterior alignment probabilities
  sstring pfold_prefix;  /// prefix for files holding posterior fold probabilities

  // Commented-out model comparison stuff:
  // sstring loadcds;
  // sstring traincds;
  // sstring savecds;
  // bool do_pairwise_model_comparison;
  // sstring schedule_filename;

  // constants
  const Alphabet& alphabet;

  // runtime vars
  // input data
  Stockholm_database align_db;
  FASTA_sequence_database seq_db;

  // null model
  vector<Prob> null_emit_prob;

  // fold envelope databases
  Fold_envelope_database foldenv_db;  // standard cache
  Fold_envelope_database multi_foldenv_db;  // multiple alignment cache

  // pairwise alignment database
  list<Pair_CFG_alignment> pair_align_db;

  // database of unaligned sequence names, and pairwise alignments of unaligned sequences
  list<const Pair_CFG_alignment*> seed_db;
  set<sstring> unaligned_seqname;

  // Probabilistic grammars.

  // HACK WARNING:
  // Currently the handling of different parameterizations is a little sloppy;
  // things like GC-content & sequence identity should be handled more systematically.
  // Log-odds ratio scores are used somewhat haphazardly at present, without much
  // consistency as to whether a single "null model", or multiple null models, are used.
  // Also, the parameterization is chosen heuristically at the primary pairwise pre-alignment stage,
  // without a clear underlying interpretation in terms of Bayesian model comparison.
  // Finally, the parameters are split over multiple PScope's, and the alternative
  // parameter sets are named using an awkward automatic naming scheme.
  // The former needlessly prevents parameters from being stored in a single large Telegraph file,
  // while the latter is just messy.

  // The Gramset object. This contains:
  //  -- the single-SCFG for pre-folding (i.e. predicting the fold envelope);
  //  -- a series of {pair-HMM,pair-SCFG} model-pairs:
  //     ... the pair-HMM is for pre-aligning (predicting the alignment envelope) and choosing which model-pair to use;
  //     ... the pair-SCFG is for simultaneous aligning-and-folding (constrained by the fold & alignment envelopes).
  Gramset gramset;

  // versions of the various models with Scores evaluated
  Pair_CFG_scores single_scfg;
  map<sstring,Pair_HMM_scores> pair_hmm;
  map<sstring,Pair_CFG_scores> pair_cfg;
  map<sstring,set<int> > pair_hmm_states;

  // main methods
  void init_opts();
  void parse_opts();
  void input_data();
  void load_gramset();
  void init_fold_params();
  void init_align_params();
  void init_cmp_params();
  void save_gramset();
  void build_pairwise_alignments();
  void build_multiple_alignments();
  Stockade build_conservative_multiple_alignment (const Pair_CFG_alignment& seed);
  Stockholm build_liberal_multiple_alignment (const Pair_CFG_alignment& seed);
  Stockholm build_motif_alignment (const Pair_CFG_alignment& seed);

  /// use AMA
  void build_multiple_alignments_ama();

  // helper methods
  void init_pair_fold_env (Fold_envelope& env, const Named_profile& np, int nfold);
  void init_multi_fold_env (Fold_envelope& env, const Named_profile& np, const Local_fold_string& local_fold);
  bool attempt_init_from_fold_string (Fold_envelope& env, const Named_profile& np);
  void sample_folds (Fold_envelope& env, const Named_profile& np, int nfold);  // on entry, env = banding envelope
  void init_align_env (Pair_envelope& pair_env, const Named_profile& npx, const Named_profile& npy, const Fold_envelope& xenv, const Fold_envelope& yenv, sstring& best_param_set, int nalign);
  void check_name_unique (const sstring& seqname);  // flushes alignment & fold caches if not
  void remove_seeds (const Alignment& multi_align);  // removes sequences in alignment from seed_db and unaligned_seqname

  // run method
  int run();

  // constructor
  Stemloc (int argc, char** argv);
};

#endif /* SLMAIN_INCLUDED */
