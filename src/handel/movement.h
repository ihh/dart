#ifndef HANDEL_MOVEMENT_INCLUDED
#define HANDEL_MOVEMENT_INCLUDED

#include "handel/ehmm.h"
#include "handel/transducer.h"
#include "handel/transducer_sexpr.h"
#include "handel/multiwaydp.h"
#include "handel/hmmoc_adapter_opts.h"

// This is a representation/dump of a Handel MCMC move
// that will (probably eventually) be one of the main Handel objects,
// along with the tree-alignment and the program environment.
//
// Essentially, handalign is a wrapper for phylocomposer
// (whose input language is remarkably similar to an old idea called Telegraph...)
//
struct Handel_movement
{
  // DATA INITIALIZED BY CALLER
  // display options
  bool quiet;
  bool composite;
  bool show_constrained_composite_paths;
  bool acyclic;
  bool normalize_peeled_profiles;
  sstring dotfile;
  sstring stockfile;

  // algorithm options
  bool simulate;
  bool viterbi;
  int nforward;  // -1 means do not fill forward; 0 means fill but don't do traceback
  bool optacc;
  bool want_expected_counts;
  bool propose_redelings_suchard_move, evaluate_redelings_suchard_inverse_move;
  bool use_centroid;
  int centroid_band_width;

  // HMMoC adapter options
  HMMoC_adapter_options hmmoc_opts;
  bool use_hmmoc_for_forward, use_hmmoc_for_viterbi;
  bool old_path_is_outside_hmmoc_band;  // can be used to reject alignment-changing MCMC moves that start outside the band

  // Transducer_SExpr_file
  // can be initialized via read_composition method, or by creating on stack (calling constructor) & assigning
  Transducer_SExpr_file composition;

  // DATA INITIALIZED BY THIS OBJECT (specifically, by dump_composition method)
  // peeled Transducer_SExpr_file object
  Transducer_SExpr_file peeled_composition;

  // log-likelihood of "old" path
  Loge old_loglike;

  // log-likelihood of peeling
  Loge peeling_loglike;

  // alignment flags
  bool peeling_required, fully_constrained, fill_forward, fill_backward, redelings_suchard, do_alignment, use_composite_as_HMM, ehmm_funcs_defined;

  // EHMM transducer objects, initialized by this object
  EHMM_transducer_funcs ehmm_funcs;
  EHMM_transducer_scores ehmm_scores;
  vector<int> states_to_eliminate;
  Eliminated_EHMM_transducer_scores elim_ehmm_scores;

  // pointers to Score_profile's in Transducer_SExpr_file objects
  vector<Score_profile*> seq_vec;

  // dynamic programming matrices & associated data
  Transducer_Viterbi_matrix vit;  // Viterbi matrix
  vector<int> vit_trace;  // Viterbi traceback
  Transducer_forward_matrix fwd;  // Forward matrix
  vector<vector<int> > fwd_trace;  // Forward traceback
  Transducer_backward_matrix back;  // Backward matrix
  PCounts pcounts;  // Baum-Welch/peeling counts
  Alignment_path optacc_path;  // optimal accuracy traceback
  Loge fwd_ll, vit_ll;

  // Redelings-Suchard kernel info
  Transducer_SExpr_file::NodePathMap redsuch_path;  // indexed by unpeeled node index
  vector<int> redsuch_trace;  // path through ehmm_funcs
  Loge redsuch_ll;
  vector<Loge> hastings_term;

  // constructor
  // default behavior: show composite transducer, but do not do any DP
  Handel_movement()
  : quiet (false),
    composite (true),
    show_constrained_composite_paths (false),
    acyclic (false),
    normalize_peeled_profiles (true),
    dotfile(),
    stockfile(),
    simulate(false),
    viterbi(false),
    nforward(-1),  // -1 means do not fill forward
    optacc(false),
    want_expected_counts(false),
    propose_redelings_suchard_move(false),
    evaluate_redelings_suchard_inverse_move(false),
    use_centroid(false),
    old_path_is_outside_hmmoc_band(false),
    old_loglike (0),
    peeling_loglike (0),
    elim_ehmm_scores (ehmm_scores, states_to_eliminate),
    fwd_ll (-InfinityLoge),
    vit_ll (-InfinityLoge)
  { }

  // methods
  void read_composition (const vector<sstring>& filenames);
  void dump_composition (ostream& out);  // does bulk of high-level phylocomposer work

  // output helpers
  void dump_score (ostream& out, const char* tag, Score sc);  // adds peeling_loglike to sc argument
  void dump_raw_loglike (ostream& out, const char* tag, Loge ll);  // does NOT add peeling_loglike to sc argument
  void dump_unpeeled_node_path_map (ostream& out, const char* tag, const Transducer_SExpr_file::NodePathMap& path, Loge path_ll);
  void dump_peeled_composite_trace (ostream& out, const char* tag, const vector<int>& path, Loge path_ll);
  void dump_and_store_forward_path (ostream& out, const vector<int>& acyclic_path);
  void dump_and_store_viterbi_path (ostream& out, const vector<int>& composite_path, Loge path_ll);

  // method to return the DP score plus the peeled score, as a log-likelihood
  Loge unpeeled_ll (Score peeled_sc);

};


#endif /* HANDEL_MOVEMENT_INCLUDED */

