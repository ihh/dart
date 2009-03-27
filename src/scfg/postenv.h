#ifndef POSTENV_INCLUDED
#define POSTENV_INCLUDED

#include "scfg/paircfgdp.h"

// structure to hold a subset of states from a Pair CFG
struct Pair_CFG_state_set
{
  // data
  const Pair_CFG_state_typing& typing;
  set<int> states;
  // constructor
  Pair_CFG_state_set (const Pair_CFG_state_typing& typing);
  // build methods
  void initialise_full();
};

// sampled single-sequence fold envelope
struct Sampled_fold_envelope : public Fold_envelope
{
  // constructor
  Sampled_fold_envelope();
  // builders
  void initialise_CYK (const Pair_CYK_matrix& cyk, int min_loop_len = 0);
  Subseq_coords_count initialise_sampled (const Pair_inside_matrix& imx, int n_samples, int min_loop_len = 0);
  Subseq_coords_count initialise_best (const Pair_CYK_KYC_matrix& cyk_kyc, int n_best, int min_loop_len = 0);
};

// sampled pairwise alignment envelope
struct Sampled_pair_envelope : Pair_envelope
{
  // Named_profile pointers used for logging only
  const Named_profile* npx;
  const Named_profile* npy;
  // Fold_envelope pointers for memory usage calculations
  const Fold_envelope* xenv;
  const Fold_envelope* yenv;
  // memory usage
  int paths_added;
  unsigned long cells_needed, cells_available;
  // constructor
  Sampled_pair_envelope();
  Sampled_pair_envelope (const Named_profile& npx, const Named_profile& npy,
			 const Fold_envelope& xenv, const Fold_envelope& yenv);
  // builders
  void initialise_sampled (const Pair_forward_DP_matrix& fwd, int n_samples, const set<int>& states);
  void initialise_best (const Pair_CYK_KYC_matrix& cyk_kyc, const Pair_HMM_scores& hmm, int n_best, const set<int>& states);
  bool add_state_path (const vector<int>& path, const Pair_HMM_scores& hmm, const set<int>& states);  // returns true if success, false if out-of-memory
};

#endif /* POSTENV_INCLUDED */
