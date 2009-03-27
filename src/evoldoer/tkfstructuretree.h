#ifndef TKF_STRUCTURE_TREE_INCLUDED
#define TKF_STRUCTURE_TREE_INCLUDED

#include "scfg/paircfg.h"
#include "handel/tkfparams.h"
#include "hsm/em_matrix.h"

struct TKFST_params
{
  // data
  EM_matrix  loop_mat;
  EM_matrix  stem_mat;
  TKF_params loop;
  TKF_params stem;
  double     stem_prob;
  // constructor
  TKFST_params();
};

struct TKFST_default_params : TKFST_params
{
  // constructor
  TKFST_default_params();
};

struct TKFST_single_CFG : Pair_CFG_scores
{
  // data
  TKFST_params& params;
  // states
  enum { L = 0, Lprime = 1, Sprime = 2, B = 3, Bprime = 4, C = 5, aL = 6, aS = 7,
	 TotalStates = 8 };
  // constructor
  TKFST_single_CFG (TKFST_params& params, bool odds_ratio);
};

struct TKFST_pair_CFG : Pair_CFG_scores
{
  // data
  TKFST_params& params;
  double time;
  // states
  enum { L1 = 0, S1 = 1, L2 = 2, L3 = 3, L4 = 4,
	 L1prime = 5, S1prime = 6, L3prime = 7, S3prime = 8, L4prime = 9, S4prime = 10,
	 B41 = 11, B32 = 12, B11 = 13, B33 = 14, B44 = 15, B11prime = 16, B33prime = 17, B44prime = 18, C11 = 19, C33 = 20, C44 = 21,
	 adL1 = 22, dL1 = 23, aL2 = 24, aL3 = 25, dL4 = 26, adS1 = 27, dS1 = 28, aS2 = 29, aS3 = 30, dS4 = 31,
	 TotalStates = 32 };
  // constructor
  TKFST_pair_CFG (TKFST_params& params, double time, bool conditional, bool odds_ratio);
};

#endif /* TKF_STRUCTURE_TREE_INCLUDED */
