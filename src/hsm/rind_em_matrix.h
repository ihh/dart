#ifndef RIND_EM_MATRIX
#define RIND_EM_MATRIX

#include "hsm/em_matrix.h"
#include "ecfg/ecfgkeywords.h"

struct RIND_EM_matrix : EM_matrix
{
  // constructor
  RIND_EM_matrix (int C, int A, int max_fork = 1, const Tree_alignment_database* align_db = 0, double timepoint_res = DEFAULT_TIMEPOINT_RES);

  // descriptor
  const char* update_policy() const { return EG_POLICY_RIND; }

  // EM algorithm(s)
  void quick_M (const Update_statistics& stats, bool intra = TRUE, bool inter = TRUE);  // the M-step previously quick_M_rind()
};

#endif /* RIND_EM_MATRIX */
