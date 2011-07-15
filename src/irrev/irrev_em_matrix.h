// pk 3/05 from em_matrix.h support assymetric matrices
// pk 4/05 cleanup - eliminate iterate_EM(), single_EM(), set_rind_params(), get_rind_counts(), quick_M_rind(),
//   param_minima(), set_params_from_vector()
// pk 4/05 perturb matrix, retry if MatrixExpEigenPrepare fails, finish elimination of Newton-Raphson params, methods
// pk 4/05 combine common em_matrix, irrev_em_matrix classes

#ifndef ASYM_EM_MATRIX
#define ASYM_EM_MATRIX

#include "hsm/em_matrix_base.h"	
#include "util/vector_output.h"
#include "ecfg/ecfgkeywords.h"

typedef EM_matrix_base::Update_statistics Update_statistics;
typedef EM_matrix_base::Column_matrix Column_matrix;

// number retries if MatrixExpEigenPrepare fails
#define MATRIXEXP_RETRIES 20
// amount to perturb matrix if MatrixExpEigenPrepare fails
#define MATRIXEXP_PERTURB 2e-14

// Irreversible EM algorithm
struct Irrev_EM_matrix : EM_matrix_base
{
  // constructor
  Irrev_EM_matrix (int C,
		   int A,
		   int max_fork = 1,
		   const Tree_alignment_database* align_db = 0,
		   double timepoint_res = DEFAULT_TIMEPOINT_RES);
  // Null constructor ,added by OW 7-15-2011
  Irrev_EM_matrix (void); 


  // descriptor
  const char* update_policy() const { return EG_POLICY_IRREV; }

  // irreversible diagonalize() method
  void diagonalize(); 	// matrix diagonalization

  // irreversible transform_symmetrised... method
  void transform_symmetrised_waits_transitions (Update_statistics& stats, bool symmetrise) const;

  // EM algorithm(s)
  Update_statistics single_quick_EM (bool intra = TRUE, bool inter = TRUE, bool infer_class_labels = FALSE);  // does partial M-step only
  void quick_M (const Update_statistics& stats, bool intra = TRUE, bool inter = TRUE);  // the M-step - irrev version
};

#endif /* ASYM_EM_MATRIX */
