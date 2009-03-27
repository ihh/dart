#ifndef PFUNC_EM_MATRIX_INCLUDED
#define PFUNC_EM_MATRIX_INCLUDED

#include "irrev/irrev_em_matrix.h"

struct EM_matrix_funcs : EM_matrix_params_template<PFunc>
{
  // initialise method
  void init_matrix_funcs (int new_C, int new_A)
  { init_matrix_template (new_C, new_A, PFunc()); }
};

struct PFunc_EM_matrix : Irrev_EM_matrix
{
  // data
  EM_matrix_funcs funcs;

  // constructor
  PFunc_EM_matrix (int C, int A, double timepoint_res = DEFAULT_TIMEPOINT_RES)
    : Irrev_EM_matrix (C, A, 1, 0, timepoint_res)
  {
    funcs.init_matrix_funcs (C, A);
  }

  // descriptor
  const char* update_policy() const { return EG_PARAMETRIC; }
};

#endif /* PFUNC_EM_MATRIX_INCLUDED */
