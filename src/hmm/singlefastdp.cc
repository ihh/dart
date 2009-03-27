#include "hmm/singlefastdp.h"
#include "util/dexception.h"

Single_fast_matrix_factory Single_fast_matrix_factory::instance = Single_fast_matrix_factory();

Single_Viterbi_interface* Single_fast_matrix_factory::new_Viterbi (const Single_HMM_scores& hmm, const Named_profile& np)
{
  return new Single_fast_Viterbi_matrix (hmm, np);
}

Single_forward_interface* Single_fast_matrix_factory::new_forward (const Single_HMM_scores& hmm, const Named_profile& np)
{
  return new Single_fast_forward_matrix (hmm, np);
}

Single_forward_backward_interface* Single_fast_matrix_factory::new_forward_backward (const Single_HMM_scores& hmm, const Named_profile& np)
{
  return new Single_fast_forward_backward_matrix (hmm, np);
}
