#ifndef TKFEM_INCLUDED
#define TKFEM_INCLUDED

#include "tkf/tkfcoeff.h"
#include "hsm/em_matrix.h"

struct TKF_EM_params
{
  // data
  EM_matrix  subst_model;
  TKF_params indel_model;

  // static data
  static const double default_insert_rate;
  static const double default_delete_rate;

  // typedefs & structs
  typedef EM_matrix::Update_statistics Subst_counts;
  struct Update_statistics
  {
    // data
    Subst_counts subst_counts;
    TKF_counts   indel_counts;
    // constructor
    Update_statistics (const TKF_EM_params& params);
    // methods
    void clear();  // resets all counts
  };

  // methods
  // constructors
  TKF_EM_params (int A, int C = 1);  // A = alphabet size, C = number of hidden alphabet classes
  TKF_EM_params (const Alphabet& alphabet, int C = 1);  // A = alphabet size, C = number of hidden alphabet classes
};

#endif /* TKFEM_INCLUDED */
