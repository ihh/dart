#ifndef SIMPLE_EHMM_INCLUDED
#define SIMPLE_EHMM_INCLUDED

#include "hsm/em_matrix.h"
#include "hmm/singlehmm.h"

// Simple evolutionary HMM (no proper modeling of gaps)
struct Simple_EHMM_scores : Single_HMM_scores
{
  // data
  vector<EM_matrix> em_matrix;  // one matrix per state; keep it simple
  // constructor
  Simple_EHMM_scores (int states);
  // I/O methods
  void read (istream& in);
  void write (ostream& out);
};

// Base class of DP matrices for simple evolutionary HMM
struct Simple_EHMM_DP_matrix_base
{
  // typedefs
  typedef EM_matrix::Column_matrix Column_matrix;
  // data
  const Simple_EHMM_scores& ehmm;
  const Tree_alignment& tree_align;
  vector<Column_matrix> colmat;
  // constructor
  Simple_EHMM_DP_matrix_base (const Simple_EHMM_scores& ehmm, const Tree_alignment& tree_align);
};

// Forward matrix
struct Simple_EHMM_forward_matrix : Simple_EHMM_DP_matrix_base
{
  Simple_EHMM_forward_matrix (const Simple_EHMM_scores& ehmm, const Tree_alignment& tree_align);
};

// Forward-backward matrix
struct Simple_EHMM_forward_backward_matrix : Simple_EHMM_DP_matrix_base
{
  Simple_EHMM_forward_backward_matrix (const Simple_EHMM_scores& ehmm, const Tree_alignment& tree_align);
};


#endif /* SIMPLE_EHMM_INCLUDED */
