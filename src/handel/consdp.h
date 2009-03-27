#ifndef CONSTRAINEDDP_INCLUDED
#define CONSTRAINEDDP_INCLUDED

#include "handel/multiwaydp.h"

// DP matrix constrained to fit a given alignment
struct Transducer_constrained_DP_matrix : Transducer_DP_base
{
  // typedefs
  typedef Alignment_path::Sequence_coords Sequence_coords;

  // Constrained alignment path description
  const Alignment_path* path;
  Alignment_path collapsed_path;  // contains only rows for observed sequences
  vector<Sequence_coords> seq_coords;

  // DP matrix
  array2d<Score> cell_sc;  // indexed by (state,column)

  // alloc method
  void alloc();

  // accessors
  int columns() const { return seq_coords.size(); }
};

// Forward DP matrix constrained to fit a given alignment
// Assumes no null states in transducer.
struct Transducer_constrained_forward_matrix : Transducer_constrained_DP_matrix, Transducer_forward_matrix_interface
{
  // fill method
  void fill();

  // sample state path
  vector<int> sample_traceback();
};

#endif /* CONSTRAINEDDP_INCLUDED */
