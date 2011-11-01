#ifndef HANDEL_INIT_ALIGN_INCLUDED
#define HANDEL_INIT_ALIGN_INCLUDED

// Code for quick initial estimation of alignment & tree in Handel

#include "handel/alitrans.h"
#include "amap/dotplot.h"
#include "amap/amap_adapter.h"
#include "hmm/pairhmm.h"
#include "tree/nj.h"
#include "tree/subdistmat.h"

struct Stockade_initializer {
  // Alignment algorithm
  // Input: Sequence_database object
  // Output: Stockade object
  Stockade align (const Sequence_database& seq_db, Pair_transducer_factory& trans_fac, double branch_len);

  // Tree algorithm
  // Input: Stockholm object
  // Output: Same Stockholm object with a Newick-format PHYLIP_tree "#=GF NH" annotation
  void build_tree (Stockholm& stock, Substitution_matrix_factory& submat_factory);

  // Combined alignment+tree algorithm
  // Input: Sequence_database object
  // Output: Stockade object whose Stockholm has a Newick-format PHYLIP_tree "#=GF NH" annotation
  Stockade align_and_build_tree (const Sequence_database& seq_db, Transducer_alignment_with_subst_model& trans_fac, double branch_len);
};

#endif /* HANDEL_INIT_ALIGN_INCLUDED */
