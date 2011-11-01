#include "handel/initalign.h"

// Alignment algorithm
// Input: Sequence_database object
// Output: Stockade object
Stockade Stockade_initializer::align (const Sequence_database& seq_db, Pair_transducer_factory& trans_fac, double branch_len)
{
  Stockade stockade;
// Create Pair_HMM_scores object from Pair_transducer_scores returned by Pair_transducer_factory::prior_pair_trans_sc()
//  (need to give Pair_transducer_scores a method that returns the equivalent Pair_HMM_scores)
//  (branch length set on command line; default is 0.5)
// Construct a Dotplot_map for the FASTA_sequence_database and Pair_HMM_scores, using a derived Dotplot class
//  (which populates dotplot by summing EmitXY-state postprobs from a Pair_forward_backward_DP_matrix)
//  (similar to PairCFG_alignment_dotplot for SCFGs)
// At this point we have the final Stockade
  THROWEXPR ("Not implemented");
  return stockade;
}

// Tree algorithm
// Input: Stockholm object
// Output: Same Stockholm object with a Newick-format PHYLIP_tree "#=GF NH" annotation
void Stockade_initializer::build_tree (Stockholm& stock, Substitution_matrix_factory& submat_factory)
{
// Initialise a Tree_alignment from the Stockholm object
// Create a Subst_dist_func_factory from the Substitution_matrix_factory
// Call Tree_alignment::estimate_tree_by_nj
// Call Tree_alignment::copy_tree_to_Stockholm
  THROWEXPR ("Not implemented");
}

// Combined alignment+tree algorithm
// Input: Sequence_database object
// Output: Stockade object whose Stockholm has a Newick-format PHYLIP_tree "#=GF NH" annotation
Stockade Stockade_initializer::align_and_build_tree (const Sequence_database& seq_db, Transducer_alignment_with_subst_model& trans_fac, double branch_len)
{
  Stockade stockade = align (seq_db, trans_fac, branch_len);
  build_tree (stockade.align, trans_fac.submat_factory());
  return stockade;
}
