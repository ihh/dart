#ifndef HANDEL_INIT_ALIGN_INCLUDED
#define HANDEL_INIT_ALIGN_INCLUDED

// Sketch of code for quick initial estimation of alignment & tree in Handel
// Input: Sequence_database object
// Output: Stockade object, containing a Stockholm object with a Newick-format PHYLIP_tree "#=GF NH" annotation

// Algorithm:
// Create Pair_HMM_scores object from Pair_transducer_scores returned by Pair_transducer_factory::prior_pair_trans_sc()
//  (need to give Pair_transducer_scores a method that returns the equivalent Pair_HMM_scores)
//  (branch length set on command line; default is 0.5)
// Construct a Dotplot_map for the FASTA_sequence_database and Pair_HMM_scores, using a derived Dotplot class
//  (which populates dotplot by summing EmitXY-state postprobs from a Pair_forward_backward_DP_matrix)
//  (similar to PairCFG_alignment_dotplot for SCFGs)
// At this point we have what will be the final Stockade; initialise a Tree_alignment from the Stockholm object
// Create a Subst_dist_func_factory from the EM_matrix_base
// Call Tree_alignment::estimate_tree_by_nj
// Call Tree_alignment::copy_tree_to_Stockholm

#include "handel/alitrans.h"
#include "amap/dotplot.h"
#include "amap/amap_adapter.h"
#include "hmm/pairhmm.h"
#include "tree/nj.h"
#include "tree/subdistmat.h"

#endif /* HANDEL_INIT_ALIGN_INCLUDED */
