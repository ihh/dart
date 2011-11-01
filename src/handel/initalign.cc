#include "handel/initalign.h"

// Alignment algorithm
Stockade Stockade_initializer::align (FASTA_sequence_database& seq_db, const Alphabet& alphabet, Pair_transducer_factory& trans_fac, double branch_len)
{
  CTAG(6,HANDEL) << "Building initial alignment of " << seq_db.size() << " sequences by sequence annealing\n";
  if (seq_db.size() < 2)
    THROWEXPR("You need at least two sequences to build an alignment");
  // convert sequences using alphabet
  seq_db.seqs2scores (alphabet);
  // Create Pair_HMM_scores object from Pair_transducer_scores returned by Pair_transducer_factory::branch_pair_trans_sc()
  Pair_transducer_scores pair_trans_sc = trans_fac.branch_pair_trans_sc (branch_len);
  Pair_HMM_scores pair_hmm = pair_trans_sc.pair_hmm (alphabet);
  // Construct a Dotplot_map for the FASTA_sequence_database and Pair_HMM_scores
  Dotplot_map dotplot_map;
  for (int i = 0; i < seq_db.index.size(); ++i)
    for (int j = i + 1; j < seq_db.index.size(); ++j)
      {
	Named_profile& np_i (*seq_db.index.profile[i]);
	Named_profile& np_j (*seq_db.index.profile[j]);
	CTAG(5,HANDEL) << "Aligning " << np_i.name << " and " << np_j.name << "\n";
	Pair_forward_backward_DP_matrix fb (pair_hmm, np_i.prof_sc, np_j.prof_sc);
	Post_pair_HMM fb_post (fb);
	Dotplot dotplot (np_i.seq, np_j.seq);
	(array2d<Prob>&) dotplot = (array2d<Prob>&) fb_post;  // ugh
	dotplot_map[i][j] = dotplot;
      }
  AMAP_adapter amap (seq_db, dotplot_map);
  Stockade stockade = amap.get_alignment();
  if (CTAGGING(5,HANDEL))
    {
      CL << "Initial alignment:\n";
      stockade.align.write_Stockholm (CL);
    }
  return stockade;
}

// Tree algorithm
void Stockade_initializer::build_tree (Stockholm& stock, Substitution_matrix_factory& submat_factory)
{
  CTAG(6,HANDEL) << "Building initial tree for " << stock.rows() << "-sequence alignment by neighbor-joining\n";
  if (stock.rows() < 2)
    THROWEXPR("You need at least two sequences to build a tree");
// Initialise a Tree_alignment from the Stockholm object
  Tree_alignment tree_align (stock, false);
// Create a Subst_dist_func_factory from the Substitution_matrix_factory
  Subst_dist_func_factory distfunc_fac (submat_factory);
// Call Tree_alignment::estimate_tree_by_nj
  tree_align.estimate_tree_by_nj (distfunc_fac);
// Call Tree_alignment::copy_tree_to_Stockholm
  tree_align.copy_tree_to_Stockholm (stock);
  // log
  if (CTAGGING(5,HANDEL))
    {
      CL << "Initial alignment, with tree:\n";
      stock.write_Stockholm (CL);
    }
}
