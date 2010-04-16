#ifndef ECFG_PLACER_INCLUDED
#define ECFG_PLACER_INCLUDED

#include "ecfg/ecfg_branch_length_em.h"

// class to compute & store ECFG_branch_state_counts for every possible attachment node of every unattached alignment row
struct ECFG_placer
{
  // typedefs
  typedef EM_matrix_base::Column_matrix Column_matrix;

  typedef map<int,ECFG_branch_state_counts> ECFG_row_attachment_counts;  // indexed by attachment node index
  typedef map<int,ECFG_row_attachment_counts> ECFG_unaligned_attachment_counts;  // indexed by unattached alignment row index

  typedef pair<int,double> Attachment_branch;
  typedef map<int,Attachment_branch> Attachment_map;  // indexed by alignment row

  // data
  ECFG_unaligned_attachment_counts attach_counts;  // attach_counts[alignment_row][tree_node][chain_idx](src_state,dest_state)
  vector<int> unattached_rows;

  ECFG_scores& ecfg;
  Stockholm& stock;
  Tree_alignment& tree_align;  // this must be constructed from stock
  double prior_param;  // P(branch_length=t) = exp(-prior_param*t)

  // constructor
  ECFG_placer (ECFG_scores& ecfg, Stockholm& stock, Tree_alignment& tree_align, double prior_param = 0.);

  // method to attach unattached alignment rows to the tree
  void attach (double tres = .01, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);

  // helper method to return the best (attachment-node,branch-length) pair for each unattached row
  // called by attach()
  Attachment_map best_attachments (double tres = .01, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);

  // helper method to populate attach_counts
  // called by best_attachments()
  // At the moment, the counts are conditioned on the CYK parse tree; probably would be better to sum over all parse trees
  void populate_counts();
};

// extension to ECFG_EM_tree_alignment_database incorporating branch length EM for entire ECFG
struct ECFG_attachable_tree_alignment_database : ECFG_EM_tree_alignment_database
{
  // constructors
  ECFG_attachable_tree_alignment_database (Sequence_database& seq_db)
    : ECFG_EM_tree_alignment_database (seq_db)
  { }

  // method to attach all unattached rows
  void attach_rows (ECFG_scores& ecfg, double prior_param = 0.,
		    double time_resolution = TINY, double time_max = DART_MAX_BRANCH_LENGTH, double time_min = 0.);
};


#endif /* ECFG_PLACER_INCLUDED */
