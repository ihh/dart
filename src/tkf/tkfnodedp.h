#ifndef TKF_NODE_DP_INCLUDED
#define TKF_NODE_DP_INCLUDED

#include "tkf/tkfdata.h"
#include "util/strsaver.h"

class TKF_node_HMM : protected Stream_saver
{
public:
  const TKF_align&                tkf;

  double                         kT;

  Phylogeny::Node                test_node;       // the node that we're re-aligning
  Phylogeny::Node                parent_node;     // the parent of that node (set this to -1 if there's no information about the parent)
  Phylogeny::Node_vector         child_node;      // the children of that node (parent/child distinction is actually unnecessary, due to reversibility)

  const TKF_seq_scores&          seq_scores;      // scores that don't depend on branch lengths
  const TKF_branch_scores*       parent_scores;   // branch-dependent scores for the parent branch
  vector<TKF_branch_scores*>     child_scores;    // branch-dependent scores for the child branches

  bool has_parent() const { return parent_node >= 0; }                              // TRUE if we have data on the parent
  int  neighbours() const { return child_node.size() + (has_parent() ? 1 : 0); }    // number of neighbour nodes on which we have data, including the parent

  bool                     is_parent (int neighbour_index) const                    // TRUE if the n'th neighbour is the parent node
    { return has_parent() && neighbour_index == 0; }
  Phylogeny::Node          neighbour_node (int neighbour_index) const               // returns the n'th neighbour node
    { return is_parent(neighbour_index) ? parent_node : child_node[neighbour_index - (has_parent() ? 1 : 0)]; }
  const TKF_branch_scores* neighbour_scores (int neighbour_index) const             // returns the scores for the n'th branch
    { return is_parent(neighbour_index) ? parent_scores : child_scores[neighbour_index - (has_parent() ? 1 : 0)]; }

  // Bit N+1 of path HMM transition label is set for an emission in neighbour sequence N.
  // Bit 0 corresponds to the test sequence.
  //
  int labels()     const { return 2 << neighbours(); }         // number of possible transition labels

  // Bits (2N,2N+1) of path HMM state indicate the state of the pair HMM for neighbour sequence N and the test sequence.
  //
  int states()     const { return 2 << (neighbours() * 2); }   // total number of states, not including end state
  
  // State transition matrix is indexed by source state and transition label
  //  i.e. label_trans[src_state][transition_label] = Dest_tscore(dest_state,transition_score)
  //
  typedef pair<int,int>          Dest_tscore;  // first = destination state, second = transition score
  vector<vector<Dest_tscore> >   label_trans;  // transition matrix indexed by source state and transition label
  vector<int>                    end_trans;    // transition scores to end state
  int                            null_state;   // index of null state (i.e. state corresponding to an emission in the test sequence only)
  int                            null_extend;  // self-loop transition score in null state
  int                            null_sum;     // 1/(1-n) in score space, where n is the null self-loop transition probability (i.e. null_extend in probability space)

  // method to return transition scores into start state
  //
  int null_start(int src_state) const { return label_trans[src_state][1].second; }

  // constructor - but the real action is in the construct() method, below
  //
  TKF_node_HMM (const TKF_align& tkf, double kT = 1);

  // call the following method AFTER the constructor, to set things up
  //
  void construct (Phylogeny::Node test, Phylogeny::Node parent, const Phylogeny::Node_vector& children);

  // output methods for logging
  //
  void show_components (ostream& o) const;
  void show_combined (ostream& o) const;
};

class TKF_node_DP : protected Stream_saver
{
 public:
  const TKF_align&              tkf;

  Phylogeny::Node              test_node;
  Phylogeny::Node              parent_node;
  Phylogeny::Node_vector       child_node;
  TKF_node_HMM                 hmm;
  Pairwise_path                align_to_subpath_map;
  Subalignment_path            path;

  vector<int> neighbour_rows;  // mapping from our rows to TKF_align rows

  // Observed data.
  //
  // column_label is an encoding of the path data structure.
  // Bit (R+1) of column_label[C] is equal to path(R,C);
  // bit 0 is always zero (denotes the test row, which is "hidden").
  //
  // column_xxx_emit_score[C] holds the substitution score for column C,
  // depending on whether the test row has a '1' (xxx=match) or a '0' (xxx=gap)
  // in that column.
  //
  vector<int>                  column_label;
  vector<int>                  column_match_emit_score;
  vector<int>                  column_gap_emit_score;

  // DP matrix:     dpm [alignment_column+1] [state]  =  score
  //
  typedef map<int,int>        DP_column;    // key = state, data = score
  typedef DP_column::iterator DP_col_iter;
  typedef pair<int,int>       DP_cell;      // first = state, second = score
  vector<DP_column>           dpm;
  int                         end_score;

  TKF_node_DP (const TKF_align& tkf, Phylogeny::Node test_node, double kT = 1);

  int            path_transition_score (const Pairwise_path& align_to_sequence_map) const;
  int            path_emit_score (const Pairwise_path& align_to_sequence_map) const;
  int            path_score (const Pairwise_path& align_to_sequence_map) const;

  void           do_max_dp();
  Pairwise_path  max_traceback() const;        // returns a Pairwise_path mapping the old alignment row to the new sequence

  void           do_sum_dp();
  Pairwise_path  sample_traceback() const;     // returns a Pairwise_path mapping the old alignment row to the new sequence

  void           show (ostream& o, vector<DP_column>* cell_traceback_count = 0) const;

 private:

  DP_cell dest_cell (const DP_cell& src_cell, int label, int emit_score) const
    {
      const TKF_node_HMM::Dest_tscore& dest_tscore = hmm.label_trans[src_cell.first][label];
      return DP_cell (dest_tscore.first, max (src_cell.second + dest_tscore.second + emit_score, -InfinityScore));
    }

  void push_max (const DP_cell& src_cell, DP_column& dest_col, int label, int emit_score) const
    {
      const TKF_node_HMM::Dest_tscore& dest_tscore = hmm.label_trans[src_cell.first][label];
      int sc = max (src_cell.second + dest_tscore.second + emit_score, -InfinityScore);
      DP_col_iter dest_cell_iter = dest_col.find (dest_tscore.first);
      if (dest_cell_iter == dest_col.end())
	dest_col.insert (DP_cell (dest_tscore.first, sc));
      else
	(*dest_cell_iter).second = max ((*dest_cell_iter).second, sc);
    }

  void push_sum (const DP_cell& src_cell, DP_column& dest_col, int label, int emit_score) const
    {
      const TKF_node_HMM::Dest_tscore& dest_tscore = hmm.label_trans[src_cell.first][label];
      int sc = max (src_cell.second + dest_tscore.second + emit_score, -InfinityScore);
      DP_col_iter dest_cell_iter = dest_col.find (dest_tscore.first);
      if (dest_cell_iter == dest_col.end())
	dest_col.insert (DP_cell (dest_tscore.first, sc));
      else
	ScorePSumAcc ((*dest_cell_iter).second, sc);
    }

  Pairwise_path convert_reversed_label_trace_to_map (const vector<int>& trace) const;
};

#endif
