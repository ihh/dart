#include "tkf/tkfnodedp.h"
#include "util/rnd.h"

TKF_node_HMM::TKF_node_HMM (const TKF_align& tkf, double kT) : tkf (tkf), kT (kT), seq_scores (tkf.seq_scores) {}

void TKF_node_HMM::construct (Phylogeny::Node test, Phylogeny::Node parent, const Phylogeny::Node_vector& children)
{
  // set up local copies of a few commonly-used variables & types
  //
  test_node = test;
  parent_node = parent;
  child_node = children;

  TKF_align& tkfd = (TKF_align&) tkf;            // cast away const
  typedef Phylogeny::Undirected_pair Upair;    // useful typedef for later on
  
  child_scores.clear();
  for (int i = 0; i < (int) child_node.size(); ++i)
    child_scores.push_back (tkfd.branch_scores[Upair(test_node,child_node[i])]);
  parent_scores = has_parent() ? tkfd.branch_scores[Upair(test_node,parent_node)] : (TKF_branch_scores*) 0;

  // the main action: set up the label_trans and end_trans structures
  //
  // remember:   label_trans[src_state][transition_label] = Dest_tscore(dest_state,transition_score)
  //
  label_trans = vector<vector<Dest_tscore> > (states(), vector<Dest_tscore>(labels(), Dest_tscore(-1,0)));
  end_trans = vector<int> (states());
  for (int src = 0; src < states(); src++)             // loop through source states
    {
      for (int label = 0; label < labels(); label++)   // loop through transition labels from source state
	{
	  // the 0'th bit of the label is TRUE if we're emitting at the test node
	  // the (n+1)'th bit of the label is TRUE if we're emitting at the n'th neighbour node
	  //
	  // we build up the compound destination state & transition score accordingly, by considering the pair HMM for each neighbour.
	  // NB we also add a factor for the equilibrium length distribution onto the transition score
	  // (how this factor is related to the transition label depends on whether we have a parent or not).
	  // this is not really necessary, but kind of nice.... maybe...
	  //
	  int dest = 0;
	  int tscore = label & (has_parent() ? 2 : 1) ? seq_scores.eqmlen : 0;
	  for (int nbr = 0; nbr < neighbours(); nbr++)      // build up compound transition score for each neighbour
	    {
	      int nbr_src = (src & 3 << nbr*2) >> nbr*2;    // pair HMM source state for this neighbour
	      
	      int test_label = label & 1;                          // emitting at test node?
	      int nbr_label = (label & 2 << nbr) / (2 << nbr);     // emitting at neighbour node?
	      
	      int nbr_dest;
	      int nbr_score;
	      if (test_label == 0 && nbr_label == 0)        // if not emitting at either neighbour or test nodes, then make no contribution to transition score
		{
		  nbr_dest = nbr_src;
		  nbr_score = 0;
		}
	      else                                          // if emitting at one or both of neighbour and test nodes, then work out the pair HMM destination state ...
		{
		  nbr_dest = is_parent(nbr)
		    ? TKF_branch_scores::make_pair_HMM_state (nbr_label, test_label)
		    : TKF_branch_scores::make_pair_HMM_state (test_label, nbr_label);
		  
		  nbr_score = neighbour_scores(nbr) -> conditional_pair_HMM_transition_score (nbr_src, nbr_dest);    // ... and add on the transition score
		}
	      
	      dest += nbr_dest << nbr*2;
	      tscore += nbr_score;
	    }
	  ScaleScore (tscore, 1/kT);
	  label_trans[src][label] = Dest_tscore (dest, tscore);
	}

      // work out the end transition scores for each state
      //
      int escore = seq_scores.not_eqmlen;
      for (int nbr = 0; nbr < neighbours(); nbr++)
	{
	  int nbr_src = (src & 3 << nbr*2) >> nbr*2;
	  escore += neighbour_scores(nbr) -> conditional_pair_HMM_transition_score (nbr_src, TKF_branch_scores::pair_HMM_end_state);
	}
      ScaleScore (end_trans[src] = escore, 1/kT);
    }

  // add on the self-looping correction for the null (insert) state
  //
  null_state = label_trans[0][1].first;                       // quick way of getting the index of the null state
  null_extend = label_trans[null_state][1].second;
  int tmp_extend = null_extend;
  ScaleScore (tmp_extend, kT);
  ScaleScore (null_sum = -Prob2Score (1 - Score2Prob (tmp_extend)), 1/kT);
  for (int label = 0; label < labels(); label++)
    label_trans[null_state][label].second += null_sum;
  label_trans[null_state][1].second = -InfinityScore;    // invalidate the self-loop transition, just in case
}

void TKF_node_HMM::show_components (ostream& o) const
{
  save_flags (o);
  left_align (o);

  o << "Node HMM info for node '" << tkf.tree.node_specifier (test_node);
  sstring::size_type maxnslen = 0;
  if (has_parent())
    {
      sstring ns = tkf.tree.node_specifier (parent_node);
      maxnslen = ns.size();
      o << "', parent '" << ns;
    }
  o << "', children";
  for_const_contents (Phylogeny::Node_vector, child_node, c)
    {
      sstring ns = tkf.tree.node_specifier (*c);
      maxnslen = max (maxnslen, ns.size());
      o << " '" << ns << "'";
    }
  o << ", kT=" << kT << "\n";
  
  o << "Individual node-to-neighbour pair HMM transition matrices (0=start, 1=c, 2=p, 3=pc, 4=end)\n";
  o << sstring (maxnslen, ' ') << "        src/dest 1 (c,lbl=01)  2 (p,lbl=10)  3 (pc,lbl=11) 4 (end)\n";

  for (int i = 0; i < neighbours(); i++)
    {
      right_align (o);
      o.width (maxnslen);
      o << tkf.tree.node_specifier (neighbour_node(i)) << ": 0,1,3 (s,c,pc)";
      left_align (o);

      const TKF_branch_scores& bs = *neighbour_scores(i);
      for (int d = 1; d <= 4; d++) { o << "  "; o.width(12); o << bs.conditional_pair_HMM_transition_score (0, d) / kT; }
      o << "\n";
      o << sstring (maxnslen, ' ') << "      2      (p)";
      for (int d = 1; d <= 4; d++) { o << "  "; o.width(12); o << bs.conditional_pair_HMM_transition_score (2, d) / kT; }
      o << "\n";
    }
  o << "Equilibrium length scores: extension score " << seq_scores.eqmlen / kT << ", termination score " << seq_scores.not_eqmlen / kT << "\n";

  restore_flags (o);
}

void TKF_node_HMM::show_combined (ostream& o) const
{
  save_flags (o);
  left_align (o);

  o << "Combined node HMM transition matrix (column title = transition label, row title = source state, entry = dest:score) for node '" << tkf.tree.node_specifier (test_node) << "'\n";
  o << "State numbers formatted in little-endian base 4 (0=start, 1=c, 2=p, 3=pc), labels in little-endian binary\n";
  o << "       ";
  for (int lbl = 0; lbl < labels(); lbl++)
    {
      o << "     ";
      for (int r = 0; r <= neighbours(); r++) o << ((lbl & 1 << r) >> r);
      o << sstring (11 - (neighbours() + 1), ' ');
    }
  o << "end\n";
  for (int src = 0; src < states(); src++)
    {
      o << sstring (7 - neighbours(), ' ');
      for (int r = 0; r < neighbours(); r++) o << ((src & 3 << 2*r) >> 2*r);
      for (int lbl = 0; lbl < labels(); lbl++)
	{
	  o << sstring (4 - neighbours(), ' ');
	  for (int r = 0; r < neighbours(); r++) o << ((label_trans[src][lbl].first & 3 << 2*r) >> 2*r);
	  o << ":";
	  o.width (10);
	  left_align (o);
	  o << label_trans[src][lbl].second << " ";
	}
      left_align (o);
      o << end_trans[src] << "\n";
    }
  o << "Null extend score = " << null_extend << ", null sum score = " << null_sum << "\n";

  restore_flags (o);
}


TKF_node_DP::TKF_node_DP (const TKF_align& tkf, Phylogeny::Node test_node, double kT)
  : tkf (tkf), test_node (test_node), hmm (tkf, kT)
{
  int test_row = tkf.node2row[test_node];
  if (test_row == -1) THROW Standard_exception ("During node DP: missing row from alignment");

  if ((parent_node = tkf.tree.parent[test_node]) != -1)
    {
      if (tkf.node2row[parent_node] == -1)
	parent_node = -1;                                      // forget about parent node if it's not in the alignment
      else
	neighbour_rows.push_back (tkf.node2row[parent_node]);
    }

  for_iterator (Phylogeny::Child_iter, c_node,
		tkf.tree.children_begin (test_node, parent_node),
		tkf.tree.children_end (test_node, parent_node))
    {
      int c_row = tkf.node2row[*c_node];
      if (c_row != -1)
	{
	  neighbour_rows.push_back (c_row);
	  child_node.push_back (*c_node);
	}
    }
  
  hmm.construct (test_node, parent_node, child_node);
  path = Subalignment_path (tkf.align.path, neighbour_rows, align_to_subpath_map);

  // set up the observed data arrays

  column_label = column_match_emit_score = column_gap_emit_score = vector<int> (path.columns());
  vector<int> seq_coords = tkf.align.path.create_seq_coords();
  int align_col = 0;
  for (int col = 0; col < path.columns(); col++)
    {
      while (align_to_subpath_map(1,align_col) == 0) tkf.align.path.inc_seq_coords (seq_coords, align_col++);      // find corresponding column in original alignment

      for (int row = 0; row < path.rows(); row++)
	column_label[col] += path(row,col) << (row+1);
      
      bool tmp = tkf.align.path(test_row,align_col);

      ((TKF_align&)tkf).align.path[test_row][align_col] = 1;
      ScaleScore (column_match_emit_score[col] = tkf.column_emit_score (align_col, seq_coords), 1/kT);

      ((TKF_align&)tkf).align.path[test_row][align_col] = 0;
      ScaleScore (column_gap_emit_score[col] = tkf.column_emit_score (align_col, seq_coords), 1/kT);

      ((TKF_align&)tkf).align.path[test_row][align_col] = tmp;

      tkf.align.path.inc_seq_coords (seq_coords, align_col++);
    }
  
  // allocate DP matrix & set boundary conditions

  dpm = vector<DP_column> (path.columns()+1, DP_column());
  dpm[0][0] = 0;
}


void TKF_node_DP::show (ostream& o, vector<DP_column>* cell_traceback_count) const
{
  save_flags (o);
  left_align (o);

  const Symbol_score_map wildcard_score_map = tkf.alphabet().flat_score_map (0);
  const char wildcard_char = tkf.alphabet().score2char (wildcard_score_map);

  o << "Node DP matrix " << (cell_traceback_count ? "after" : "before") << " traceback (kT=" << hmm.kT << "):\n";
  o << "Col Aln Label Match      Gap        State:score";
  if (cell_traceback_count) o << " (or state*score for cells on traceback path)";
  o << "\n";
  vector<int> seq_coords = tkf.align.path.create_seq_coords();
  for (int col = 0; col < path.columns(); col++)
    {
      o.width(4);
      o << col;
      
      for (int row = 0; row < path.rows(); row++)
	{
	  const int tkf_row = neighbour_rows[row];
	  o << (path(row,col)
		? (tkf.align.prof[tkf_row]
		   ? tkf.alphabet().int2char((*tkf.align.prof[tkf_row]).consensus(seq_coords[tkf_row]))
		   : wildcard_char)
		: '-');
	}
      tkf.align.path.inc_seq_coords (seq_coords, col);
      for (int spc = path.rows(); spc < 4; spc++) o << ' ';
      
      o << 'x';
      for (int r = 1; r <= path.rows(); r++) o << ((column_label[col] & 1 << r) >> r);
      o << sstring (6 - (path.rows() + 1), ' ');
      
      o.width(10);
      o << column_match_emit_score[col] << ' ';
      
      o.width(10);
      o << column_gap_emit_score[col];
      
      for_const_contents (DP_column, dpm[col+1], d)
	{
	  o << " ";
	  for (int r = 0; r < path.rows(); r++) o << (((*d).first & 3 << 2*r) >> 2*r);
	  o << (cell_traceback_count ? ((*cell_traceback_count)[col+1][(*d).first] ? "*" : ":") : ":");
	  o.width(10);
	  o << (*d).second;
	}
      
      o << '\n';
    }
  o << "End score: " << end_score << "\n";

  restore_flags (o);
}

int TKF_node_DP::path_transition_score (const Pairwise_path& align_to_sequence_map) const
{
  Pairwise_path tmp_path = align_to_subpath_map;
  tmp_path.swap_parent_child();
  tmp_path *= align_to_sequence_map;
  Pairwise_path& subpath_to_sequence_map = tmp_path;

  int score = 0;
  int state = 0;
  int subpath_col = -1;
  for (int col = 0; col < subpath_to_sequence_map.columns(); ++col)
    {
      if (subpath_to_sequence_map.parent() [col])
	{
	  int label = column_label [++subpath_col] + subpath_to_sequence_map.child() [col];
	  ScorePMulAcc (score, hmm.label_trans[state][label].second);
	  state = hmm.label_trans[state][label].first;
	}
      else
	if (subpath_to_sequence_map.child() [col])
	  {
	    ScorePMulAcc (score, hmm.label_trans[state][1].second);
	    state = hmm.label_trans[state][1].first;
	  }
    }
  ScorePMulAcc (score, hmm.end_trans[state]);
  return score;
}

int TKF_node_DP::path_emit_score (const Pairwise_path& align_to_sequence_map) const
{
  Pairwise_path tmp_path = align_to_subpath_map;
  tmp_path.swap_parent_child();
  tmp_path *= align_to_sequence_map;
  Pairwise_path& subpath_to_sequence_map = tmp_path;

  int score = 0;
  int subpath_col = -1;
  for (int col = 0; col < subpath_to_sequence_map.columns(); ++col)
    if (subpath_to_sequence_map.parent() [col])
      ScorePMulAcc (score, subpath_to_sequence_map.child() [col] ? column_match_emit_score [++subpath_col] : column_gap_emit_score [++subpath_col]);

  return score;
}

int TKF_node_DP::path_score (const Pairwise_path& align_to_sequence_map) const
{
  return ScorePMul (path_transition_score (align_to_sequence_map), path_emit_score (align_to_sequence_map));
}

void TKF_node_DP::do_max_dp()
{
  // DP is done as a "push"
  //
  for (int col = 0; col < path.columns(); col++)
    {
      const DP_column& src_col  = dpm[col];
      DP_column&       dest_col = dpm[col+1];
      
      const int gap_label   = column_label[col];
      const int match_label = gap_label | 1;

      const int gap_emit     = column_gap_emit_score[col];
      const int match_emit   = column_match_emit_score[col];

      // each state has two possible outgoing transitions for each column,
      // corresponding to whether the test sequence has a match or a gap
      // in that column

      for_const_contents (DP_column, src_col, src_cell_iter)
	{
	  push_max (*src_cell_iter, dest_col, gap_label,   gap_emit);
	  push_max (*src_cell_iter, dest_col, match_label, match_emit);
	}
    }

  end_score = -InfinityScore;
  const DP_column& src_col = dpm.back();
  for_const_contents (DP_column, src_col, src_cell_iter)
    end_score = max (end_score, (*src_cell_iter).second + hmm.end_trans[(*src_cell_iter).first]);

  if (CTAGGING(1,TKFNODEDP)) hmm.show_components (CL);
  if (CTAGGING(-1,TKFNODEDP)) hmm.show_combined (CL);
  if (CTAGGING(-2,TKFNODEDP)) show (CL);
}

Pairwise_path TKF_node_DP::max_traceback() const
{
  vector<int>       label_trace;
  vector<DP_column> cell_traceback_count (dpm.size(), DP_column());

  // find last cell on traceback path
  //
  DP_cell current_cell;
  for_const_contents (DP_column, dpm.back(), last_col_iter)
    if (max (-InfinityScore, (*last_col_iter).second + hmm.end_trans[(*last_col_iter).first]) == end_score)
      {
	current_cell = *last_col_iter;
	break;
      }
  cell_traceback_count.back()[current_cell.first] = 1;

  // do traceback
  //
  for (int col = path.columns() - 1; col >= 0; col--)
    {
      const DP_column& src_col  = dpm[col];
      
      const int gap_label   = column_label[col];
      const int match_label = gap_label | 1;
      
      const int gap_emit     = column_gap_emit_score[col];
      const int match_emit   = column_match_emit_score[col];
      
      bool traceback_worked = 0;
      for_const_contents (DP_column, src_col, src_cell_iter)
	{
	  if (dest_cell (*src_cell_iter, gap_label, gap_emit) == current_cell)          // gap-label transition
	    {
	      current_cell = *src_cell_iter;
	      label_trace.push_back (gap_label);
	      traceback_worked = 1;
	      break;
	    }
	  if (dest_cell (*src_cell_iter, match_label, match_emit) == current_cell)      // match-label transition
	    {
	      current_cell = *src_cell_iter;
	      label_trace.push_back (match_label);
	      traceback_worked = 1;
	      break;
	    }
	}
      cell_traceback_count[col][current_cell.first] = 1;
      if (!traceback_worked) THROW Standard_exception ("Traceback failed");
    }
  
  if (CTAGGING(2,TKFNODEDP)) show (CL, &cell_traceback_count);

  return convert_reversed_label_trace_to_map (label_trace);
};


void TKF_node_DP::do_sum_dp()
{
  // DP is done as a "push"
  //
  for (int col = 0; col < path.columns(); col++)
    {
      const DP_column& src_col  = dpm[col];
      DP_column&       dest_col = dpm[col+1];
      
      const int gap_label   = column_label[col];
      const int match_label = gap_label | 1;

      const int gap_emit     = column_gap_emit_score[col];
      const int match_emit   = column_match_emit_score[col];

      // each state has two possible outgoing transitions for each column,
      // corresponding to whether the test sequence has a match or a gap
      // in that column

      for_const_contents (DP_column, src_col, src_cell_iter)
	{
	  push_sum (*src_cell_iter, dest_col, gap_label,   gap_emit);
	  push_sum (*src_cell_iter, dest_col, match_label, match_emit);
	}

      // do the null transitions..
      //
      int null_score = -InfinityScore;
      for_const_contents (DP_column, dest_col, dest_col_iter)
	ScorePSumAcc (null_score, (*dest_col_iter).second + hmm.null_start ((*dest_col_iter).first));
      dest_col.insert (DP_cell (hmm.null_state, null_score));
    }

  end_score = -InfinityScore;
  const DP_column& src_col = dpm.back();
  for_const_contents (DP_column, src_col, src_cell_iter)
    end_score = max (end_score, (*src_cell_iter).second + hmm.end_trans[(*src_cell_iter).first]);

  if (CTAGGING(1,TKFNODEDP)) hmm.show_components (CL);
  if (CTAGGING(-1,TKFNODEDP)) hmm.show_combined (CL);
  if (CTAGGING(-2,TKFNODEDP)) show (CL);
}


Pairwise_path TKF_node_DP::sample_traceback() const
{
  vector<int>       label_trace;
  vector<DP_column> cell_traceback_count (dpm.size(), DP_column());

  // find last cell on traceback path
  //
  vector<DP_cell> last_cell;
  vector<int>     last_sc;
  for_const_contents (DP_column, dpm.back(), last_col_iter)
    {
      last_cell.push_back (*last_col_iter);
      last_sc.push_back ((*last_col_iter).second + hmm.end_trans[(*last_col_iter).first]);
    }
  DP_cell current_cell = last_cell [Rnd::choose (Score2ProbVecNorm (last_sc))];
  cell_traceback_count.back()[current_cell.first] = 1;

  // do traceback
  //
  int col = path.columns() - 1;
  while (col >= 0)
    {
      const int current_state = current_cell.first;
      const int src_col_index = current_state == hmm.null_state ? col+1 : col;

      const DP_column& dest_col = dpm[col+1];
      
      // check if we're in a null insertion; if so, do a geometrically distributed number of self-loops
      //
      if (current_state == hmm.null_state)
	{
	  double null_extend_weight = Score2Prob (hmm.null_extend);
	  double null_stop_weight   = Score2Prob (-hmm.null_sum);
	  double null_extend_prob   = null_extend_weight / (null_extend_weight + null_stop_weight);
	  while (Rnd::decide(null_extend_prob))
	    {
	      label_trace.push_back(1);
	      ++cell_traceback_count[src_col_index][current_cell.first];
	    }
	}

      // find previous cell on traceback path
      //
      vector<DP_cell> src_cell;
      vector<int>     src_label;
      vector<int>     src_sc;
      DP_cell         tmp_cell;
      
      if (current_state == hmm.null_state)             // if we're in the null state, look for null-label transitions from this column
	{
	  for_const_contents (DP_column, dest_col, src_cell_iter)
	    if ((*src_cell_iter).first != current_state)             // already done the self-looping, so don't do it again
	      {
		tmp_cell = dest_cell (*src_cell_iter, 1, 0);
		if (tmp_cell.first == current_state)
		  {
		    src_cell.push_back (*src_cell_iter);
		    src_label.push_back (1);
		    src_sc.push_back (tmp_cell.second);
		  }
	      }
	}
      else                                             // not in the null state, so look for gap-label and match-label transitions from previous column
	{
	  const DP_column& src_col  = dpm[col];
	  
	  const int gap_label   = column_label[col];
	  const int match_label = gap_label | 1;
	  
	  const int gap_emit     = column_gap_emit_score[col];
	  const int match_emit   = column_match_emit_score[col];

	  for_const_contents (DP_column, src_col, src_cell_iter)
	    {
	      tmp_cell = dest_cell (*src_cell_iter, gap_label, gap_emit);                       // gap-label transition
	      if (tmp_cell.first == current_state)
		{
		  src_cell.push_back (*src_cell_iter);
		  src_label.push_back (gap_label);
		  src_sc.push_back (tmp_cell.second);
		}
	      tmp_cell = dest_cell (*src_cell_iter, match_label, match_emit);                   // match-label transition
	      if (tmp_cell.first == current_state)
		{
		  src_cell.push_back (*src_cell_iter);
		  src_label.push_back (match_label);
		  src_sc.push_back (tmp_cell.second);
		}
	    }
	}
      if (src_cell.size() == 0) THROW Standard_exception ("Traceback failed");

      int i = Rnd::choose (Score2ProbVecNorm (src_sc));
      current_cell = src_cell[i];
      label_trace.push_back (src_label[i]);
      cell_traceback_count[src_col_index][current_cell.first] = 1;
      
      if (current_state != hmm.null_state) --col;
    }
  
  if (CTAGGING(2,TKFNODEDP)) show (CL, &cell_traceback_count);

  Pairwise_path result = convert_reversed_label_trace_to_map (label_trace);
  if (CTAGGING(5,TKFNODEDP))
    CL << "Sampled node-HMM path posterior log-probability = " << Score2Bits (path_score(result) - end_score) << " bits\n";
  return result;
}


Pairwise_path TKF_node_DP::convert_reversed_label_trace_to_map (const vector<int>& trace) const
{
  Pairwise_path subpath_to_new_row_map;
  for (int i = trace.size()-1; i >= 0; i--)
    subpath_to_new_row_map.append_column (trace[i] & hmm.labels()-2 ? 1 : 0, trace[i] & 1);
  return align_to_subpath_map * subpath_to_new_row_map;
}
