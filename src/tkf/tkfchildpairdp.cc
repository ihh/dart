#include "tkf/tkfchildpairdp.h"
#include "util/rnd.h"

TKF_child_pair_HMM::TKF_child_pair_HMM (const TKF_align& tkf, int parent_node, double kT) :
  tkf (tkf),

  kT (kT),

  parent_node (parent_node),
  x_node (tkf.tree.children (parent_node, tkf.tree.parent[parent_node]) [0]),
  y_node (tkf.tree.children (parent_node, tkf.tree.parent[parent_node]) [1]),

  seq_scores (tkf.seq_scores),
  x_branch_scores (*((TKF_align&)tkf).branch_scores[Phylogeny::Undirected_pair (parent_node, x_node)]),
  y_branch_scores (*((TKF_align&)tkf).branch_scores[Phylogeny::Undirected_pair (parent_node, y_node)]),
  xy_substitution_matrix (Prob2ScoreArray2d (tkf.params.submat_factory.create_joint_substitution_matrix (((TKF_align&)tkf).tree.branch_length (parent_node, x_node) + ((TKF_align&)tkf).tree.branch_length (parent_node, y_node)))),
  
  transition_score (6, 6, -InfinityScore),
  start_score (6, -InfinityScore),
  end_score (6, -InfinityScore),
  transition_uses_null (6, 6, 0),
  start_transition_uses_null (6, 0),

  state_da (6, 0),
  state_dx (6, 0),
  state_dy (6, 0)
{
  vector<int> score_to_loop (6);
  vector<int> score_from_loop (6);

  score_to_loop[0] = x_branch_scores.not_descen + y_branch_scores.not_descen;
  score_to_loop[1] = x_branch_scores.not_descen + y_branch_scores.not_orphan;
  score_to_loop[2] = x_branch_scores.not_orphan + y_branch_scores.not_descen;
  score_to_loop[3] = x_branch_scores.not_descen + y_branch_scores.not_descen;
  score_to_loop[4] = x_branch_scores.not_orphan + y_branch_scores.not_descen;
  score_to_loop[5] = x_branch_scores.not_descen;

  score_from_loop[0] = x_branch_scores.ancest     + y_branch_scores.ancest;
  score_from_loop[1] = x_branch_scores.ancest     + y_branch_scores.not_ancest;
  score_from_loop[2] = x_branch_scores.not_ancest + y_branch_scores.ancest;
  score_from_loop[3] = -InfinityScore;
  score_from_loop[4] = x_branch_scores.not_ancest + y_branch_scores.not_ancest + y_branch_scores.orphan;
  score_from_loop[5] = x_branch_scores.not_ancest + y_branch_scores.not_ancest + y_branch_scores.not_orphan + x_branch_scores.orphan;
  
  for (int dest = 0; dest < 6; dest++)
    for (int src = 0; src < 6; src++)
      transition_score(src,dest) = score_to_loop[src] + seq_scores.eqmlen + score_from_loop[dest];

  for (int src = 0; src < 6; src++) end_score[src] = score_to_loop[src] + seq_scores.not_eqmlen;
  for (int dest = 0; dest < 6; dest++) start_score[dest] = seq_scores.eqmlen + score_from_loop[dest];

  start_transition_uses_null[4] = 1;
  start_transition_uses_null[5] = 1;

  transition_score(0,3) = y_branch_scores.descen;
  transition_score(1,3) = y_branch_scores.orphan;
  transition_score(3,3) = y_branch_scores.descen;

  transition_uses_null (0, 4) = 1;
  transition_uses_null (1, 4) = 1;
  transition_uses_null (3, 4) = 1;
  transition_uses_null (5, 4) = 1;

  try_direct_transition (2, 4, y_branch_scores.descen);
  try_direct_transition (4, 4, y_branch_scores.descen);

  try_direct_transition (0, 5, x_branch_scores.descen + y_branch_scores.not_descen);
  try_direct_transition (1, 5, x_branch_scores.descen + y_branch_scores.not_orphan);
  try_direct_transition (2, 5, x_branch_scores.orphan + y_branch_scores.not_descen);
  try_direct_transition (3, 5, x_branch_scores.descen + y_branch_scores.not_descen);
  try_direct_transition (4, 5, x_branch_scores.orphan + y_branch_scores.not_descen);
  try_direct_transition (5, 5, x_branch_scores.descen);

  state_da[0] = state_da[1] = state_da[2] = 1;
  state_dx[0] = state_dx[1] = state_dx[5] = 1;
  state_dy[0] = state_dy[2] = state_dy[3] = state_dy[4] = 1;

  // multiply all scores by 1/kT

  for (int src = 0; src < 6; src++)
    {
      ScaleScore (start_score[src], 1/kT);      
      ScaleScore (end_score[src], 1/kT);
      for (int dest = 0; dest < 6; dest++)
	ScaleScore (transition_score(src,dest), 1/kT);
    }  
}

void TKF_child_pair_HMM::try_direct_transition (int src, int dest, int score)
{
  if (score > transition_score(src,dest))
    {
      transition_score(src,dest) = score;
      transition_uses_null(src,dest) = 0;
    }
  else
    transition_uses_null(src,dest) = 1;
}

void TKF_child_pair_HMM::show (ostream& o) const
{
  save_flags (o);
  left_align (o);

  o << "Child-pair HMM for node " << tkf.tree.node_specifier (parent_node);
  o << ", children " << tkf.tree.node_specifier (x_node) << " and " << tkf.tree.node_specifier (y_node) << ":\n";
  o << "       AXY        AX-        A-Y        --Y1       --Y2       -X-        End\n";
  o << "Start ";
  for (int s = 0; s < 6; s++) { o << " "; o.width(10); o << start_score[s]; }
  o << " -\n";
  for (int s = 0; s < 6; s++)
    {
      o << (state_da[s] ? "A" : "-");
      o << (state_dx[s] ? "X" : "-");
      o << (state_dy[s] ? "Y" : "-");
      if (s == 3 || s == 4) o << s-2; else o << " ";
      o << "  ";
      for (int d = 0; d < 6; d++) { o << " "; o.width(10); o << ((TKF_child_pair_HMM*)this)->transition_score(s,d); }
      o << " "; o.width(10); o << end_score[s] << "\n";
    }
  
  restore_flags (o);
}

  

TKF_child_pair_DP::TKF_child_pair_DP (const TKF_align& tkf, int parent_node, const Score_profile& x_prof, const Score_profile& y_prof, double kT) :
  tkf (tkf),

  parent_node (parent_node),
  x_node (tkf.tree.children (parent_node, tkf.tree.parent[parent_node]) [0]),
  y_node (tkf.tree.children (parent_node, tkf.tree.parent[parent_node]) [1]),

  x_prof (x_prof),
  y_prof (y_prof),

  hmm (tkf, parent_node),

  x_gap_score (x_prof.size(), -InfinityScore),
  y_gap_score (y_prof.size(), -InfinityScore),

  dpm (x_prof.size() + 1, y_prof.size() + 1)

{
  // check that parent sequence isn't already specified
  //
  if (tkf.align.prof[tkf.node2row[parent_node]] != 0)
    THROW Standard_exception ("Attempt to infer an already specified sequence");

  // pre-calculate unpaired emission scores for x & y sequences
  //
  for (int i = 0; i < (int) x_prof.size(); i++)
    for_const_contents (Symbol_score_map, x_prof[i], ssx)
      x_gap_score[i] = ScorePSum (x_gap_score[i], hmm.seq_scores.prior[(*ssx).first] + (*ssx).second);

  for (int i = 0; i < (int) y_prof.size(); i++)
    for_const_contents (Symbol_score_map, y_prof[i], ssy)
      y_gap_score[i] = ScorePSum (y_gap_score[i], hmm.seq_scores.prior[(*ssy).first] + (*ssy).second);
}

void TKF_child_pair_DP::do_max_dp()
{
  // boundary conditions
  //
  for (int s = 0; s < 6; s++) dpm(0,0)[s] = -InfinityScore;
  dpm(0,0)[0] = hmm.start_score[0] - hmm.transition_score(0,0);
  
  if (x_prof.size())
    {
      for (int s = 0; s < 6; s++) dpm(1,0)[s] = -InfinityScore;
      dpm(1,0)[1] = hmm.start_score[1] + x_gap_score[0];
      dpm(1,0)[5] = hmm.start_score[5] + x_gap_score[0];
      for (int x = 1; x < (int) x_prof.size(); x++)
	{
	  dpm(x+1,0)[1] = max (x_gap_score[x] + max (dpm(x,0)[1] + hmm.transition_score(1,1),
						     dpm(x,0)[5] + hmm.transition_score(5,1)),
			       -InfinityScore);
	  dpm(x+1,0)[5] = max (x_gap_score[x] + max (dpm(x,0)[1] + hmm.transition_score(1,5),
						     dpm(x,0)[5] + hmm.transition_score(5,5)),
			       -InfinityScore);
	  dpm(x+1,0)[0] = dpm(x+1,0)[2] = dpm(x+1,0)[3] = dpm(x+1,0)[4] = -InfinityScore;
	}
    }

  if (y_prof.size())
    {
      for (int s = 0; s < 6; s++) dpm(0,1)[s] = -InfinityScore;
      dpm(0,1)[2] = max (hmm.start_score[2] + y_gap_score[0], -InfinityScore);
      dpm(0,1)[3] = max (hmm.start_score[3] + y_gap_score[0], -InfinityScore);
      dpm(0,1)[4] = max (hmm.start_score[4] + y_gap_score[0], -InfinityScore);
      for (int y = 1; y < (int) y_prof.size(); y++)
	{
	  dpm(0,y+1)[2] = max (y_gap_score[y] + max (max (dpm(0,y)[2] + hmm.transition_score(2,2),
							  dpm(0,y)[3] + hmm.transition_score(3,2)),
						     dpm(0,y)[4] + hmm.transition_score(4,2)),
			       -InfinityScore);
	  dpm(0,y+1)[3] = max (y_gap_score[y] + max (max (dpm(0,y)[2] + hmm.transition_score(2,3),
							  dpm(0,y)[3] + hmm.transition_score(3,3)),
						     dpm(0,y)[4] + hmm.transition_score(4,3)),
			       -InfinityScore);
	  dpm(0,y+1)[4] = max (y_gap_score[y] + max (max (dpm(0,y)[2] + hmm.transition_score(2,4),
							  dpm(0,y)[3] + hmm.transition_score(3,4)),
						     dpm(0,y)[4] + hmm.transition_score(4,4)),
			       -InfinityScore);
	  dpm(0,y+1)[0] = dpm(0,y+1)[1] = dpm(0,y+1)[5] = -InfinityScore;
	}
    }

  // main recursion - done as a pull
  //
  for (int y = 1; y <= (int) y_prof.size(); y++)
    for (int x = 1; x <= (int) x_prof.size(); x++)
      {
	DP_cell& cur_cell = dpm(x,y);
	for (int dest = 0; dest < 6; dest++)
	  {
	    DP_cell& prev_cell = dpm (x - hmm.state_dx[dest], y - hmm.state_dy[dest]);
	    int score = -InfinityScore;
	    for (int src = 0; src < 6; src++) score = max (score, prev_cell[src] + hmm.transition_score (src, dest));
	    cur_cell[dest] = max (score + state_emit_score (x, y, dest), -InfinityScore);
	  }
      }
  
  // find end score
  //
  end_score = -InfinityScore;
  DP_cell& end_cell = dpm (x_prof.size(), y_prof.size());
  for (int src = 0; src < 6; src++) end_score = max (end_score, end_cell[src] + hmm.end_score[src]);

  if (CTAGGING(1,TKF_HMM)) hmm.show (CL);
  if (CTAGGING(-2,TKF_PRETRACE)) show (CL);
}


Alignment_path TKF_child_pair_DP::max_traceback (Score_profile& parent_profile, bool parent_wild) const
{
  // parent_profile[A][X] =
  //  sum of likelihoods of all subtrees rooted at parent that have symbol X at position A of parent.
  //
  parent_profile.clear();
  Alignment_path axy (3);
  Symbol_score_map log1 = tkf.alphabet().flat_score_map(0);

  array2d<int> cell_is_on_traceback (x_prof.size() + 1, y_prof.size() + 1, 0);

  int x = x_prof.size();
  int y = y_prof.size();
  int state;
  for (state = 5; state >= 0; --state) if (max (((TKF_child_pair_DP*)this)->dpm(x,y)[state] + hmm.end_score[state], -InfinityScore) == end_score) break;
  if (state < 0) THROW Standard_exception ("Traceback error");
  int path_score = hmm.end_score[state];

  while (x != 0 || y != 0)
    {
      cell_is_on_traceback(x,y) |= 1 << state;
      int emit_sc = state_emit_score (x, y, state);
      ScorePMulAcc (path_score, emit_sc);
      int score_without_emit = ((TKF_child_pair_DP*)this)->dpm(x,y)[state] - emit_sc;

      // calculate next position in parent profile
      //
      if (hmm.state_da[state])
	{
	  parent_profile.push_back (log1);
	  for_contents (Symbol_score_map, parent_profile.back(), ssa)
	    {
	      if (hmm.state_dx[state])
		{
		  int x_score = -InfinityScore;
		  for_const_contents (Symbol_score_map, x_prof[x-1], ssx)
		    x_score = ScorePSum (x_score, ((array2d<int>&)hmm.x_branch_scores.cond_submat) ((*ssa).first, (*ssx).first) + (*ssx).second);
		  ScorePMulAcc ((*ssa).second, x_score);
		}
	      if (hmm.state_dy[state])
		{
		  int y_score = -InfinityScore;
		  for_const_contents (Symbol_score_map, y_prof[y-1], ssy)
		    y_score = ScorePSum (y_score, ((array2d<int>&)hmm.y_branch_scores.cond_submat) ((*ssa).first, (*ssy).first) + (*ssy).second);
		  ScorePMulAcc ((*ssa).second, y_score);
		}
	    }
	}

      // continue with traceback
      //
      x -= hmm.state_dx[state];
      y -= hmm.state_dy[state];

      int prev_st;
      for (prev_st = 5; prev_st >= 0; --prev_st)
	if (max (((TKF_child_pair_DP*)this)->dpm(x,y)[prev_st] + ((TKF_child_pair_HMM&)hmm).transition_score (prev_st, state), -InfinityScore) == score_without_emit) break;
      if (prev_st < 0) THROW Standard_exception ("Traceback error");
      ScorePMulAcc (path_score, hmm.transition_score (prev_st, state));

      axy.row(0).push_back (hmm.state_da[state]);
      axy.row(1).push_back (hmm.state_dx[state]);
      axy.row(2).push_back (hmm.state_dy[state]);

      if (((TKF_child_pair_HMM&)hmm).transition_uses_null (prev_st, state))
	{
	  axy.row(0).push_back (1);
	  axy.row(1).push_back (0);
	  axy.row(2).push_back (0);
	}
      
      state = prev_st;
    }
  
  if (((TKF_child_pair_HMM&)hmm).start_transition_uses_null[state])
    {
      axy.row(0).push_back (1);
      axy.row(1).push_back (0);
      axy.row(2).push_back (0);
    }
  ScorePMulAcc (path_score, hmm.start_score[state]);

  cell_is_on_traceback(0,0) |= 1 << state;
  if (CTAGGING(2,TKF_TRACE)) show (CL, &cell_is_on_traceback);

  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < axy.columns(); ++c)
      {
	// STL's swap() doesn't work here with gcc 4.0, hence neither does reverse()... hey ho...
	const bool old_val = axy[r][c];
	axy[r][c] = axy[r][axy.columns()-1-c];
	axy[r][axy.columns()-1-c] = old_val;
      }

  reverse (parent_profile.begin(), parent_profile.end());

  if (CTAGGING(5,TKF_PATH_POST)) CL << "Optimal child pair-HMM path posterior log-probability = " << Score2Bits (path_score - end_score) << " bits\n";

  if (!parent_wild)
    {
      parent_profile.normalise();
      Score_profile pp = Score_profile (parent_profile.consensus_dsq());
      if (CLOGGING(5)) CL << "Optimal parent sequence has posterior log-probability = " << Score2Bits (parent_profile.inner_product (pp)) << " bits\n";
      parent_profile.swap (pp);
    }

  return axy;
}


void TKF_child_pair_DP::do_sum_dp()
{
  // boundary conditions
  //
  for (int s = 0; s < 6; s++) dpm(0,0)[s] = -InfinityScore;
  dpm(0,0)[0] = hmm.start_score[0] - hmm.transition_score(0,0);
  
  if (x_prof.size())
    {
      for (int s = 0; s < 6; s++) dpm(1,0)[s] = -InfinityScore;
      dpm(1,0)[1] = hmm.start_score[1] + x_gap_score[0];
      dpm(1,0)[5] = hmm.start_score[5] + x_gap_score[0];
      for (int x = 1; x < (int) x_prof.size(); x++)
	{
	  dpm(x+1,0)[1] = max (x_gap_score[x] + ScorePSum (dpm(x,0)[1] + hmm.transition_score(1,1),
								    dpm(x,0)[5] + hmm.transition_score(5,1)),
			       -InfinityScore);
	  dpm(x+1,0)[5] = max (x_gap_score[x] + ScorePSum (dpm(x,0)[1] + hmm.transition_score(1,5),
								    dpm(x,0)[5] + hmm.transition_score(5,5)),
			       -InfinityScore);
	  dpm(x+1,0)[0] = dpm(x+1,0)[2] = dpm(x+1,0)[3] = dpm(x+1,0)[4] = -InfinityScore;
	}
    }

  if (y_prof.size())
    {
      for (int s = 0; s < 6; s++) dpm(0,1)[s] = -InfinityScore;
      dpm(0,1)[2] = max (hmm.start_score[2] + y_gap_score[0], -InfinityScore);
      dpm(0,1)[3] = max (hmm.start_score[3] + y_gap_score[0], -InfinityScore);
      dpm(0,1)[4] = max (hmm.start_score[4] + y_gap_score[0], -InfinityScore);
      for (int y = 1; y < (int) y_prof.size(); y++)
	{
	  dpm(0,y+1)[2] = max (y_gap_score[y] + ScorePSum (ScorePSum (dpm(0,y)[2] + hmm.transition_score(2,2),
											dpm(0,y)[3] + hmm.transition_score(3,2)),
								    dpm(0,y)[4] + hmm.transition_score(4,2)),
			       -InfinityScore);
	  dpm(0,y+1)[3] = max (y_gap_score[y] + ScorePSum (ScorePSum (dpm(0,y)[2] + hmm.transition_score(2,3),
											dpm(0,y)[3] + hmm.transition_score(3,3)),
								    dpm(0,y)[4] + hmm.transition_score(4,3)),
			       -InfinityScore);
	  dpm(0,y+1)[4] = max (y_gap_score[y] + ScorePSum (ScorePSum (dpm(0,y)[2] + hmm.transition_score(2,4),
											dpm(0,y)[3] + hmm.transition_score(3,4)),
								    dpm(0,y)[4] + hmm.transition_score(4,4)),
			       -InfinityScore);
	  dpm(0,y+1)[0] = dpm(0,y+1)[1] = dpm(0,y+1)[5] = -InfinityScore;
	}
    }

  // main recursion - done as a pull
  //
  for (int y = 1; y <= (int) y_prof.size(); y++)
    for (int x = 1; x <= (int) x_prof.size(); x++)
      {
	DP_cell& cur_cell = dpm(x,y);
	for (int dest = 0; dest < 6; dest++)
	  {
	    DP_cell& prev_cell = dpm (x - hmm.state_dx[dest], y - hmm.state_dy[dest]);
	    int score = -InfinityScore;
	    for (int src = 0; src < 6; src++) score = ScorePSum (score, prev_cell[src] + hmm.transition_score (src, dest));
	    cur_cell[dest] = max (score + state_emit_score (x, y, dest), -InfinityScore);
	  }
      }
  
  // find end score
  //
  end_score = -InfinityScore;
  DP_cell& end_cell = dpm (x_prof.size(), y_prof.size());
  for (int src = 0; src < 6; src++) end_score = ScorePSum (end_score, end_cell[src] + hmm.end_score[src]);

  if (CTAGGING(1,TKF_HMM)) hmm.show (CL);
  if (CTAGGING(-2,TKF_PRETRACE)) show (CL);
}


Alignment_path TKF_child_pair_DP::sample_traceback (Score_profile& parent_profile, bool parent_wild) const
{
  // parent_profile[A][X] =
  //  sum of likelihoods of all subtrees rooted at parent that have symbol X at position A of parent.
  //
  parent_profile.clear();
  Alignment_path axy (3);
  Symbol_score_map log1 = tkf.alphabet().flat_score_map(0);

  array2d<int> cell_is_on_traceback (x_prof.size() + 1, y_prof.size() + 1, 0);

  int x = x_prof.size();
  int y = y_prof.size();
  vector<int> sc;
  for (int s = 0; s < 6; s++) sc.push_back (((TKF_child_pair_DP*)this)->dpm(x,y)[s] + hmm.end_score[s]);
  int state = Rnd::choose (Score2ProbVecNorm (sc));
  int path_score = hmm.end_score[state];

  while (x != 0 || y != 0)
    {
      cell_is_on_traceback(x,y) |= 1 << state;
      int emit_sc = state_emit_score (x, y, state);
      ScorePMulAcc (path_score, emit_sc);

      // calculate next position in parent profile
      //
      if (hmm.state_da[state])
	{
	  parent_profile.push_back (log1);
	  for_contents (Symbol_score_map, parent_profile.back(), ssa)
	    {
	      if (hmm.state_dx[state])
		{
		  int x_score = -InfinityScore;
		  for_const_contents (Symbol_score_map, x_prof[x-1], ssx)
		    x_score = ScorePSum (x_score, ((array2d<int>&)hmm.x_branch_scores.cond_submat) ((*ssa).first, (*ssx).first) + (*ssx).second);
		  ScorePMulAcc ((*ssa).second, x_score);
		}
	      if (hmm.state_dy[state])
		{
		  int y_score = -InfinityScore;
		  for_const_contents (Symbol_score_map, y_prof[y-1], ssy)
		    y_score = ScorePSum (y_score, ((array2d<int>&)hmm.y_branch_scores.cond_submat) ((*ssa).first, (*ssy).first) + (*ssy).second);
		  ScorePMulAcc ((*ssa).second, y_score);
		}
	    }
	}

      // continue with traceback
      //
      x -= hmm.state_dx[state];
      y -= hmm.state_dy[state];

      sc.clear();
      for (int s = 0; s < 6; s++)
	sc.push_back (((TKF_child_pair_DP*)this)->dpm(x,y)[s] + hmm.transition_score (s, state));
      int prev_st = Rnd::choose (Score2ProbVecNorm (sc));
      ScorePMulAcc (path_score, hmm.transition_score (prev_st, state));

      axy.row(0).push_back (hmm.state_da[state]);
      axy.row(1).push_back (hmm.state_dx[state]);
      axy.row(2).push_back (hmm.state_dy[state]);

      if (((TKF_child_pair_HMM&)hmm).transition_uses_null (prev_st, state))
	{
	  parent_profile.push_back (log1);
	  axy.row(0).push_back (1);
	  axy.row(1).push_back (0);
	  axy.row(2).push_back (0);
	}
      
      state = prev_st;
    }
  
  if (((TKF_child_pair_HMM&)hmm).start_transition_uses_null[state])
    {
      parent_profile.push_back (log1);
      axy.row(0).push_back (1);
      axy.row(1).push_back (0);
      axy.row(2).push_back (0);
    }
  ScorePMulAcc (path_score, hmm.start_score[state]);

  cell_is_on_traceback(0,0) |= 1 << state;
  if (CTAGGING(2,TKF_TRACE)) show (CL, &cell_is_on_traceback);

  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < axy.columns(); ++c)
      {
	// STL's swap() doesn't work here with gcc 4.0, hence neither does reverse()... hey ho...
	const bool old_val = axy[r][c];
	axy[r][c] = axy[r][axy.columns()-1-c];
	axy[r][axy.columns()-1-c] = old_val;
      }

  reverse (parent_profile.begin(), parent_profile.end());

  if (CTAGGING(5,TKF_PATH_POST)) CL << "Sampled child pair-HMM path posterior log-probability = " << Score2Bits (path_score - end_score) << " bits\n";

  if (!parent_wild)
    {
      parent_profile.normalise();
      Score_profile pp = Score_profile (parent_profile.sample_dsq (hmm.kT));
      if (CLOGGING(5)) CL << "Sampled parent sequence has posterior log-probability = " << Score2Bits (parent_profile.inner_product (pp)) << " bits\n";
      parent_profile.swap (pp);
    }

  return axy;
}


void TKF_child_pair_DP::show (ostream& o, array2d<int>* cell_is_on_traceback) const
{
  save_flags (o);
  left_align(o);

  o << "Child-pair DP matrix " << (cell_is_on_traceback ? "after" : "before") << " traceback for node " << tkf.tree.node_specifier (parent_node) << ", children " << tkf.tree.node_specifier (x_node) << " and " << tkf.tree.node_specifier (y_node) << "\n";
  o << "Matrix scores are listed for states AXY, AX-, A-Y, --Y1(from AXY/AX-), --Y2(from A-Y), -X-";
  if (cell_is_on_traceback) o << " ('*' before score indicates cell is on traceback path)";
  o << "\n";
  o << "  State";
  const Alphabet& alphabet = tkf.alphabet();
  for (int x = 0; x <= (int) x_prof.size(); x++)
    { o << " "; o.width(10); o << (x > 0 ? alphabet.int2char (x_prof.consensus (x-1)) : '*') << " "; }
  o << "\n";
  for (int y = 0; y <= (int) y_prof.size(); y++)
    for (int s = 0; s < 6; s++)
      {
	o << (s==0 ? (y > 0 ? alphabet.int2char (y_prof.consensus (y-1)) : '*') : ' ') << " ";
	o << (hmm.state_da[s] ? "A" : "-");
	o << (hmm.state_dx[s] ? "X" : "-");
	o << (hmm.state_dy[s] ? "Y" : "-");
	if (s == 3 || s == 4) o << s-2; else o << " ";
	for (int x = 0; x <= (int) x_prof.size(); x++)
	  {
	    o << " ";
	    o << ((cell_is_on_traceback ? ((*cell_is_on_traceback)(x,y) & 1<<s) : 0) ? "*" : " ");
	    o.width(10);
	    o << ((TKF_child_pair_DP*)this)->dpm(x,y)[s];
	  }
	o << "\n";
      }
  o << "End score is " << end_score << "\n";
  restore_flags (o);
}

