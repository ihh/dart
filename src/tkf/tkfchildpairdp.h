#ifndef TKF_CHILD_PAIR_DP_INCLUDED
#define TKF_CHILD_PAIR_DP_INCLUDED

#include "tkf/tkfdata.h"
#include "util/strsaver.h"

class TKF_child_pair_HMM : protected Stream_saver
{
 public:
  const TKF_align&          tkf;

  double                   kT;

  int                      parent_node;
  int                      x_node;
  int                      y_node;

  const TKF_seq_scores&    seq_scores;
  const TKF_branch_scores& x_branch_scores;
  const TKF_branch_scores& y_branch_scores;
  array2d<int>             xy_substitution_matrix;

  // XY pair HMM states (A=ancestor):
  // 0=AXY, 1=AX-, 2=A-Y, 3=--Y (from 0 or 1), 4=--Y (from 2), 5=-X-
  // Null state A-- is folded into (some) transitions.
  //
  array2d<int>             transition_score;       // XY pair HMM transition matrix
  vector<int>              start_score;
  vector<int>              end_score;
  array2d<int>             transition_uses_null;
  vector<int>              start_transition_uses_null;

  vector<int>              state_da;
  vector<int>              state_dx;
  vector<int>              state_dy;

  TKF_child_pair_HMM (const TKF_align& data, int parent_node, double kT = 1);

  void show (ostream& o) const;

 private:
  void try_direct_transition (int src, int dest, int score);  // either adopts proposed direct-transition score, or uses default null-transition score & sets transition_uses_null flag
};


struct TKF_child_pair_DP : Stream_saver
{
  const TKF_align&       tkf;

  int                   parent_node;
  int                   x_node;
  int                   y_node;

  const Score_profile&  x_prof;
  const Score_profile&  y_prof;

  TKF_child_pair_HMM    hmm;

  vector<int>           x_gap_score;
  vector<int>           y_gap_score;

  struct DP_cell { int score[6]; int& operator[] (int i) { return score[i]; } };
  array2d<DP_cell>      dpm;
  int                   end_score;

  TKF_child_pair_DP (const TKF_align& data, int parent_node, const Score_profile& x_prof, const Score_profile& y_prof, double kT = 1);

  void            do_max_dp();
  Alignment_path  max_traceback (Score_profile& parent_profile, bool parent_wild = 1) const;
  
  void            do_sum_dp();
  Alignment_path  sample_traceback (Score_profile& parent_profile, bool parent_wild = 1) const;    // NB this currently isn't done quite right: should sample hidden transitions through null state.
  
  void show (ostream& o, array2d<int>* cell_is_on_traceback = 0) const;

  int state_emit_score (int x_index, int y_index, int state) const
    {
      if (state == 1 || state == 5) return x_gap_score [x_index - 1];
      if (state == 2 || state == 3 || state == 4) return y_gap_score [y_index - 1];
      int score = -InfinityScore;
      for_const_contents (Symbol_score_map, x_prof[x_index-1], ssx)
	{
	  const int xsym = (*ssx).first;
	  int xscore = -InfinityScore;
	  for_const_contents (Symbol_score_map, y_prof[y_index-1], ssy)
	    xscore = ScorePSum (xscore, ((TKF_child_pair_HMM&)hmm).xy_substitution_matrix (xsym, (*ssy).first) + (*ssy).second);
	  score = ScorePSum (score, xscore + (*ssx).second);
	}
      return score;
    }
};

#endif
