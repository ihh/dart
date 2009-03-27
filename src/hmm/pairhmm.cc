#include "hmm/pairhmm.h"
#include "util/rnd.h"

// fraction by which forward & backward scores are allowed to differ, before a warning is issued
#define MAX_PERMISSIBLE_FORWARD_BACKWARD_DIFFERENCE_ERROR 0.01

void Pair_HMM_scores::scale_all_scores (double beta)
{
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	ScaleScore (transition(src,dest), beta);
      for_contents (vector<int>, single_emit[src], sc)
	ScaleScore (*sc, beta);
      for_contents (array2d<int>, pair_emit[src], sc)
	ScaleScore (*sc, beta);
      ScaleScore (start[src], beta);
      ScaleScore (end[src], beta);
    }
  ScaleScore (start_to_end(), beta);
}

Regexp Pair_HMM_scores::header_line_regexp ("^[ \t]*#");
Regexp Pair_HMM_scores::nonempty_line_regexp ("[^ \t]");

void Pair_HMM_scores::read_BLAST_matrix (istream& in, int state)
{
  if (state_type[state] != EmitXY)
    THROWEXPR ("Attempt to read a BLAST matrix into state " << state << ", which has type " << state_type[state]);

  array2d<Score>& state_pair_emit = pair_emit[state];
  for (int i = 0; i < state_pair_emit.xsize(); ++i)
    for (int j = 0; j < state_pair_emit.ysize(); ++j)
      state_pair_emit(i,j) = -InfinityScore;

  sstring s;
  vector<int> col;
  int row = 0;
  while (in && !in.eof())
    {
      s.getline (in);
      if (!nonempty_line_regexp.Match (s.c_str()))
	continue;
      if (header_line_regexp.Match (s.c_str()))
	continue;
      s.chomp();

      if (row == 0)
	{
	  vector<sstring> alph_char = s.split();
	  col = vector<int> (alph_char.size());
	  if (alphabet)
	    for (int i = 0; i < (int) alph_char.size(); ++i)
	      {
		if (alph_char[i].size() != 1)
		  THROWEXPR ("Bad alphabet symbol '" << alph_char[i] << "' in BLAST matrix");
		const int j = alphabet->char2int_deg (alph_char[i][0]);
		col[i] = (j < 0 || j >= alphabet->size()) ? -1 : j;
	      }
	  else
	    for (int i = 0; i < (int) alph_char.size(); ++i)
	      col[i] = (i >= state_pair_emit.xsize()) ? -1 : i;
	}
      else
	{
	  vector<sstring> row_data = s.split();
	  if (row_data.size() != col.size() + 1)
	    THROWEXPR ("Bad number of elements in row of BLAST matrix '" << s << "'");
	  int i = row - 1;
	  if (alphabet)
	    {
	      if (row_data[0].size() != 1)
		THROWEXPR ("Bad alphabet symbol '" << row_data[0] << "' in BLAST matrix");
	      const int k = alphabet->char2int_deg (row_data[0][0]);
	      i = (k < 0 || k >= alphabet->size()) ? -1 : k;
	    }
	  if (i >= 0)
	    for (int j = 0; j < (int) col.size(); ++j)
	      if (col[j] >= 0)
		state_pair_emit(i,col[j]) = row_data[j+1].to_int();
	}

      ++row;
    }
}

Score Pair_HMM_scores::pairwise_path_score (const Score_profile& xseq, const Score_profile& yseq, const Pairwise_path& path)
{
  // current implementation is compact but inefficient:
  // creates an alignment envelope allowing only one alignment path,
  // then creates Forward DP matrix for pair HMM
  // thus time & memory complexity is O(L^2) rather than O(L)
  array2d<int> pair_env (xseq.size() + 1, yseq.size() + 1, (int) 0);
  int x = 0, y = 0;
  for (int i = 0;; ++i)
    {
      pair_env(x,y) = 1;
      if (i == path.columns())
	break;
      if (path(0,i)) ++x;
      if (path(1,i)) ++y;
    }
  Pair_forward_DP_matrix fwd (*this, xseq, yseq, pair_env);
  return fwd.final_score;
}

vector<int> Pair_HMM_scores::find_optimal_state_path (const Pairwise_path& path) const
{
  vector<int> state_path;
  if (path.columns())
    {
      array2d<int> score (path.columns(), states());
      array2d<int> ok (path.columns(), states());
      for (int s = 0; s < states(); s++)
	if (dx(s)==path(0,0) && dy(s)==path(1,0))
	  {	  
	    ok(0,s) = 1;
	    score(0,s) = start[s];
	  }
	else
	  ok(0,s) = 0;
      
      for (int x = 1; x < path.columns(); x++)
	for (int s = 0; s < states(); s++)
	  if (dx(s)==path(0,x) && dy(s)==path(1,x))
	    {
	      ok(x,s) = 1;
	      int cell_score = -InfinityScore;
	      for (int t = 0; t < states(); t++)
		if (ok(x-1,t))
		  cell_score = max (cell_score, score(x-1,t) + transition(t,s));
	      score(x,s) = cell_score;
	    }
	  else
	    ok(x,s) = 0;

      const int last_col = path.columns() - 1;
      int bt = -1;
      int bt_score = -InfinityScore;
      for (int s = 0; s < states(); s++)
	if (ok(last_col,s) && (score(last_col,s) + end[s] > bt_score || bt < 0))
	  bt_score = score(last_col,s) + end[bt=s];
      if (bt < 0) THROW Standard_exception ("No state path compatible with alignment");

      state_path.push_back (bt);
      bt_score -= end[bt];
      for (int x = path.columns()-2; x >= 0; --x)
	{
	  int bt_prev = -1;
	  for (int s = 0; s < states(); s++)
	    if (ok(x,s) && max(-InfinityScore,score(x,s) + transition(s,bt)) == bt_score)
	      bt_prev = s;
	  if (bt_prev < 0) THROW Standard_exception ("Traceback failed");
	  state_path.push_back (bt = bt_prev);
	  bt_score = score(x,bt_prev);
	}

      reverse (state_path.begin(), state_path.end());
    }
  return state_path;
}

Pairwise_path Pair_HMM_scores::convert_state_path_to_alignment (const vector<int>& state_path) const
{
  Pairwise_path align;
  for_const_contents (vector<int>, state_path, s) align.append_column (dx(*s), dy(*s));
  return align;
}

int Pair_HMM_scores::path_transition_score (const vector<int>& state_path) const
{
  if (state_path.size() == 0) return start_to_end();
  int sc = start[state_path[0]];
  for (int i = 1; i < (int) state_path.size(); ++i) ScorePMulAcc (sc, transition (state_path[i-1], state_path[i]));
  ScorePMulAcc (sc, end[state_path.back()]);
  return sc;
}

int Pair_HMM_scores::path_emit_score (const vector<int>& state_path, const Score_profile& xseq, const Score_profile& yseq) const
{
  Pair_emit_score_calculator calculator (*this, xseq, yseq);
  int xpos = 0;
  int ypos = 0;
  int sc = 0;
  for_const_contents (vector<int>, state_path, s)
    {
      xpos += dx(*s);
      ypos += dy(*s);
      ScorePMulAcc (sc, calculator.calc_emit_score (*s, xpos, ypos));
    }
  return sc;
}

int Pair_HMM_scores::path_score (const vector<int>& state_path, const Score_profile& xseq, const Score_profile& yseq) const
{
  return ScorePMul (path_transition_score (state_path), path_emit_score (state_path, xseq, yseq));
}


Pair_HMM_derivatives::Pair_HMM_derivatives (int states, const Alphabet* alphabet) : Pair_HMM<double>(states,0.,alphabet), wrt ("derivatives") { }
Pair_HMM_derivatives::Pair_HMM_derivatives (int states, const char* x, const Alphabet* alphabet) : Pair_HMM<double>(states,0.,alphabet), wrt()
{
  wrt << "derivatives w.r.t " << x;
}


vector<int> Pair_HMM_scores::null_states_unsorted() const
{
  vector<int> null;
  null.reserve (states());
  for (int s = 0; s < states(); ++s) if (state_type[s] == Null) null.push_back(s);
  return null;
}

vector<int> Pair_HMM_scores::emit_states() const
{
  vector<int> result;
  for (int s = 0; s < states(); ++s) if (state_type[s] != Null) result.push_back(s);
  return result;
}

vector<int> Pair_HMM_scores::null_states() const
{
  const vector<int> unsorted = null_states_unsorted();
  return Transition_methods::topological_sort (*this, unsorted);
}

void Pair_HMM_scores::eliminate_null_states()
{
  const Concrete_transition_probs old_tp = Transition_methods::score2prob (*this);
  const Concrete_transition_probs new_tp = Transition_methods::eliminate (old_tp, null_states_unsorted());
  assign_transition_matrix (Transition_methods::prob2score (new_tp));
}

Pair_HMM_mask::Pair_HMM_mask (const Pair_HMM_scores& hmm, bool init_flag) : Pair_HMM<int> (hmm.states(), hmm.alphabet)
{
  for (int s = 0; s < states(); s++)
    {
      state_type[s] = (State_type) hmm.state_type[s];
      switch (state_type[s])
	{
	case Null:
	  break;
	case EmitXY:
	  pair_emit[s] = array2d<int> (hmm.pair_emit[s].xsize(), hmm.pair_emit[s].ysize(), init_flag);
	  break;
	case EmitX:
	case EmitY:
	  single_emit[s] = vector<int> (hmm.single_emit[s].size(), init_flag);
	  break;
	default:
	  THROW Standard_exception ("Unknown state type!");
	  break;
	}
    }
  for (int s = 0; s < states(); ++s)
    {
      for (int d = 0; d < states(); ++d)
	transition(s,d) = init_flag;
      start[s] = end[s] = init_flag;
    }
  start_to_end() = init_flag;
}

void Pair_HMM_mask::set_incoming_transitions (int state) { for (int s = 0; s < states(); ++s) transition(state,s) = 1; }
void Pair_HMM_mask::clear_incoming_transitions (int state) { for (int s = 0; s < states(); ++s) transition(state,s) = 0; }
void Pair_HMM_mask::set_outgoing_transitions (int state) { for (int s = 0; s < states(); ++s) transition(s,state) = 1; }
void Pair_HMM_mask::clear_outgoing_transitions (int state) { for (int s = 0; s < states(); ++s) transition(s,state) = 0; }
void Pair_HMM_mask::set_emit (int state)
{
  for_contents (array2d<int>, pair_emit[state], i) *i = 1;
  for_contents (vector<int>, single_emit[state], i) *i = 1;
}
void Pair_HMM_mask::clear_emit (int state)
{
  for_contents (array2d<int>, pair_emit[state], i) *i = 0;
  for_contents (vector<int>, single_emit[state], i) *i = 0;
}
void Pair_HMM_mask::set_all_start() { for (int s = 0; s < states(); ++s) start[s] = 1; start_to_end() = 1; }
void Pair_HMM_mask::clear_all_start() { for (int s = 0; s < states(); ++s) start[s] = 0; start_to_end() = 0; }
void Pair_HMM_mask::set_all_end() { for (int s = 0; s < states(); ++s) end[s] = 1; start_to_end() = 1; }
void Pair_HMM_mask::clear_all_end() { for (int s = 0; s < states(); ++s) end[s] = 0; start_to_end() = 0; }
void Pair_HMM_mask::set_all_transitions() { for (int s = 0; s < states(); ++s) set_outgoing_transitions(s); }
void Pair_HMM_mask::clear_all_transitions() { for (int s = 0; s < states(); ++s) clear_outgoing_transitions(s); }
void Pair_HMM_mask::set_all_emit() { for (int s = 0; s < states(); ++s) set_emit(s); }
void Pair_HMM_mask::clear_all_emit() { for (int s = 0; s < states(); ++s) clear_emit(s); }

Pair_HMM_counts::Pair_HMM_counts (const Pair_HMM_scores& hmm) : Pair_HMM<double> (hmm.states(), hmm.alphabet)
{
  for (int s = 0; s < states(); s++)
    {
      state_type[s] = (State_type) hmm.state_type[s];
      single_emit[s] = vector<double> (hmm.single_emit[s].size(), (double) 0);
      pair_emit[s] = array2d<double> (hmm.pair_emit[s].xsize(), hmm.pair_emit[s].ysize(), (double) 0);
      start[s] = end[s] = 0.;
      for (int d = 0; d < states(); ++d)
	transition(s,d) = 0.;
    }
  start_to_end() = 0;
  log_likelihood = 0;
}

Score Pair_HMM_counts::add_counts_from_unaligned_sequences (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq)
{
  Pair_forward_backward_DP_matrix fb (hmm, xseq, yseq);
  add_counts (fb.counts);
  return fb.forward_score;
}

Score Pair_HMM_counts::add_counts_from_aligned_sequences (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const Pairwise_path& path)
{
  array2d<int> pair_env (xseq.size() + 1, yseq.size() + 1, (int) 0);
  int x = 0, y = 0;
  for (int i = 0;; ++i)
    {
      pair_env(x,y) = 1;
      if (i == path.columns())
	break;
      if (path(0,i)) ++x;
      if (path(1,i)) ++y;
    }
  Pair_forward_backward_DP_matrix fb (hmm, xseq, yseq, pair_env);
  add_counts (fb.counts);
  return fb.forward_score;
}

Score Pair_HMM_counts::add_counts_from_state_path (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const vector<int>& state_path)
{
  Score score = 0;
  Score indel_score = 0;
  Score subst_score = 0;
  if (state_path.size() == 0)
    {
      score = indel_score = hmm.start_to_end();
      start_to_end() += 1;
    }
  else
    {
      Pair_emit_score_calculator calculator (hmm, xseq, yseq);
      int x = 1;
      int y = 1;
      int s = End;
      for (int i = 0; i < (int) state_path.size(); i++)
	{
	  s = state_path[i];
	  if (i == 0)
	    {
	      ScorePMulAcc (score, hmm.start[s]);
	      ScorePMulAcc (indel_score, hmm.start[s]);
	      start[s] += 1;
	    }
	  else
	    {
	      ScorePMulAcc (score, hmm.transition(state_path[i-1],s));
	      ScorePMulAcc (indel_score, hmm.transition(state_path[i-1],s));
	      transition(state_path[i-1],s) += 1;
	    }
	  const Score emit_sc = calculator.calc_emit_score(s,x,y);
	  ScorePMulAcc (score, emit_sc);
	  ScorePMulAcc (subst_score, emit_sc);
	  calculator.accum_emit_score (-emit_sc, *this, s, x, y);
	  x += hmm.dx(s);
	  y += hmm.dy(s);
	}
      ScorePMulAcc (score, hmm.end[s]);
      ScorePMulAcc (indel_score, hmm.end[s]);
      end[s] += 1;
    }

  if (CLOGGING(1)) { CL << "Pair-HMM counts from state path:\n"; show(CL); }

  log_likelihood += Score2Nats(score);
  return score;
}

void Pair_HMM_counts::add_counts (const Pair_HMM_counts& counts)
{
  if (!same_dimensions (counts)) THROW Standard_exception ("While adding pair-HMM counts: state spaces don't match");
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest) += counts.transition(src,dest);

      switch (state_type[src])
	{
	case Null:
	  break;
	case EmitXY:
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    for (int j = 0; j < pair_emit[src].ysize(); ++j)
	      pair_emit[src] (i, j) += counts.pair_emit[src] (i, j);
	  break;
	case EmitX:
	case EmitY:
	  for (int i = 0; i < (int) single_emit[src].size(); ++i)
	    single_emit[src][i] += counts.single_emit[src][i];
	  break;
	default:
	  THROW Standard_exception ("Unknown state type!");
	  break;
	}
      
      start[src] += counts.start[src];
      end[src] += counts.end[src];
    }
  start_to_end() += counts.start_to_end();
  log_likelihood += counts.log_likelihood;
}

void Pair_HMM_counts::subtract_counts (const Pair_HMM_counts& counts)
{
  if (!same_dimensions (counts)) THROW Standard_exception ("While adding pair-HMM counts: state spaces don't match");
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest) -= counts.transition(src,dest);

      switch (state_type[src])
	{
	case Null:
	  break;
	case EmitXY:
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    for (int j = 0; j < pair_emit[src].ysize(); ++j)
	      pair_emit[src] (i, j) -= counts.pair_emit[src] (i, j);
	  break;
	case EmitX:
	case EmitY:
	  for (int i = 0; i < (int) single_emit[src].size(); ++i)
	    single_emit[src][i] -= counts.single_emit[src][i];
	  break;
	default:
	  THROW Standard_exception ("Unknown state type!");
	  break;
	}
      
      start[src] -= counts.start[src];
      end[src] -= counts.end[src];
    }
  start_to_end() -= counts.start_to_end();
  log_likelihood -= counts.log_likelihood;
}

void Pair_HMM_counts::update_HMM_scores (Pair_HMM_scores& hmm, const Pair_HMM_mask& mask, bool sample, double kT) const
{
  if (!same_dimensions (hmm)) THROW Standard_exception ("While updating pair-HMM scores: state spaces of counts & scores structures don't match");
  if (!same_dimensions (mask)) THROW Standard_exception ("While updating pair-HMM scores: state spaces of counts & update mask structures don't match");
  int score_pr_product;
  vector<double> counts;
  vector<double> prob;
  int k;
  for (int src = 0; src < states(); ++src)
    {
      score_pr_product = -InfinityScore;
      counts.clear();
      for (int dest = 0; dest < states(); dest++)
	if (mask.transition(src,dest))
	  {
	    ScorePSumAcc (score_pr_product, hmm.transition(src,dest));
	    counts.push_back (transition(src,dest) / kT);
	  }
      prob = sample ? Rnd::sample_dirichlet(counts) : Rnd::dirichlet_max(counts);
      k = 0;
      for (int dest = 0; dest < states(); dest++)
	if (mask.transition(src,dest))
	  hmm.transition(src,dest) = Prob2Score (prob[k++]) + score_pr_product;
      
      switch (state_type[src])
	{
	case Null:
	  break;
	case EmitXY:
	  score_pr_product = -InfinityScore;
	  counts.clear();
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    for (int j = 0; j < pair_emit[src].ysize(); ++j)
	      if (mask.pair_emit[src](i,j))
		{
		  ScorePSumAcc (score_pr_product, hmm.pair_emit[src](i,j));
		  counts.push_back (pair_emit[src](i,j));
		}
	  prob = sample ? Rnd::sample_dirichlet(counts) : Rnd::dirichlet_max(counts);
	  k = 0;
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    for (int j = 0; j < pair_emit[src].ysize(); ++j)
	      if (mask.pair_emit[src](i,j))
		hmm.pair_emit[src](i,j) = Prob2Score (prob[k++]) + score_pr_product;
	  break;
	case EmitX:
	case EmitY:
	  score_pr_product = -InfinityScore;
	  counts.clear();
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    if (mask.single_emit[src][i])
	      {
		ScorePSumAcc (score_pr_product, hmm.single_emit[src][i]);
		counts.push_back (single_emit[src][i]);
	      }
	  prob = sample ? Rnd::sample_dirichlet(counts) : Rnd::dirichlet_max(counts);
	  k = 0;
	  for (int i = 0; i < (int) single_emit[src].size(); ++i)
	    if (mask.single_emit[src][i])
	      hmm.single_emit[src][i] = Prob2Score (prob[k++]) + score_pr_product;
	  break;
	default:
	  THROW Standard_exception ("Unknown state type!");
	  break;
	}
    }
  
  score_pr_product = -InfinityScore;
  counts.clear();
  for (int s = 0; s < states(); ++s)
    if (mask.start[s])
      {
	ScorePSumAcc (score_pr_product, hmm.start[s]);
	counts.push_back (start[s]);
      }
  if (mask.start_to_end())
    {
      ScorePSumAcc (score_pr_product, hmm.start_to_end());
      counts.push_back (start_to_end());
    }
  prob = sample ? Rnd::sample_dirichlet(counts) : Rnd::dirichlet_max(counts);
  k = 0;
  for (int s = 0; s < states(); ++s)
    if (mask.start[s])
      hmm.start[s] = Prob2Score (prob[k++]) + score_pr_product;
  if (mask.start_to_end())
    hmm.start_to_end() = Prob2Score (prob[k++]) + score_pr_product;
  
  score_pr_product = -InfinityScore;
  counts.clear();
  for (int s = 0; s < states(); ++s)
    if (mask.end[s])
      {
	ScorePSumAcc (score_pr_product, hmm.end[s]);
	counts.push_back (end[s]);
      }
  prob = sample ? Rnd::sample_dirichlet(counts) : Rnd::dirichlet_max(counts);
  k = 0;
  for (int s = 0; s < states(); ++s)
    if (mask.end[s])
      hmm.end[s] = Prob2Score (prob[k++]) + score_pr_product;
}

double Pair_HMM_counts::dloglike_dx (const Pair_HMM_scores& hmm, const Pair_HMM_derivatives& deriv) const
{
  //    dlog(P)/dx   =   sum_{T} dlog(P)/dT(partial) * dT/dX
  //                 =   (1/P) * sum_{T} dP/dT * dT/dX
  //                 =   sum_{T} n_T * dT/dX / T
  //
  //  where   n_T  =  expected T count  =  dP/dT * T / P  =  dlog(P)/dlog(T)
  //
  double sum = 0;
  for (int s = 0; s < states(); s++)
    {
      sum += start[s] * deriv.start[s] * Score2Prob(-hmm.start[s]);
      sum += end[s] * deriv.end[s] * Score2Prob(-hmm.end[s]);
      for (int t = 0; t < states(); t++)
	sum += transition(s,t) * deriv.transition(s,t) * Score2Prob(-hmm.transition(s,t));
      switch (hmm.state_type[s])
	{
	case Null:
	  break;
	case EmitX:
	case EmitY:
	  for (int i = 0; i < (int) single_emit[s].size(); i++)
	    sum += single_emit[s][i] * deriv.single_emit[s][i] * Score2Prob(-hmm.single_emit[s][i]);
	  break;
	case EmitXY:
	  for (int i = 0; i < pair_emit[s].xsize(); i++)
	    for (int j = 0; j < pair_emit[s].ysize(); j++)
	      sum += pair_emit[s](i,j) * deriv.pair_emit[s](i,j) * Score2Prob(-hmm.pair_emit[s](i,j));
	  break;
	default:
	  CLOG(6) << "Unknown state type!\n";
	  break;
	}
    }
  sum += start_to_end() * deriv.start_to_end() * Score2Prob(-hmm.start_to_end());
  return sum;
}

Pair_emit_score_calculator::Pair_emit_score_calculator (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq) :
  hmm(hmm), xseq(xseq), yseq(yseq), xseq_evolved (hmm.states(), (Score_profile*) 0)
{
  for (int s = 0; s < hmm.states(); ++s)
    if (hmm.state_type[s] == Pair_HMM_scores::EmitXY)
      {
	xseq_evolved[s] = new Score_profile;
	*xseq_evolved[s] = xseq.evolve (hmm.pair_emit[s]);
      }
}

Pair_emit_score_calculator::~Pair_emit_score_calculator()
{
  for_contents (vector<Score_profile*>, xseq_evolved, x) if (*x != 0) delete *x;
}

void Pair_DP_matrix_base::show (ostream& o) const
{
  save_flags (o);
  
  o << "                 0:*";
  for (int x = 1; x < xsize; x++)
    {
      o.width(9);
      right_align (o);
      o << x << ':' << score2char(xseq[x-1]);
    }
  o << "\n";
  for (int y = 0; y < ysize; y++)
    {
      for (int s = 0; s < states; s++)
	{
	  if (s == 0)
	    {
	      o.width(4);
	      right_align (o);
	      o << y << ':' << (y == 0 ? '*' : score2char(yseq[y-1])) << ' ';
	    }
	  else
	    o << "       ";
	  left_align (o);
	  o.width(2);
	  o << s;
	  for (int x = 0; x < xsize; x++)
	    {
	      o << " ";
	      right_align (o);
	      o.width(10);
	      ShowScore (cell[s](x,y), o);
	    }
	  o << "\n";
	}
    }
  o << "Final score: ";
  ShowScore (final_score, o);
  o << "\n";

  restore_flags (o);
}

void Pair_DP_matrix_base::show_sparse (ostream& o) const
{
  save_flags (o);
  
  o << "                 0:*";
  for (int x = 1; x < xsize; x++)
    {
      o.width(9);
      right_align (o);
      o << x << ':' << score2char(xseq[x-1]);
    }
  o << "\n";

  for (int y = 0; y < ysize; y++)
    for (int x = 0; x < xsize; x++)
      if (pair_env(x,y))
	{
	  o << "x=" << x;
	  if (x > 0)
	    o << '(' << score2char(xseq[x-1]) << ')';
	  o << ", y=" << y;
	  if (y > 0)
	    o << '(' << score2char(yseq[y-1]) << ')';
	  o << ':';
	  for (int s = 0; s < states; s++)
	    {
	      o << ' ' << hmm.state_name[s] << '(';
	      ShowScore (cell[s](x,y), o);
	      o << ')';
	    }
	  o << '\n';
	}
  o << "Final score (x=" << xsize-1 << ", y=" << ysize-1 << "): ";
  ShowScore (final_score, o);
  o << "\n";

  restore_flags (o);
}

sstring Pair_DP_matrix_base::hmm_dump() const
{
  // until sstring::operator<<(int) bug is fixed, we use the following hack
  ostream& dump = CLOGERR;
  //  sstring dump;
  dump << "Pair HMM:\n";
  hmm.show (dump);
  dump << "Pair HMM DP matrix:\n";
  show (dump);
  //  return dump;
  return sstring();
}

Pair_Viterbi_DP_matrix::Pair_Viterbi_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq)
  : Pair_DP_matrix_base (hmm, xseq, yseq)
{
  fill();
}

void Pair_Viterbi_DP_matrix::fill()
{
  // check pair envelope is sane
  if (!pair_env(0,0) || !pair_env(xsize-1,ysize-1))
    {
      CLOGERR << "Warning -- pair envelope excludes start or end states\n";
      return;
    }

  // get null states
  const vector<int> null_states = hmm.null_states();
  const vector<int> emit_states = hmm.emit_states();

  // initialise top-left cell
  for_const_contents (vector<int>, emit_states, s)
    cell[*s](0,0) = -InfinityScore;
  for_const_contents (vector<int>, null_states, s)
    cell[*s] (0,0) = calc_score (*s, 0, 0, hmm.start[*s]);

  if (xsize == 1 && ysize == 1)
    {
      final_score = hmm.start_to_end();
      for_const_contents (vector<int>, null_states, s)
	final_score = max (final_score, ScorePMul (cell[*s](0,0), hmm.end[*s]));
    }
  else
    {
      // do top and left boundary conditions
      //
      for (int x = 1; x < xsize; x++)
	if (pair_env(x,0))
	  for (int s = 0; s < states; s++)
	    switch (hmm.state_type[s])
	      {
	      case Pair_HMM_scores::EmitX:
		cell[s] (x,0) = calc_score (s, x, 0, x==1 ? hmm.start[s] : -InfinityScore);
		break;
	      case Pair_HMM_scores::Null:
	      case Pair_HMM_scores::EmitY:
	      case Pair_HMM_scores::EmitXY:
		cell[s] (x,0) = -InfinityScore;
		break;
	      default:
		CLOG(6) << "Unknown state type!\n";
		cell[s] (x,0) = -InfinityScore;
		break;
	      }
      
      for (int y = 1; y < ysize; y++)
	if (pair_env(0,y))
	  for (int s = 0; s < states; s++)
	    switch (hmm.state_type[s])
	      {
	      case Pair_HMM_scores::EmitY:
		cell[s] (0,y) = calc_score (s, 0, y, y==1 ? hmm.start[s] : -InfinityScore);
		break;
	      case Pair_HMM_scores::Null:
	      case Pair_HMM_scores::EmitX:
	      case Pair_HMM_scores::EmitXY:
		cell[s] (0,y) = -InfinityScore;
		break;
	      default:
		CLOG(6) << "Unknown state type!\n";
		cell[s] (0,y) = -InfinityScore;
		break;
	      }
      
      // fill the matrix
      //
      for (int x = 1; x < xsize; x++)
	for (int y = 1; y < ysize; y++)
	  if (pair_env(x,y))
	    {
	      for_const_contents (vector<int>, emit_states, s)
		cell[*s] (x,y) = calc_score (*s, x, y, x==1 && y==1 && hmm.state_type[*s] == Pair_HMM_scores::EmitXY ? hmm.start[*s] : -InfinityScore);
	      for_const_contents (vector<int>, null_states, s)
		cell[*s] (x,y) = calc_score (*s, x, y, -InfinityScore);
	    }
      
      // get the final score
      //
      for (int s = 0; s < states; s++)
	final_score = max (final_score, ScorePMul (cell[s] (xsize-1, ysize-1), hmm.end[s]));
    }
  
  if (CTAGGING(-1,VITERBI_HMM)) { CL << "Viterbi HMM:\n"; hmm.show(CL); }
  if (CTAGGING(-1,DP_MATRIX VITERBI_MATRIX)) { CL << "Viterbi matrix:\n"; show(CL); }
  if (CTAGGING(-2,DP_MATRIX_SPARSE VITERBI_MATRIX_SPARSE)) { CL << "Viterbi matrix (sparse):\n"; show_sparse(CL); }
  CTAG(4,DP_SCORE VITERBI_SCORE) << "Viterbi log-likelihood is " << Score2Bits(final_score) << " bits\n";
}


vector<int> Pair_Viterbi_DP_matrix::optimal_state_path() const
{
  if (final_score == -InfinityScore)
    THROWEXPR ("Final score is -infinity; traceback likely to break");

  vector<int> state_path;

  if (xsize > 1 || ysize > 1 || final_score != hmm.start_to_end())
    {
      int x = xsize - 1;
      int y = ysize - 1;
      int current_state = -1;
      for (int s = 0; s < hmm.states(); s++)
	if (ScorePMul (cell[s](x,y), hmm.end[s]) == final_score) { current_state = s; break; }
      if (current_state < 0)
	THROWEXPR (hmm_dump() <<  "Pair HMM: traceback failed at first hurdle");
      
      while (1)
	{
	  state_path.push_back (current_state);
	  const Pair_HMM_scores::State_type t = hmm.state_type [current_state];

	  if (((t & Pair_HMM_scores::EmitX) && x == 0) || ((t & Pair_HMM_scores::EmitY) && y == 0))
	    THROWEXPR (hmm_dump() << "Pair HMM: traceback out of bounds at x=" << x << ", y=" << y << ", state=" << current_state);

	  const Score emit_sc = calc_emit_score (current_state, x, y);
	  const Score current_score = cell[current_state] (x, y);
	  x -= hmm.dx (current_state);
	  y -= hmm.dy (current_state);
	  if (x == 0 && y == 0) break;
	  
	  int prev_state = -1;
	  for (int s = 0; s < hmm.states(); s++)
	    if (ScorePMul3 (cell[s](x,y), hmm.transition (s, current_state), emit_sc) == current_score) { prev_state = s; break; }
	  if (prev_state < 0)
	    THROWEXPR (hmm_dump() << "Pair HMM: traceback failed at x=" << x << ", y=" << y << ", state=" << current_state);
	  current_state = prev_state;
	}

      reverse (state_path.begin(), state_path.end());
    }

  return state_path;
}


Pair_forward_DP_matrix::Pair_forward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq)
  : Pair_DP_matrix_base (hmm, xseq, yseq)
{
  fill();
}

Pair_forward_DP_matrix::Pair_forward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const array2d<int>& pair_env)
  : Pair_DP_matrix_base (hmm, xseq, yseq, pair_env)
{
  fill();
}

void Pair_forward_DP_matrix::fill()
{
  // check pair envelope is sane
  if (!pair_env(0,0) || !pair_env(xsize-1,ysize-1))
    {
      CLOGERR << "Warning -- pair envelope excludes start or end states\n";
      return;
    }

  // get null states
  const vector<int> null_states = hmm.null_states();
  const vector<int> emit_states = hmm.emit_states();

  // initialise top-left cell
  for_const_contents (vector<int>, emit_states, s)
    cell[*s](0,0) = -InfinityScore;
  for_const_contents (vector<int>, null_states, s)
    cell[*s] (0,0) = calc_score (*s, 0, 0, hmm.start[*s]);
  
  if (xsize == 1 && ysize == 1)
    {
      final_score = hmm.start_to_end();
      for_const_contents (vector<int>, null_states, s)
	ScorePSumAcc (final_score, ScorePMul (cell[*s](0,0), hmm.end[*s]));
    }
  else
    {
      // do top and left boundary conditions
      //
      for (int x = 1; x < xsize; x++)
	if (pair_env(x,0))
	  for (int s = 0; s < states; s++)
	    switch (hmm.state_type[s])
	      {
	      case Pair_HMM_scores::Null:
	      case Pair_HMM_scores::EmitX:
		cell[s] (x,0) = calc_score (s, x, 0, x==1 ? hmm.start[s] : -InfinityScore);
		break;
	      case Pair_HMM_scores::EmitY:
	      case Pair_HMM_scores::EmitXY:
		cell[s] (x,0) = -InfinityScore;
		break;
	      default:
		CLOG(6) << "Unknown state type!\n";
		cell[s] (x,0) = -InfinityScore;
		break;
	      }

      for (int y = 1; y < ysize; y++)
	if (pair_env(0,y))
	  for (int s = 0; s < states; s++)
	    switch (hmm.state_type[s])
	      {
	      case Pair_HMM_scores::Null:
	      case Pair_HMM_scores::EmitY:
		cell[s] (0,y) = calc_score (s, 0, y, y==1 ? hmm.start[s] : -InfinityScore);
		break;
	      case Pair_HMM_scores::EmitX:
	      case Pair_HMM_scores::EmitXY:
		cell[s] (0,y) = -InfinityScore;
		break;
	      default:
		CLOG(6) << "Unknown state type!\n";
		cell[s] (0,y) = -InfinityScore;
		break;
	      }

      // fill the matrix
      //
      for (int x = 1; x < xsize; x++)
	for (int y = 1; y < ysize; y++)
	  if (pair_env(x,y))
	    {
	      for_const_contents (vector<int>, emit_states, s)
		cell[*s] (x,y) = calc_score (*s, x, y, x==1 && y==1 && hmm.state_type[*s] == Pair_HMM_scores::EmitXY ? hmm.start[*s] : -InfinityScore);
	      for_const_contents (vector<int>, null_states, s)
		cell[*s] (x,y) = calc_score (*s, x, y, -InfinityScore);
	    }

      // get the final score
      //
      for (int s = 0; s < states; s++)
	ScorePSumAcc (final_score, ScorePMul (cell[s] (xsize-1, ysize-1), hmm.end[s]));
    }

  if (CTAGGING(-1,FORWARD_HMM)) { CL << "Forward HMM:\n"; hmm.show(CL); }
  if (CTAGGING(-1,DP_MATRIX FORWARD_MATRIX)) { CL << "Forward matrix:\n"; show(CL); }
  if (CTAGGING(-2,DP_MATRIX_SPARSE FORWARD_MATRIX_SPARSE)) { CL << "Forward matrix (sparse):\n"; show_sparse(CL); }
  CTAG(4,DP_SCORE FORWARD_SCORE) << "Forward log-likelihood is " << Score2Bits(final_score) << " bits\n";
}


vector<int> Pair_forward_DP_matrix::sample_state_path() const
{
  vector<int> state_path;

  if (xsize > 1 || ysize > 1)
    {
      vector<int> sc (hmm.states());
      int x = xsize - 1;
      int y = ysize - 1;
      for (int s = 0; s < hmm.states(); s++)
	sc[s] = ScorePMul (cell[s](x,y), hmm.end[s]);
      int current_state = Rnd::choose (Score2ProbVecNorm (sc));
      
      while (1)
	{
	  state_path.push_back (current_state);
	  x -= hmm.dx (current_state);
	  y -= hmm.dy (current_state);
	  if (x == 0 && y == 0) break;
	  
	  for (int s = 0; s < hmm.states(); s++)
	    sc[s] = ScorePMul ((x < hmm.dx(s) || y < hmm.dy(s)) ? -InfinityScore : cell[s](x,y), hmm.transition(s,current_state));
	  vector<double> weight = Score2ProbVecNorm (sc);
	  int dx, dy;
	  bool reject_state;
	  do
	    {
	      current_state = Rnd::choose(weight);
	      dx = hmm.dx(current_state);
	      dy = hmm.dy(current_state);
	      reject_state = x < dx || y < dy ? TRUE : !pair_env(x-dx,y-dy);
	      if (reject_state)
		weight[current_state] = 0.;
	    }
	  while (reject_state);
	}
      
      reverse (state_path.begin(), state_path.end());
    }
  if (CTAGGING(3,TRACEBACK)) CL << "Sampled pair-HMM path posterior log-probability = " << Score2Bits (hmm.path_score (state_path, xseq, yseq) - final_score) << " bits\n";
  return state_path;
}


Pair_forward_backward_DP_matrix::Pair_forward_backward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq) :
  Pair_DP_matrix_base (hmm, xseq, yseq),
  forward_matrix (hmm, xseq, yseq),
  counts (hmm),
  forward_cell (forward_matrix.cell),
  forward_score (forward_matrix.final_score)
{
  fill();
}

Pair_forward_backward_DP_matrix::Pair_forward_backward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const array2d<int>& pair_env) :
  Pair_DP_matrix_base (hmm, xseq, yseq, pair_env),
  forward_matrix (hmm, xseq, yseq, pair_env),
  counts (hmm),
  forward_cell (forward_matrix.cell),
  forward_score (forward_matrix.final_score)
{
  fill();
}

void Pair_forward_backward_DP_matrix::fill()
{
  // check pair envelope is sane
  if (!pair_env(0,0) || !pair_env(xsize-1,ysize-1))
    {
      CLOGERR << "Warning -- pair envelope excludes start or end states\n";
      return;
    }

  // make const references for null_states & all types of emit_states
  vector<vector<int> > states_by_type (4);
  vector<int>& null_states (states_by_type[0]);
  null_states = hmm.null_states();
  reverse (null_states.begin(), null_states.end());  // we want the null states sorted in *reverse* topological order
  for (int t = 1; t < 4; ++t)
    for (int s = 0; s < hmm.states(); ++s)
      if (hmm.state_type[s] == t)
	states_by_type[t].push_back (s);
  
  // set log-likelihood
  counts.log_likelihood = Score2Nats(forward_score);

  if (xsize == 1 && ysize == 1)
    {
      counts.start_to_end() = Score2Prob (ScorePMul (hmm.start_to_end(), -forward_score));
      for_const_contents (vector<int>, null_states, s)
	{
	  const Score cell_sc = back_calc_score (*s, 0, 0);
	  cell[*s] (0,0) = cell_sc;
	  counts.start[*s] += Score2Prob (ScorePMul3 (cell_sc, hmm.start[*s], -forward_score));
	}
    }
  else
    {
      // fill the matrix
      //
      for (int x = xsize-1; x >= 0; x--)
	for (int y = ysize-1; y >= 0; y--)
	  if (pair_env(x,y))
	    for (int t = 0; t < 4; ++t)   // NB state type Null is filled first
	      {
		const int dx = hmm.type_dx ((Pair_HMM_scores::State_type) t);
		const int dy = hmm.type_dy ((Pair_HMM_scores::State_type) t);
		if (x >= dx && y >= dy)
		  if (pair_env(x-dx,y-dy))
		    for_const_contents (vector<int>, states_by_type[t], s)
		      cell[*s] (x,y) = back_calc_score (*s, x, y);
	      }
      
      // do start transitions
      //
      for (int s = 0; s < states; s++)
	{
	  const int x = hmm.dx(s);
	  const int y = hmm.dy(s);
	  if (x < xsize && y < ysize ? pair_env(x,y) : FALSE)
	    {
	      int state_start_sc = ScorePMul (cell[s] (x, y), hmm.start[s]);
	      counts.start[s] += Score2Prob (ScorePMul (state_start_sc, -forward_score));
	      ScorePSumAcc (final_score, state_start_sc);
	    }
	}

      // warn if over 1% difference in forward & backward scores
      if (abs (final_score - forward_score) > (int) abs (MAX_PERMISSIBLE_FORWARD_BACKWARD_DIFFERENCE_ERROR * (double) forward_score))
	CLOGERR << "Warning: forward score = " << forward_score << ", backward score = " << final_score << "\n";

      // dump matrices to log
      if (CTAGGING(-2,BACKWARD_HMM)) { CL << "Backward HMM:\n"; hmm.show(CL); }
      if (CTAGGING(-1,DP_MATRIX BACKWARD_MATRIX)) { CL << "Backward matrix:\n"; show(CL); }
      if (CTAGGING(-2,DP_MATRIX_SPARSE BACKWARD_MATRIX_SPARSE)) { CL << "Backward matrix (sparse):\n"; show_sparse(CL); }
      CTAG(3,DP_SCORE BACKWARD_SCORE) << "Backward log-likelihood is " << Score2Bits(final_score) << " bits\n";
    }

  if (CTAGGING(4,DP_COUNTS)) { CL << "Pair-HMM counts from forward-backward matrix:\n"; counts.show(CL); }
}

Score Pair_forward_backward_DP_matrix::post_state_score (int state, int x, int y) const
{
  const Score src_sc = forward_cell [state] (x, y);
  Score sc = -InfinityScore;
  for (int dest_st = 0; dest_st < states; dest_st++)
    {
      const int dest_x = x + hmm.dx(dest_st);
      const int dest_y = y + hmm.dy(dest_st);
      if (dest_x < xsize && dest_y < ysize)
	{
	  const Score dest_sc = cell[dest_st] (dest_x, dest_y);
	  const Score trans_sc = hmm.transition (state, dest_st);
	  ScorePSumAcc (sc, ScorePMul (trans_sc, dest_sc));
	}
    }
  if (x == xsize - 1 && y == ysize - 1)
    ScorePSumAcc (sc, hmm.end[state]);
  ScorePMulAcc (sc, ScorePMul (src_sc, -forward_score));
  return sc;
}

void Pair_forward_backward_DP_matrix::dump_forward_cell (ostream& o) const
{
  for (int s = 0; s < states; ++s)
    {
      o << "State " << s << ":\n";
      o << forward_cell[s];
    }
}
