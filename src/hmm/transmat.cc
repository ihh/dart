#include <algorithm>
#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "hmm/singlehmm.h"
#include "hmm/singledp.h"
#include "util/rnd.h"

vector<int> Transition_methods::make_global_path (const vector<int>& local_path)
{
  vector<int> global_path;
  global_path.reserve (local_path.size() + 2);
  global_path.push_back (Start);
  global_path.insert (global_path.end(), local_path.begin(), local_path.end());
  global_path.push_back (End);
  return global_path;
}

vector<int> Transition_methods::make_local_path (const vector<int>& global_path)
{
  if (global_path.size() < 2 || global_path.front() != Start || global_path.back() != End) THROWEXPR ("Not a global path");
  return vector<int> (global_path.begin() + 1, global_path.end() - 1);
}

vector<int> Transition_methods::sample_local_path (const Transition_scores& ts, int start_state, const set<int>& end_states)
{
  vector<int> state_path;
  int s = start_state;
  vector<Score> t (ts.tm_states() + 1);
  while (1)
    {
      state_path.push_back (s);
      if (s == End || end_states.find(s) != end_states.end()) break;
      for (int d = 0; d < ts.tm_states(); d++) t[d] = ts.transition (s, d);
      t.back() = ts.transition (s, End);
      s = Rnd::choose (Score2ProbVecNorm (t));
      if (s == (int) t.size() - 1) s = End;
    }
  return state_path;
}

vector<int> Transition_methods::sample_global_path (const Transition_scores& ts)
{
  set<int> dummy_end_states;
  return sample_local_path (ts, Start, dummy_end_states);
}

vector<int> Transition_methods::consensus_local_path (const Transition_scores& ts, int start_state, const set<int>& end_states)
{
  vector<int> state_path;
  int s = start_state;
  vector<int> t (ts.tm_states() + 1);
  vector<int> pos (ts.tm_states(), (int) -1);
  while (1)
    {
      state_path.push_back (s);
      if (s == End || end_states.find(s) != end_states.end()) break;
      if (s >= 0) pos[s] = state_path.size() - 1;
      for (int d = 0; d < ts.tm_states(); d++)
	{
	  const Score tr = ts.transition (s, d);
	  if (pos[d] >= 0 && tr > -InfinityScore)
	    t[d] = -InfinityScore + state_path.size() - pos[d];
	  else
	    t[d] = tr;
	}
      t.back() = ts.transition (s, End);
      s = max_element (t.begin(), t.end()) - t.begin();
      if (s == (int) t.size() - 1) s = End;
    }
  if (CTAGGING(2,CONSENSUS_PATH)) { CL << "Consensus state path is (" << state_path << ")\n"; }
  return state_path;
}

vector<int> Transition_methods::consensus_global_path (const Transition_scores& ts)
{
  set<int> dummy_end_states;
  return consensus_local_path (ts, Start, dummy_end_states);
}

Score Transition_methods::path_transition_score (const Transition_scores& ts, const vector<int>& state_path)
{
  Score sc = 0;
  for (int i = 1; i < (int) state_path.size(); ++i) ScorePMulAcc (sc, ts.transition (state_path[i-1], state_path[i]));
  return sc;
}

void Transition_methods::add_transition_counts_from_path (Transition_counts& tc, const vector<int>& state_path)
{
  for (int i = 1; i < (int) state_path.size(); i++)
    tc.transition (state_path[i-1], state_path[i]) += 1.0;
}

Concrete_transition_probs Transition_methods::eliminate (const Transition_probs& tp,
							 const vector<int>& null_states,
							 bool zero2null,
							 Concrete_transition_probs* s_mx)
{
  // print log message
  if (CTAGGING(0,TRANSMAT_ELIMINATE TRANSMAT_PRE_ELIMINATE))
    {
      CL << "[prior to eliminating states " << null_states << "] ";
      tp.show_transitions (CL);
    }

  // get size
  const int states = tp.tm_states();
  const int size = tp.tm_index (states);  // number of matrix entries = number of states + 2 (Start,End)

  // make vector of non-null (emit) states, including Start & End (unless specified null!)
  vector<bool> isNull (size, FALSE);
  for_const_contents (vector<int>, null_states, ns)
    isNull[tp.tm_index (*ns)] = TRUE;
  vector<int> emit_states;
  for (int s = 0; s < size; ++s)
    if (!isNull[s])
      emit_states.push_back (tp.tm_inverse (s));

  // make a,b,c,d so M=a+b+c+d, a=emit->emit, b=emit->null, c=null->emit, d=null->null
  // see (Holmes, ISMB 2003)
  Matrix a (size, size);
  Matrix b (size, size);
  Matrix c (size, size);
  Matrix d (size, size);
  a = 0.;
  b = 0.;
  c = 0.;
  d = 0.;
  for_const_contents (vector<int>, emit_states, i)
    {
      for_const_contents (vector<int>, emit_states, j)
	if (*i != End && *j != Start)
	  a(tp.tm_index(*i)+1,tp.tm_index(*j)+1) = tp.transition(*i,*j);
      for_const_contents (vector<int>, null_states, j)
	if (*i != End && *j != Start)
	  b(tp.tm_index(*i)+1,tp.tm_index(*j)+1) = tp.transition(*i,*j);
    }
  for_const_contents (vector<int>, null_states, i)
    {
      for_const_contents (vector<int>, emit_states, j)
	if (*i != End && *j != Start)
	  c(tp.tm_index(*i)+1,tp.tm_index(*j)+1) = tp.transition(*i,*j);
      for_const_contents (vector<int>, null_states, j)
	if (*i != End && *j != Start)
	  d(tp.tm_index(*i)+1,tp.tm_index(*j)+1) = tp.transition(*i,*j);
    }

  // make identity matrix I
  DiagonalMatrix I (size);
  I = 1.;
  // find e = (I - d)  and f = e^{-1}
  const Matrix e = I - d;
  const Matrix f = e.i();
  // find q matrix:
  // q = a + b * f * c
  //   = a + sum_{n=0}^infty (b * d^n * c)
  const Matrix bf = b * f;
  const Matrix q = a + bf * c;

  // log
  if (CTAGGING(-3,TRANSMAT_ELIMINATE TRANSMAT_ELIMINATE_INVERT))
    {
      const Matrix fc = f * c;
      CL << "emit_states: " << emit_states << "\nnull_states: " << null_states
	 << "\na:\n" << a << "b:\n" << b << "c:\n" << c << "d:\n" << d
	 << "e=(I-d):\n" << e << "f=(I-d)^{-1}:\n" << f
	 << "bf:\n" <<  bf<< "fc:\n" << fc
	 << "q=a+bfc:\n" << q;
    }

  // create new Transition_probs & fill it
  // emit->emit transitions
  Concrete_transition_probs new_tp (states, 0.);
  for_const_contents (vector<int>, emit_states, i)
    if (*i != End)
      for_const_contents (vector<int>, emit_states, j)
	if (*j != Start)
	  new_tp.transition(*i,*j) = q(tp.tm_index(*i)+1,tp.tm_index(*j)+1);

  // The transitions to null states in new_tp are zero by default.
  // If zero2null is false, then set the emit->null transitions to the sum-over-null-path probabilities of
  //  sum_n (b * d^n) = b * (I - d)^{-1} = bf
  // and the null->null transitions to those of
  //  sum_n (d^n) = (I - d)^{-1} = f
  // (Matrix bf is equivalent to matrix r in Holmes, 2003)
  // These transitions are useful for sampling eliminated paths during stochastic traceback.
  // Returning these probabilities in new_tp is a bit hacky,
  // but "safe" vis-a-vis null state elimination,
  // since there are no complete paths through the null states.
  if (!zero2null)
    {
      // emit->null
      for_const_contents (vector<int>, emit_states, i)
	if (*i != End)
	  for_const_contents (vector<int>, null_states, j)
	    if (*j != Start)
	      new_tp.transition(*i,*j) = bf(tp.tm_index(*i)+1,tp.tm_index(*j)+1);

      // null->null
      for_const_contents (vector<int>, null_states, i)
	if (*i != End)
	  for_const_contents (vector<int>, null_states, j)
	    if (*j != Start)
	      new_tp.transition(*i,*j) = f(tp.tm_index(*i)+1,tp.tm_index(*j)+1);
    }

  // null->emit transitions
  // For other applications (such as counting eliminated transitions,
  // or retaining some null states as bifurcation targets in an SCFG),
  // it's useful to have the null->emit sum-over paths probabilities
  //  sum_n (d^n * c) = (I - d)^{-1} * c = fc
  // If a pointer was specified, store the transitions in s_mx (the s matrix of Holmes 2003).
  if (s_mx)
    {
      const Matrix fc = f * c;
      *s_mx = Concrete_transition_probs (states, 0.);
      for_const_contents (vector<int>, null_states, i)
	if (*i != End)
	  for_const_contents (vector<int>, emit_states, j)
	    if (*j != Start)
	      s_mx->transition(*i,*j) = fc(tp.tm_index(*i)+1,tp.tm_index(*j)+1);
    }

  // print log message
  if (CTAGGING(0,TRANSMAT_ELIMINATE TRANSMAT_POST_ELIMINATE))
    {
      CL << "[after eliminating states " << null_states << "] ";
      new_tp.show_transitions (CL);
    }

  // return
  return new_tp;
}

vector<int> Transition_methods::sample_eliminated (const Transition_probs& tp_orig,
						   const Transition_probs& tp_elim,
						   const vector<int>& null_states,
						   const vector<int>& elim_path,
						   bool choose_ML_path)
{
  // get incoming_states
  const Concrete_transition_scores ts_orig = prob2score (tp_orig);
  const vector<vector<int> > incoming = selected_incoming_states (ts_orig, null_states);

  // get incoming from end state
  vector<int> incoming_from_end;
  for_const_contents (vector<int>, null_states, s)
    if (*s != End)
      if (ts_orig.transition (*s, End) > -InfinityScore)
	incoming_from_end.push_back (*s);

  // sample path
  vector<int> path;
  if (elim_path.size())
    {
      int current = elim_path.back();
      path.push_back (current);
      for (int i = elim_path.size() - 2; i >= 0; --i)
	{
	  const int previous = elim_path[i];
	  bool finished = false;
	  do
	    {
	      vector<int> candidate_state (current == End ? incoming_from_end : incoming[current]);
	      vector<Prob> candidate_prob;
	      for_const_contents (vector<int>, candidate_state, null_state)
		candidate_prob.push_back (tp_elim.transition (previous, *null_state) * tp_orig.transition (*null_state, current));

	      // incoming[] and incoming_from_end only include null states, whereas 'previous' is an emit state: we need to account for the possibility of a direct transition
	      candidate_state.push_back (previous);
	      candidate_prob.push_back (tp_orig.transition (previous, current));

	      const int candidate_state_index
		= choose_ML_path

		? (max_element (candidate_prob.begin(),
				candidate_prob.end())
		   - candidate_prob.begin())

		: Rnd::choose (candidate_prob);

	      current = candidate_state[candidate_state_index];
	      path.push_back (current);
	      if (current == previous)
		{
		  const bool prev_is_null = find (null_states.begin(), null_states.end(), current) != null_states.end();
		  finished = prev_is_null ? Rnd::decide (1 / tp_elim.transition (previous, previous)) : true;
		}
	    }
	  while (!finished);
	}

      reverse (path.begin(), path.end());
    }
  return path;
}

Concrete_transition_counts Transition_methods::count_eliminated (const Transition_probs& tp_orig,
								 const Transition_probs& tp_elim,
								 const Transition_probs& s_mx,
								 const vector<int>& null_states,
								 const vector<int>& emit_states,
								 const Transition_counts& elim_counts)
{
  // In notation of (Holmes, 2003):
  // tp_orig = a + b + c + d
  // tp_elim = q + r + (1-d)^{-1} = q + (1+b)(1-d)^{-1}
  //    s_mx = s

  // dest_states = emit_states \cup { End }
  vector<int> dest_states (emit_states);
  dest_states.push_back (Grammar_state_enum::End);

  // src_states = all states except End
  vector<int> src_states (1, Grammar_state_enum::Start);
  for (int s = 0; s < elim_counts.tm_states(); ++s)
    src_states.push_back (s);

  // now do it
  Concrete_transition_counts counts (elim_counts.tm_states(), 0.);
  for_const_contents (vector<int>, src_states, i)
    for_const_contents (vector<int>, dest_states, j)
    {
      const Prob c_ij = elim_counts.transition(*i,*j);
      if (c_ij > 0.)
	{
	  const Prob q_ij = tp_elim.transition(*i,*j);
	  const Prob pp = c_ij / q_ij;

	  counts.transition(*i,*j) += pp * tp_orig.transition(*i,*j);

	  for_const_contents (vector<int>, null_states, k)
	    {
	      counts.transition(*i,*k) += pp * tp_orig.transition(*i,*k) * s_mx.transition(*k,*j);
	      counts.transition(*k,*j) += pp * tp_elim.transition(*i,*k) * tp_orig.transition(*k,*j);

	      for_const_contents (vector<int>, null_states, l)
		counts.transition(*k,*l) += pp * tp_elim.transition(*i,*k) * tp_orig.transition(*k,*l) * s_mx.transition(*l,*j);
	    }
	}
    }

  return counts;
}

Concrete_transition_probs Transition_methods::score2prob (const Transition_scores& ts)
{
  const int size = ts.tm_states();
  Concrete_transition_probs tp (size, 0.);
  for (int i = 0; i < size; ++i)
    {
      tp.transition(Start,i) = Score2Prob (ts.transition(Start,i));
      tp.transition(i,End) = Score2Prob (ts.transition(i,End));
      for (int j = 0; j < size; ++j)
	tp.transition(i,j) = Score2Prob (ts.transition(i,j));
    }
  tp.transition(Start,End) = Score2Prob (ts.transition(Start,End));
  return tp;
}

Concrete_transition_scores Transition_methods::prob2score (const Transition_probs& tp)
{
  const int size = tp.tm_states();
  Concrete_transition_scores ts (size, -InfinityScore);
  for (int i = 0; i < size; ++i)
    {
      ts.transition(Start,i) = Prob2Score (tp.transition(Start,i));
      ts.transition(i,End) = Prob2Score (tp.transition(i,End));
      for (int j = 0; j < size; ++j)
	ts.transition(i,j) = Prob2Score (tp.transition(i,j));
    }
  ts.transition(Start,End) = Prob2Score (tp.transition(Start,End));
  return ts;
}
