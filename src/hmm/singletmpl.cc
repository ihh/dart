#include <algorithm>
#include "hmm/singlehmm.h"
#include "hmm/singledp.h"
#include "util/rnd.h"
#include "util/vector_output.h"

Single_HMM_scores::Single_HMM_scores() : Single_meta_HMM<Score>() { }

Single_HMM_scores::Single_HMM_scores (int states, const Alphabet& alphabet)
  : Single_meta_HMM<Score> (states, alphabet, -InfinityScore)
{ }

Single_HMM_scores& Single_HMM_scores::operator= (const Single_HMM_scores& s)
{
  Single_HMM<int>::operator= (s);
  metascore_idx = s.metascore_idx;
  return *this;
}

vector<int> Single_HMM_scores::sample_state_path() const
{
  const vector<int> global_path = Transition_methods::sample_global_path (*this);
  return Transition_methods::make_local_path (global_path);
}

void Single_HMM_scores::sample_sequence (const vector<int>& state_path, Digitized_biosequence& dsq) const
{
  vector<vector<Prob> > emit_prob;
  for (int s = 0; s < states(); ++s) emit_prob.push_back (Score2ProbVecNorm (emit[s]));
  dsq.clear();
  dsq.reserve (state_path.size());  // hopefully there won't be so many null states that this is mega-wasteful
  for (int i = 0; i < (int) state_path.size(); ++i)
    dsq.push_back (Rnd::choose (emit_prob [state_path[i]]));
}

vector<int> Single_HMM_scores::consensus_state_path() const
{
  const vector<int> global_path = Transition_methods::consensus_global_path (*this);
  return Transition_methods::make_local_path (global_path);
}

void Single_HMM_scores::consensus_sequence (const vector<int>& state_path, Digitized_biosequence& dsq) const
{
  dsq.clear();
  dsq.reserve (state_path.size());  // hopefully there won't be so many null states that this is mega-wasteful
  for (int i = 0; i < (int) state_path.size(); ++i)
    dsq.push_back (max_element (emit[state_path[i]].begin(), emit[state_path[i]].end()) - emit[state_path[i]].begin());
}

void Single_HMM_scores::make_profile (const vector<int>& state_path, Score_profile& prof_sc) const
{
  prof_sc.clear();
  prof_sc.reserve (state_path.size());  // hopefully there won't be so many null states that this is mega-wasteful
  for_const_contents (vector<int>, state_path, state)
    if (state_type[*state] == Emit)
      {
	Symbol_score_map sc;
	for (int sym = 0; sym < alphabet().size(); ++sym) sc[sym] = emit[*state][sym];
	prof_sc.push_back (sc);
      }
}

vector<int> Single_HMM_scores::null_states_unsorted() const
{
  vector<int> null;
  null.reserve (states());
  for (int s = 0; s < states(); ++s) if (state_type[s] == Null) null.push_back(s);
  return null;
}

vector<int> Single_HMM_scores::emit_states() const
{
  vector<int> result;
  for (int s = 0; s < states(); ++s) if (state_type[s] == Emit) result.push_back(s);
  return result;
}

vector<vector<int> > Single_HMM_scores::incoming_states() const
{
  return Transition_methods::incoming_states (*this);
}

vector<vector<int> > Single_HMM_scores::selected_outgoing_states (const vector<int>& selection) const
{
  return Transition_methods::selected_outgoing_states (*this, selection);
}

vector<int> Single_HMM_scores::null_states() const
{
  const vector<int> unsorted = null_states_unsorted();
  return Transition_methods::topological_sort (*this, unsorted);
}

Score Single_HMM_scores::path_transition_score (const vector<int>& state_path) const
{
  const vector<int> global_path = Transition_methods::make_global_path (state_path);
  return Transition_methods::path_transition_score (*this, global_path);
}

Score Single_HMM_scores::path_emit_score (const vector<int>& state_path, const Meta_profile& mp) const
{
  Single_emit_calculator calculator (*this, mp);
  int seqpos = 0;
  Score sc = 0;
  for_const_contents (vector<int>, state_path, s)
    if (state_type[*s] == Emit)
      ScorePMulAcc (sc, calculator.calc_emit_score (*s, seqpos++));
  return sc;
}

Score Single_HMM_scores::path_score (const vector<int>& state_path, const Meta_profile& mp) const
{
  return ScorePMul (path_transition_score (state_path), path_emit_score (state_path, mp));
}

Score Single_HMM_scores::path_emit_score_dsq (const vector<int>& state_path, const Named_profile& np) const
{
  Single_fast_emit_calculator calculator (*this, np);
  int seqpos = 0;
  Score sc = 0;
  for_const_contents (vector<int>, state_path, s)
    if (state_type[*s] == Emit)
      ScorePMulAcc (sc, calculator.calc_emit_score (*s, seqpos++));
  return sc;
}

Score Single_HMM_scores::path_score_dsq (const vector<int>& state_path, const Named_profile& np) const
{
  return ScorePMul (path_transition_score (state_path), path_emit_score_dsq (state_path, np));
}

void Single_HMM_scores::scale_all_scores (double beta)
{
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest) = (int) (beta * (double) transition(src,dest));
      for_contents (vector<int>, emit[src], sc)
	*sc = (int) (beta * (double) *sc);
      start[src] = (int) (beta * (double) start[src]);
      end[src]   = (int) (beta * (double) end[src]);
    }
  ScaleScore (start_to_end(), beta);
}

Single_HMM_derivatives::Single_HMM_derivatives (int states, const Alphabet& alphabet) : Single_HMM<double>(states,alphabet), wrt ("derivatives") { }
Single_HMM_derivatives::Single_HMM_derivatives (int states, const char* x, const Alphabet& alphabet) : Single_HMM<double>(states,alphabet), wrt()
{
  wrt << "derivatives w.r.t " << x;
}

Single_HMM_mask::Single_HMM_mask (const Single_HMM_scores& hmm, bool init_flag) :
  Single_HMM<int> (hmm.states(), hmm.alphabet(), init_flag)
{ }

void Single_HMM_mask::set_incoming_transitions (int state)
{ start[state] = 1; for (int s = 0; s < states(); ++s) transition(state,s) = 1; }

void Single_HMM_mask::clear_incoming_transitions (int state)
{ start[state] = 0; for (int s = 0; s < states(); ++s) transition(state,s) = 0; }

void Single_HMM_mask::set_outgoing_transitions (int state)
{ end[state] = 1; for (int s = 0; s < states(); ++s) transition(s,state) = 1; }

void Single_HMM_mask::clear_outgoing_transitions (int state)
{ end[state] = 0; for (int s = 0; s < states(); ++s) transition(s,state) = 0; }

void Single_HMM_mask::set_emit (int state) { for_contents (vector<int>, emit[state], i) *i = 1; }
void Single_HMM_mask::clear_emit (int state) { for_contents (vector<int>, emit[state], i) *i = 0; }
void Single_HMM_mask::set_all_start() { for (int s = 0; s < states(); ++s) start[s] = 1; start_to_end() = 1; }
void Single_HMM_mask::clear_all_start() { for (int s = 0; s < states(); ++s) start[s] = 0; start_to_end() = 0; }
void Single_HMM_mask::set_all_end() { for (int s = 0; s < states(); ++s) end[s] = 1; start_to_end() = 1; }
void Single_HMM_mask::clear_all_end() { for (int s = 0; s < states(); ++s) end[s] = 0; start_to_end() = 0; }
void Single_HMM_mask::set_all_transitions() { set_all_start(); set_all_end(); for (int s = 0; s < states(); ++s) set_outgoing_transitions(s); }
void Single_HMM_mask::clear_all_transitions() { clear_all_start(); clear_all_end(); for (int s = 0; s < states(); ++s) clear_outgoing_transitions(s); }
void Single_HMM_mask::set_all_emit() { for (int s = 0; s < states(); ++s) set_emit(s); }
void Single_HMM_mask::clear_all_emit() { for (int s = 0; s < states(); ++s) clear_emit(s); }


Single_HMM_counts::Single_HMM_counts (const Single_HMM_scores& hmm) :
  Single_HMM<Prob> (hmm),
  log_likelihood (0.0)
{
  reset (0.0);
}

void Single_HMM_counts::add_counts (const Single_HMM_counts& counts)
{
  if (!same_dimensions (counts)) THROW Standard_exception ("While adding pair-HMM counts: state spaces don't match");
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest) += counts.transition(src,dest);
      for (int i = 0; i < (int) emit[src].size(); ++i)
	emit[src][i] += counts.emit[src][i];
      start[src] += counts.start[src];
      end[src] += counts.end[src];
    }
  start_to_end() += counts.start_to_end();
  log_likelihood += counts.log_likelihood;
}

void Single_HMM_counts::subtract_counts (const Single_HMM_counts& counts)
{
  if (!same_dimensions (counts)) THROW Standard_exception ("While subtracting pair-HMM counts: state spaces don't match");
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest) -= counts.transition(src,dest);
      for (int i = 0; i < (int) emit[src].size(); ++i)
	emit[src][i] -= counts.emit[src][i];
      start[src] -= counts.start[src];
      end[src] -= counts.end[src];
    }
  start_to_end() -= counts.start_to_end();
  log_likelihood -= counts.log_likelihood;
}

void Single_HMM_counts::update_HMM_scores (Single_HMM_scores& hmm, const Single_HMM_mask& mask, bool sample, double kT) const
{
  if (!same_dimensions (hmm)) THROW Standard_exception ("While updating single-HMM scores: state spaces of counts & scores structures don't match");
  if (!same_dimensions (mask)) THROW Standard_exception ("While updating single-HMM scores: state spaces of counts & update mask structures don't match");
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
      
      score_pr_product = -InfinityScore;
      counts.clear();
      for (int i = 0; i < (int) emit[src].size(); ++i)
	if (mask.emit[src][i])
	  {
	    ScorePSumAcc (score_pr_product, hmm.emit[src][i]);
	    counts.push_back (emit[src][i]);
	  }
      prob = sample ? Rnd::sample_dirichlet(counts) : Rnd::dirichlet_max(counts);
      k = 0;
      for (int i = 0; i < (int) emit[src].size(); ++i)
	if (mask.emit[src][i])
	  hmm.emit[src][i] = Prob2Score (prob[k++]) + score_pr_product;
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

double Single_HMM_counts::dloglike_dx (const Single_HMM_scores& hmm, const Single_HMM_derivatives& deriv) const
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
      for (int i = 0; i < (int) emit[s].size(); i++)
	sum += emit[s][i] * deriv.emit[s][i] * Score2Prob(-hmm.emit[s][i]);
    }
  sum += start_to_end() * deriv.start_to_end() * Score2Prob(-hmm.start_to_end());
  return sum;
}

void Single_HMM_counts::add_counts_from_state_path (const Single_HMM_scores& hmm, const Named_profile& np, const vector<int>& state_path, vector<Metaprob>& metacounts)
{
  const vector<int> global_path = Transition_methods::make_global_path (state_path);
  Transition_methods::add_transition_counts_from_path (*this, global_path);

  log_likelihood += Score2Nats (hmm.path_score_dsq (state_path, np));  // this involves making the global path again. oh well

  Single_fast_emit_calculator calculator (hmm, np);
  int pos = 0;
  for_const_contents (vector<int>, state_path, s)
    if (state_type[*s] != Null)
      calculator.accum_emit_score (0, *this, metacounts, *s, pos++);

  if (CLOGGING(1)) { CL << "Single-HMM counts from state path:\n"; show(CL); }
}

void Single_HMM_counts::add_counts_from_state_path (const Single_HMM_scores& hmm, const Named_profile& np, const vector<int>& state_path)
{
  vector<Metaprob> dummy_metacounts (hmm.max_metascore_idx() + 1, Metaprob (np.dsq.size(), (Prob) 0));
  add_counts_from_state_path (hmm, np, state_path, dummy_metacounts);
}

void Single_HMM_counts::write (ostream& out) const
{
  Single_HMM<Prob>::write (out);
  out << log_likelihood << "\n";
}

void Single_HMM_counts::read (istream& in)
{
  Single_HMM<Prob>::read (in);
  in >> log_likelihood;
}
