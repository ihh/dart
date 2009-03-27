#include "empath/multihit.h"

Sequence_hitter::Sequence_hitter (Local_trainer& trainer)
  : trainer (trainer),
    hmm (trainer.local_hmm),
    prior (trainer.prior),
    pscore (*prior.pscores),
    dp_factory (trainer.dp_factory),
    pcount (pscore)
{ }

void Sequence_hitter::seed (const Named_profile& np, int pos)
{
  trainer.local_seed (np, pos);
}

void Sequence_hitter::initialise()
{
  hmm_scores = hmm.eval_hmm_scores (pscore);
  pcount.clear();
}

Loge Sequence_hitter::calc_loglike (const Named_profile& np)
{
  Single_forward_interface* fwd = dp_factory.new_forward (hmm_scores, np);
  const Loge loglike = Score2Nats (fwd->get_forward_score());
  delete fwd;
  return loglike;
}

void Sequence_hitter::send_update (ostream& out, const Named_profile& np, Prob weight)
{
  Single_forward_backward_interface* fb = dp_factory.new_forward_backward (hmm_scores, np);
  const Single_HMM_counts seq_hmm_counts = fb->get_expected_counts();
  seq_hmm_counts.write (out);
  out << weight << "\n";
  delete fb;
}

void Sequence_hitter::receive_update (istream& in)
{
  Single_HMM_counts seq_hmm_counts (hmm_scores);
  seq_hmm_counts.read (in);
  Prob weight;
  in >> weight;
  hmm.inc_var_counts (pcount, pscore, seq_hmm_counts, weight);
}

void Sequence_hitter::optimise()
{
  prior.optimise (pcount);
}

vector<double> Sequence_hitter::get_params()
{
  vector<double> params;
  for (int g = 0; g < pscore.groups(); ++g)
    for (int p = 0; p < pscore.group_size(g); ++p)
      params.push_back ((double) pscore.group[g][p]);
  return params;
}

void Sequence_hitter::set_params (const vector<double>& params)
{
  int i = 0;
  for (int g = 0; g < pscore.groups(); ++g)
    for (int p = 0; p < pscore.group_size(g); ++p)
      {
	if (i >= (int) params.size()) THROWEXPR ("Param vector too small");
	pscore.group[g][p] = (Score) params[i++];
      }
  if (i < (int) params.size()) THROWEXPR ("Param vector too big");
}

void Sequence_hitter::display (ostream& out)
{
  trainer.display_vars (out);
}

Multihitter::Multihitter (const Transition_funcs& trans, Hitter& hitter, PScores& pscore, int max_fork, double threshold) :
  Piper (max_fork),
  db_hmm (trans.tm_states() * 2 + 1, Dummy_alphabet),
  hit_meta_idx (0),
  dummy_sym (0),
  threshold (threshold),
  hitter (hitter),
  pscore (pscore),
  prior (pscore)
{
  db_hmm.transition (Start,              start_miss_state()) = 1.0;
  db_hmm.transition (start_miss_state(), start_miss_state()) = 1.0;
  db_hmm.emit[start_miss_state()][dummy_sym] = 1.0;
  for (int s = 0; s < trans.tm_states(); ++s)
    {
      db_hmm.transition (Start,              hit_state(s))  = trans.transition (Start, s);
      db_hmm.transition (start_miss_state(), hit_state(s))  = trans.transition (Start, s);
      db_hmm.transition (hit_state(s),       End)           = trans.transition (s,     End);
      db_hmm.transition (miss_state(s),      End)           = trans.transition (s,     End);
      db_hmm.transition (hit_state(s),       miss_state(s)) = 1.0;
      db_hmm.transition (miss_state(s),      miss_state(s)) = 1.0;
      for (int t = 0; t < trans.tm_states(); ++t)
	{
	  db_hmm.transition (hit_state(s),  hit_state(t)) = trans.transition (s, t);
	  db_hmm.transition (miss_state(s), hit_state(t)) = trans.transition (s, t);
	}
      db_hmm.metascore_idx[hit_state(s)].push_back (hit_meta_idx);
      db_hmm.emit[hit_state(s)][dummy_sym] = 1.0;
      db_hmm.emit[miss_state(s)][dummy_sym] = 1.0;
    }
  db_hmm.transition (Start,              End)  = trans.transition (Start, End);
  db_hmm.transition (start_miss_state(), End)  = trans.transition (Start, End);
}

Named_profile Multihitter::make_hit_profile (const vector<Named_profile*>& db)
{
  const int dummy_sym = 0;
  // initialise Hitter
  hitter.initialise();
  // get score for each sequence, create dummy sequence with forward-likelihood Metascore
  Named_profile dummy_np;
  // allocate space
  dummy_np.dsq = Digitized_biosequence (db.size(), dummy_sym);
  dummy_np.meta_sc.push_back (Metascore (db.size()));
  // decide whether or not to fork
  if (max_fork < 2)
    for (int i = 0; i < (int) db.size(); ++i)
      dummy_np.meta_sc[hit_meta_idx][i] = Nats2Score (hitter.calc_loglike (*db[i]));
  else
    {
      // open some pipes
      open_pipes();
      // fill the profile, forkily
      for (int i = 0; i < (int) db.size(); i += max_fork)
	{
	  const int j_end = min (i + max_fork, (int) db.size());
	  for (int j = i; j < j_end; ++j)
	    if (fork_child(j-i))  // child
	      {
		ostream_cfile pipe_out (write_fd(j-i));
		pipe_out.out << Nats2Score (hitter.calc_loglike (*db[j])) << "\n";;
		pipe_out.close();
		exit(0);
	      }
	  
	  // parent
	  for (int j = i; j < j_end; ++j)
	    {
	      istream_cfile pipe_in (read_fd(j-i));
	      pipe_in.in >> dummy_np.meta_sc[hit_meta_idx][j];
	      wait_for_child(j-i);  // wait for each child to exit
	    }
	}
      // close the pipes
      close_pipes();
    }
  // return
  return dummy_np;
}

Loge Multihitter::multihit_forward (const vector<Named_profile*>& db)
{
  // get dummy sequence
  const Named_profile dummy_np = make_hit_profile (db);
  // do forward on dummy sequence
  Single_HMM_scores db_hmm_scores = db_hmm.eval_hmm_scores (pscore);
  Single_fast_forward_matrix db_fwd (db_hmm_scores, dummy_np);
  return Score2Nats (db_fwd.get_forward_score());
}

Loge Multihitter::multihit_EM (const vector<Named_profile*>& db)
{
  // get dummy sequence
  const Named_profile dummy_np = make_hit_profile (db);
  // do forward-backward on dummy sequence
  Single_HMM_scores db_hmm_scores = db_hmm.eval_hmm_scores (pscore);
  Single_fast_forward_backward_matrix db_fb (db_hmm_scores, dummy_np);
  // get metacounts (i.e. posterior expectation of number of times each seq in db is a hit)
  const Metaprob& post_hit_count = db_fb.get_expected_metacounts() [hit_meta_idx];
  // decide which sequences to look at
  vector<int> strong_hits;
  for (int i = 0; i < (int) db.size(); ++i)
    {
      const Loge seq_log_like = Score2Nats (dummy_np.meta_sc[hit_meta_idx][i]);
      if (CTAGGING(2,MULTIHIT)) CL << "Sequence '" << db[i]->name << "' log-likelihood " << Nats2Bits(seq_log_like) << " bits, expected hit count " << post_hit_count[i] << "\n";
      if (post_hit_count[i] >= threshold / (double) db.size())
	strong_hits.push_back (i);
    }
  // sort strong_hits by increasing posterior hit count, and print the top hits to the logfile
  Schwartzian<double> by_post_hit_count (post_hit_count);
  sort (strong_hits.begin(), strong_hits.end(), by_post_hit_count);
  if (CTAGGING(3,MULTIHIT MULTIHIT_TOP))
    {
      const int n_top = min (db_hmm.states()-1, (int) strong_hits.size());
      CL << "Top " << n_top << " posterior expected hit-counts for sequences:";
      for (int n = 0; n < n_top; ++n)
	{
	  const int i = strong_hits[strong_hits.size() - 1 - n];
	  CL << ' ' << db[i]->name << '=' << post_hit_count[i];
	}
      CL << "\n";
    }
  // decide whether or not to fork
  if (max_fork < 2)
    {
      open_pipes();  // set up a dummy pipe
      istream_cfile pipe_in (read_fd(0));
      ostream_cfile pipe_out (write_fd(0));
      for (int i = 0; i < (int) strong_hits.size(); ++i)
	{
	  hitter.send_update (pipe_out.out, *db[strong_hits[i]], post_hit_count[strong_hits[i]]);
	  hitter.receive_update (pipe_in.in);
	}
      close_pipes();
      }
  else
    {
      // open a bunch of pipes
      open_pipes();
      // add all the chosen sequence counts, weighted by posterior hit probability
      for (int i = 0; i < (int) strong_hits.size(); i += max_fork)
	{
	  const int j_end = min (i + max_fork, (int) strong_hits.size());
	  for (int j = i; j < j_end; ++j)
	    if (fork_child(j-i))  // child
	      {
		ostream_cfile pipe_out (write_fd (j-i));
		hitter.send_update (pipe_out.out, *db[strong_hits[j]], post_hit_count[strong_hits[j]]);
		pipe_out.close();
		exit(0);
	      }
	  
	  // parent
	  for (int j = i; j < j_end; ++j)
	    {
	      istream_cfile pipe_in (read_fd (j-i));
	      hitter.receive_update (pipe_in.in);
	      wait_for_child(j-i);  // wait for each child to exit
	    }
	}
      // close the pipes
      close_pipes();
    }
  // optimise the Hitter
  hitter.optimise();
  // get counts for multi-hit transition functions and give them to our personal Transition_funcs prior
  PCounts pcount (pscore);
  const Single_HMM_counts& trans_counts = db_fb.get_expected_counts();
  db_hmm.inc_var_counts (pcount, pscore, trans_counts, 1.0);
  prior.optimise (pcount);
  // return forward score
  return Score2Nats (db_fb.get_forward_score());
}

void Multihitter::iterate_multihit_EM (const vector<Named_profile*>& db, double min_inc, int forgive)
{
  vector<double> best_params = hitter.get_params();
  Loge best_loglike = -InfinityLoge;
  int dec_count = 0;
  for (int em_iter = 1; 1; ++em_iter)  // EM outer loop
    {
      const vector<double> old_params = hitter.get_params();
      const Loge loglike = multihit_EM (db);  // do one round of EM
      CTAG(6,MULTIHIT_EM) << "EM iteration #" << em_iter << ": score " << Nats2Bits(loglike) << " bits\n";
      if (CTAGGING(5,MULTIHIT_EM)) hitter.display (CL);
      const double inc = (loglike - best_loglike) / abs(best_loglike);
      if (inc > 0)
	{
	  best_loglike = loglike;  // record hiscores
	  best_params = old_params;
	  dec_count = 0;
	}
      if (inc < min_inc)
	if (++dec_count > forgive)
	  {
	    CTAG(6,MULTIHIT_EM) << "Likelihood increase negligible; stopping EM\n";
	    break;
	  }
    }
  hitter.set_params (best_params);
}

Range_transitions::Range_transitions (int range_min, int range_max) : Transition_funcs_concrete (range_max)
{
  const double start_prob = 1.0 / (double) (range_max + 1 - range_min);
  for (int i = 0; i < range_max; ++i)
    {
      if (i <= range_max - range_min)
	transition (Start, i) = start_prob;
      if (i < range_max - 1)
	transition (i, i + 1) = 1;
    }
  if (range_min <= 0)
    transition (Start, End) = start_prob;
  transition (range_max - 1, End) = 1;
}
