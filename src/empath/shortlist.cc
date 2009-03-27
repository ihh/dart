#include <algorithm>
#include "empath/shortlist.h"
#include "util/vector_output.h"
#include "seq/suffix.h"

Shortlist::Shortlist (const FASTA_sequence_database& seq_db, int max_fork)
  : Piper (max_fork),
    seq_db (seq_db),
    index (seq_db.index),
    by_loglike (loglike),
    diagnostic_flag (CTAGGING(1,SHORTLIST_DIAGNOSTIC))
{
  alloc_loglike();
}

Shortlist::~Shortlist()
{ }

void Shortlist::alloc_loglike()
{
  loglike = vector<vector<Loge> > (index.size());
  for (int seq = 0; seq < index.size(); ++seq)
    loglike[seq] = vector<Loge> (index.profile[seq]->size());
}

void Shortlist::free_loglike()
{
  loglike.clear();
}

void Shortlist::seed_everywhere (int motif_len)
{
  seed.clear();
  for (int seq = 0; seq < index.size(); ++seq)
    for (int pos = 0; pos < ((int) index.profile[seq]->size()) - motif_len; ++pos)
      seed.push_back (Seq_pos (seq, pos));
}

void Shortlist::seed_locally_repeated (int motif_len, int word_len, int min_reps)
{
  seed.clear();
  for (int seq = 0; seq < index.size(); ++seq)
    {
      const Digitized_biosequence& dsq = index.profile[seq]->dsq;
      Suffix_tree suff_tree (FALSE);
      suff_tree.add_suffices (dsq, word_len);
      for (int pos = 0; pos < ((int) dsq.size()) - motif_len; ++pos)
	for (int word_start = pos; word_start < pos + motif_len - word_len; ++word_start)
	  if (suff_tree.suffix_count (&dsq[word_start], &dsq[word_start + word_len]) >= min_reps)
	    {
	      seed.push_back (Seq_pos (seq, pos));
	      break;
	    }
    }
}

void Shortlist::copy_seeds (const Shortlist& source, int N)
{
  if (N == 0 || N > (int) source.seed.size())
    N = source.seed.size();
  seed = vector<Seq_pos> (source.seed.begin(), source.seed.begin() + N);
}

void Shortlist::full_sort()
{
  reset_cache();
  CTAG(3,SHORTLIST) << "Computing likelihoods for seed list\n";
  const int one_percent = max ((int) (seed.size() / 100), (int) 1);
  if (max_fork < 2)
    for (int i = 0; i < (int) seed.size(); ++i)
      {
	const Loge new_loglike = compute_loglike (seed[i]);
	if (diagnostic_flag) run_diagnostic (seed[i], new_loglike, FALSE);
	loglike[seed[i].seq][seed[i].pos] = new_loglike;
	if (i % one_percent == 0)
	  CTAG(3,SHORTLIST) << "Processed " << i+1 << "/" << seed.size() << " seeds\n";
      }
  else
    {
      open_pipes();
      for (int i = 0; i < (int) seed.size(); i += max_fork)
	{
	  const int j_end = min (i + max_fork, (int) seed.size());
	  for (int j = i; j < j_end; ++j)
	    if (fork_child(j-i))  // child
	      {
		const Loge new_loglike = compute_loglike (seed[j]);
		if (diagnostic_flag) run_diagnostic (seed[j], new_loglike, FALSE);
		ostream_cfile pipe_out (write_fd(j-i));
		pipe_out.out << new_loglike << "\n";
		pipe_out.close();
		exit(0);
	      }
	  
	  for (int j = i; j < j_end; ++j)
	    {
	      istream_cfile pipe_in (read_fd(j-i));
	      pipe_in.in >> loglike[seed[j].seq][seed[j].pos];
	      wait_for_child(j-i);  // wait for each child to exit
	    }

	  if (i % one_percent == 0)
	    CTAG(3,SHORTLIST) << "Processed " << i+1 << "/" << seed.size() << " seeds\n";
	}
      close_pipes();
    }
  CTAG(3,SHORTLIST) << "Processed all seeds; sorting\n";
  sort (seed.begin(), seed.end(), by_loglike);
}

void Shortlist::partial_sort (int N, Loge min_loglike)
{
  reset_cache();
  CTAG(3,SHORTLIST) << "Recomputing top " << N << " likelihoods for seed list\n";
  set<Seq_pos> fresh;  // set of refreshed seeds
  while (1)
    {
      int n;
      for (n = 0; n < N; ++n) if (fresh.find(seed[n]) == fresh.end()) break;  // find first seed that needs recomputing
      if (N > 0 && n == N) break;  // stop if first N are fresh
      const vector<Seq_pos>::iterator old_iter = seed.begin() + n;
      const Seq_pos update = *old_iter;
      const Loge old_loglike = loglike[update.seq][update.pos];
      if (old_loglike < min_loglike) break;
      const Loge new_loglike = compute_loglike (update);
      const double TINY_INCREASE = .01;
      if (new_loglike > old_loglike + TINY_INCREASE)
	{
	  CLOGERR << "WARNING: Seed likelihood at seq '" << index.profile[update.seq]->name << "' pos " << update.pos << " increased from " << old_loglike << " to " << new_loglike << "\n";
	  run_diagnostic (update, new_loglike, TRUE);
	}
      else if (diagnostic_flag)
	run_diagnostic (update, new_loglike, FALSE);
      loglike[update.seq][update.pos] = new_loglike;
      vector<Seq_pos>::iterator update_iter = lower_bound (old_iter + 1, seed.end(), update, by_loglike);
      if (update_iter > old_iter + 1)
	{
	  copy (old_iter + 1, update_iter, old_iter);
	  *--update_iter = update;
	}
      fresh.insert (update);
    }
  // consistency check commented out, because is_sorted() is an SGI-specific STL extension
  //  if (!is_sorted (seed.begin(), seed.end(), by_loglike))
  //    THROWEXPR ("Shortlist unsorted following partial refresh");
}

int Shortlist::n_above (Loge min_loglike) const
{
  // create a dummy Seq_pos with log-likelihood min_loglike
  ((Shortlist&) *this).loglike.push_back (vector<Loge> (1, min_loglike));
  Seq_pos dummy_seqpos (((int) loglike.size()) - 1, 0);
  // find last place dummy Seq_pos could go in seed list
  const vector<Seq_pos>::const_iterator seed_iter = upper_bound (seed.begin(), seed.end(), dummy_seqpos, by_loglike);
  // delete dummy Seq_pos from loglike, and return
  ((Shortlist&) *this).loglike.pop_back();
  return seed_iter - seed.begin();
}

void Shortlist::reset_cache() { }

void Shortlist::show (ostream& o) const
{
  typedef map <int, set<int> > Pos_by_seq;
  Pos_by_seq pos_by_seq;
  for_const_contents (vector<Seq_pos>, seed, sp) pos_by_seq[sp->seq].insert (sp->pos);
  for_const_contents (Pos_by_seq, pos_by_seq, sp)
    {
      o << "Sequence '" << index.name[sp->first] << "' " << name << " scores (bits)";
      for_const_contents (set<int>, sp->second, pos)
	o << " " << *pos << "=" << Nats2Bits (loglike[sp->first][*pos]);
      o << "\n";
    }
}

void Shortlist::run_diagnostic (const Seq_pos& seqpos, Loge loglike, bool be_verbose) { }

Seed_drill::Seed_drill() : first_pass (TRUE) { }

void Seed_drill::add_fixed_length_shortlist (Shortlist* s, int s_len)
{
  shortlist.push_back (s);
  len.push_back (s_len);
  min_loglike.push_back (-InfinityLoge);
}

void Seed_drill::add_variable_length_shortlist (Shortlist* s, Loge s_loglike)
{
  shortlist.push_back (s);
  len.push_back (0);
  min_loglike.push_back (s_loglike);
}

void Seed_drill::refresh()
{
  const int n_lists = shortlist.size();
  for (int n = 0; n < n_lists; ++n)
    {
      Shortlist& sl = *shortlist[n];
      const int s_len = (len[n] <= 0 || len[n] > (int) sl.seed.size()) ? ((int) sl.seed.size()) : len[n];
      const Loge s_loglike = min_loglike[n];

      CTAG(6,DRILL) << "Sorting " << sl.name << " seed list\n";

      if (n == 0 && !first_pass)
	sl.partial_sort (n_lists == 1 ? 1 : s_len, s_loglike);
      else
	sl.full_sort();

      if (CTAGGING(4,DRILL)) sl.show (CL);

      if (n < n_lists - 1)
	{
	  const int n_copy = min (s_len, sl.n_above (s_loglike));
	  shortlist[n+1]->copy_seeds (sl, n_copy);
	}
    }
  first_pass = FALSE;
}

const Seq_pos& Seed_drill::best_seed() const
{
  if (shortlist.back()->seed.size() == 0) THROWEXPR ("No valid seeds");
  return shortlist.back()->seed.front();
}

void Seed_drill::delete_shortlists()
{
  for_contents (vector<Shortlist*>, shortlist, s)
    delete *s;
  clear();
}

void Seed_drill::clear()
{
  shortlist.clear();
  len.clear();
  min_loglike.clear();
}

Loge Self_similarity_shortlist::compute_loglike (const Seq_pos& seqpos)
{
  const Named_profile& np = *index.profile[seqpos.seq];
  const int dsq_size = np.dsq.size();
  // check if word is in cache
  const int word_pos = min (seqpos.pos + offset, dsq_size - motif_len);
  const int first_pos = first_instance[seqpos.seq][word_pos];
  const Seq_pos first_seqpos (seqpos.seq, first_pos);
  if (fresh.find (first_seqpos) != fresh.end())
    return Score2Nats (fresh[first_seqpos]);
  // initialise
  Digitized_biosequence::const_iterator dsq_iter = np.dsq.begin();
  Metascore::const_iterator mask_iter = np.meta_sc[mask_metascore_idx].begin();
  Digitized_biosequence::const_iterator target_iter = dsq_iter + word_pos;
  Score fwd_loop_sc = 0;  // current forward score
  Score rev_loop_sc = 0;  // current reverse score
  vector<Score> fwd_start_sc (motif_len);  // rolling forward start score
  vector<Score> rev_start_sc (motif_len);  // rolling reverse start score
  // make revcomp of target
  if (revcomp)
    for (int i = 0; i < motif_len; ++i)
      revcomp_target[i] = alphabet.complement (target_iter[motif_len - 1 - i]);
  Digitized_biosequence::const_iterator revcomp_target_iter = revcomp_target.begin();
  // precalculate max achievable score from each motif position in forward & reverse
  vector<Score> max_fwd_sc (motif_len);
  vector<Score> max_rev_sc (motif_len);
  Score tmp_sc = 0; // do forward...
  for (int i = motif_len - 1; i >= 0; --i) max_fwd_sc[i] = (tmp_sc += max_sym_sc[target_iter[i]]);
  if (revcomp)  // do reverse...
    {
      tmp_sc = 0;
      for (int i = motif_len - 1; i >= 0; --i) max_rev_sc[i] = (tmp_sc += max_sym_sc[revcomp_target_iter[i]]);
    }
  Score max_metascore = -InfinityScore;
  if (dp_logging || long_dp_logging)
    {
      CTAG(10,SELFSIM_DP SELFSIM_LONG_DP);
      CL << "Sequence '" << np.name << "' seed position " << seqpos.pos << "\n";
      CL << "Max forward score from each pos = (" << max_fwd_sc << ")\n";
      if (revcomp)
	CL << "Max reverse score from each pos = (" << max_rev_sc << ")\n";
      CL << "Forward";
      if (revcomp) CL << "/reverse";
      CL << " scores at each endpos:";
    }
  // do DP
  const Score trunc_sc = -InfinityScore / motif_len;
  for (int endpos = 0; endpos < dsq_size; ++endpos)
    {
      max_metascore = max (max_metascore, mask_iter[endpos]);
      const int ptr = endpos % motif_len;
      Score fwd_sc = fwd_start_sc[ptr];
      Score rev_sc = rev_start_sc[ptr];
      const int startpos = endpos - motif_len + 1;
      if (startpos >= 0)
	{
	  Digitized_biosequence::const_iterator query_iter = dsq_iter + startpos;
	  Metascore::const_iterator query_mask_iter = mask_iter + startpos;

	  for (int i = 0; i < motif_len; ++i)
	    {
	      const int j = motif_len - i;
	      const Score fwd_ceiling_sc = (max_metascore <= trunc_sc) ? InfinityScore : (fwd_loop_sc - j * max_metascore);
	      const Score rev_ceiling_sc = (max_metascore <= trunc_sc) ? InfinityScore : (rev_loop_sc - j * max_metascore);
	      const bool fwd_hopeless = fwd_sc + max_fwd_sc[i] <= fwd_ceiling_sc;
	      const bool rev_hopeless = !revcomp || (rev_sc + max_rev_sc[i] <= rev_ceiling_sc);
	      if (fwd_hopeless && rev_hopeless)
		{
		  fwd_sc = rev_sc = -InfinityScore;
		  break;
		}
	      fwd_sc += submat (target_iter[i], query_iter[i]) + query_mask_iter[i];
	      if (revcomp)
		rev_sc += submat (revcomp_target_iter[i], query_iter[i]) + query_mask_iter[i];
	    }
	  const bool improved = fwd_sc > fwd_loop_sc || rev_sc > rev_loop_sc;
	  fwd_start_sc[ptr] = fwd_loop_sc = max (fwd_loop_sc, fwd_sc);
	  rev_start_sc[ptr] = rev_loop_sc = max (rev_loop_sc, rev_sc);
	  if (long_dp_logging || (dp_logging && improved))
	    {
	      CL << ' ' << endpos << '=';
	      ShowScore(fwd_loop_sc,CL);
	      if (revcomp) { CL << '/'; ShowScore(rev_loop_sc,CL); }
	    }
	}
    }
  const Score best_sc = max (fwd_loop_sc, rev_loop_sc);
  fresh[first_seqpos] = best_sc;  // insert into cache
  if (dp_logging)
    {
      CL << "\tFinal score ";
      ShowScore(best_sc,CL);
      CL << '\n';
    }
  return Score2Nats (best_sc);
}

Score Self_similarity_shortlist::compute_match_sc (const int* target, const int* query, const Score* query_mask) const
{
  Score sc = 0;
  for (int i = 0; i < motif_len; ++i)
    ScorePMulAcc (sc, ScorePMul (submat (target[i], query[i]), query_mask[i]));
  return max (sc, 0);
}

Loge Self_similarity_shortlist::compute_loglike_slow (const Seq_pos& seqpos)
{
  const Named_profile& np = *index.profile[seqpos.seq];
  const int dsq_size = np.dsq.size();
  const int word_pos = min (seqpos.pos + offset, dsq_size - motif_len);
  // initialise
  Digitized_biosequence::const_iterator dsq_iter = np.dsq.begin();
  Metascore::const_iterator mask_iter = np.meta_sc[mask_metascore_idx].begin();
  Digitized_biosequence::const_iterator target_iter = dsq_iter + word_pos;
  vector<Score> fwd_loop_sc (dsq_size + 1);
  vector<Score> rev_loop_sc (dsq_size + 1);
  fwd_loop_sc[0] = rev_loop_sc[0] = 0;
  // make revcomp of target
  if (revcomp)
    for (int i = 0; i < motif_len; ++i)
      revcomp_target[i] = alphabet.complement (target_iter[motif_len - 1 - i]);
  Digitized_biosequence::const_iterator revcomp_target_iter = revcomp_target.begin();
  // do DP
  for (int endpos = 0; endpos < dsq_size; ++endpos)
    {
      Score fwd_sc = fwd_loop_sc[endpos];
      Score rev_sc = rev_loop_sc[endpos];
      const int startpos = endpos - motif_len + 1;

      if (startpos >= 0)
	{
	  fwd_sc = max (fwd_sc, ScorePMul (fwd_loop_sc[startpos], compute_match_sc (&*target_iter, &*(dsq_iter + startpos), &*(mask_iter + startpos))));
	  rev_sc = max (rev_sc, ScorePMul (rev_loop_sc[startpos], compute_match_sc (&*revcomp_target_iter, &*(dsq_iter + startpos), &*(mask_iter + startpos))));
	}
      fwd_loop_sc[endpos+1] = fwd_sc;
      rev_loop_sc[endpos+1] = rev_sc;
    }
  const Score best_sc = max (fwd_loop_sc.back(), rev_loop_sc.back());
  if (dp_logging || long_dp_logging)
    {
      CTAG(10,SELFSIM_DP SELFSIM_LONG_DP);
      show_dp_vec (fwd_loop_sc, rev_loop_sc, CL);
    }
  return Score2Nats (best_sc);
}

void Self_similarity_shortlist::show_dp_vec (const vector<Score>& fwd, const vector<Score>& rev, ostream& out) const
{
  out << "Forward";
  if (rev.size()) out << "/reverse";
  out << " scores at each endpos:";
  Score best_fwd = 0;
  Score best_rev = 0;
  for (int i = 1; i < (int) fwd.size(); ++i)
    if (long_dp_logging || fwd[i] > best_fwd || (rev.size() ? (rev[i] > best_rev) : FALSE))
      {
	out << ' ' << i-1 << '=';
	ShowScore(fwd[i],out);
	best_fwd = max (best_fwd, fwd[i]);
	if (rev.size())
	  {
	    out << '/';
	    ShowScore(rev[i],out);
	    best_rev = max (best_rev, rev[i]);
	  }
      }
  out << " Final score ";
  ShowScore(max(best_fwd,best_rev),out);
  out << '\n';
}

void Self_similarity_shortlist::run_diagnostic (const Seq_pos& seqpos, Loge loglike, bool be_verbose)
{
  const Loge slow_loglike = compute_loglike_slow (seqpos);
  const Named_profile& np = *index.profile[seqpos.seq];
  const sstring& seqname = np.name;
  const bool found_error = Nats2Score(loglike) != Nats2Score(slow_loglike);  // convert to Scores to avoid floating-point errors
  if (found_error)
    CLOGERR << "WARNING: at seq '" << seqname << "' pos " << seqpos.pos << ": fast DP gives log-likelihood " << loglike << ", slow DP gives " << slow_loglike << "\n";
  if (found_error || be_verbose)
    {
      if (!(dp_logging || long_dp_logging))
	{
	  const bool tmp_dl = dp_logging;
	  map<Seq_pos,Score> tmp_fresh;
	  tmp_fresh.swap (fresh);  // make the cache disappear for a moment
	  dp_logging = TRUE;
	  CLOGERR << "Repeating fast DP with logging turned on\n";
	  compute_loglike (seqpos);
	  CLOGERR << "Repeating slow DP with logging turned on\n";
	  compute_loglike_slow (seqpos);
	  dp_logging = tmp_dl;
	  tmp_fresh.swap (fresh);
	}
      CLOGERR << "Mask for sequence '" << seqname << "': (" << np.meta_sc[mask_metascore_idx] << ")\n";
    }
}

Self_similarity_shortlist::Self_similarity_shortlist (const FASTA_sequence_database& seq_db,
						      int max_fork,
						      const Alphabet& alphabet,
						      int motif_len,
						      int offset,
						      Prob p_match,
						      bool revcomp,
						      int mask_metascore_idx)
  : Shortlist (seq_db, max_fork),
    first_instance (seq_db.size()),
    dp_logging (CTAGGING(1,SELFSIM_DP)),
    long_dp_logging (CTAGGING(1,SELFSIM_LONG_DP)),
    alphabet (alphabet),
    motif_len (motif_len),
    offset (offset),
    mask_metascore_idx (mask_metascore_idx),
    revcomp (revcomp),
    submat (alphabet.size(), alphabet.size()),
    max_sym_sc (alphabet.size(), -InfinityScore),
    revcomp_target (motif_len)
{
  name = "self-similarity";
  for (int i = 0; i < alphabet.size(); ++i)
    for (int j = 0; j < alphabet.size(); ++j)
      {
	submat(i,j) = Prob2Score (i == j ? p_match : (1 - p_match));
	max_sym_sc[i] = max (max_sym_sc[i], submat(i,j));
      }
  for (int seq = 0; seq < seq_db.size(); ++seq)
    {
      const Digitized_biosequence& dsq = seq_db.index.profile[seq]->dsq;
      first_instance[seq] = vector<int> (dsq.size());
      map<vector<int>,int> first_pos;
      for (int pos = 0; pos < ((int) dsq.size()) - motif_len; ++pos)
	{
	  const vector<int> word (dsq.begin() + pos, dsq.begin() + pos + motif_len);
	  if (first_pos.find (word) == first_pos.end())
	    first_pos[word] = pos;
	  first_instance[seq][pos] = first_pos[word];
	}
    }
  if (CTAGGING(2,SHORTLIST))
    {
      CL << "Substitution matrix for " << name << " shortlist:\n" << submat;
      CL << "Max score for each symbol: (" << max_sym_sc << ")\n";
    }
}

void Self_similarity_shortlist::reset_cache()
{
  fresh.clear();
}

Loge Local_shortlist::compute_loglike (const Seq_pos& seqpos)
{
  const Named_profile& np = *index.profile[seqpos.seq];
  trainer.local_seed (np, seqpos.pos);
  trainer.local_EM (np);  // do one round of EM
  return trainer.local_forward (np);  // return forward score
}

Local_shortlist::Local_shortlist (const FASTA_sequence_database& seq_db, int max_fork, Local_trainer& trainer)
  : Shortlist (seq_db, max_fork),
    trainer (trainer)
{
  name = "single-sequence-modeling";
}

Loge Multihit_shortlist::compute_loglike (const Seq_pos& seqpos)
{
  model.hitter.seed (*index.profile[seqpos.seq], seqpos.pos);
  const vector<Named_profile*>* targets = 0;
  if (similar_np.size())
    targets = &similar_np[seqpos.seq];
  else
    targets = &index.profile;
  model.multihit_EM (*targets);  // do one round of EM
  return model.multihit_forward (*targets);  // return forward score
}

Multihit_shortlist::Multihit_shortlist (const FASTA_sequence_database& seq_db, Multihitter& model)
  : Shortlist (seq_db, 1),
    model (model)
{
  name = "database-modeling";
}

Multihit_shortlist::Multihit_shortlist (const FASTA_sequence_database& seq_db,
					Multihitter& model,
					Dist_func& dist,
					int n_similar)
  : Shortlist (seq_db, 1),
    model (model),
    similar_np (seq_db.size(), vector<Named_profile*> (n_similar, (Named_profile*) 0))
{
  n_similar = min (n_similar, seq_db.size());
  name << "database-nearest-" << n_similar << "-modeling";
  for (int i = 0; i < seq_db.size(); ++i)
    {
      const Nearest sim_idx (dist, n_similar, seq_db.size(), i);
      for (int j = 0; j < n_similar; ++j)
	similar_np[i][j] = index.profile[sim_idx[j]];
    }
}
