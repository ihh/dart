#include "empath/trirep.h"
#include "util/vector_output.h"

bool Trirep_set::contains (Digitized_biosequence& dsq)
{
  for (int len = 1; len <= (int) dsq.size(); ++len)
    {
      if (dsq.size() % len != 0) continue;
      // check for periodicity
      bool periodic = 1;
      for (int i = 0; i < (int) dsq.size(); i += len)
	for (int j = 0; j < len; ++j)
	  if (dsq[i+j] != dsq[j])
	    {
	      periodic = 0;
	      break;
	    }
      if (!periodic) continue;
      // check for all cyclic permutations
      Digitized_biosequence subseq (len);
      for (int perm = 0; perm < len; ++perm)
	{
	  for (int i = 0; i < len; ++i)
	    subseq[i] = dsq[(i+perm) % len];
	  if (reps.find(subseq) != reps.end())
	    return 1;
	}
    }
  return 0;
}

void Trirep_set::add_degenerates()
{
  sstring degen = "rymksw";
  for (int i = 0; i < (int) degen.size(); ++i)
    {
      Digitized_biosequence dsq;
      dsq.push_back (alphabet.char2int_deg (degen[i]));
      reps.insert (dsq);
    }
}

Trirep_set::Trirep_set (int max_motif_len, int max_palindrome_len, const Alphabet& alphabet) : max_len (max_len), alphabet (alphabet)
{
  for (int len = 1; len <= max (max_motif_len, max_palindrome_len); ++len)
    {
      if (len > max_motif_len && len % 2 != 0) continue;
      Digitized_biosequence dsq (len, (int) 0);
      while (1)
	{
	  bool accept = 1;
	  if (len > max_motif_len)  // check for palindrome
	    for (int i = 0; i < len / 2; ++i)
	      if (dsq[i] != alphabet.complement (dsq[len-1-i]))
		{ accept = 0; break; }
	  if (accept && !contains(dsq)) reps.insert (dsq);
	  int i = 0;
	  while (i < len)
	    {
	      if (++dsq[i] < alphabet.size()) break;
	      dsq[i++] = 0;
	    }
	  if (i == len) break;
	}
    }
}

ostream& operator<< (ostream& o, const Trirep_set& trirep_set)
{
  Biosequence seq;
  for_const_contents (set<Digitized_biosequence>, trirep_set.reps, rep)
    {
      trirep_set.alphabet.dsq2seq (*rep, seq);
      o << seq << "\n";
    }
  return o;
}

int Trirep_model::n_states (const Trirep_set& trirep_set, int units)
{
  int n = 0;
  for_const_contents (set<Digitized_biosequence>, trirep_set.reps, rep)
    n += (rep->size() + 2) * units;
  return n;
}

Trirep_model::Repdata::Repdata (const Digitized_biosequence& dsq, const Alphabet& alph)
{
  alph.dsq2weight (dsq, rep_profile);
}

Trirep_model::Trirep_model (Trirep_set& trirep_set, PScores& pscore, int min_reps, int shape, double p_mismatch, double p_ins, double p_del)
  : Trainable (n_states (trirep_set, shape + min_reps), trirep_set.alphabet, pscore),
    trirep_set (trirep_set),
    min_reps (min_reps),
    shape (shape)
{
  // pseudocounts/parameters
  const double pcount_loopfwd = .1;     // pseudocount for loopforward
  const double pcount_loopback = .1;    // pseudocount for loopback
  const double p_rep_extend = .5;       // probability of extending repeat (for which_rep[])
  const double pcount_which_rep = 100;  // number of pseudocounts to share between which_rep[]'s

  // set up match, insert, delete PGroups
  match = pscore.new_boolean_group ("Microsatellite match");
  pscore[match.NO] = Prob2Score (p_mismatch);
  pscore[match.YES] = Prob2Score (1 - p_mismatch);

  ins = pscore.new_boolean_group ("Microsatellite stutter");
  pscore[ins.NO] = Prob2Score (1 - p_ins);
  pscore[ins.YES] = Prob2Score (p_ins);

  del = pscore.new_boolean_group ("Microsatellite delete");
  pscore[del.NO] = Prob2Score (1 - p_del);
  pscore[del.YES] = Prob2Score (p_del);

  // set up Repdata's
  Biosequence seq;  // variable holding text representation of current rep
  int next_state = 0;  // next free state
  map<int,int> len_count;  // number of repeats with each particular length
  for_const_contents (set<Digitized_biosequence>, trirep_set.reps, rep)
    {
      const int len = rep->size();
      ++len_count[len];
      trirep_set.alphabet.dsq2seq (*rep, seq);

      // initialise Repdata
      Repdata rd (*rep, trirep_set.alphabet);

      // assign start states
      for (int i = 0; i < shape + min_reps; ++i)
	{
	  rd.start_state.push_back (next_state);
	  next_state += len + 2;
	}

      // create start transition group
      sstring start_group_name;
      start_group_name << "Poly-'" << seq << "' start";
      rd.start = pscore.new_group (shape, start_group_name.c_str());

      // set Dirichlet prior for start
      vector<Prob> start_pseudocount (shape, 1.0);
      Dirichlet_mixture start_mixture (start_pseudocount);
      prior.assign (rd.start, start_mixture);

      // create loop transition group
      sstring loop_group_name;
      loop_group_name << "Poly-'" << seq << "' loop";
      rd.loop = pscore.new_boolean_group (loop_group_name.c_str());
      
      // set Dirichlet prior for loop
      vector<Prob> loop_pseudocount (2);
      loop_pseudocount[0] = pcount_loopfwd;
      loop_pseudocount[1] = pcount_loopback;
      Dirichlet_mixture loop_mixture (loop_pseudocount);
      prior.assign (rd.loop, loop_mixture);

      // save rd
      repdata.push_back (rd);
      repname.push_back (seq);
    }

  // create which_rep group, set Dirichlet prior
  which_rep = pscore.new_group (repdata.size(), "Repeat selector", repname);
  vector<Prob> rep_pseudocount (repdata.size());
  for (int i = 0; i < (int) repdata.size(); ++i)
    {
      const int len = repdata[i].rep_profile.size();
      rep_pseudocount[i] = pcount_which_rep * pow (p_rep_extend, len) * (1 - p_rep_extend) / (double) len_count[len];
    }
  Dirichlet_mixture rep_mixture (rep_pseudocount);
  prior.assign (which_rep, rep_mixture);
  
  // assign PFunc's
  for (int nrep = 0; nrep < (int) repdata.size(); ++nrep)
    {
      const Repdata& rd = repdata[nrep];
      const int len = rd.rep_profile.size();
      
      // visit each loop unit
      for (int unit = 0; unit < shape + min_reps; ++unit)
	{
	  const int start = rd.start_state[unit];
	  const int end = start + len + 1;

	  // unit start and end states are null
	  hmm.state_type[start] = hmm.state_type[end] = Single_PHMM::Null;

	  // loopforward transition (mandatory for last min_reps units)
	  if (unit < shape)
	    hmm.transition (end, rd.start_state[unit + 1]) = rd.loop.NO;
	  else if (unit < shape + min_reps - 1)
	    hmm.transition (end, rd.start_state[unit + 1]) = 1;
	  
	  // loopback transition
	  if (unit < shape)
	    hmm.transition (end, start) = rd.loop.YES;

	  for (int pos = 0; pos <= len; ++pos)
	    {
	      const int state = start + pos;

	      // emissions
	      if (pos > 0)
		{
		  hmm.state_type[state] = Single_PHMM::Emit;
		  for (int sym = 0; sym < hmm.alphabet().size(); ++sym)
		    hmm.emit[state][sym] = match.NO * null_emit[sym];
		  for_const_contents (Symbol_weight_map, rd.rep_profile[pos-1], sw)
		    hmm.emit[state][sw->first] += match.YES * sw->second;
		}
	      
	      // stutter
	      PFunc no_indel = 1;
	      if (pos > 0 && len > 1)
		{
		  no_indel = ins.NO;
		  hmm.transition (state, state) = ins.YES;
		}

	      // end transition (last loop unit only)
	      if (pos > 0 && unit == shape + min_reps - 1)
		hmm.transition (pos < len ? state : state + 1, End) = no_indel / (double) len;

	      // deletions
	      const int last_allowed_dest = pos > 0 ? end : end - 1;  // can't delete entire motif
	      for (int dest = state + 1; dest < last_allowed_dest; ++dest)
		{
		  hmm.transition (state, dest) = no_indel * del.NO;
		  no_indel *= del.YES;
		}
	      hmm.transition (state, last_allowed_dest) = no_indel;
	  
	      // start transition
	      if (pos > 0 && unit < shape)
		hmm.transition (Start, state) = which_rep[nrep] * rd.start[unit] / (double) len;
	    }
	}
    }
  // add mask metascores
  mask_metascore_idx = 0;
  for (int i = 0; i < hmm.states(); ++i)
    hmm.metascore_idx[i].push_back (mask_metascore_idx);
}

void Trirep_model::write_repdist (ostream& out) const
{
  for_const_contents (vector<Repdata>, repdata, rd)
    {
      Biosequence seq;
      hmm.alphabet().weight2seq (rd->rep_profile, seq);

      out << seq;
      out << " mat " << Score2Prob (pscore[match.YES]);
      out << " ins " << Score2Prob (pscore[ins.YES]);

      out << " start (";
      for (int i = 0; i < rd->start.group_size; ++i) out << " " << Score2Prob (pscore[rd->start[i]]);
      out << " ) loop " << Score2Prob (pscore[rd->loop.YES]);
      out << "\n";
    }
}
