#include "empath/ungapped.h"
#include "util/math_fn.h"

Ungapped_model::Ungapped_model (int min_size, int max_size, bool revcomp, PScores& pscore, bool mask) :
  Trainable ((revcomp ? 2 : 1) * max_size + StateOffset, DNA_alphabet, pscore),
  min_size (min_size),
  max_size (max_size),
  revcomp (revcomp)
{
  // declare the PGroups
  for (int i = 0; i < max_size; ++i)
    {
      sstring gname;
      gname << "Match emit #" << i+1;
      match_emit.push_back (new_emit_group (gname.c_str()));
    }
  const int max_skip = max_size - min_size;
  skip_start = pscore.new_group (max_skip + 1, "Start skip");
  skip_end = pscore.new_group (max_skip + 1, "End skip");
  // set up the emissions
  for (int pos = 0; pos < max_size; ++pos)
    for (int sym = 0; sym < hmm.alphabet().size(); ++sym)
      {
	hmm.emit[StateOffset + pos][sym] = match_emit[pos][sym];
	if (revcomp)
	  hmm.emit[StateOffset + 2*max_size - pos - 1][hmm.alphabet().complement(sym)] = match_emit[pos][sym];
      }
  // set up the transitions
  for (int pos = 0; pos < max_size; ++pos)
    {
      if (pos <= max_skip)
	{
	  hmm.transition (FwdStart, StateOffset + pos) = skip_start[pos];
	  if (revcomp)
	    hmm.transition (RevStart, StateOffset + max_size + pos) = skip_end[pos];
	}
      if (pos < max_size - 1)
	{
	  hmm.transition (StateOffset + pos, StateOffset + pos + 1) = 1;
	  if (revcomp)
	    hmm.transition (StateOffset + max_size + pos, StateOffset + max_size + pos + 1) = 1;
	}
      if (pos >= min_size - 1)
	{
	  const int skip = max_size - 1 - pos;
	  hmm.transition (StateOffset + pos, FwdEnd) = skip_end[skip];
	  if (revcomp)
	    hmm.transition (StateOffset + max_size + pos, RevEnd) = skip_start[skip];
	}
    }
  // set up the prior for transitions into spacer block (emit prior is handled by Trainable)
  const Prob skip_pseudocount = 1.0e6;  // nice big number forces a hard decision; but not TOO big (or else rounding errors)
  const double n_cpt = .5 * (max_size + 1 - min_size) * (max_size + 2 - min_size);  // no. of components
  vector<int> group_size (2, (int) skip_start.group_size);
  Dirichlet_mixture skip_mixture (group_size, (int) n_cpt);
  int cpt = 0;
  for (int s = 0; s <= max_skip; ++s)
    {
      vector<Prob> skip_start_pseudocount (skip_start.group_size, (Prob) 0);
      skip_start_pseudocount[s] = skip_pseudocount;
      for (int e = 0; e + s <= max_skip; ++e)
	{
	  vector<Prob> skip_end_pseudocount (skip_end.group_size, (Prob) 0);
	  skip_end_pseudocount[e] = skip_pseudocount;
	  skip_mixture.log_cpt_prior[cpt] = -Prob2Nats(n_cpt);
	  skip_mixture.alpha[cpt][0] = skip_start_pseudocount;
	  skip_mixture.alpha[cpt][1] = skip_end_pseudocount;
	  ++cpt;
	}
    }
  if (cpt != n_cpt) THROWEXPR ("fuckup!");
  // do the assignment
  Multigroup skip_multi = skip_start * skip_end;
  prior.assign (skip_multi, skip_mixture);
  // set mask metascore
  int next_free_meta = 0;
  if (mask)
    {
      mask_metascore_idx = next_free_meta++;
      for (int pos = StateOffset; pos < hmm.states(); ++pos)
	hmm.metascore_idx[pos].push_back (mask_metascore_idx);
    }
  // record reverse states
  if (revcomp)
    for (int s = max_size + StateOffset; s < 2 * max_size + StateOffset; ++s)
      reverse_states.push_back(s);
  // set up the seed path
  for (int pos = StateOffset; pos < max_size + StateOffset; ++pos)
    seed_path.push_back (pos);
  // attach PScope for Telegraph
  hmm.set_pscope (pscore);
  // print HMM out to the logfiles
  if (CTAGGING(2,UNGAPPED_DISPLAY)) hmm.show(CL);
}

Fixed_model::Fixed_model (int size, bool revcomp, PScores& pscore, bool mask)
  : Ungapped_model (size, size, revcomp, pscore, mask)
{ }

void Fixed_model::get_local_scores (const Named_profile& np, Metascore& result, bool viterbi)
{
  // make local copy of match scores & other params
  const int motiflen = size();
  const int seqlen = np.size();
  vector<vector<Score>::const_iterator> match_iter (motiflen);
  for (int i = 0; i < motiflen; ++i) match_iter[i] = pscore[match_emit[i]].begin();
  const Alphabet& alphabet = hmm.alphabet();
  const Digitized_biosequence::const_iterator fwd_iter = np.dsq.begin();
  const bool using_metascores = mask_metascore_idx >= 0;
  Metascore::const_iterator mask_iter;
  if (using_metascores) mask_iter = np.meta_sc[mask_metascore_idx].begin();
  // make complement sequence (NB not reversed)
  Digitized_biosequence comp;
  if (revcomp)
    {
      comp.reserve (seqlen);
      for (int i = 0; i < seqlen; ++i) comp[i] = alphabet.complement (fwd_iter[i]);
    }
  const Digitized_biosequence::const_iterator rev_iter = comp.begin();
  // initialise result vector, fill first few cells with -inf
  result = Metascore (seqlen);
  for (int endpos = 0; endpos < min (seqlen, motiflen-1); ++endpos)
    result[endpos] = -InfinityScore;
  // fill remaining cells
  for (int startpos = 0; startpos < seqlen - motiflen; ++startpos)
    {
      // get forward-strand score
      Score sc = 0;
      for (int i = 0; i < motiflen; ++i)
	ScorePMulAcc (sc, match_iter[i][fwd_iter[startpos+i]]);
      // get reverse-strand score
      if (revcomp)
	{
	  Score rev = 0;
	  for (int i = 0; i < motiflen; ++i)
	    ScorePMulAcc (rev, match_iter[motiflen-i-1][rev_iter[startpos+i]]);
	  sc = viterbi ? max (sc, rev) : ScorePMul (sc, rev);
	}
      // add mask
      if (using_metascores)
	for (int i = 0; i < motiflen; ++i)
	  ScorePMulAcc (sc, mask_iter[startpos+i]);
      // store
      result[startpos+motiflen-1] = sc;
    }
}
