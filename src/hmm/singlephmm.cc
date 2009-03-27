#include "hmm/singlephmm.h"
#include "util/vector_output.h"

Single_PHMM::Single_PHMM (int states, const Alphabet& alphabet) :
  Single_meta_HMM<PFunc> (states, alphabet, PFunc(0.0)),
  group_suffix (0)
{ }

Single_PHMM::Single_PHMM (int states, const Alphabet& alphabet, const vector<vector<sstring> >& group_suffix) :
  Single_meta_HMM<PFunc> (states, alphabet, PFunc(0.0)),
  group_suffix (&group_suffix)
{ }

void Single_PHMM::set_pscope (const PScope& ps)
{
  pscope = &ps;
}

void Single_PHMM::set_group_suffix (const vector<vector<sstring> >& g_suffix)
{
  group_suffix = &g_suffix;
}

Single_HMM_scores Single_PHMM::eval_hmm_scores_gs (const PScores& var_scores) const
{
  ((Single_PHMM*) this) -> set_group_suffix (var_scores.group_suffix);
  return eval_hmm_scores (var_scores);
}

Single_HMM_scores Single_PHMM::eval_hmm_scores (const PScores& var_scores) const
{
  if (CTAGGING(1,HMM_PSCORES))
    {
      CL << "PScores:\n";
      var_scores.show (CL);
    }
  if (CTAGGING(3,HMM_PFUNC))
    {
      CL << "Evaluating Single_HMM_scores from PFunc's:\n";
      show(CL);
      if (CTAGGING(2,HMM_PFUNC)) var_scores.show(CL);
    }
  Single_HMM_scores hmm_scores (*this);
  hmm_scores.name = name;
  hmm_scores.state_name = state_name;
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	hmm_scores.transition(src,dest) = transition(src,dest).eval_sc(var_scores);
      for (int i = 0; i < (int) emit[src].size(); ++i)
	hmm_scores.emit[src][i] = emit[src][i].eval_sc(var_scores);
      hmm_scores.start[src] = start[src].eval_sc(var_scores);
      hmm_scores.end[src] = end[src].eval_sc(var_scores);
    }
  hmm_scores.start_to_end() = start_to_end().eval_sc(var_scores);
  if (CTAGGING(-1,HMM_SCORES))
    {
      CL << "Evaluated the following HMM scores from PFunc's:\n";
      hmm_scores.show (CL);
    }
  return hmm_scores;
}

void Single_PHMM::inc_var_counts (PCounts& var_counts,
				  const PScores& var_scores,
				  const Single_HMM_counts& hmm_counts,
				  const Prob model_count) const
{
  assert_same_dimensions (hmm_counts);
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest).inc_var_counts (var_counts, var_scores, hmm_counts.transition(src,dest) * model_count);
      for (int i = 0; i < (int) emit[src].size(); ++i)
	emit[src][i].inc_var_counts (var_counts, var_scores, hmm_counts.emit[src][i] * model_count);
      start[src].inc_var_counts (var_counts, var_scores, hmm_counts.start[src] * model_count);
      end[src].inc_var_counts (var_counts, var_scores, hmm_counts.end[src] * model_count);
    }
  start_to_end().inc_var_counts (var_counts, var_scores, hmm_counts.start_to_end() * model_count);
}

Single_local_PHMM::Single_local_PHMM (const Single_PHMM& global,
				      const Alphabet_group& null_emit,
				      const Boolean_group& null_extend,
				      PScope& pscope,
				      int mask_metascore_idx,
				      bool factor_in_null_model) :
  Single_PHMM (global.states() + GlobalOffset, global.alphabet()),
  global (global),
  null_emit (null_emit),
  null_extend (null_extend),
  mask_metascore_idx (mask_metascore_idx),
  factor_in_null_model (factor_in_null_model)
{
  set_pscope (pscope);
  update();
}

void Single_local_PHMM::update()
{
  // debugging output
  if (CTAGGING(-2,LOCAL_PHMM))
    {
      CL << "Building local model from the following global model:\n";
      global.show (CL);
    }
  // test that we've got the right number of states
  const int gs = global.states();
  if (states() != local_state(gs))
    THROWEXPR ("Global model has " << gs << " states, local model has " << states() << " != " << gs << " + " << GlobalOffset);
  // clear all transitions, emit vectors & metascore indices
  reset (PFunc(0.0));
  reset_meta();
  // set padding state_types
  state_type[LeftPad] = state_type[RightPad] = Emit;
  state_type[LocalStart] = state_type[LocalEnd] = Null;
  // set padding state emit "probabilities" to 1.0 (actually these are odds ratios to null model)
  for (int sym = 0; sym < alphabet().size(); ++sym)
    {
      emit[LeftPad][sym]  = 1.0;
      emit[RightPad][sym] = 1.0;
    }
  // set transition "probabilities". again, these are odds ratios, with null_extend factored out
  // left padding transitions
  transition (Start,    LocalStart) = 1;
  transition (Start,    End)        = 1;
  transition (Start,    LeftPad)    = 1;
  transition (LeftPad,  LeftPad)    = 1;
  transition (LeftPad,  LocalStart) = 1;
  transition (LeftPad,  End)        = 1;
  // right padding transitions
  transition (LocalEnd, End)        = 1;
  transition (LocalEnd, RightPad)   = 1;
  transition (RightPad, RightPad)   = 1;
  transition (RightPad, End)        = 1;
  // copy the global model
  for (int s = 0; s < gs; ++s)
    {
      int ls = local_state(s);
      // state type
      const int type = state_type[ls] = global.state_type[s];
      // start & end transitions
      transition (LocalStart, ls) = global.transition (Start, s);
      transition (ls, LocalEnd)   = global.transition (s, End);
      // other transitions
      for (int t = 0; t < gs; ++t)
	transition (local_state(t), ls) = global.transition(t,s);
      // emit profiles. NB null model factored out
      if (type == Emit)
	for (int sym = 0; sym < alphabet().size(); ++sym)
	  {
	    emit[ls][sym] = global.emit[s][sym];
	    if (factor_in_null_model)
	      emit[ls][sym] /= null_emit[sym] * null_extend.YES;
	  }
      // metascore indices
      metascore_idx[ls] = global.metascore_idx[s];
    }
}

void Single_local_PHMM::allow_multiple_hits()
{
  transition (LocalEnd, LocalStart) = 1;
  transition (RightPad, LocalStart) = 1;
}

void Single_local_PHMM::search (const PScores& pscores,
				const Sequence_database& db,
				GFF_list& results,
				const vector<int>& reverse_states,
				Single_matrix_factory& matrix_factory,
				const char* source,
				const char* feature,
				const char* pattern) const
{
  Single_HMM_scores hmm_scores = eval_hmm_scores (pscores);
  for_const_contents (Sequence_database, db, np)
    {
      if (CLOGGING(6)) CL << "Searching sequence '" << (*np).name << "'\n";
      const Single_Viterbi_interface* viterbi_matrix = matrix_factory.new_Viterbi (hmm_scores, *np);
      const vector<int> viterbi_path = viterbi_matrix->optimal_state_path();
      int seq_pos = 0;

      GFF gff;
      gff.seqname = (*np).name;
      gff.source = source;
      gff.feature = feature;
      gff.frame = GFF::NoFrame;

      bool in_match = 0;
      for (int path_pos = 0; path_pos < (int) viterbi_path.size(); ++path_pos)
	{
	  const bool was_in_match = in_match;
	  const int state = viterbi_path[path_pos];
	  const GFF::Strand strand = find (reverse_states.begin(), reverse_states.end(), global_state(state)) == reverse_states.end() ? GFF::PlusStrand : GFF::MinusStrand;
	  const bool state_masked = find (hmm_scores.metascore_idx[state].begin(), hmm_scores.metascore_idx[state].end(), mask_metascore_idx) != hmm_scores.metascore_idx[state].end();

	  // check for end of match
	  if (in_match && (!state_masked || gff.strand != strand))
	    {
	      gff.group.clear();
	      if (strlen (pattern))
		gff.group << '/' << pattern << "/ ";
	      if (gff.strand == GFF::PlusStrand)
		for (int i = gff.start; i <= gff.end; ++i)
		  gff.group << np->seq[i-1];
	      else
		for (int i = gff.end; i >= gff.start; --i)
		  gff.group << alphabet().complement_char (np->seq[i-1]);
	      
	      results.push_back (gff);
	      in_match = 0;
	    }

	  // check for start of match
	  if (!in_match && state_masked)
	    {
	      // initialise GFF. note that score is set up to include null_extend.NO
	      gff.strand = strand;
	      gff.start = gff.end = seq_pos + 1;
	      gff.score = 0;
	      in_match = 1;
	    }

	  if (in_match)
	    {
	      // calculate transition score
	      if (was_in_match)
		gff.score += Score2Bits (hmm_scores.transition (viterbi_path[path_pos-1], state));
	      // calculate emit score (uses Score_profile)
	      if (state_type[state] == Emit)
		{
		  gff.end = seq_pos + 1;  // update GFF end coord if this is an emit state
		  Score emit_sc = hmm_scores.emit[state][np->dsq[seq_pos]];
		  for_const_contents (vector<int>, hmm_scores.metascore_idx[state], meta_idx)
		    ScorePMulAcc (emit_sc, np->meta_sc[*meta_idx][seq_pos]);
		  gff.score += Score2Bits (emit_sc);
		}
	    }
	  // update sequence pointer
	  if (state_type[state] == Emit) ++seq_pos;
	}
      if (in_match) CLOGERR << "Doh! faulty path through local model";
      delete viterbi_matrix;
    }
}

void Single_local_PHMM::reset_mask (Named_profile& np) const
{
  if ((int) np.meta_sc.size() < mask_metascore_idx)
    CLOGERR << "Warning -- sequence '" << np.name << "' is missing metascores\n";
  while ((int) np.meta_sc.size() < mask_metascore_idx + 1)
    np.meta_sc.push_back (Metascore());
  np.meta_sc[mask_metascore_idx] = Metascore (np.size(), (Score) 0);
}

void Single_local_PHMM::reset_masks (Sequence_database& db) const
{
  for_contents (Sequence_database, db, np) reset_mask (*np);
}

void Single_local_PHMM::mask (Named_profile& np, const vector<Metaprob>& expected_metacounts) const
{
  const Metaprob& in_motif_prob = expected_metacounts[mask_metascore_idx];
  Metascore delta (in_motif_prob.size());
  
  for (int pos = 0; pos < (int) in_motif_prob.size(); ++pos)
    delta[pos] = Prob2Score (1.0 - in_motif_prob[pos]);
  
  if (CTAGGING(4,LOCALHMM_MASK)) CL << "Masking '" << np.name << "' with (" << delta << ")\n";
  np.add_metascores (mask_metascore_idx, delta);
}

void Single_local_PHMM::get_mask_gff (GFF_list& gff_list, const vector<Metaprob>& expected_metacounts, Prob min_prob, const char* seqname, const char* seq, const char* source, const char* feature) const
{
  gff_list.acquire_mask (expected_metacounts[mask_metascore_idx], min_prob, seqname, seq, source, feature, 0);
}

