#include <numeric>

#include "handel/transducer.h"
#include "handel/transducer_sexpr.h"

#define FAKE_ALPHABET_PREFIX "SYM"
#define PVAR_TRANSITION_CHAR '>'

Score Pair_transducer_scores::pairwise_path_score (const Pairwise_path& path)
{
  const bool logging = CTAGGING(2,PAIR_TRANSDUCER_PATH_SCORE);
  if (logging)
    {
      CL << "Calculating pairwise alignment score by summing over pair transducer state paths\n";

      CL << "Transducer:\n";
      show (CL);

      CL << "Alignment path:\n";
      vector<sstring> row_name;
      row_name.push_back (sstring ("In"));
      row_name.push_back (sstring ("Out"));
      path.show (CL, row_name);
    }

  Score result_sc = -InfinityScore;
  const int cols = path.columns();

  if (cols == 0)
    {
      result_sc = start_to_end();
      for (int s = 0; s < states(); ++s)
	if (state_type[s] == TransducerWaitType)  // loop over wait states
	  ScorePSumAcc (result_sc, ScorePMul (transition (Start, s),
					      transition (s, End)));
    }
  else
    {
      // get null states
      vector<int> null;
      for (int s = 0; s < states(); ++s)
	if (state_type[s] == TransducerWaitType)
	  null.push_back (s);
      const vector<int> null_sorted = Transition_methods::topological_sort (*this, null);
      const vector<vector<int> > incoming = Transition_methods::incoming_states (*this);

      // create DP matrix
      if (logging)
	CTAG(3,PAIR_TRANSDUCER_PATH_SCORE) << "Allocating pairwise transducer path DP matrix: "
					   << path.columns() + 1 << " columns * " << states() << " states\n";
      array2d<Score> cell (path.columns() + 1, states(), -InfinityScore);
      for (int col = 0; col <= cols; ++col)
	{
	  if (logging)
	    CTAG(2,PAIR_TRANSDUCER_PATH_SCORE) << "Column " << col << ':';

	  // make ordered list of destination states to fill
	  vector<int> col_state (null_sorted.rbegin(), null_sorted.rend());  // reverse order of null states for col_state
	  if (col > 0)
	    {
	      const unsigned int col_type = (path(0,col-1) ? 1 : 0) | (path(1,col-1) ? 2 : 0);  // 0==TransducerWaitType, 1==TransducerDeleteType, 2==TransducerInsertType, 3==TransducerMatchType
	      if (col_type)
		for (int dest = 0; dest < states(); ++dest)
		  if (state_type[dest] == col_type)
		    col_state.push_back (dest);
	    }

	  // fill dest states in order
	  for_const_reverse_contents (vector<int>, col_state, dest)  // traverse col_state in reverse
	    {
	      // get address of cell
	      Score& cell_sc = cell (col, *dest);
	      // do start transitions
	      if (col == (state_type[*dest] ? 1 : 0))  // start goes to column #0 for null states, #1 for emit states
		ScorePSumAcc (cell_sc, transition (Start, *dest));
	      // loop over incoming transitions
	      for_const_contents (vector<int>, incoming[*dest], src)
		{
		  const Score src_sc = cell (state_type[*dest] ? (col-1) : col, *src);
		  const Score trans_sc = transition (*src, *dest);
		  ScorePSumAcc (cell_sc, ScorePMul (src_sc, trans_sc));
		}

	      // log
	      if (logging)
		{
		  CL << ' ' << *dest << '[';
		  ShowScore (cell_sc, CL);
		  CL << ']';
		}
	    }
	  if (logging)
	    CL << '\n';
	}

      // transitions to End state
      result_sc = -InfinityScore;
      for (int src = 0; src < states(); ++src)
	ScorePSumAcc (result_sc, ScorePMul (cell (cols, src), transition (src, End)));

      if (logging)
	{
	  CL << "End: ";
	  ShowScore (result_sc, CL);
	  CL << '\n';
	}
    }

  return result_sc;
}

Score Pair_transducer_scores::effective_trans_score (int src, int dest) const
{
  if (dest >= 0 && state_type[dest] == TransducerInsertType)
    return transition (src, dest);

  Score sc = -InfinityScore;
  for (int s = 0; s < states(); ++s)
    if (state_type[s] == TransducerWaitType)
      ScorePSumAcc (sc, ScorePMul (transition (src, s), transition (s, dest)));
  return sc;
}

void pfuncify_transition (PScores& pscores, PFunc& pfunc, Score sc, const sstring& src_name, const sstring& dest_name, const char* prefix)
{
  if (sc > -InfinityScore)
    {
      sstring trans_name;
      trans_name << prefix << src_name << PVAR_TRANSITION_CHAR << dest_name;
      const PGroup pg = pscores.new_group ("dummy_group_name", trans_name.c_str());
      const PVar pv (pg.group_idx, 0);
      pscores[pv] = sc;
      pfunc = pv;
    }
}

Pair_transducer_funcs::Pair_transducer_funcs (const Pair_transducer_scores& pair_trans_sc, PScores& pscores, const vector<sstring>& alphabet, const char* pvar_prefix)
  : Pair_transducer<PFunc> (pair_trans_sc, PFunc())
{
  for (int src = 0; src < states(); src++)
    {
      // transitions
      for (int dest = 0; dest < states(); dest++)
	pfuncify_transition (pscores, transition(src,dest), pair_trans_sc.transition(src,dest), state_name[src], state_name[dest], pvar_prefix);
      pfuncify_transition (pscores, start[src], pair_trans_sc.start[src], start_name, state_name[src], pvar_prefix);
      pfuncify_transition (pscores, end[src], pair_trans_sc.end[src], state_name[src], end_name, pvar_prefix);

      // emissions
      Alphabet_group pg;
      sstring pg_name;
      pg_name << pvar_prefix << (emit_label[src].size() ? emit_label[src] : state_name[src]);

      switch (state_type[src])
	{
	case TransducerMatchType:
	  emit_label[src] = pg_name;
	  pg = pscores.new_alphabet_group (alphabet, 2, pg_name.c_str(), false);
	  for (int xsym = 0; xsym < alphabet_size; ++xsym)
	    for (int ysym = 0; ysym < alphabet_size; ++ysym)
	      {
		vector<int> word (2);
		word[0] = xsym;
		word[1] = ysym;
		const PVar pv = pg[word];
		pscores[pv] = ((Pair_transducer_scores&)pair_trans_sc).match_val(src,xsym,ysym);  // cast away const
		match_val(src,xsym,ysym) = pv;
	      }
	  break;

	case TransducerDeleteType:
	  emit_label[src] = pg_name;
	  pg = pscores.new_alphabet_group (alphabet, 1, pg_name.c_str(), false);
	  for (int xsym = 0; xsym < alphabet_size; ++xsym)
	    {
	      const PVar pv = pg[xsym];
	      pscores[pv] = ((Pair_transducer_scores&)pair_trans_sc).delete_val(src,xsym);  // cast away const
	      delete_val(src,xsym) = pv;
	    }
	  break;

	case TransducerInsertType:
	  emit_label[src] = pg_name;
	  pg = pscores.new_alphabet_group (alphabet, 1, pg_name.c_str(), false);
	  for (int ysym = 0; ysym < alphabet_size; ++ysym)
	    {
	      const PVar pv = pg[ysym];
	      pscores[pv] = ((Pair_transducer_scores&)pair_trans_sc).insert_val(src,ysym);  // cast away const
	      insert_val(src,ysym) = pv;
	    }
	  break;

	case TransducerWaitType:
	  break;

	default:
	  THROWEXPR ("Unexpected state type");
	  break;
	}

      if (pg.group_idx >= 0)
	Transducer_SExpr_file::munge_group_suffix (pscores, pg, pg_name);
    }
  // start->end transition
  pfuncify_transition (pscores, start_to_end(), pair_trans_sc.start_to_end(), start_name, end_name, pvar_prefix);
}

Eliminated_EHMM_transducer_scores::Eliminated_EHMM_transducer_scores (const EHMM_transducer_scores& ehmm_trans_sc_ref,
								      const vector<int>& states_to_eliminate_ref)
  : EHMM_transducer_scores (ehmm_trans_sc_ref.states() - (int) states_to_eliminate_ref.size()),
    ehmm_trans_sc (&ehmm_trans_sc_ref),
    states_to_eliminate (&states_to_eliminate_ref),
    loopy_probs (Transition_methods::score2prob (ehmm_trans_sc_ref)),
    loop_exit(0,0.),  // will be re-initialized by Transition_methods::eliminate
    elim_probs (Transition_methods::eliminate (loopy_probs, states_to_eliminate_ref, false, &loop_exit)),
    elim_scores (Transition_methods::prob2score (elim_probs))
{
  // figure out states to keep
  const set<int> elim_set (states_to_eliminate->begin(), states_to_eliminate->end());
  for (int s = 0; s < ehmm_trans_sc->states(); ++s)
    if (elim_set.find (s) == elim_set.end())
      states_to_keep.push_back (s);
  if ((int) states_to_keep.size() != states())
    THROWEXPR ("Oops, miscalculated number of states to keep");

  // copy over elements from the original EHMM, one by one
  start_to_end() = elim_scores.start_to_end();
  for (int i = Grammar_state_enum::End; i <= Grammar_state_enum::Start; ++i)
    state2estate[i] = ((EHMM_transducer_scores*)ehmm_trans_sc)->state2estate[i];
  for (int i = 0; i < (int) states_to_keep.size(); ++i)
    {
      const int ei = states_to_keep[i];

      state_type[i] = ehmm_trans_sc->state_type[ei];
      state_name[i] = ehmm_trans_sc->state_name[ei];
      node_record[i] = ((map<int,sstring>&) ehmm_trans_sc->node_record)[ei];  // cast away const
      node_html[i] = ((map<int,sstring>&) ehmm_trans_sc->node_html)[ei];  // cast away const

      branch_trans_states[i] = ehmm_trans_sc->branch_trans_states[ei];
      emission[i] = ehmm_trans_sc->emission[ei];
      absorption[i] = ehmm_trans_sc->absorption[ei];
      state2estate[i] = ((EHMM_transducer_scores*)ehmm_trans_sc)->state2estate[ei];

      start[i] = elim_scores.start[ei];
      end[i] = elim_scores.end[ei];
      for (int j = 0; j < (int) states_to_keep.size(); ++j)
	{
	  const int ej = states_to_keep[j];
	  transition (i, j) = elim_scores.transition (ei, ej);
	}
    }

  start_name = ehmm_trans_sc->start_name;
  end_name = ehmm_trans_sc->end_name;

  etree = ehmm_trans_sc->etree;
  branch_transducer = ehmm_trans_sc->branch_transducer;
  alphabet_size = ehmm_trans_sc->alphabet_size;

  tape_name = ehmm_trans_sc->tape_name;

  // print some log messages
  if (states_to_eliminate->size()
      && CTAGGING(5,TRANSDUCER TRANSDUCER_ELIM TRANSDUCER_ELIM_STATES))
    {
      CL << "Eliminating the following \"null\" states:\n";
      for_const_contents (vector<int>, *states_to_eliminate, s)
	CL << ' ' << ehmm_trans_sc->get_state_name(*s) << '\n';
    }


  if (CTAGGING(0,TRANSDUCER TRANSDUCER_PRE_ELIM))
    {
      CL << "EHMM transducer scores, prior to null state elimination... ";
      ehmm_trans_sc->show (CL);
    }

  if (CTAGGING(-1,TRANSDUCER))
    {
      CL << "EHMM transducer probabilities, prior to null state elimination... ";
      loopy_probs.show_transitions (CL);
    }

  if (CTAGGING(-1,TRANSDUCER))
    {
      CL << "EHMM transducer probabilities, after null state elimination... ";
      elim_probs.show_transitions (CL);
    }

  if (CTAGGING(0,TRANSDUCER))
    {
      CL << "EHMM transducer scores, after null state elimination... ";
      elim_scores.show_transitions (CL);
    }

  if (CTAGGING(-2,TRANSDUCER TRANSDUCER_ELIM))
    {
      CL << "EHMM transducer scores, with null states eliminated & trimmed... ";
      show (CL);
    }

  if (CTAGGING(-2,TRANSDUCER TRANSDUCER_LOOP_EXIT))
    {
      CL << "loop_exit matrix following EHMM transducer null state elimination... ";
      loop_exit.show_transitions (CL);
    }
}

vector<int> Eliminated_EHMM_transducer_scores::sample_eliminated (const vector<int>& path, bool choose_ML_path)
{
  vector<int> orig_path;
  for_const_contents (vector<int>, path, s)
    orig_path.push_back (*s < 0 ? *s : states_to_keep[*s]);

  return Transition_methods::sample_eliminated (loopy_probs, elim_probs, *states_to_eliminate, orig_path, choose_ML_path);
}

Transducer_counts Eliminated_EHMM_transducer_scores::count_eliminated (const Transducer_counts& transition_counts) const
{
  if (CTAGGING(3,TRANSDUCER_COUNTS))
    {
      CL << "Printing composite transducer state names on entry to count_eliminated method...\n";
      for (int s = 0; s < ehmm_trans_sc->states(); ++s)
	CL << s << ": " << ehmm_trans_sc->state_name[s] << '\n';
      CL << "Eliminated states: " << *states_to_eliminate << '\n';
    }

  // re-map the state indices onto the original, loopy matrix
  Concrete_transition_counts elim_counts (loopy_probs.tm_states(), 0.);
  for (int i = 0; i < (int) states_to_keep.size(); ++i)
    {
      const int ei = states_to_keep[i];
      elim_counts.transition (Start, ei) = transition_counts.transition (Start, i);
      elim_counts.transition (ei, End) = transition_counts.transition (i, End);
      for (int j = 0; j < (int) states_to_keep.size(); ++j)
	{
	  const int ej = states_to_keep[j];
	  elim_counts.transition (ei, ej) = transition_counts.transition (i, j);
	}
    }
  elim_counts.transition (Start, End) = transition_counts.transition (Start, End);

  // log
  if (CTAGGING(3,TRANSDUCER_COUNTS))
    {
      CL << "Correcting the following counts for null state elimination:\n";
      elim_counts.show_transitions (CL);
    }

  // do the null state correction
  const Concrete_transition_counts corrected_counts = Transition_methods::count_eliminated (loopy_probs, elim_probs, loop_exit, *states_to_eliminate, states_to_keep, elim_counts);
  if (CTAGGING(3,TRANSDUCER_COUNTS))
    {
      CL << "Counts, corrected for null state elimination:\n";
      corrected_counts.show_transitions (CL);
    }

  // copy to Transducer_counts and return
  Transducer_counts loopy_counts (*ehmm_trans_sc, 0.);
  loopy_counts.assign_transition_matrix (corrected_counts);
  return loopy_counts;
}

void Transducer_methods::inc_emit_label_counts (const Pair_transducer<Prob>& trans_counts,
						const Pair_transducer<PFunc>& trans_funcs,
						PCounts& var_counts,
						const PScores& var_scores,
						const Prob weight)
{
  if (CTAGGING(3,TRANSDUCER_COUNTS))
    {
      CL << "Adding the following Pair_transducer_counts to PCounts:\n";
      trans_counts.show_pair_emit (CL);
    }

  // emissions
  for (int s = 0; s < trans_counts.states(); s++)
    switch (trans_counts.state_type[s])
      {
      case TransducerMatchType:
	for (int xsym = 0; xsym < trans_counts.alphabet_size; ++xsym)
	  for (int ysym = 0; ysym < trans_counts.alphabet_size; ++ysym)
	    ((Pair_transducer_funcs&)trans_funcs).match_val(s,xsym,ysym).inc_var_counts (var_counts, var_scores, ((Pair_transducer_counts&)trans_counts).match_val(s,xsym,ysym) * weight);
	break;

      case TransducerDeleteType:
	for (int xsym = 0; xsym < trans_counts.alphabet_size; ++xsym)
	  ((Pair_transducer_funcs&)trans_funcs).delete_val(s,xsym).inc_var_counts (var_counts, var_scores, ((Pair_transducer_counts&)trans_counts).delete_val(s,xsym) * weight);
	break;

      case TransducerInsertType:
	for (int ysym = 0; ysym < trans_counts.alphabet_size; ++ysym)
	  ((Pair_transducer_funcs&)trans_funcs).insert_val(s,ysym).inc_var_counts (var_counts, var_scores, ((Pair_transducer_counts&)trans_counts).insert_val(s,ysym) * weight);
	break;

      case TransducerWaitType:
	break;

      default:
	THROWEXPR ("Unexpected state type");
	break;
      }
}

void Transducer_methods::inc_transition_counts (const Transducer<Prob>& trans_count,
						const Transducer<PFunc>& trans_func,
						PCounts& var_counts,
						const PScores& var_scores,
						const Prob weight)
{
  if (CTAGGING(3,TRANSDUCER_COUNTS))
    {
      CL << "Adding the following transition matrix counts to PCounts:\n";
      const Concrete_transition_counts concrete_trans_count (trans_count);
      concrete_trans_count.show_transitions (CL);
    }

  for (int src = 0; src < trans_count.states(); src++)
    {
      for (int dest = 0; dest < trans_count.states(); dest++)
	trans_func.transition(src,dest).inc_var_counts (var_counts, var_scores, trans_count.transition(src,dest) * weight);
      trans_func.start[src].inc_var_counts (var_counts, var_scores, trans_count.start[src] * weight);
      trans_func.end[src].inc_var_counts (var_counts, var_scores, trans_count.end[src] * weight);
    }
  trans_func.start_to_end().inc_var_counts (var_counts, var_scores, trans_count.start_to_end() * weight);
}

void Pair_transducer_scores::sample (const Digitized_biosequence& parent_dsq, vector<int>& parent_child_path, Digitized_biosequence& child_dsq, int max_tries)
{
  int n_tries;
  for (n_tries = 0; true; ++n_tries)
    {
      parent_child_path.clear();
      child_dsq.clear();

      if (n_tries >= max_tries)
	{
	  CLOGERR << "Tried " << n_tries << " times to feed parent sequence through transducer, unsuccessfully; giving up.\n";
	  break;
	}

      int current_state = Grammar_state_enum::Start;
      int next_parent_pos = 0;
      bool failed = false;
      while (!failed)
	{
	  parent_child_path.push_back (current_state);
	  if (current_state == Grammar_state_enum::End)
	    break;

	  const State_type t = current_state == Grammar_state_enum::Start ? TransducerStartType : state_type[current_state];

	  // handle emissions
	  switch (t)
	    {

	    case TransducerInsertType:
	      {
		vector<Prob> ins_prob;
		Prob total = 0.;
		for (int child_sym = 0; child_sym < alphabet_size; ++child_sym)
		  {
		    const Prob psym = Score2Prob (insert_val (current_state, child_sym));
		    ins_prob.push_back (psym);
		    total += psym;
		  }
		if (total > 0.)
		  child_dsq.push_back (Rnd::choose (ins_prob));
		else
		  failed = true;
		break;
	      }

	    case TransducerMatchType:
	      {
		const int parent_sym = parent_dsq[next_parent_pos++];
		vector<Prob> match_prob;
		Prob total = 0.;
		for (int child_sym = 0; child_sym < alphabet_size; ++child_sym)
		  {
		    const Prob psym = Score2Prob (match_val (current_state, parent_sym, child_sym));
		    match_prob.push_back (psym);
		    total += psym;
		  }
		if (total > 0.)
		  child_dsq.push_back (Rnd::choose (match_prob));
		else
		  failed = true;
		break;
	      }

	    case TransducerDeleteType:
	      {
		const int parent_sym = parent_dsq[next_parent_pos++];
		failed = !Rnd::decide (Score2Prob (delete_val (current_state, parent_sym)));
		break;
	      }

	    default:
	      break;

	    }

	  // bail now if failed
	  if (failed)
	    continue;

	  // handle transitions
	  if (t == TransducerWaitType && next_parent_pos == (int) parent_dsq.size())
	    {
	      failed = !Rnd::decide (Score2Prob (transition (current_state, Grammar_state_enum::End)));
	      current_state = Grammar_state_enum::End;
	    }
	  else
	    {
	      vector<Prob> next_state_prob;
	      Prob total = 0.;
	      for (int next_state = 0; next_state < states(); ++next_state)
		{
		  const Prob ptrans = Score2Prob (transition (current_state, next_state));
		  next_state_prob.push_back (ptrans);
		  total += ptrans;
		}

	      if (total > 0.)
		current_state = Rnd::choose (next_state_prob);
	      else
		failed = true;
	    }
	}

      if (failed)
	CLOGERR << "Failed to sample path through transducer (attempt #" << n_tries+1 << ")\n";
      else
	break;
    }
}

Pair_HMM_scores Pair_transducer_scores::pair_hmm (const Alphabet& alph) const
{
  Pair_HMM_scores hmm (states(), &alph);
  hmm.assign_transition_matrix (*this);
  for (int s = 0; s < states(); ++s)
    switch (state_type[s])
      {
      case TransducerMatchType:
	hmm.init_emit (s, Pair_HMM_state_type_enum::EmitXY, -InfinityScore);
	for (int x = 0; x < alphabet_size; ++x)
	  for (int y = 0; y < alphabet_size; ++y)
	    hmm.pair_emit[s](x,y) = ((Pair_transducer_scores&)*this).match_val(s,x,y);
	break;

      case TransducerDeleteType:
	hmm.init_emit (s, Pair_HMM_state_type_enum::EmitX, -InfinityScore);
	for (int x = 0; x < alphabet_size; ++x)
	  hmm.single_emit[s][x] = ((Pair_transducer_scores&)*this).delete_val(s,x);
	break;

      case TransducerInsertType:
	hmm.init_emit (s, Pair_HMM_state_type_enum::EmitY, -InfinityScore);
	for (int y = 0; y < alphabet_size; ++y)
	  hmm.single_emit[s][y] = ((Pair_transducer_scores&)*this).insert_val(s,y);
	break;

      case TransducerWaitType:
	hmm.state_type[s] = Pair_HMM_state_type_enum::Null;
      default:
	break;
      }

  return hmm;
}
