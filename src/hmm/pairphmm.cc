#include "hmm/pairphmm.h"
#include "util/vector_output.h"

Pair_PHMM::Pair_PHMM (Pair_HMM<PFunc>& hmm)
  : Pair_HMM<PFunc> (hmm),
    group_suffix (0)
{ }    

Pair_PHMM::Pair_PHMM (int states, const Alphabet& alphabet) :
  Pair_HMM<PFunc> (states, PFunc(0.0), &alphabet),
  group_suffix (0)
{ }

Pair_PHMM::Pair_PHMM (int states, const Alphabet& alphabet, const vector<vector<sstring> >& group_suffix) :
  Pair_HMM<PFunc> (states, PFunc(0.0), &alphabet),
  group_suffix (&group_suffix)
{ }

void Pair_PHMM::set_pscope (const PScope& ps)
{
  pscope = &ps;
}

void Pair_PHMM::set_group_suffix (const vector<vector<sstring> >& g_suffix)
{
  group_suffix = &g_suffix;
}

Pair_HMM_scores Pair_PHMM::eval_hmm_scores_gs (const PScores& var_scores) const
{
  ((Pair_PHMM*) this) -> set_group_suffix (var_scores.group_suffix);
  return eval_hmm_scores (var_scores);
}

Pair_HMM_scores Pair_PHMM::eval_hmm_scores (const PScores& var_scores) const
{
  if (CTAGGING(1,HMM_PSCORES))
    {
      CL << "PScores:\n";
      var_scores.show (CL);
    }
  if (CTAGGING(3,HMM_PFUNC))
    {
      CL << "Evaluating Pair_HMM_scores from PFunc's:\n";
      show(CL);
      if (CTAGGING(2,HMM_PFUNC)) var_scores.show(CL);
    }
  Pair_HMM_scores hmm_scores (states(), alphabet);
  hmm_scores.name = name;
  hmm_scores.state_name = state_name;
  for (int src = 0; src < states(); src++)
    {
      hmm_scores.state_type[src] = state_type[src];
      for (int dest = 0; dest < states(); dest++)
	hmm_scores.transition(src,dest) = transition(src,dest).eval_sc(var_scores);
      switch (state_type[src])
	{
	case Null:
	  break;
	case EmitX:
	case EmitY:
	  hmm_scores.single_emit[src] = vector<Score> (single_emit[src].size());
	  for (int i = 0; i < (int) single_emit[src].size(); ++i)
	    hmm_scores.single_emit[src][i] = single_emit[src][i].eval_sc(var_scores);
	  break;
	case EmitXY:
	  hmm_scores.pair_emit[src] = array2d<Score> (pair_emit[src].xsize(), pair_emit[src].ysize());
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    for (int j = 0; j < pair_emit[src].ysize(); ++j)
	      hmm_scores.pair_emit[src](i,j) = pair_emit[src](i,j).eval_sc(var_scores);
	  break;
	default:
	  THROWEXPR ("State type " << state_type[src] << " undefined\n");
	  break;
	}
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

void Pair_PHMM::inc_var_counts (PCounts& var_counts,
				const PScores& var_scores,
				const Pair_HMM_counts& hmm_counts,
				const Prob model_count) const
{
  if (!same_dimensions (hmm_counts))
    THROWEXPR ("Pair_PHMM and Pair_HMM_counts: not same dimensions");
  if (CTAGGING(3,HMM_PCOUNTS))
    {
      CL << "Adding the following Pair_HMM_counts to PCounts:\n";
      hmm_counts.show(CL);
    }
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest).inc_var_counts (var_counts, var_scores, hmm_counts.transition(src,dest) * model_count);
      switch (state_type[src])
	{
	case Null:
	  break;
	case EmitX:
	case EmitY:
	  for (int i = 0; i < (int) single_emit[src].size(); ++i)
	    single_emit[src][i].inc_var_counts (var_counts, var_scores, hmm_counts.single_emit[src][i] * model_count);
	  break;
	case EmitXY:
	  for (int i = 0; i < pair_emit[src].xsize(); ++i)
	    for (int j = 0; j < pair_emit[src].ysize(); ++j)
	      pair_emit[src](i,j).inc_var_counts (var_counts, var_scores, hmm_counts.pair_emit[src](i,j) * model_count);
	  break;
	default:
	  THROWEXPR ("State type " << state_type[src] << " undefined\n");
	  break;
	}
      start[src].inc_var_counts (var_counts, var_scores, hmm_counts.start[src] * model_count);
      end[src].inc_var_counts (var_counts, var_scores, hmm_counts.end[src] * model_count);
    }
  start_to_end().inc_var_counts (var_counts, var_scores, hmm_counts.start_to_end() * model_count);
}

void Pair_PHMM::factor_in_null_model (const Alphabet_group& null_emit)
{
  if (alphabet == 0) THROWEXPR ("Null alphabet in Pair_PHMM");
  const int A = alphabet->size();
  for (int s = 0; s < states(); ++s)
    {
      const State_type t = state_type[s];
      if (t == Null) continue;
      const int A_xl = t & EmitX ? A : 1;
      const int A_yl = t & EmitY ? A : 1;
      for (int xl = 0; xl < A_xl; ++xl)
	for (int yl = 0; yl < A_yl; ++yl)
	  {
	    PFunc& e =
	      t == EmitX ? single_emit[s][xl] :
	      (t == EmitY ? single_emit[s][yl] :
	       pair_emit[s](xl,yl));
	    
	    if (t & EmitX) e /= null_emit[xl];
	    if (t & EmitY) e /= null_emit[yl];
	  }
    }
}
