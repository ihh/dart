#include "stemloc/sldefaults.h"
#include "indiegram/igdefaults.h"
#include "indiegram/tripletscfg.h"
#include "indiegram/tkfst.h"

void Indiegram_defaults::init_superstem_params (Telegraph_PScores_adaptor& fold_tgio)
{
  Stemloc_defaults::init_superstem_params (fold_tgio);
}

Triplet_SCFG Indiegram_defaults::init_tkfst_scfg (double t, double u, double v)
{
  // create a TKFST_Triplet_SCFG object and return it
  TKFST_Triplet_SCFG scfg (t, u, v);
  return scfg;
}

// taken straight from testtripletdp.cc
Triplet_SCFG Indiegram_defaults::init_test_scfg()
{

  // Create a toy grammar.
  // It has "stem" emit states which emit only WC-paired nucleotides as well as "loop" emit states.
  int num_null_states = 1;
  int num_emit_states = 63;
  int num_bifurc_states = 1;
  int num_states = num_null_states + num_emit_states + num_bifurc_states;

  Triplet_SCFG scfg (num_states);
  scfg.transition_scores.start_to_end() = 9999;

  // initialize states: 1 state for each state type
  for (int i = 0; i < num_states; ++i) {

    // initialize state
    const State_type t = static_cast<State_type> (i);
    scfg.init_emit (i, t); // 1 state for each state type, so let 'state' = 'type'

    // set transition scores unfavorably (we'll make some favorable later)
    scfg.transition_scores.start[i] = -999;
    scfg.transition_scores.end[i] = -999;
    for (int j = 0; j < num_states; ++j)
      scfg.transition_scores.transition (i, j) = -999;

    // if an emit state, initialize scores for emissions
    if (SCFG_state_typing::is_emit_type (t))
      {
	for (int xlc = 0; xlc < SCFG_alphabet_size; ++xlc)
	  for (int xrc = 0; xrc < SCFG_alphabet_size; ++xrc)
	    for (int ylc = 0; ylc < SCFG_alphabet_size; ++ylc)
	      for (int yrc = 0; yrc < SCFG_alphabet_size; ++yrc)
		for (int zlc = 0; zlc < SCFG_alphabet_size; ++zlc)
		  for (int zrc = 0; zrc < SCFG_alphabet_size; ++zrc)
		    {
		      Score emit_score = 0;
		      // Now we need to assign a score to each possible emission from every emit state
		      // We have designed our test sequences and folds so that paired emissions are always WC-basepaired,
		      // and so it's sufficient to assign scores only to basepaired paired emissions.
		      // Note that this does properly assign scores to single (unpaired) emissions as well.

		      // if paired emission
		      const bool xpaired = (t & SCFG_state_typing::EmitXL) && (t & SCFG_state_typing::EmitXR);
		      const bool ypaired = (t & SCFG_state_typing::EmitYL) && (t & SCFG_state_typing::EmitYR);
		      const bool zpaired = (t & SCFG_state_typing::EmitZL) && (t & SCFG_state_typing::EmitZR);
		      // if unpaired emissions (to left or right) to > 1 sequence
		      const bool xyl = (t & SCFG_state_typing::EmitXL) && (t & SCFG_state_typing::EmitYL);
		      const bool xzl = (t & SCFG_state_typing::EmitXL) && (t & SCFG_state_typing::EmitZL);
		      const bool yzl = (t & SCFG_state_typing::EmitYL) && (t & SCFG_state_typing::EmitZL);
		      const bool xyr = (t & SCFG_state_typing::EmitXR) && (t & SCFG_state_typing::EmitYR);
		      const bool xzr = (t & SCFG_state_typing::EmitXR) && (t & SCFG_state_typing::EmitZR);
		      const bool yzr = (t & SCFG_state_typing::EmitYR) && (t & SCFG_state_typing::EmitZR);

		      // disallowed paired emissions which aren't WC-paired
		      bool notwc = 0;
		      if (xpaired && (xlc+xrc != 3)) notwc = 1;
		      if (ypaired && (ylc+yrc != 3)) notwc = 1;
		      if (zpaired && (zlc+zrc != 3)) notwc = 1;
		      if (notwc)
			emit_score = -InfinityScore;

		      // disallow mismatches for paired emissions
		      bool mismatch = 0;
		      if ((xpaired && ypaired) && (xlc != ylc || xrc != yrc)) mismatch = 1;
		      if ((xpaired && zpaired) && (xlc != zlc || xrc != zrc)) mismatch = 1;
		      if ((ypaired && zpaired) && (ylc != zlc || yrc != zrc)) mismatch = 1;

		      // disallow mismatches for single emissions
		      if (xyl && !(xpaired || ypaired) && (xlc != ylc)) mismatch = 1;
		      if (xzl && !(xpaired || zpaired) && (xlc != zlc)) mismatch = 1;
		      if (yzl && !(ypaired || zpaired) && (ylc != zlc)) mismatch = 1;
		      if (xyr && !(xpaired || ypaired) && (xrc != yrc)) mismatch = 1;
		      if (xzr && !(xpaired || zpaired) && (xrc != zrc)) mismatch = 1;
		      if (yzr && !(ypaired || zpaired) && (yrc != zrc)) mismatch = 1;

		      if (mismatch)
			emit_score = -InfinityScore;

		      // randomized emission scores for the rest
		      else
			{
			  // penalty for mismatches
			  int penalty = 0;
			  if (xyl && (xlc != ylc)) ++penalty;
			  if (xzl && (xlc != zlc)) ++penalty;
			  if (yzl && (ylc != zlc)) ++penalty;
			  if (xyr && (xrc != yrc)) ++penalty;
			  if (xzr && (xrc != zrc)) ++penalty;
			  if (yzr && (yrc != zrc)) ++penalty;
			  // prefer paired emissions over unpaired emissions
			  if (xpaired || ypaired || zpaired)
			    emit_score = -100 * penalty;
			  else
			    emit_score = -1000 * penalty;
			}

		      // now set the score
		      const int emit_idx = 
			xlc * scfg.emit_xl_mul (t) + 
			xrc * scfg.emit_xr_mul (t) + 
			ylc * scfg.emit_yl_mul (t) + 
			yrc * scfg.emit_yr_mul (t) + 
			zlc * scfg.emit_zl_mul (t) + 
			zrc * scfg.emit_zr_mul (t);
		      scfg.emit[i][emit_idx] = emit_score;
		    }
      }

  }
  // deal with the bifurcation state
  scfg.init_bifurc (64, 0, 0);

  // set some transition scores by hand to favor the parse we want
  // (see state_array[], defined below)
  scfg.transition_scores.start[63] = 0; // root
  scfg.transition_scores.transition (63,63) = 0;
  scfg.transition_scores.transition (63,32) = 0;
  scfg.transition_scores.transition (32,32) = 0;
  scfg.transition_scores.transition (32,64) = 0;
  scfg.transition_scores.transition (0,55) = 0; // root->left
  scfg.transition_scores.transition (55,55) = 0;
  scfg.transition_scores.transition (55,63) = 0;
  scfg.transition_scores.end[63] = 0;
  scfg.transition_scores.transition (0,12) = 0; // root->right
  scfg.transition_scores.transition (12,12) = 0;
  scfg.transition_scores.end[12] = 0;

  return scfg;

}
