#include "stemloc/tetrastem.h"

Tetra_stem::Tetra_stem (Score pair_score, Score mispair_score, int min_loop_len, int max_loop_len,
			bool allow_GNRA, bool allow_UNCG, bool allow_XUUY,
			bool allow_UYU, bool allow_ARU)
  : Pair_CFG_scores (LoopStart + max_loop_len)
{
  // GNRA tetraloops
  if (allow_GNRA)
    {
      init_emit (GNRA1, EmitXLR, -InfinityScore);
      init_emit (GNRA2, EmitXLR, -InfinityScore);
      emit_by_chars (GNRA1, "ga") = 0;
      for (int i = 0; i < 4; ++i)
	{
	  sstring na = "xa";
	  sstring ng = "xg";
	  na[0] = ng[0] = RNA_alphabet.int2char(i);
	  emit_by_chars (GNRA2, na.c_str()) = 0;
	  emit_by_chars (GNRA2, ng.c_str()) = 0;
	}
      transition (Stem, GNRA1) = 0;
      transition (GNRA1, GNRA2) = 0;
      transition (GNRA2, End) = 0;
    }
  else
    {
      init_emit (GNRA1, Null, -InfinityScore);
      init_emit (GNRA2, Null, -InfinityScore);
    }
  // UNCG tetraloops
  if (allow_UNCG)
    {
      init_emit (UNCG1, EmitXLR, -InfinityScore);
      init_emit (UNCG2, EmitXLR, -InfinityScore);
      emit_by_chars (UNCG1, "ug") = 0;
      for (int i = 0; i < 4; ++i)
	{
	  sstring nc = "xc";
	  nc[0] = RNA_alphabet.int2char(i);
	  emit_by_chars (UNCG2, nc.c_str()) = 0;
	}
      transition (Stem, UNCG1) = 0;
      transition (UNCG1, UNCG2) = 0;
      transition (UNCG2, End) = 0;
    }
  else
    {
      init_emit (UNCG1, Null, -InfinityScore);
      init_emit (UNCG2, Null, -InfinityScore);
    }
  // XUUY tetraloops
  if (allow_XUUY)
    {
      init_emit (XUUY1, EmitXLR, -InfinityScore);
      init_emit (XUUY2, EmitXLR, -InfinityScore);
      emit_by_chars (XUUY1, "ua") = 0;
      emit_by_chars (XUUY1, "cg") = 0;
      emit_by_chars (XUUY1, "au") = 0;
      emit_by_chars (XUUY2, "uu") = 0;
      transition (Stem, XUUY1) = 0;
      transition (XUUY1, XUUY2) = 0;
      transition (XUUY2, End) = 0;
    }
  else
    {
      init_emit (XUUY1, Null, -InfinityScore);
      init_emit (XUUY2, Null, -InfinityScore);
    }
  // UYU triloops
  if (allow_UYU)
    {
      init_emit (UYU1, EmitXLR, -InfinityScore);
      init_emit (UYU2, EmitXL, -InfinityScore);
      emit_by_chars (UYU1, "uu") = 0;
      emit_by_chars (UYU2, "c") = 0;
      emit_by_chars (UYU2, "u") = 0;
      transition (Stem, UYU1) = 0;
      transition (UYU1, UYU2) = 0;
      transition (UYU2, End) = 0;
    }
  else
    {
      init_emit (UYU1, Null, -InfinityScore);
      init_emit (UYU2, Null, -InfinityScore);
    }
  // ARU triloops
  if (allow_ARU)
    {
      init_emit (ARU1, EmitXLR, -InfinityScore);
      init_emit (ARU2, EmitXL, -InfinityScore);
      emit_by_chars (ARU1, "au") = 0;
      emit_by_chars (ARU2, "a") = 0;
      emit_by_chars (ARU2, "g") = 0;
      transition (Stem, ARU1) = 0;
      transition (ARU1, ARU2) = 0;
      transition (ARU2, End) = 0;
    }
  else
    {
      init_emit (ARU1, Null, -InfinityScore);
      init_emit (ARU2, Null, -InfinityScore);
    }
  // generic loops
  if (max_loop_len >= min_loop_len)
    {
      const int LoopEnd = LoopStart + max_loop_len - 1;
      for (int s = LoopStart; s <= LoopEnd; ++s)
	init_emit (s, EmitXL, 0);
      for (int s = LoopStart; s < LoopEnd; ++s)
	transition (s, s+1) = 0;
      transition (Stem, LoopStart) = 0;
      for (int s = LoopStart + max(min_loop_len-1,0); s <= LoopEnd; ++s)
	transition (s, End) = 0;
    }
  else
    for (int s = LoopStart; s < LoopStart + max_loop_len; ++s)
      init_emit (s, Null, -InfinityScore);
  // stem
  init_emit (Stem, EmitXLR, mispair_score);
  emit_by_chars (Stem, "au") = pair_score;
  emit_by_chars (Stem, "cg") = pair_score;
  emit_by_chars (Stem, "gc") = pair_score;
  emit_by_chars (Stem, "ua") = pair_score;
  emit_by_chars (Stem, "gu") = pair_score;
  emit_by_chars (Stem, "ug") = pair_score;
  // set transitions
  transition (Start, Stem) = 0;
  transition (Stem, Stem) = 0;
  // set the list of paired states
  paired_states.insert (Stem);
}
