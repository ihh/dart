#ifndef TETRASTEM_INCLUDED
#define TETRASTEM_INCLUDED

#include "scfg/paircfg.h"

struct Tetra_stem : public Pair_CFG_scores
{
  // state enums
  enum { GNRA1 = 0, GNRA2 = 1,
	 UNCG1 = 2, UNCG2 = 3,
	 XUUY1 = 4, XUUY2 = 5,
	 UYU1 = 6, UYU2 = 7,
	 ARU1 = 8, ARU2 = 9,
	 Stem = 10, LoopStart = 11 };
  // set of paired states
  set<int> paired_states;
  // constructor
  // NB if max_loop_len < 5, then only named tetraloops & triloops will be allowed
  Tetra_stem (Score pair_score, Score mispair_score, int min_loop_len, int max_loop_len,
	      bool allow_GNRA = TRUE, bool allow_UNCG = TRUE, bool allow_XUUY = TRUE,
	      bool allow_UYU = TRUE, bool allow_ARU = TRUE);
};

#endif /* TETRASTEM_INCLUDED */
