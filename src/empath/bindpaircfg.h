#ifndef BIND_PAIR_CFG_INCLUDED
#define BIND_PAIR_CFG_INCLUDED

#include "scfg/paircfg.h"
#include "scfg/energy21.h"

struct Bind_pair_CFG : Pair_CFG_scores, Duplex_lookup
{
  enum { xPad5 = 0, yPad3 = 1, xPad3 = 2, yPad5 = 3,
	 xyDuplex = 4, xyIntBeg = 10, xyIntEnd = 16, xyIntLoop = 22,
	 TotalStates = 23 };
  set<int> mask_states;
  Bind_pair_CFG (const RNA_energy& energy);
};

#endif /* BIND_PAIR_CFG_INCLUDED */
