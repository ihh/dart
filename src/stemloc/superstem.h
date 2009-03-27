#ifndef SUPERSTEM_INCLUDED
#define SUPERSTEM_INCLUDED

#include "scfg/pairpcfg.h"
#include "seq/dirichlet.h"

class Super_stem : public Odds_PCFG
{
public:
  // the PScores object
  PScores& pscore;

  // state indices
  enum { Loop = 0, LBulge = 1, RBulge = 2, LRBulge = 3,
	 StemStart = 4, Stem = 5,  // 16 states
	 StemEnd = 21, Multi = 22, LMulti = 23, RMulti = 24,
	 TotalStates = 25 };

  // parameters
  Boolean_group start_in_stem;
  Boolean_group loop_extend;
  Boolean_group bulge_extend;
  Boolean_group stem_extend;
  Boolean_group multi_extend;
  Boolean_group multi_bulge_open;
  PGroup        post_stem;   // StemEnd -> Loop | (LBulge|RBulge) | LRBulge | Multi

  // emit parameters
  Alphabet_group         nuc;            // 4 vars
  Alphabet_group         opening_dinuc;  // 16 vars
  vector<Alphabet_group> stacked_dinuc;  // 16 elements, each has 16 vars

  // constructor
  Super_stem (PScores& pscore);

  // method to generate prior
  Dirichlet_prior default_prior() const;
};

#endif /* SUPERSTEM_INCLUDED */
