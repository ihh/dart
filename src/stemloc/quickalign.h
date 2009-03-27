#ifndef QUICKALIGN_INCLUDED
#define QUICKALIGN_INCLUDED

#include "hmm/pairphmm.h"
#include "seq/dirichlet.h"

// simple, trainable pair HMM for finding local homology using CFG_alphabet
// similar to the "ABC" model -- tries to maintain diagonals
class Quick_align : public Pair_PHMM
{
public:
  // the PScores object
  PScores* pscore;

  // local flag
  bool local;

  // state indices
  enum { LeftPadX = 0, LeftPadY = 1,
	 Match = 2, InsX = 3, InsY = 4,
	 Gap = 5, GapX = 6, GapY = 7,
	 RightPadX = 8, RightPadY = 9,
	 TotalStates = 10 };

  // global/local flag
  Boolean_group is_global;

  // indel parameters
  PGroup mat_pg;  // Match -> Match | Ins | Gap
  PGroup ins_pg;  // InsX -> Match | InsX | InsY
  PGroup gap_pg;  // Gap -> Gap | GapX | Match
  PGroup gapx_pg;  // GapX -> GapX | Match

  // mutable pgroups
  set<int> mutable_pgroups;

  // emit parameters
  const Alphabet_group*  null_emit;     // 4 vars
  vector<Alphabet_group> pair_nuc;      // 4*4 vars

  // constructor
  Quick_align (PScores& pscore, const Alphabet_group& null_emit);

  // method to return set of non-pad states
  set<int> non_pad_states() const;

  // method to generate prior
  Dirichlet_prior default_prior() const;
};

#endif /* QUICKALIGN_INCLUDED */
