#ifndef SUPERPAIR_INCLUDED
#define SUPERPAIR_INCLUDED

#include "scfg/pairpcfg.h"
#include "seq/dirichlet.h"

class Super_pair : public Odds_PCFG
{
public:
  // the PScores object
  PScores* pscore;

  // state indices
  enum { Loop = 0, LoopMatch = 1, LoopIns = 2, LoopDel = 3,
	 LBulge = 4, LBulgeMatch = 5, LBulgeIns = 6, LBulgeDel = 7,
	 RBulge = 8, RBulgeMatch = 9, RBulgeIns = 10, RBulgeDel = 11,
	 LRBulge = 12, LRBulgeMatch = 13, LRBulgeIns = 14, LRBulgeDel = 15,
	 Stem = 16, StemMatch = 17, StemIns = 18, StemDel = 19, StemEnd = 20,
	 Multi = 21, LMulti = 22, RMulti = 23,
	 TotalStates = 24 };

  // parameters
  Boolean_group start_in_stem;
  Boolean_group loop_extend;
  Boolean_group loop_gap_open;
  Boolean_group loop_gap_extend;
  Boolean_group loop_gap_swap;  // e.g. LoopIns->LoopDel
  Boolean_group stem_extend;
  Boolean_group stem_gap_open;
  Boolean_group stem_gap_extend;
  Boolean_group stem_gap_swap;  // e.g. StemIns->StemDel
  Boolean_group multi_extend;
  Boolean_group multi_bulge_open;
  PGroup        post_stem;   // StemEnd -> Loop | (LBulge|RBulge) | LRBulge | Multi

  // emit parameters
  Alphabet_group         single_nuc;    // 4 vars
  vector<Alphabet_group> pair_nuc;      // 4*4 vars
  Alphabet_group         single_dinuc;  // 16 vars
  vector<Alphabet_group> pair_dinuc;    // 16*16 vars

  // constructors
  Super_pair (PScores* pscore);
  Super_pair();

  // builder method for initialising single-stranded regions
  void init_unpaired (int start, int match, int ins, int del, int end, State_type xemit, State_type yemit);

  // method to generate prior
  Dirichlet_prior default_prior() const;
};

#endif /* SUPERPAIR_INCLUDED */
