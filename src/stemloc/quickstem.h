#ifndef QUICKSTEM_INCLUDED
#define QUICKSTEM_INCLUDED

#include "scfg/pairpcfg.h"
#include "seq/dirichlet.h"

class Quick_stem : public Pair_PCFG
{
public:
  // the PScores object
  PScores& pscore;

  // state indices
  enum { loop_states = 4 };
  enum { PrePad = 0, PadL = 1, StemBif = 2, PreStem = 3,
	 StemLR = 4, SymL = 5, SymR = 6, AsymL = 7, AsymR = 8,
	 Loop = 9 };

  // secondary structure parameters
  Boolean_group pad_extend_pg;
  Boolean_group stem_prior_pg;
  Boolean_group stem_extend_pg;
  Boolean_group stem_bifurc_pg;
  Boolean_group sym_bulge_open_pg;
  Boolean_group sym_bulge_extend_pg;
  Boolean_group asym_bulge_open_pg;
  Boolean_group asym_bulge_extend_pg;
  Boolean_group loop_extend_pg;
  PGroup loop_start_pos;

  // emit parameters
  const Alphabet_group& null_emit;
  Alphabet_group        dinuc;  // 16 vars

  // state accessors
  static inline int loop (int i) { return Loop + i; }

  // constructor
  Quick_stem (PScores& pscore, const Alphabet_group& null_emit);

  // method to generate prior
  Dirichlet_prior default_prior() const;
};

#endif /* QUICKSTEM_INCLUDED */
