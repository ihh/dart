#ifndef STEMPAIR_INCLUDED
#define STEMPAIR_INCLUDED

#include "scfg/pairpcfg.h"
#include "seq/dirichlet.h"

class Stem_pair : public Pair_PCFG
{
public:
  // the PScores object
  PScores& pscore;

  // state indices
  enum { StemMatLR = 0, StemMatL = 1, StemMatR = 2,
	 StemDelLR = 3, StemDelL = 4, StemDelR = 5,
	 StemInsLR = 6, StemInsL = 7, StemInsR = 8,
	 LoopMat = 9, LoopDel = 10, LoopIns = 11,
         Bif = 12, BifStart = 13,
         StatesWithBif = 14, StatesWithoutBif = 12 };

  // underlying secondary structure parameters
  PGroup stem_LR_pg;   // LR   -> L | LR | B
  PGroup stem_L_pg;    // L    -> L | R | LR
  PGroup stem_pg;      // stem -> stem | loop
  PGroup loop_pg;      // loop -> loop | end

  // indel parameters
  PGroup stem_mat_pg;  // smat -> smat | sins
  PGroup stem_ins_pg;  // sins -> smat | sins
  PGroup loop_mat_pg;  // lmat -> lmat | lins
  PGroup loop_ins_pg;  // lins -> lmat | lins

  // emit parameters
  Alphabet_group         single_nuc;    // 4 vars
  vector<Alphabet_group> pair_nuc;      // 4*4 vars
  Alphabet_group         single_dinuc;  // 16 vars
  vector<Alphabet_group> pair_dinuc;    // 16*16 vars

  // null model
  Boolean_group          null_extend;
  Alphabet_group         null_emit;     // 4 vars

  // parameter accessors
  PVar stem_LR_to_L()     const { return stem_LR_pg[0]; }
  PVar stem_LR_to_LR()    const { return stem_LR_pg[1]; }
  PVar stem_LR_to_B()     const { return stem_LR_pg[2]; }

  PVar stem_L_to_L()      const { return stem_L_pg[0]; }
  PVar stem_L_to_R()      const { return stem_L_pg[1]; }
  PVar stem_L_to_LR()     const { return stem_L_pg[2]; }

  PVar stem_mat_to_mat()  const { return stem_mat_pg[0]; }
  PVar stem_mat_to_ins()  const { return stem_mat_pg[1]; }

  PVar stem_ins_to_mat()  const { return stem_ins_pg[0]; }
  PVar stem_ins_to_ins()  const { return stem_ins_pg[1]; }

  PVar loop_mat_to_mat()  const { return loop_mat_pg[0]; }
  PVar loop_mat_to_ins()  const { return loop_mat_pg[1]; }

  PVar loop_ins_to_mat()  const { return loop_ins_pg[0]; }
  PVar loop_ins_to_ins()  const { return loop_ins_pg[1]; }

  PVar stem_to_stem()     const { return stem_pg[0]; }
  PVar stem_to_loop()     const { return stem_pg[1]; }

  PVar loop_to_loop()     const { return loop_pg[0]; }
  PVar loop_to_end()      const { return loop_pg[1]; }

  // constructor
  Stem_pair (PScores& pscore, bool use_bif = FALSE);

  // method to generate prior
  Dirichlet_prior default_prior() const;
};

#endif /* STEMPAIR_INCLUDED */
