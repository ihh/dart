#ifndef P53_INCLUDED
#define P53_INCLUDED

#include "hmm/singlehmm.h"
#include "empath/trainer.h"
#include "seq/dirichlet.h"
#include "seq/gff.h"

// A P53 binding site motif: double palindrome with spacer

struct P53_model : Trainable
{
  int unit_sz;   // 5 for P53 binding site; actual size is 4*unit_sz + spacer
  int max_spacer;

  vector<Alphabet_group> match_emit;  // the "unit" motif
  PGroup                 space_len;   // the amount of space between the two palindromes

  // constructor
  P53_model (int unit_sz, int max_spacer, PScores& pscore);

  // methods to return state indices
  int match_state (int unit, int pos) const { return unit * unit_sz + (unit >= 2 ? max_spacer : 0) + (unit & 1 ? unit_sz-1-pos : pos); }
  int spacer_state (int pos) const { return 2 * unit_sz + pos; }
};

#endif
