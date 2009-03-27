#ifndef SEQLOGO_INCLUDED
#define SEQLOGO_INCLUDED

#include "seq/biosequence.h"

struct Sequence_logo
{
  vector<vector<Symbol_weight> > sort_prof_w;
  const Alphabet& alphabet;

  int text_height;

  Sequence_logo (const Score_profile& prof_sc, const Alphabet& alphabet);

  void show (ostream& o) const;
};

#endif
