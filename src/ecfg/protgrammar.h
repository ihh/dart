#ifndef PROTEIN_GRAMMAR_INCLUDED
#define PROTEIN_GRAMMAR_INCLUDED

#include "ecfg/pfold.h"

struct Protein_grammar : ECFG_scores
{
  Protein_grammar (int zones, int hidden = 1);
};

#endif /* PROTEIN_GRAMMAR_INCLUDED */
