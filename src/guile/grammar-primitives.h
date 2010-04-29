#ifndef GRAMMAR_PRIMITIVES_INCLUDED
#define GRAMMAR_PRIMITIVES_INCLUDED

#include "util/svisitor.h"
#include "seq/stockholm.h"

struct ECFG_Scheme_evaluator : SExpr_Scheme_evaluator
{
  // data
  const Stockholm* stock;
  // constructors
  ECFG_Scheme_evaluator (const Stockholm* stock = 0);
};

#endif /* GRAMMAR_PRIMITIVES_INCLUDED */
