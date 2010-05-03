#ifndef GUILE_ECFG_INCLUDED
#define GUILE_ECFG_INCLUDED

#include "util/svisitor.h"
#include "seq/stockholm.h"
#include "ecfg/ecfg.h"

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include "guile/guile-defs.h"
#endif /* GUILE_INCLUDED */

// ECFG
#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
SCM ecfg_to_scm (const ECFG_scores& ecfg, const ECFG_counts* counts = 0);
#endif /* GUILE_INCLUDED */

// xrate functions
#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
void init_xrate_primitives (void);
#endif /* GUILE_INCLUDED */

// ECFG_Scheme_evaluator
struct ECFG_Scheme_evaluator : SExpr_Scheme_evaluator
{
  // data
  const Stockholm* stock;
  // constructors
  ECFG_Scheme_evaluator (const Stockholm* stock = 0);
};

#endif /* GUILE_ECFG_INCLUDED */
