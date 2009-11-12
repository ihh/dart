#ifndef GUILE_DEFS_INCLUDED
#define GUILE_DEFS_INCLUDED

#include <libguile.h>
#include "seq/stockholm.h"
#include "util/sexpr.h"
#include "ecfg/ecfg.h"

// xrate functions
void init_xrate_primitives (void);

// SExpr
SExpr* scm_to_new_sexpr (SCM scm);
SCM string_to_scm (const char* s);
SCM sexpr_to_scm (SExpr* sexpr);

// ECFG
SCM ecfg_to_scm (const ECFG_scores& ecfg, const ECFG_counts* counts = 0);


#endif /* GUILE_DEFS_INCLUDED */
