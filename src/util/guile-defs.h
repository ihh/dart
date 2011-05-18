#ifndef GUILE_DEFS_INCLUDED
#define GUILE_DEFS_INCLUDED

#include <libguile.h>
#include "util/sexpr.h"

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED

// SExpr
SExpr* scm_to_new_sexpr (SCM scm);
SCM string_to_scm (const char* s);
SCM sexpr_to_scm (SExpr* sexpr);

#endif /* GUILE_INCLUDED */

#endif /* GUILE_DEFS_INCLUDED */
