#ifndef GUILE_DEFS_INCLUDED
#define GUILE_DEFS_INCLUDED

#include <libguile.h>
#include "util/sexpr.h"

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED

// SExpr
SExpr* scm_to_parent_sexpr (SCM scm);  // returns the "parent SExpr" (a list containing exactly one element, since we don't require enclosing parentheses, ugh)
SExpr scm_to_sexpr (SCM scm);  // calls scm_to_parent_sexpr, returns first element of list
sstring scm_to_string (SCM scm);
SCM string_to_scm (const char* s);
SCM sexpr_to_scm (SExpr* sexpr);

#endif /* GUILE_INCLUDED */

#endif /* GUILE_DEFS_INCLUDED */
