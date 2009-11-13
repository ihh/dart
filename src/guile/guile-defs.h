#ifndef GUILE_DEFS_INCLUDED
#define GUILE_DEFS_INCLUDED

#include <libguile.h>
#include "util/sexpr.h"

// SExpr
SExpr* scm_to_new_sexpr (SCM scm);
SCM string_to_scm (const char* s);
SCM sexpr_to_scm (SExpr* sexpr);


#endif /* GUILE_DEFS_INCLUDED */
