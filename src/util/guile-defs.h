#ifndef GUILE_DEFS_INCLUDED
#define GUILE_DEFS_INCLUDED

#include <libguile.h>
#include "util/sexpr.h"

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED

// SExpr
SExpr* scm_to_parent_sexpr (SCM scm);  // returns the "parent SExpr" (a list containing exactly one element, since we don't require enclosing parentheses, ugh)
SExpr* scm_to_new_sexpr (SCM scm);  // calls scm_to_parent_sexpr, returns first element of list (as a new SExpr object; caller must delete)
SExpr scm_to_sexpr (SCM scm);  // calls scm_to_parent_sexpr, returns first element of list (as an SExpr object on the stack)
sstring scm_to_string (SCM scm);  // quotes included
sstring scm_to_string_unquoted (SCM scm);  // strings only, quotes not included
SCM string_to_scm (const char* s);
SCM string_to_scm (const sstring& s);
SCM vector_to_scm (const vector<sstring>& sv);
SCM sexpr_to_scm (SExpr* sexpr);  // NOTE: this does **not** currently convert lists of the form (A . B) to SCM pairs!

// logging
#define GUILE_LOG_DIRECTIVE "dart-log"
SCM scm_dart_log_directive (SCM scm);
void init_dart_primitives (void);

#endif /* GUILE_INCLUDED */
#endif /* GUILE_DEFS_INCLUDED */
