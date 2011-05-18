#include "util/guile-defs.h"

SExpr* scm_to_new_sexpr (SCM scm)
{
  // four guile API calls to get an SCM as a char* string? feel like I'm doing something the hard way here
  const char *s = scm_to_locale_string (scm_object_to_string (scm, scm_variable_ref (scm_c_lookup ("write"))));
  sstring str (s);
  SExpr* sexpr = new SExpr (str.begin(), str.end());
  free((void*) s);
  return sexpr;
}

SCM string_to_scm (const char* s)
{
  sstring str;
  str << "(quote " << s << ")";
  SCM scm = scm_c_eval_string(str.c_str());
  return scm;
}

SCM sexpr_to_scm (SExpr* sexpr)
{
  sstring str;
  str << *sexpr;
  SCM scm = string_to_scm(str.c_str());
  return scm;
}
