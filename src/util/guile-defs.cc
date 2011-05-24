#include "util/guile-defs.h"
#include "util/logfile.h"

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED

sstring scm_to_string (SCM scm)
{
  const char *s = scm_to_locale_string (scm_object_to_string (scm, scm_variable_ref (scm_c_lookup ("write"))));
  sstring str (s);
  free((void*) s);
  return str;
}

SExpr* scm_to_parent_sexpr (SCM scm)
{
  // four guile API calls to get an SCM as a char* string? feel like I'm doing something the hard way here
  sstring str = scm_to_string (scm);
  SExpr* sexpr = new SExpr (str.begin(), str.end());
  return sexpr;
}

SExpr scm_to_sexpr (SCM scm)
{
  SExpr return_sexpr;
  SExpr *result_sexpr = scm_to_parent_sexpr (scm);
  if (!result_sexpr)
    THROWEXPR ("In scm_to_sexpr: scm_to_parent_sexpr returned null");
  CTAG(3,GUILE) << "In scm_to_sexpr: " << result_sexpr->to_string() << '\n';
  if (!(result_sexpr->is_list() && result_sexpr->child.size() == 1))
    THROWEXPR ("In scm_to_sexpr: scm_to_parent_sexpr should return a list with one element");
  return_sexpr.swap ((*result_sexpr)[0]);
  delete result_sexpr;
  return return_sexpr;
}

SExpr* scm_to_new_sexpr (SCM scm)
{
  SExpr* new_sexpr = new SExpr();
  SExpr result_sexpr = scm_to_sexpr (scm);
  new_sexpr->swap (result_sexpr);
  return new_sexpr;
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
  str << sexpr->to_parenthesized_string();
  SCM scm = string_to_scm(str.c_str());
  return scm;
}

SCM scm_log_directive (SCM scm)
{
  const sstring log_dir (scm_to_string_unquoted (scm));
  const bool ok = clog_stream.clog_directive (log_dir);
  return ok ? SCM_BOOL_T : SCM_BOOL_F;
}

sstring scm_to_string_unquoted (SCM scm)
{
  THROWASSERT (scm_is_string(scm));
  char* cstr = scm_to_locale_string (scm);
  const sstring str (cstr);
  free (cstr);
  return str;
}

void init_dart_primitives (void)
{
  scm_c_define_gsubr (GUILE_LOG_DIRECTIVE, 1, 0, 0, (SCM (*)()) scm_log_directive);
}

#endif /* GUILE_INCLUDED */
