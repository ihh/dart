#include "util/guile-defs.h"
#include "util/logfile.h"
#include "util/guile-keywords.h"
#include "util/discrete_gamma.h"

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
  CTAG(1,GUILE) << "In scm_to_parent_sexpr: " << str << '\n';
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

SCM string_to_scm (const sstring& s)
{
  return string_to_scm (s.c_str());
}

SCM vector_to_scm (const vector<sstring>& sv)
{
  SCM list_scm = scm_list_n (SCM_UNDEFINED);  // empty list
  for_const_reverse_contents (vector<sstring>, sv, s)
    list_scm = scm_cons (scm_from_locale_string (s->c_str()), list_scm);
  return list_scm;
}

SCM vector_to_scm (const vector<double>& dv)
{
  SCM list_scm = scm_list_n (SCM_UNDEFINED);  // empty list
  for_const_reverse_contents (vector<double>, dv, d)
    list_scm = scm_cons (scm_from_double(*d), list_scm);
  return list_scm;
}

SCM sexpr_to_scm (SExpr* sexpr)
{
  sstring str;
  str << sexpr->to_parenthesized_string();
  SCM scm = string_to_scm(str.c_str());
  return scm;
}

sstring scm_to_string_unquoted (SCM scm)
{
  THROWASSERT (scm_is_string(scm));
  char* cstr = scm_to_locale_string (scm);
  const sstring str (cstr);
  free (cstr);
  return str;
}

SCM scm_log_directive (SCM scm)
{
  const sstring log_dir (scm_to_string_unquoted (scm));
  const bool ok = clog_stream.clog_directive (log_dir);
  return ok ? SCM_BOOL_T : SCM_BOOL_F;
}

SCM scm_discrete_gamma (SCM alpha_scm, SCM beta_scm, SCM K_scm, bool medians)
{
  THROWASSERT (scm_is_real(alpha_scm));
  THROWASSERT (scm_is_real(beta_scm));
  THROWASSERT (scm_is_integer(K_scm));
  const double alpha = scm_to_double (alpha_scm);
  const double beta = scm_to_double (beta_scm);
  const int K = scm_to_int (K_scm);
  THROWASSERT (alpha > 0);
  THROWASSERT (beta > 0);
  THROWASSERT (K > 0);
  Discrete_gamma dg (alpha, beta, K, medians);
  return vector_to_scm (dg);
}

SCM scm_discrete_gamma_medians (SCM alpha_scm, SCM beta_scm, SCM K_scm)
{
  return scm_discrete_gamma (alpha_scm, beta_scm, K_scm, true);
}

SCM scm_discrete_gamma_means (SCM alpha_scm, SCM beta_scm, SCM K_scm)
{
  return scm_discrete_gamma (alpha_scm, beta_scm, K_scm, false);
}

SCM scm_incomplete_gamma_inverse (SCM prob_scm, SCM alpha_scm, SCM beta_scm)
{
  THROWASSERT (scm_is_real(alpha_scm));
  THROWASSERT (scm_is_real(beta_scm));
  THROWASSERT (scm_is_real(prob_scm));
  const double alpha = scm_to_double (alpha_scm);
  const double beta = scm_to_double (beta_scm);
  const double prob = scm_to_double (prob_scm);
  THROWASSERT (alpha > 0);
  THROWASSERT (beta > 0);
  THROWASSERT (prob >= 0);
  return scm_from_double (POINTGAMMA (prob, alpha, beta));
}

SCM scm_incomplete_gamma (SCM x_scm, SCM alpha_scm, SCM beta_scm)
{
  THROWASSERT (scm_is_real(alpha_scm));
  THROWASSERT (scm_is_real(beta_scm));
  THROWASSERT (scm_is_real(x_scm));
  const double alpha = scm_to_double (alpha_scm);
  const double beta = scm_to_double (beta_scm);
  const double x = scm_to_double (x_scm);
  THROWASSERT (alpha > 0);
  THROWASSERT (beta > 0);
  THROWASSERT (x >= 0);
  return scm_from_double (IncompleteGamma (x*beta, alpha, LnGamma (alpha)));
}

SCM scm_ln_gamma (SCM alpha_scm)
{
  THROWASSERT (scm_is_real(alpha_scm));
  const double alpha = scm_to_double (alpha_scm);
  THROWASSERT (alpha > 0);
  return scm_from_double (LnGamma (alpha));
}

SCM scm_gamma_density (SCM x_scm, SCM alpha_scm, SCM beta_scm)
{
  THROWASSERT (scm_is_real(alpha_scm));
  THROWASSERT (scm_is_real(beta_scm));
  THROWASSERT (scm_is_real(x_scm));
  const double alpha = scm_to_double (alpha_scm);
  const double beta = scm_to_double (beta_scm);
  const double x = scm_to_double (x_scm);
  THROWASSERT (alpha > 0);
  THROWASSERT (beta > 0);
  THROWASSERT (x >= 0);
  return scm_from_double (GammaDensity (x, alpha, beta));
}

void init_dart_primitives (void)
{
  scm_c_define_gsubr (GUILE_LOG_DIRECTIVE, 1, 0, 0, (scm_t_subr) scm_log_directive);

  scm_c_define_gsubr (GUILE_DISCRETE_GAMMA_MEDIANS, 3, 0, 0, (scm_t_subr) scm_discrete_gamma_medians);
  scm_c_define_gsubr (GUILE_DISCRETE_GAMMA_MEANS, 3, 0, 0, (scm_t_subr) scm_discrete_gamma_means);

  scm_c_define_gsubr (GUILE_LOG_GAMMA, 1, 0, 0, (scm_t_subr) scm_ln_gamma);
  scm_c_define_gsubr (GUILE_GAMMA_DENSITY, 3, 0, 0, (scm_t_subr) scm_gamma_density);
  scm_c_define_gsubr (GUILE_INCOMPLETE_GAMMA, 3, 0, 0, (scm_t_subr) scm_incomplete_gamma);
  scm_c_define_gsubr (GUILE_INCOMPLETE_GAMMA_INVERSE, 3, 0, 0, (scm_t_subr) scm_incomplete_gamma_inverse);
}

#endif /* GUILE_INCLUDED */
