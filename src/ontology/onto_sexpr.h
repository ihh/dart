#ifndef ONTO_SEXPR_INCLUDED
#define ONTO_SEXPR_INCLUDED

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include <libguile.h>
#else
#error Terminatrix requires guile! Please install guile from http://www.gnu.org/software/guile/ and re-run the configure script.
#endif

#include "ecfg/ecfgsexpr.h"
#include "util/svisitor.h"

struct Terminatrix
{
  // primary object interface to Scheme
  SExpr_Scheme_evaluator& scheme;

  // SExpr's
  SExpr init_sexpr, model_sexpr, tree_db_sexpr, knowledge_sexpr;

  // model
  const Alphabet dummy_alph;
  list<Alphabet> alph_list;
  Alphabet_dictionary alph_dict;
  PScores pscores;
  PCounts pcounts, var_counts;
  set<int> mutable_pgroups;
  ECFG_matrix_set matrix_set;

  // results
  Update_statistics stats;

  // status variables
  bool got_counts;  // true if model has been trained

  // constructor
  Terminatrix (SExpr_Scheme_evaluator& scheme);  // the SExpr_Scheme_evaluator is assumed to have been initialized

  // helpers
  void eval_funcs();

  // accessors
  ECFG_chain& chain();
  EM_matrix_base& rate_matrix();
};

struct Terminatrix_builder : ECFG_builder
{
  // load
  static void init_terminatrix (Terminatrix& terminatrix, SExpr& terminatrix_sexpr);
  // save
  static void terminatrix2stream (ostream& out, Terminatrix& terminatrix);

  // helpers
  static void init_terminatrix_model (Terminatrix& terminatrix, SExpr& terminatrix_model_sexpr);
  static void init_terminatrix_member_sexpr (SExpr& member_sexpr, SExpr& parent_sexpr, const char* tag);

  static void terminatrix_model2stream (ostream& out, Terminatrix& terminatrix);
  static void terminatrix_member_sexpr2stream (ostream& out, SExpr& member_sexpr);
};


#endif /* ONTO_SEXPR_INCLUDED */
