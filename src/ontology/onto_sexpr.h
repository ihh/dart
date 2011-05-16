#ifndef ONTO_SEXPR_INCLUDED
#define ONTO_SEXPR_INCLUDED

#include "ecfg/ecfgsexpr.h"

struct Terminatrix
{
  Alphabet alph;
  PScores pscores;
  PCounts pcounts, var_counts;
  set<int> mutable_pgroups;
  ECFG_matrix_set matrix_set;
  Update_statistics stats;

  bool got_counts;  // true if model has been trained

  // constructor
  Terminatrix();

  // accessors
  ECFG_chain& chain();
  EM_matrix_base& rate_matrix();
};

struct Terminatrix_builder : ECFG_builder
{
  // load
  static void init_terminatrix (Terminatrix& terminatrix, SExpr& terminatrix_sexpr);
  // save
  static void terminatrix2stream (ostream& out, const Terminatrix& terminatrix);
};


#endif /* ONTO_SEXPR_INCLUDED */
