#ifndef ONTO_SEXPR_INCLUDED
#define ONTO_SEXPR_INCLUDED

#include "ecfg/ecfgsexpr.h"

struct Terminatrix
{
  Alphabet alph;
  PScores pscores;
  PCounts pcounts;
  ECFG_matrix_set matrix_set;
  set<int> mutable_pgroups;

  // constructor
  Terminatrix();

  // accessors
  EM_matrix_base& rate_matrix();
};

struct Terminatrix_builder : ECFG_builder
{
  // initialise a chain, a multi-char alphabet, and a parameter set
  static void init_terminatrix (Terminatrix& terminatrix, SExpr& terminatrix_sexpr);
};


#endif /* ONTO_SEXPR_INCLUDED */
