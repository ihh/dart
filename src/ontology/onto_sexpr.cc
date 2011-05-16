#include "ontology/onto_sexpr.h"
#include "seq/pkeywords.h"

Terminatrix::Terminatrix()
  : alph("Uninitialized",1),
    matrix_set(alph)
{ }

EM_matrix_base& Terminatrix::rate_matrix()
{
  if (matrix_set.chain.size() != 1 || matrix_set.chain[0].matrix == NULL)
    THROWEXPR ("Oops -- bad rate matrix in Terminatrix");
  return *matrix_set.chain[0].matrix;
}

void Terminatrix_builder::init_terminatrix (Terminatrix& term, SExpr& terminatrix_sexpr)
{
  SymPVar sym2pvar;
  SymIndex term2chain;
  // initialise PScores
  init_pgroups_and_rates (term.pscores, sym2pvar, terminatrix_sexpr, &term.mutable_pgroups);
  // initialise PCounts
  init_pseudocounts (term.pcounts, term.pscores, sym2pvar, terminatrix_sexpr);
  // initialise alphabet
  init_alphabet (term.alph, terminatrix_sexpr.find_or_die (PK_ALPHABET), true);
  // initialise chain
  init_chain (term.matrix_set, term2chain, sym2pvar, terminatrix_sexpr.find_or_die (EG_CHAIN), DEFAULT_TIMEPOINT_RES, true);
}

