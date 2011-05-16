#include "ontology/onto_sexpr.h"
#include "ontology/onto_keywords.h"
#include "seq/pkeywords.h"

Terminatrix::Terminatrix()
  : alph("Uninitialized",1),
    matrix_set(alph),
    got_counts(false)
{ }

EM_matrix_base& Terminatrix::rate_matrix()
{
  if (matrix_set.chain.size() != 1 || matrix_set.chain[0].matrix == NULL)
    THROWEXPR ("Oops -- no rate matrix in Terminatrix");
  return *matrix_set.chain[0].matrix;
}

ECFG_chain& Terminatrix::chain()
{
  if (matrix_set.chain.size() != 1)
    THROWEXPR ("Oops -- no rate matrix in Terminatrix");
  return matrix_set.chain[0];
}

void Terminatrix::eval_funcs()
{
  matrix_set.eval_funcs (pscores);
}

void Terminatrix_builder::init_terminatrix (Terminatrix& term, SExpr& sexpr)
{
  SExpr terminatrix_sexpr = sexpr.find_or_die (TERMINATRIX);
  SymPVar sym2pvar;
  SymIndex term2chain;
  // initialise PScores
  init_pgroups_and_rates (term.pscores, sym2pvar, terminatrix_sexpr, &term.mutable_pgroups);
  // initialise PCounts
  init_pseudocounts (term.pcounts, term.pscores, sym2pvar, terminatrix_sexpr);
  term.var_counts = PCounts (term.pscores);
  // initialise alphabet
  init_alphabet (term.alph, terminatrix_sexpr.find_or_die (PK_ALPHABET), true);
  // initialise chain
  init_chain (term.matrix_set, term2chain, sym2pvar, terminatrix_sexpr.find_or_die (EG_CHAIN), DEFAULT_TIMEPOINT_RES, true);
}

void Terminatrix_builder::terminatrix2stream (ostream& out, const Terminatrix& term)
{
  out << '(' << TERMINATRIX << '\n';
  pscores2stream (out, term.pscores, term.mutable_pgroups, &term.var_counts);
  pcounts2stream (out, term.pcounts, EG_PSEUDOCOUNTS, (const PCounts*) 0, true, false);
  if (term.got_counts)
    pcounts2stream (out, term.var_counts, EG_EXPECTED_COUNTS, &term.pcounts, true, true);
  chain2stream (out, term.alph, term.pscores, ((Terminatrix&)term).chain(), term.got_counts ? &term.stats : NULL);
  alphabet2stream (out, term.alph);
  out << ')';
}
