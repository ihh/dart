#include "ontology/onto_sexpr.h"
#include "ontology/onto_keywords.h"
#include "seq/pkeywords.h"

Terminatrix::Terminatrix()
  : dummy_alph(),
    matrix_set(dummy_alph),
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
  SExpr model_sexpr = terminatrix_sexpr.find_or_die (TERMINATRIX_MODEL);
  SymPVar sym2pvar;
  SymIndex term2chain;
  // initialise PScores
  init_pgroups_and_rates (term.pscores, sym2pvar, model_sexpr, &term.mutable_pgroups);
  // initialise PCounts
  init_pseudocounts (term.pcounts, term.pscores, sym2pvar, model_sexpr);
  term.var_counts = PCounts (term.pscores);
  // initialise alphabet(s)
  const vector<SExpr*> alph_sexprs = model_sexpr.find_all (PK_ALPHABET);
  if (alph_sexprs.size() == 0)
    THROWEXPR ("You need at least one alphabet");
  for_const_contents (vector<SExpr*>, alph_sexprs, alph_sexpr)
    {
      term.alph_list.push_back (term.dummy_alph);
      Alphabet& alph (term.alph_list.back());
      init_alphabet (alph, **alph_sexpr, true);
      term.alph_dict.add (alph);
    }
  // initialise chain
  init_chain (term.matrix_set, term.alph_dict, term.alph_list.front(), term2chain, sym2pvar, model_sexpr.find_or_die (EG_CHAIN), DEFAULT_TIMEPOINT_RES, true);
}

void Terminatrix_builder::terminatrix2stream (ostream& out, const Terminatrix& term)
{
  out << '(' << TERMINATRIX << '\n';
  out << " (" << TERMINATRIX_MODEL << '\n';
  pscores2stream (out, term.pscores, term.mutable_pgroups, &term.var_counts);
  pcounts2stream (out, term.pcounts, EG_PSEUDOCOUNTS, (const PCounts*) 0, true, false);
  if (term.got_counts)
    pcounts2stream (out, term.var_counts, EG_EXPECTED_COUNTS, &term.pcounts, true, true);
  chain2stream (out, term.pscores, ((Terminatrix&)term).chain(), term.alph_dict, term.got_counts ? &term.stats : NULL);
  for_const_contents (list<Alphabet>, term.alph_list, alph)
    alphabet2stream (out, *alph);
  out << "))\n";
}
