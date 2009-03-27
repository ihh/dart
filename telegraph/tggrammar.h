#ifndef TELEGRAPH_GRAMMAR_INCLUDED
#define TELEGRAPH_GRAMMAR_INCLUDED


/* Structure representing the contents of a Telegraph grammar file */
  typedef struct tg_grammar_str
  {
    /* Identifier table */
    tg_type_table     idtab;

    /* Declarations: tokens, nonterminals, terminals, parameters */
    tg_alph_vec       alph;
    tg_term_vec       term;
    tg_symbol_vec     nonterm;
    tg_prob_vec       prob;

    /* List of production rules */
    tg_rule_vec       rule;

    /* Production rule probabilities & intermediate expressions */
    tg_expression_vec expr;

    /* Results */
    tg_func_vec       func;

    /* Parses */
    tg_pexpr_vec      pexpr;

    /* Assignments of values to parameters */
    tg_expression_vec assign;

  } tg_grammar;

  /* Grammar constructor, destructor */
  tg_grammar* tg_new_grammar();
  void tg_delete_grammar (tg_grammar* tg);

  /* ID table lookup */
  /* Returns type of ID, and stores value in *value_ret */
  tg_type tg_lookup_type (tg_grammar* tg, const char* name, void** value_ret);
  tg_idtab_index tg_set_type (tg_grammar* tg, const char* name,
			      tg_type type, void* value);
  const char* tg_get_name (tg_grammar* tg, tg_idtab_index id);



#endif /* TELEGRAPH_GRAMMAR_INCLUDED */
