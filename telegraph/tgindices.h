#ifndef TELEGRAPH_INDICES_INCLUDED
#define TELEGRAPH_INDICES_INCLUDED

/* Indices into tg_grammar arrays */
typedef int tg_prob_index;   /* prob params */
typedef int tg_alph_index;   /* token alphabets */
typedef int tg_term_index;   /* terminal alphabets */
typedef int tg_expr_index;   /* expressions */
typedef int tg_rule_index;   /* rules */
typedef int tg_idtab_index;  /* ID table (must be signed: -1 = unknown)  */

/* Indices into symbol alphabets */
typedef int tg_sym_index;  /* terminal & nonterminal symbols */

/* Indices of symbols on the RHS side of production rules */
typedef int tg_rhs_index;  /* RHS symbol position */

#endif /* TELEGRAPH_INDICES_INCLUDED */
