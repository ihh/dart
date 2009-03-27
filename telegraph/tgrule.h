#ifndef TELEGRAPH_RULE_INCLUDED
#define TELEGRAPH_RULE_INCLUDED

  /* Production rules */
  typedef struct tg_rule_str
  {
    tg_rule_index     index; /* index of this rule */
    tg_sym*           lhs;   /* index of nonterminal on LHS */
    tg_identifier_vec rhs;   /* array of (tg_idtab_entry*)'s on RHS */
    tg_expr*          expr;  /* index of rule probability expression */
  } tg_rule;

typedef tg_vec tg_rule_vec;  /* elements have type (tg_rule*) */

#endif /* TELEGRAPH_RULE_INCLUDED */
