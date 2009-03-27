#ifndef TELEGRAPH_GRAMMAR_BUILDERS_INCLUDED
#define TELEGRAPH_GRAMMAR_BUILDERS_INCLUDED


  /* Expression builders */
tg_expr* tg_make_binary_expr (tg_grammar* tg, char op, void* op1, void* op2);
tg_expr* tg_make_const_expr (tg_grammar* tg, double r);
tg_expr* tg_make_prob_expr (tg_grammar* tg, tg_prob* prob,
			    tg_vec subscript);
tg_expr* tg_make_func_expr (tg_grammar* tg, tg_func* func);
tg_expr* tg_make_psum_expr (tg_grammar* tg, tg_pexpr* pexpr);

/* Token alphabet builders */
tg_alph* tg_new_alphabet (tg_grammar* tg);
void tg_add_token (tg_grammar* tg, tg_alph* alph, const char* token);
void tg_declare_alphabet (tg_grammar* tg, tg_alph* alph, const char* name);

/* Terminal builder */
tg_term* tg_declare_terminal (tg_grammar* tg, tg_alph* alph,
			      const char* name);

tg_sym* tg_get_terminal (tg_grammar* tg, tg_term* term, tg_sym* tok);

/* Nonterminal builder */
void tg_declare_nonterminal (tg_grammar* tg, const char* name);

/* Probability parameter builder */
tg_prob* tg_declare_prob (tg_grammar* tg, tg_vec dim, const char* name);

/* Probability parameter assignments */
void tg_declare_func (tg_grammar* tg, const char* name, tg_expr* rhs);
void tg_assign_expr (tg_grammar* tg, tg_expr* lhs, tg_expr* rhs);
void tg_assign_prob (tg_grammar* tg, tg_expr* lhs, double val);
void tg_assign_log2prob (tg_grammar* tg, tg_expr* lhs, double log2_val);

/* Production rules */
void tg_declare_rule (tg_grammar* tg,
		      tg_sym* lhs, tg_vec rhs, tg_expr* expr);

/* Method to check prob subscripts in an expression for validity */
int tg_subscripts_valid (tg_grammar* tg, tg_vec rhs, tg_expr* expr);

/* Method to create parse tree nodes */
tg_parse_node* tg_new_parse_node (tg_grammar* grammar, tg_sym* label,
				  tg_parse_vec child);

/* Method to create parse expressions */
tg_pexpr* tg_new_pexpr (tg_grammar* tg, int op, void* operand);

/* Method to declare pset's */
void tg_declare_pset (tg_grammar* grammar, tg_pexpr* pexpr,
		      const char* name);

/* Methods to build parse expressions */
tg_pexpr* tg_make_dynamic_pexpr (tg_grammar* tg, tg_sequence_vec seqs,
				 tg_sym_index start);

/* Method to check that a given result maps to a psum */
void tg_assert_is_psum (tg_grammar* grammar, tg_expr* result_expr);


#endif /* TELEGRAPH_GRAMMAR_BUILDERS_INCLUDED */
