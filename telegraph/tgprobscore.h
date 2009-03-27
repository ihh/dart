/* Lookup table mapping all prob's to scores */

/*typedef tg_vector(tg_score*) tg_prob_score_table;*/
typedef tg_vector(tg_score_array) tg_prob_score_table;

/* Functions to manipulate the prob score table */
tg_prob_score_table tg_alloc_prob_score_table (tg_grammar* tg, FILE* file);
void tg_free_prob_score_table (tg_grammar* tg, tg_prob_score_table prob_tab, FILE* file);
int tg_map_subscript2index (tg_subscript_vec subscript_vec, tg_prob* prob);
int tg_map_bound_subscript2index (tg_subscript_vec subscript_vec, tg_prob* prob, tg_unsigned_vec dollar_binding);
void tg_write_prob_score_table(tg_grammar* tg, tg_prob_score_table prob_tab);
void tg_print_prob_score_table(tg_grammar* tg, tg_prob_score_table prob_tab, FILE* file);

/* Fn to evaluate the rule expression, also used for evaluating assignment bindings */
tg_score tg_eval_expr (tg_grammar* tg, tg_expr* expr, tg_prob_score_table prob_tab, FILE* file, tg_unsigned_vec dollar_binding);
void tg_eval_rules(tg_grammar* tg, tg_prob_score_table prob_tab, FILE* file);
tg_unsigned_vec get_dollar_binding(tg_grammar* tg, tg_rule* rule);

