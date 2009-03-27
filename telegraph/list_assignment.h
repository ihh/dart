/* Fn for list assignment, a modification to the grammar */


typedef tg_vector(double*) tg_assign_vec;
typedef tg_vector(int) tg_int_vector;

tg_assign_vec tg_return_assign_vec(double prob_const);
void tg_list_assign_prob(tg_grammar* tg, tg_prob* prob, tg_assign_vec prob_vec);
