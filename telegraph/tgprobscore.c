#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "tgcompiler.h"
#include "gram.tab.h"

/* Suppose we have the following Telegraph grammar:

	token DNA { A C G T };
	prob paramZero[5], paramOne[DNA], paramTwo[DNA][DNA][3];

Let "tg" be a pointer to the tg_grammar struct:

	tg_grammar* tg;

Then

	tg_size(tg->prob) == 3

	tg_get_name (tg, tg->prob[0]->id) == "paramZero"
	tg_get_name (tg, tg->prob[1]->id) == "paramOne"
	tg_get_name (tg, tg->prob[2]->id) == "paramTwo"                                                         
 
To initialise the tg_prob_score_table we would need the following:	

	tg_prob_score_table x = (tg_prob_score_table) tg_vec_init (3, (void*) 0);
	x[0] = malloc (5 * sizeof(tg_score));  // paramZero
	x[1] = malloc (4 * sizeof(tg_score));  // paramOne
	x[2] = malloc (4 * 4 * 3 * sizeof(tg_score));  // paramTwo
 
then to access e.g. paramZero[3], we need x[0][3]

Let's define an MSB-first indexing convention for accessing the
elements of multidimensional prob's.
(MSB = "Most Significant Byte" - should strictly be "Most Significant Digit")
Thus,
 
	paramTwo[C][G][0] == paramTwo[1][2][0] == x[2][1*4*3 + 2*3 + 0]
 
...formally, if dim[n] is the dimension of the n'th square bracket,
so that (for paramTwo) dim[0]=4, dim[1]=4, dim[2]=3,

then the multidimensional index [x][y][z] translates
to the integer index (x*dim[1]*dim[2] + y*dim[2] + z)
 
...or, more generally, if idx[n] is the index in the n'th square bracket
(so x==idx[0], y==idx[1], z==idx[2])
then the single-dimensional integer index is defined to be
	\sum_{i=0}^{N-1} (idx[i] * \prod_{j=i+1}^{N-1} dim[j])
 
*/
/*typedef tg_vector(tg_subscript*) tg_subscript_vec;*/
/* Lookup table mapping all prob's to scores */
/*	typedef tg_vector(tg_score*) tg_prob_score_table;*/

tg_prob_score_table tg_alloc_prob_score_table (tg_grammar* tg, FILE* file){
	tg_prob_score_table prob_tab = (tg_prob_score_table) tg_vec_init (tg_size(tg->prob), (void*) 0);
	/* for each of the probs, we need to malloc enough space based on the dimensions of the prob params */
	int i, j, dim_size;
	
	for (i=0; i < tg_size(tg->prob); i++) {
		/* Calculate the total size to malloc using dim vector of each probability */
		
		dim_size=1;	
		for (j=0; j<(tg_size(tg->prob[i]->dim)); j++) {
			dim_size *= tg->prob[i]->dim[j];
		}	
		prob_tab[i] = (tg_score_array) tg_vec_init (dim_size, (void*) 0);
 	}
	return prob_tab;
}
		 


void tg_free_prob_score_table (tg_grammar* tg, tg_prob_score_table prob_tab, FILE* file){
 	int i;
	
	for (i=0; i<tg_size(prob_tab); i++) {
		if (tg_size(prob_tab[i])>1){
			tg_vec_delete((tg_vec*)&prob_tab[i]);  
		}				
	}
}
	
/* We want the single integer index of the subscript */
int tg_map_subscript2index (tg_subscript_vec subscript_vec, tg_prob* prob){
	/* roughly, subscript_vec is "idx", and prob contains "dim" - see above 
		\sum_{i=0}^{N-1} (idx[i] * \prod_{j=i+1}^{N-1} dim[j])*/
	int i,j, index=0;
	for (i=0; i<tg_size(subscript_vec); i++) {
		int dimlen = 1;
		for (j=i+1; j<tg_size(prob->dim); j++){
			dimlen *= prob->dim[j];
		} 
		if (!subscript_vec[i]->dollar) 
			/* This part is wrong */
			/* index += (subscript_vec[i]->index.dollar_subscript) * dimlen;
		else 
		*/
			index += (subscript_vec[i]->index.int_subscript) * dimlen; 
	}
	return index;	
}

 
	/* To fill out the prob score table, we want to go through the assignment
	and evaluate them. */
void tg_write_prob_score_table(tg_grammar* tg, tg_prob_score_table prob_tab) { 
	int i;
 	tg_score score;   

 	/*temp*/
 	tg_unsigned_vec empty_dollar_binding;
 	
	if (tg_size(tg->assign)) {   
	   for (i = 0; i < tg_size(tg->assign); ++i) {
			score = tg_eval_expr (tg, tg->assign[i], prob_tab, stdout, empty_dollar_binding);
      }
 	}
}


void tg_print_prob_score_table(tg_grammar* tg, tg_prob_score_table prob_tab, FILE* file) {
	fprintf(file, "************* Prob Score Table *************\n");
	int i,j;
	for (i=0; i<tg_size(prob_tab); i++) {
		fprintf(file, "* prob_tab[%d]", i);	 
 		for (j=0;j<tg_size(prob_tab[i]);j++){
			fprintf(file, "%d ", prob_tab[i][j]);
		}
		fprintf(file, "\n");	 
	}
	fprintf(file, "********************************************\n");
}

 /* Here we want to evaluate the expression and return the evaluated score 
	Do we need to update the score tab here or after this function is called? */
tg_score tg_eval_expr (tg_grammar* tg, tg_expr* expr, tg_prob_score_table prob_tab, FILE* file, tg_unsigned_vec dollar_binding){
	int index1, index2;
  tg_expr* expr1;
  tg_expr* expr2;
  tg_score score1;
  tg_score score2;
  tg_prob* prob;
  tg_vec vec;
  tg_score score;
 
	switch (expr->op) {
	  case tg_Multiplication:	
	    expr1 = (tg_expr*) expr->operand1;
	    expr2 = (tg_expr*) expr->operand2;
	    score1 = tg_eval_expr (tg, expr1, prob_tab, file, dollar_binding);
	    score2 = tg_eval_expr (tg, expr2, prob_tab, file, dollar_binding);
	    return score1 * score2;
	  	break;
 	  
	  case tg_Division:
	    expr1 = (tg_expr*) expr->operand1;
	    expr2 = (tg_expr*) expr->operand2;
	    score1 = tg_eval_expr (tg, expr1, prob_tab, file, dollar_binding);
	    score2 = tg_eval_expr (tg, expr2, prob_tab, file, dollar_binding);
	    return score1 / score2;
	   break;
 
		case tg_Addition:
			expr1 = (tg_expr*) expr->operand1;
			expr2 = (tg_expr*) expr->operand2;
			score1 = tg_eval_expr (tg, expr1, prob_tab, file, dollar_binding);
		   score2 = tg_eval_expr (tg, expr2, prob_tab,file, dollar_binding);
		   return score1 + score2;
    	break;
 
  case tg_ProbAssign:
    expr1 = (tg_expr*) expr->operand1;
    expr2 = (tg_expr*) expr->operand2;
 	
	if (expr1->op == tg_ProbBinding) {
		/* first we need to find the index 'a' in prob_tab[a][b] */
		prob = (tg_prob*) expr1->operand1;
		index1 = prob->index;
		index2 = tg_map_subscript2index ((tg_subscript_vec)expr1->operand2, prob);
 
	   /* If we already have a constant on the rhs, then add the assignment to the assign vector */
    	if (expr2->op == tg_Constant && (score2 = -(tg_score) expr2->operand1) > 0) { 	 
			prob_tab[index1][index2] = score2;
    	}  
    	else {
    	/* This is the case that we don't have valid constant on the rhs of the probBinding*/
		  	score2 = tg_eval_expr (tg, expr2, prob_tab, file, dollar_binding);
	 		prob_tab[index1][index2] = score2;
	 	}
 			return score2;
	}
 		/* else, in the case we don't have  prob binding on the lhs 
			Here,  we should be handling FUNC's */
	else {
		return score2;
	}
    break;

 	
    /* "[" 
    	 of the form p[A][G]
       We are given the tg_prob in operand1 
       and the tg_vec in operand2*/
  case tg_ProbBinding:
    prob = (tg_prob*) expr->operand1;
    vec = (tg_vec) expr->operand2;
    index1 = prob->index;
    /*printf("In ProbBinding\n");*/
	
    if (tg_size(vec)) {
	    if (((tg_subscript*) vec[0])->dollar){
	    	/*printf("We found a dollar subscript vector\n");*/
			index2 = tg_map_bound_subscript2index ((tg_subscript_vec)vec, prob, dollar_binding);
		 	/*printf("After map bound subscript2 index\n");*/
		  	/*printf("index2 = %d\n", index2);	*/
			} 
		 else {
		 	/*printf("We DIDNT find a dollar subscript vector\n");*/
		 	/*printf("Before map subscript to index\n");*/
		 	index2 = tg_map_subscript2index ((tg_subscript_vec)vec, prob);
		 	/*printf("After map subscript to index\n");	 		*/
		 	}
	 	 } 	
 	 else
 	 		index2 = 0; 
 	 		
		/*printf("right before returning in probbinding\n"); */
  		return prob_tab[index1][index2];
    break;
 
  case tg_Constant:
    score = (tg_score) expr->operand1;
    /*printf("tg_Constant: %d\n", score);*/
    return score;
    break;
  
  default:
	return 0;
    break;
  }
}	

 
/* with dollar symbols: something like
	tg_score tg__bound_expr (tg_grammar* tg, tg_expr* expr, tg_prob_score_table score_tab,
										  tg_unsigned_vec dollar_binding);

so dollar_binding[0] == value of $1
   dollar_binding[1] == value of $2, etc.
*/

/* I think we can bypass this and go directly to tg map bound subscript to index 
	when evaluating prob binding, unless, that is, we want to make another operator
	specifically for dollar bindings */
tg_score tg__bound_expr (tg_grammar* tg, tg_expr* expr, tg_prob_score_table prob_tab, tg_unsigned_vec dollar_binding){
	tg_score score;
	tg_prob* prob;
	tg_vec vec;
	
	/* If we have a ProbBinding, it will be bound */
	if (expr->op == tg_ProbBinding) {
	  		prob = (tg_prob*) expr->operand1;
    	vec = (tg_vec) expr->operand2;
 	/* This should be the same */	
    int index1 = prob->index;
    /* This will not depend on the subscript vec alone anymore */
	 int index2 = tg_map_bound_subscript2index ((tg_subscript_vec)vec, prob, dollar_binding);
    return prob_tab[index1][index2];
	}
	 
	return score;
}


int tg_map_bound_subscript2index (tg_subscript_vec subscript_vec, tg_prob* prob, tg_unsigned_vec dollar_binding){
	/* First, we need to get the values of the subscripts $1, $3, etc... */
	/* We use these values to calculate the index */
	int i,j, index=0;
	for (i=0; i<tg_size(subscript_vec); i++) {
		int dimlen = 1;
		for (j=i+1; j<tg_size(prob->dim); j++){
			dimlen *= prob->dim[j];
		} 
		if (subscript_vec[i]->dollar) 
			index += dollar_binding[subscript_vec[i]->index.dollar_subscript] *dimlen;
			/* We're assuming that dollar_binding returns the integer number of the subscript*/
	}
	return index;	
}

/* probably just a temp fn until we can find out how we are going to change 
dollar binding */
tg_unsigned_vec get_dollar_binding(tg_grammar* tg, tg_rule* rule) {
	/* the size of dollar binding should be the same as the num of rhs in the rule */
	tg_unsigned_vec dollar_binding = (tg_unsigned_vec) tg_vec_init (tg_size(rule->rhs), (void*) 0);
	return dollar_binding;
}

void tg_eval_rules(tg_grammar* tg, tg_prob_score_table prob_tab, FILE* file){
	tg_rule* tg_rule_temp;
	tg_expr* tg_expr_temp; 
	tg_score score;
	int i;
	tg_unsigned_vec dollar_binding;
	
	if (tg_size(tg->rule)) {
		/* for each rule */  
		for (i=0; i<tg_size(tg->rule); i++) {
      	/* pointer to tg_rule */
      	tg_rule_temp = (tg_rule*)(tg->rule[i]);      
			tg_expr_temp = (tg_expr*)(tg_rule_temp->expr);
			/* calling this fn will most likely be temp 
				This will probably all be folded into a dynamic algorithm program later */
			dollar_binding = get_dollar_binding(tg, tg_rule_temp);		
			score = tg_eval_expr (tg, tg_expr_temp, prob_tab, file, dollar_binding);
			/* Debug */
			fprintf(file, "Rule %d has prob = %d \n", i, score);
			/* Q: Where do we store the prob? */ 
			/* Right now I guess we'll just put it tg_constant form and tuck it back into the expr */
			/* The thing with the fn used below is that we hafta give it a prob, which will be 
			converted to a score */
		
			tg_expr_temp = tg_make_const_expr (tg, score);
 		}
	}
}
