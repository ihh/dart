#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "tgcompiler.h"
#include "gram.tab.h"

/* An easy way to handle this is to make an assign for each 
element in the list.  We should check that the number in the list
is correct first */		

tg_assign_vec tg_return_assign_vec(double prob_const){
	double* prob;
	*prob = prob_const;
	tg_assign_vec vec =  (tg_assign_vec)tg_vec_init(1, (void*) prob);
	return vec;
}


void tg_list_assign_prob(tg_grammar* tg, tg_prob* prob, tg_assign_vec prob_vec){
	/* cycle through each of the elements of the list, making each an assign along the way, 
	append it to the (tg_vec*) &tg->assign by  tg_vec_append ((tg_vec*) &tg->assign, (void*) assign);
	For  tg_expr* assign;  assign = tg_make_binary_expr (tg, tg_ProbAssign, lhs, rhs);
	lhs should be tg_prob
	rhs should be a tg_constant
 	We should keep track of which position in the list we are and make a tg_prob for the lhs.
	*/

	int i, j=0, k, m, n, p, q, r;
	int prob_vec_counter=0;
	int dim_size=1;
	tg_expr* lhs;
	tg_subscript* sub;
	
	tg_int_vector int_subs = (tg_int_vector) tg_vec_init(tg_size(prob->dim), (void*) 0); 
	
		/* In the case we have something like gapOpen = 0.8 */
	 if (tg_size(prob->dim) == 0) {
		tg_vec sub_vec = tg_vec_new (0);
		lhs = tg_make_prob_expr (tg, prob, sub_vec);		
	 	tg_assign_prob (tg, lhs, *prob_vec[0]);
	 }
	 
	/* First do a check on the number in the probability list */ 
	for (p=0; p<tg_size(prob->dim); p++) {
		dim_size *= prob->dim[p];
	}
	
	if (tg_size(prob_vec) != dim_size) 
			printf("There aren't enough parameters in the probability list assignment\n");	
	
	/* account for the first value ex. probbo[0][0][0] */
	tg_subscript_vec sub_vec = (tg_subscript_vec) tg_vec_init(tg_size(prob->dim), (void*) 0);
	for (n=0; n<tg_size(int_subs); n++) {
			/*printf("int_subs[%d]=%d\n", n, int_subs[n]);*/
			sub = tg_new_subscript (0, int_subs[n]);		
			sub_vec[n] = sub;		
	}	
	lhs = tg_make_prob_expr (tg, prob, (tg_vec) sub_vec);
	/*	printf("prob_vec_counter=%d\n", prob_vec_counter);
		printf("Assigning log2prob, *prob_vec[prob_vec_counter]=%f\n", *prob_vec[prob_vec_counter]);
	*/
	tg_assign_prob(tg, lhs ,*prob_vec[prob_vec_counter]);
	prob_vec_counter++;
	
	
	/* For each of the dimensions in prob, ex. 5 iterations for prob[DNA][DNA}[DNA][DNA][DNA]*/
	for (i=tg_size(prob->dim)-1; i>=0; i--) {
			/*printf("Going through outer loop\n");
			printf("i=%d\n",i);*/
		/* For each of the possible indices in each dim */
		for (j=1; j<prob->dim[i]; j++){
			/*printf("Going through middle loop\n");
			printf("j=%d\n",j);*/
			
			int_subs[i]=j; 
			/* we need to have this clause for the last dimension
			cuz this particular case
			would not go to the inner loop below */ 
			
			if (i==(tg_size(prob->dim)-1)) {
				tg_subscript_vec sub_vec = (tg_subscript_vec) tg_vec_init(tg_size(prob->dim), (void*) 0);
				for (n=0; n<tg_size(int_subs); n++) {
						/*printf("int_subs[%d]=%d\n", n, int_subs[n]);*/
						sub = tg_new_subscript (0, int_subs[n]);		
						sub_vec[n] = sub;
				}	
				lhs = tg_make_prob_expr (tg, prob, (tg_vec) sub_vec);
				/*printf("prob_vec_counter=%d\n", prob_vec_counter);
				printf("Assigning log2prob, *prob_vec[prob_vec_counter]=%f\n", *prob_vec[prob_vec_counter]);*/
				tg_assign_prob(tg, lhs ,*prob_vec[prob_vec_counter]);
				prob_vec_counter++;
			}
		
			/* This is for when we have more than one dimension, and 
			we need to cycle through the dimensions to the right of 
			the one used above */
			for (p=tg_size(prob->dim)-2; p>i; p--)
				int_subs[p]=0;
				
			for (k=tg_size(prob->dim)-1; k>i; k--) {
				/*printf("in the second loop \n");
				printf("k=%d\n", k);*/
				for (m=0; m<prob->dim[k]; m++) {
					/* update int_subs */
					/*printf("in the second loop: INNer Loop \n");
					printf("m=%d\n", m);*/
					int_subs[k]=m;
					
					/* If we have the case that the index being changed is the last dim*/
					if (k==tg_size(prob->dim)-1) {
							tg_subscript_vec sub_vec = (tg_subscript_vec) tg_vec_init(tg_size(prob->dim), (void*) 0);
							/* make a subscript vector */
							for (n=0; n<tg_size(int_subs); n++) {
								sub = tg_new_subscript (0, int_subs[n]);
								/*printf("in the second loop: int_subs[%d]=%d\n", n, int_subs[n]);	*/
								sub_vec[n] = sub;
						 	}	
							/* Now that we have a subscript vector, we can
							make a prob expr */
							lhs = tg_make_prob_expr (tg, prob , (tg_vec) sub_vec);
							/*printf("prob_vec_counter=%d\n", prob_vec_counter);
							printf("Assigning prob, *prob_vec[prob_vec_counter]=%f\n", *prob_vec[prob_vec_counter]);*/
							tg_assign_prob(tg, lhs ,*prob_vec[prob_vec_counter]);
							prob_vec_counter++;		
					}
					/* Otherwise, we're going to need to update everything to the right as well */
					else if (m!=0) {						 
						for (q=tg_size(prob->dim)-1; q>k; q--) {
							/*printf("the third loop \n");
							printf("q=%d\n", q);*/
							for (r=0; r<prob->dim[q]; r++) {
								int_subs[q]=r;		
								tg_subscript_vec sub_vec = (tg_subscript_vec) tg_vec_init(tg_size(prob->dim), (void*) 0);
								/* make a subscript vector */
								for (n=0; n<tg_size(int_subs); n++) {
									sub = tg_new_subscript (0, int_subs[n]);
									/*printf("in the second loop: int_subs[%d]=%d\n", n, int_subs[n]);	*/
									sub_vec[n] = sub;
								} 	
								/* Now that we have a subscript vector, we can
								make a prob expr */
								lhs = tg_make_prob_expr (tg, prob , (tg_vec) sub_vec);
								/*printf("prob_vec_counter=%d\n", prob_vec_counter);
								printf("Assigning prob, *prob_vec[prob_vec_counter]=%f\n", *prob_vec[prob_vec_counter]);*/
								tg_assign_prob(tg, lhs ,*prob_vec[prob_vec_counter]);
								prob_vec_counter++;		
							}
						}		
					}
				}
			}
		}	
	}
}








