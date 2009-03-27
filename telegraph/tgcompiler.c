#include <stdio.h>
#include "tgcompiler.h"
#include "tgcheck.h"
#include "gram.tab.h"


extern void yyerror (char*);

void tg_compile (tg_grammar* tg) {
  tg_write_grammar (tg, stdout);

  printf ("\n/* Results of grammar tests:\n");

  if (tg_is_Chomsky(tg))
    printf("Yay,Chomsky \n");
  else
    printf("Boo,not Chomsky \n");

  if (tg_is_RNA_Normal(tg,stdout))
    printf("Yay, RNA normal form\n");
  else
    printf("Boo, not RNA normal form\n");

	tg_prob_score_table prob_tab = tg_alloc_prob_score_table(tg, stdout);
/*	printf("TGCOMPILER: After allocating prob score table\n");*/
	tg_write_prob_score_table(tg, prob_tab); 
/*	printf("TGCOMPILER: After writing score table\n");*/
	tg_print_prob_score_table(tg, prob_tab, stdout); 
/* printf("TGCOMPILER: After Printing Score Table\n");*/
	tg_eval_rules(tg, prob_tab, stdout);
/*	printf("TGCOMPILER: After Evaluating the rules\n");*/
  	tg_free_prob_score_table (tg, prob_tab, stdout); 
/*	printf("TGCOMPILER: After Freeing Score Table\n");*/
	printf ("*/\n");
}
