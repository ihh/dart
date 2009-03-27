/***********************************************************
 * tgcheck.c
 * 16 March 2005
 * 
 * Checks if the grammar is in Chomsky Normal Form
 * or RNA Normal Form
 ***********************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "tgcompiler.h"
#include "tgcheck.h"
#include "gram.tab.h"

/* Check for grammar format */
int tg_is_Chomsky (tg_grammar* tg){

  tg_rule* tg_rule_temp;
  tg_sym* lhs_temp;
  tg_idtab_entry* rhs_temp0;
  tg_idtab_entry* rhs_temp1;
  int i;
  int rhs_count; 

  if (tg_size(tg->rule)) {
    /* for each rule */  
    for (i=0; i<tg_size(tg->rule); i++) {
      /* pointer to tg_rule */
      tg_rule_temp = (tg_rule*)(tg->rule[i]);      
      lhs_temp = tg_rule_temp -> lhs;    /* pointer to the tg_sym index of the LHS nonterminal */
      /* LHS of rules should be non-terminal   */
      
      /* RHS of rules should be terminals followed by a nonterminal */
    
      rhs_count = tg_size(tg_rule_temp->rhs);
      rhs_temp0 = (tg_idtab_entry*) tg_rule_temp->rhs[0];
      rhs_temp1 = (tg_idtab_entry*) tg_rule_temp->rhs[1];
  
      if (rhs_count == 2){           /* should be two nonterminals */
	if (rhs_temp0->type != NONTERMID || rhs_temp1->type != NONTERMID)
	  return 0; 
      }
      else if (rhs_count == 1 ){     /* should be one terminal */
	if (rhs_temp0->type != TERMID && rhs_temp0->type !=TERMALPHID)       
	  return 0; 
      }
      else                           /* return 0 if rhs count is >2 or 0 */
	return 0;
    } /* end for this rule */
  } 
  /* Yes, this is Chomsky Normal form  at this point */
  return 1;
}
  

int tg_is_RNA_Normal (tg_grammar* tg, FILE* file ){

  tg_rule* tg_rule_temp;

  int rhs_count = tg_size(tg_rule_temp->rhs);
  int rhs_index = 0;   
  tg_idtab_entry* rhs_temp;
  tg_idtab_entry* rhs_temp2;
  tg_alph_index a_i;
  tg_alph_index a_temp;


  int i;
  int n; 
  int NumNonTerms;
  int r; 
  int old_rhs_index;

  if (tg_size(tg->rule)) {

    /* for each rule */  
    for (n=0; n<tg_size(tg->rule); n++) {
      NumNonTerms=0;
      /* pointer to tg_rule */
      tg_rule_temp = (tg_rule*)(tg->rule[n]);      
      /* commented out by IH - doesn't seem necessary any more? */
      /* fprintf(file, "Rule number...%d\n", n); */
      rhs_count = tg_size(tg_rule_temp->rhs);
      rhs_index = 0;   
      
      /*
	For bifurcation, check if more than one NTERM, if yes and rhs_count=2 then okk
	else, return 0
      */

      for (r=0; r<rhs_count; r++) {
	if (((tg_idtab_entry*) (tg_rule_temp->rhs[r]))->type == NONTERMID) 
	  NumNonTerms++;
      }

      if (NumNonTerms==2 && rhs_count!=2) 
	return 0;

      if (NumNonTerms>2)
	return 0;

      for (rhs_index=0; 
	   rhs_index < rhs_count && 
	     (rhs_temp=((tg_idtab_entry*)(tg_rule_temp->rhs[rhs_index])))->type != NONTERMID;
	   rhs_index++){
	if (rhs_temp->type == TERMID) {
	  a_i = ((tg_sym*)(rhs_temp->value))->alph;
	  if (rhs_index > 0) {
	    for (i =rhs_index-1; i>=0; i--) {
		  rhs_temp2 = (tg_idtab_entry*)(tg_rule_temp->rhs[i]);

		  if (rhs_temp2->type == TERMID) {
		    a_temp = ((tg_sym*)(rhs_temp2->value))->alph;
		    if (a_i  == a_temp) 
		      return 0;
		  }
		  else if (rhs_temp2->type == TERMALPHID){
		    a_temp = ((tg_alph*)(rhs_temp2->value))->alph;
		    if (a_i  == ((tg_alph*)(rhs_temp2->value))->alph) 
		      return 0;
		  }
		}
	      }
	    }

	    if (rhs_temp->type == TERMALPHID) {
	      a_i = ((tg_alph*)(rhs_temp->value))->alph;
	      if (rhs_index > 0) {
		for (i=rhs_index-1; i>=0; i--) 
		  {
		    rhs_temp2 = (tg_idtab_entry*)(tg_rule_temp->rhs[i]);
		    if (rhs_temp2->type == TERMID){
		      a_temp = ((tg_sym*)(rhs_temp2->value))->alph;
		      if (a_i  == ((tg_sym*)(rhs_temp2->value))->alph) 
			return 0;
		    }
		    
		    else if (rhs_temp2->type == TERMALPHID){
		      a_temp = ((tg_alph*)(rhs_temp2->value))->alph;
		      if (a_i  == ((tg_alph*)(rhs_temp2->value))->alph) 
			return 0;
		    }
		  }
	      }
	      }
	  }

      /* at this point we have run through all the beginning terminals
	 if we don't find a NONTERMINAL, then this is incorrect
      */

      if (rhs_count == rhs_index) {
	return 0;
      }

      /* skip the NONTERMINAL 
       and start checking the latter half of the rhs
      */

      old_rhs_index = rhs_index;

      for (rhs_index+=1; rhs_index<rhs_count; rhs_index++)
	  {
	    rhs_temp=(tg_idtab_entry*)(tg_rule_temp->rhs[rhs_index]);
	    if (rhs_temp->type == TERMID) {
	      a_i = ((tg_sym*)(rhs_temp->value))->alph;
	      if (rhs_index>0) { 
		for (i=rhs_index-1; i>old_rhs_index; i--) 
		{
		  rhs_temp2 = (tg_idtab_entry*)(tg_rule_temp->rhs[i]);
 		  if (rhs_temp2->type == TERMID) {
		    if (a_i  == ((tg_sym*)(rhs_temp2->value))->alph) 
		      return 0;
		  }
		  else if (rhs_temp2->type == TERMALPHID){
		    if (a_i  == ((tg_alph*)(rhs_temp2->value))->alph) 
		      return 0;
		  }
		}
	      }
	    }

	    if (rhs_temp->type == TERMALPHID) {
	      a_i = ((tg_alph*)(rhs_temp->value))->alph;
	      if (rhs_index > 0) {
		for (i=rhs_index-1; i>old_rhs_index; i--) 
		  {
		    rhs_temp2 =(tg_idtab_entry*)(tg_rule_temp->rhs[i]);
		  if (rhs_temp2->type == TERMID){
		    if (a_i  == ((tg_sym*)(rhs_temp2->value))->alph) 
		      return 0;
		  }
		  
		  else if (rhs_temp2->type == TERMALPHID){
		    if (a_i  == ((tg_alph*)(rhs_temp2->value))->alph) 
		      return 0;
		  }
		}
	      }
	      }
	  }
    }
    /* at this point we will have a grammar in RNA normal form */
  }
  return 1;
}
