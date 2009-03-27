/***********************************************************
 * tgcheck.h
 * 16 March 2005
 * 
 * Checks if the grammar is in Chomsky Normal Form,
 * RNA Normal Form or right-regular form
 ***********************************************************/

int tg_is_Chomsky (tg_grammar* tg);
int tg_is_RNA_Normal (tg_grammar* tg, FILE* file);

int tg_is_right_regular (tg_grammar* tg);
