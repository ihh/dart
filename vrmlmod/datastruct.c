#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "datastruct.h"

tg_vec tg_vec_new (int initAlloc) {
    tg_vec vec;
    vec = ((void**) malloc ((initAlloc + 2) * sizeof(void*))) + 2;
    tg_alloc(vec) = initAlloc;
    tg_size(vec) = 0;
    return vec;
}

tg_vec tg_vec_init (int initSize, void* val) {
    tg_vec vec;
    int i;
    vec = tg_vec_new (initSize);
    for (i = 0; i < initSize; ++i)
	tg_vec_append (&vec, val);
    return vec;
}

void tg_vec_delete (tg_vec* vec) {
    if (*vec)
	free (*vec - 2);
    *vec = 0;
}

int tg_vec_append (tg_vec* vec, void* value) {
    unsigned int new_alloc;
    tg_vec newvec;
    int i;

    if (++tg_size(*vec) >= tg_alloc(*vec)) {
	new_alloc = tg_alloc(*vec) * 2;
	if (new_alloc == 0)
	    new_alloc = 1;
	newvec = tg_vec_new (new_alloc);
	tg_size(newvec) = tg_size(*vec);
	for (i = 0; i < tg_size(*vec) - 1; ++i)
	    newvec[i] = (*vec)[i];
	tg_vec_delete (vec);
	*vec = newvec;
    }
    (*vec)[tg_size(*vec)-1] = value;

    return tg_size(*vec) - 1;
}

void tg_vec_clear (tg_vec vec)
{
  tg_size(vec) = 0;
}

void tg_vec_reserve (tg_vec* vec,
		     int alloc)
{
    tg_vec newvec;
    int i;

    if (alloc > tg_alloc(*vec)) {
	newvec = tg_vec_new (alloc);
	tg_size(newvec) = tg_size(*vec);
	for (i = 0; i < tg_size(*vec); ++i)
	    newvec[i] = (*vec)[i];
	tg_vec_delete (vec);
	*vec = newvec;
    }
}

int tg_vecs_equal (tg_vec vec1, tg_vec vec2)
{
    int i;
    if (tg_size(vec1) != tg_size(vec2))
	return 0;
    for (i = 0; i < tg_size(vec1); ++i)
	if (vec1[i] != vec2[i])
	    return 0;
    return 1;
}


/*SS_queue functions*/
ss_elmt* new_ss_elmt(char c, int i) {
   ss_elmt* ret = malloc(sizeof(ss_elmt));
   ret->c = c;
   ret->index = i;
   ret->under = NULL;
   return ret;
}

ss_stack* new_ss_stack(){
   ss_stack* ret = malloc(sizeof(ss_stack));
   ret->top = NULL;
   return ret;
}

ss_elmt* pop(ss_stack* s){
   ss_elmt* ret;
   if (!isempty(s)) {
      ret = s->top;
      s->top = s->top->under;
	  ret->under = NULL; /*protection against data struct corruption*/
	  return ret;
   } else {
      return NULL;
   }
}
   
void push(ss_stack* s, ss_elmt* e){
   e->under = s->top;
   s->top = e;
}

int isempty(ss_stack* s) {
   return (s->top ? 0 : 1);
}

