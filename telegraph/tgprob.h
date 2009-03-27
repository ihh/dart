#ifndef TELEGRAPH_PARAMETER_INCLUDED
#define TELEGRAPH_PARAMETER_INCLUDED

/* Telegraph arrays */
typedef struct tg_array_str
{
  tg_idtab_index  id;    /* index of identifier */
  tg_prob_index   index; /* index of this parameter */
  tg_unsigned_vec dim;   /* dimension of each subscript */
  tg_alph_vec     label; /* label of each subscript, or 0 if unlabeled */
} tg_array;

/* A 'prob' is represented as an array */
typedef tg_array tg_prob;
typedef tg_vector(tg_prob*) tg_prob_vec;

/* Probability array dimensions */
typedef struct tg_dim_str
{
  tg_alph*     label; /* null if subscript lacks an alphabet label */
  unsigned int size;  /* the dimension of this subscript */
} tg_dim;

/* Probability subscripts (bound) */
typedef unsigned int tg_const_subscript;
typedef tg_vector(tg_const_subscript) tg_const_subscript_vec;

/* Probability subscripts (unbound)
   Can include both constant integer subscripts, e.g. X[3],
   and roaming "dollar" subscripts, e.g. X[$2].
*/
typedef struct tg_subscript_str
{
  /* flag; true for $1, $2 etc */
  char dollar;

  /* the subscript, or N-1 for $N */
  union tg_subscript_union {
    tg_const_subscript int_subscript;     /* used if dollar==0 */
    tg_rhs_index       dollar_subscript;  /* used if dollar!=0 */
  } index;

} tg_subscript;

typedef tg_vector(tg_subscript*) tg_subscript_vec;

/* Named functions */
typedef tg_prob tg_func;
typedef tg_vector(tg_func*) tg_func_vec;

/* Subscript and dimension builders */
tg_dim* tg_new_dim (tg_alph* alph, int size);
tg_subscript* tg_new_subscript (char dollar, tg_expr_index index);

#endif /* TELEGRAPH_PARAMETER_INCLUDED */
