#ifndef TELEGRAPH_EXPRESSION_INCLUDED
#define TELEGRAPH_EXPRESSION_INCLUDED

/* Simple expressions */
#define tg_Addition       '+'
#define tg_Multiplication '*'
#define tg_Division       '/'
#define tg_ProbBinding    '['
#define tg_Constant       '9'

/* Script expressions */
#define tg_ProbAssign     '='
#define tg_ParseAssign    ':'
#define tg_Func           'f'
#define tg_PSum           'p'
#define tg_Derivative     'd'

  /* Opcode table:
 op1 type     | op2 type         | Pseudo-syntax       | op
--------------+------------------+---------------------+----------------
 tg_expr*     | tg_expr*         | op1 + op2           | tg_Addition
 tg_expr*     | tg_expr*         | op1 * op2           | tg_Multiplication
 tg_expr*     | tg_expr*         | op1 / op2           | tg_Division
 tg_prob*     | tg_subscript_vec | op1[A][G][T][G]     | tg_ProbBinding
 tg_result*   |                  | op1                 | tg_Result
 tg_score     |                  |  ... => op1         | tg_Constant
 tg_expr*     | tg_expr*         | op1 = op2           | tg_ProbAssign
 tg_idtab_index tg_pexpr*        | op1 = op2           | tg_ParseAssign
 tg_pexpr*    |                  | psum(op1)           | tg_PSum
 tg_expr*     | tg_expr*         | dlog(op1)/dlog(op2) | tg_Derivative
  */

  typedef struct tg_expr_str
  {
    char op;
    struct tg_expr_str *operand1, *operand2;
  } tg_expr;

typedef tg_vector(tg_expr*) tg_expression_vec;

/* Recursive function to return the max '$N' value for an expr */
int tg_max_dollar (tg_expr* expr);  /* returns 0 if expr is dollar-free */

#endif /* TELEGRAPH_EXPRESSION_INCLUDED */
