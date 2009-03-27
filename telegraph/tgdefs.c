#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "tgcompiler.h"
#include "gram.tab.h"

extern void yyerror (char*);

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

tg_dim* tg_new_dim (tg_alph* alph, int size) {
    tg_dim* dim;
    dim = malloc (sizeof (tg_dim));
    dim->label = alph;
    dim->size = size;
    return dim;
}

tg_subscript* tg_new_subscript (char dollar, int index) {
    tg_subscript* sub;
    sub = malloc (sizeof (tg_subscript));
    sub->dollar = dollar;
    sub->index.int_subscript = index;
    return sub;
}

tg_grammar* tg_new_grammar() {
    tg_grammar* tg;
    tg = malloc (sizeof (tg_grammar));
    tg->alph = (tg_alph_vec) tg_vec_new (0);
    tg->term = (tg_term_vec) tg_vec_new (0);
    tg->nonterm = (tg_symbol_vec) tg_vec_new (0);
    tg->prob = (tg_prob_vec) tg_vec_new (0);
    tg->expr = (tg_expression_vec) tg_vec_new (0);
    tg->rule = (tg_rule_vec) tg_vec_new (0);
    tg->idtab = (tg_type_table) tg_vec_new (0);
    tg->func = (tg_func_vec) tg_vec_new (0);
    tg->assign = (tg_expression_vec) tg_vec_new (0);
    tg->pexpr = (tg_pexpr_vec) tg_vec_new (0);
    return tg;
}

void tg_delete_grammar (tg_grammar* tg) {
    /* MEMORY LEAK -- should free the objects here */
    tg_vec_delete ((tg_vec*) &tg->alph);
    tg_vec_delete ((tg_vec*) &tg->term);
    tg_vec_delete ((tg_vec*) &tg->nonterm);
    tg_vec_delete ((tg_vec*) &tg->prob);
    tg_vec_delete ((tg_vec*) &tg->expr);
    tg_vec_delete ((tg_vec*) &tg->rule);
    tg_vec_delete ((tg_vec*) &tg->idtab);
    tg_vec_delete ((tg_vec*) &tg->func);
    tg_vec_delete ((tg_vec*) &tg->assign);
    tg_vec_delete ((tg_vec*) &tg->pexpr);
    free (tg);
}

tg_type tg_lookup_type (tg_grammar* tg, const char* name, void** value_ret) {
    int i;
    tg_idtab_entry* entry;
    for (i = 0; i < tg_size(tg->idtab); ++i) {
	entry = (tg_idtab_entry*) tg->idtab[i];
	if (strcmp (name, entry->name) == 0) {
	    *value_ret = entry->value;
	    return entry->type;
	}
    }
    *value_ret = (void*) strdup (name);
    return IDENTIFIER;
}

tg_idtab_index tg_set_type (tg_grammar* tg, const char* name,
			    tg_type type, void* value) {
    tg_idtab_entry* entry;

    /* fprintf (stderr, "Setting type of %s to %d\n", name, type); */

    entry = malloc (sizeof (tg_idtab_entry));
    entry->name = name;
    entry->type = type;
    entry->value = value;
    return tg_vec_append ((tg_vec*) &tg->idtab, (void*) entry);
}

const char* tg_get_name (tg_grammar* tg, tg_idtab_index id) {
  tg_idtab_entry* entry;
  entry = (tg_idtab_entry*) tg->idtab[id];
  return entry->name;
}

tg_expr* tg_make_binary_expr (tg_grammar* tg, char op, void* op1, void* op2) {
    tg_expr* expr;
    expr = malloc (sizeof (tg_expr));
    expr->op = op;
    expr->operand1 = op1;
    expr->operand2 = op2;
    return expr;
}

tg_expr* tg_make_const_expr (tg_grammar* tg, double r) {
    tg_score score;
    score = tg_prob2score (r);
    return tg_make_binary_expr (tg, tg_Constant, (void*) score, (void*) 0);
}

tg_expr* tg_make_prob_expr (tg_grammar* tg, tg_prob* prob, tg_vec subscript) {
    return tg_make_binary_expr (tg, tg_ProbBinding, (void*) prob, (void*) subscript);
}

tg_expr* tg_make_func_expr (tg_grammar* tg, tg_func* func) {
    return tg_make_binary_expr (tg, tg_Func, (void*) func, (void*) 0);
}

tg_expr* tg_make_psum_expr (tg_grammar* tg, tg_pexpr* pexpr) {
    return tg_make_binary_expr (tg, tg_PSum, (void*) pexpr, (void*) 0);
}

tg_alph* tg_new_alphabet (tg_grammar* tg) {
    tg_alph* alph;
    alph = malloc (sizeof (tg_alph));
    alph->token = (tg_symbol_vec) tg_vec_new (0);
    alph->alph = tg_vec_append ((tg_vec*) &tg->alph, (void*) alph);
    return alph;
}

void tg_add_token (tg_grammar* tg, tg_alph* alph, const char* name) {
    tg_sym* sym;
    int i, j, match;
    const char* prevstr;

    /* Check that the terminal string isn't a prefix/suffix of any others */
    for (i = 0; i < tg_size(alph->token); ++i) {
      prevstr = tg_get_name (tg, alph->token[i]->id);
      match = 1;
      for (j = 0; match && j < strlen(prevstr) && j < strlen(name); ++j)
	  if (name[i] != prevstr[j])
	      match = 0;
      if (match) {
	  yyerror ("Token prefix overlap\n");
	  fprintf (stderr, "Tokens \"%s\" and \"%s\" clash\n",
		   prevstr, name);
	  exit(1);
      }
    }

    sym = malloc (sizeof (tg_sym));
    sym->sym = tg_vec_append ((tg_vec*) &alph->token, (void*) sym);
    sym->alph = alph->alph;
    sym->id = tg_set_type (tg, name, TOKENID, (void*) sym);
}

void tg_declare_alphabet (tg_grammar* tg, tg_alph* alph, const char* name) {
    alph->id = tg_set_type (tg, name, TOKALPHID, (void*) alph);
}

tg_term* tg_declare_terminal (tg_grammar* tg, tg_alph* alph,
			      const char* name) {
    tg_term* term;
    tg_sym* sym;
    const char *tokid;
    char* fake_id;
    int i;

    term = malloc (sizeof (tg_term));
    term->id = tg_set_type (tg, name, TERMALPHID, (void*) term);
    term->term = tg_vec_append ((tg_vec*) &tg->term, (void*) term);
    term->tok = alph->alph;

    term->termsym = (tg_symbol_vec) tg_vec_new (0);
    for (i = 0; i < tg_size (alph->token); ++i) {
	sym = malloc (sizeof (tg_sym));
	sym->alph = term->term;
	sym->sym = i;

	/* hack: put fake entry for terminal into symbol table */
	tokid = tg_get_name (tg, alph->token[i]->id);
	fake_id = malloc (strlen(name) + strlen(tokid) + 3);
	sprintf (fake_id, "%s[%s]", name, tokid);
	sym->id = tg_set_type (tg, fake_id, TERMID, (void*) sym);

	tg_vec_append ((tg_vec*) &term->termsym, (void*) sym);
    }

    return term;
}

tg_sym* tg_get_terminal (tg_grammar* tg, tg_term* term, tg_sym* tok) {
    if (tok->alph != term->tok) {
	yyerror ("Terminal index does not match token alphabet");
	exit(1);
    }
    return term->termsym[tok->sym];
}

void tg_declare_nonterminal (tg_grammar* tg, const char* name) {
    tg_sym* sym;
    sym = malloc (sizeof (tg_sym));
    sym->sym = tg_vec_append ((tg_vec*) &tg->nonterm, (void*) sym);
    sym->alph = -1;
    sym->id = tg_set_type (tg, name, NONTERMID, (void*) sym);
}

tg_prob* tg_declare_prob (tg_grammar* tg, tg_vec dim, const char* name) {
    int i;
    tg_prob* prob;
    tg_dim* subscript_dim;

    prob = malloc (sizeof (tg_prob));
    prob->index = tg_vec_append ((tg_vec*) &tg->prob, (void*) prob);
    prob->id = tg_set_type (tg, name, PROBID, (void*) prob);

    prob->dim = (tg_unsigned_vec) tg_vec_new (tg_size(dim));
    prob->label = (tg_alph_vec) tg_vec_new (tg_size(dim));
    for (i = 0; i < tg_size(dim); ++i) {
	subscript_dim = (tg_dim*) dim[i];
	tg_vec_append ((tg_vec*) &prob->dim, (void*) subscript_dim->size);
	tg_vec_append ((tg_vec*) &prob->label, (void*) subscript_dim->label);
    }

    return prob;
}

void tg_declare_func (tg_grammar* tg, const char* name, tg_expr* rhs) {
    tg_func* result;
    tg_expr *lhs, *assign;

    result = malloc (sizeof (tg_func));
    result->index = tg_vec_append ((tg_vec*) &tg->func, (void*) result);
    result->id = tg_set_type (tg, name, FUNCID, (void*) result);

    result->dim = (tg_unsigned_vec) 0;
    result->label = (tg_alph_vec) 0;

    lhs = tg_make_func_expr (tg, result);
    assign = tg_make_binary_expr (tg, tg_ProbAssign, lhs, rhs);
    tg_vec_append ((tg_vec*) &tg->assign, assign);
}

void tg_assign_expr (tg_grammar* tg, tg_expr* lhs, tg_expr* rhs) {
  tg_expr* assign;
  assign = tg_make_binary_expr (tg, tg_ProbAssign, lhs, rhs);
  tg_vec_append ((tg_vec*) &tg->assign, (void*) assign);
}

void tg_assign_prob (tg_grammar* tg, tg_expr* lhs, double val) {
  tg_expr* rhs;
  rhs = tg_make_const_expr (tg, val);
  tg_assign_expr (tg, lhs, rhs);
}

void tg_assign_log2prob (tg_grammar* tg, tg_expr* lhs, double log2_val) {
  tg_assign_prob (tg, lhs, pow (2., log2_val));
}

int tg_subscripts_valid (tg_grammar* tg, tg_vec rhs, tg_expr* expr) {
  tg_prob* prob;
  tg_vec vec;
  tg_subscript* sub;
  unsigned int dim;
  tg_idtab_entry* entry;
  tg_alph* alph;
  tg_term* term;
  int i, valid;

  valid = 1;
  switch (expr->op) {
  case tg_Addition:
  case tg_Multiplication:
  case tg_Division:
    valid = valid && tg_subscripts_valid (tg, rhs, expr->operand1);
    valid = valid && tg_subscripts_valid (tg, rhs, expr->operand2);
    break;

  case tg_ProbBinding:
    prob = (tg_prob*) expr->operand1;
    vec = (tg_vec) expr->operand2;
    valid = valid && (tg_size(vec) == tg_size(prob->dim));

    /* check each subscript */
    for (i = 0; valid && i < tg_size(vec); ++i) {
      dim = (unsigned int) prob->dim[i];
      sub = (tg_subscript*) vec[i];
      if (sub->dollar) {
	  valid = valid && rhs;
	valid = valid && (sub->index.dollar_subscript >= 0
			  && sub->index.dollar_subscript < tg_size(rhs));
	if (valid) {
	  entry = (tg_idtab_entry*) rhs[sub->index.dollar_subscript];
	  valid = valid && (entry->type == TERMALPHID);
	  if (valid) {
	      term = (tg_term*) entry->value;
	      alph = tg->alph[term->tok];
	      valid = valid && (tg_size(alph->token) == dim);
	  }
	}
      } else /* sub->dollar */
	valid = valid && (sub->index.int_subscript >= 0
			  && sub->index.int_subscript < dim);
    }
    break;
  default:
    break;
  }

  return valid;
}

void tg_declare_rule (tg_grammar* tg,
		      tg_sym* lhs, tg_vec rhs, tg_expr* expr)
{
  tg_rule* rule;

  if (!tg_subscripts_valid (tg, rhs, expr)) {
    yyerror ("Bad prob subscript in rule probability expression");
    exit(1);
  }

  rule = malloc (sizeof (tg_rule));
  rule->index = tg_size (tg->rule);
  rule->lhs = lhs;
  rule->rhs = (tg_identifier_vec) rhs;
  rule->expr = expr;

  tg_vec_append ((tg_vec*) &tg->rule, (void*) rule);
}


void tg_write_grammar (tg_grammar* tg, FILE* file) {
  int i, j, indent, lastindent;
  tg_alph* alph;
  tg_sym* sym;
  tg_prob* prob;
  unsigned int subscript_dim;
  tg_alph* subscript_label;
  tg_rule* rule;
  tg_sym* lastlhs;
  tg_idtab_entry* rhs_entry;

  /* Symbols */
  fprintf (file, "/* Symbols */\n");

  /* token */
  if (tg_size (tg->alph)) {
      fprintf (file, "token");
      for (i = 0; i < tg_size(tg->alph); ++i) {
	  alph = (tg_alph*) tg->alph[i];
	  fprintf (file, "%s%s {",
		   i == 0 ? " " : ", ",
		   tg_get_name (tg, alph->id));
	  for (j = 0; j < tg_size(alph->token); ++j) {
	      sym = (tg_sym*) alph->token[j];
	      fprintf (file, " %s", tg_get_name (tg, sym->id));
	  }
	  fprintf (file, " };\n");
      }
  }

  /* term */
  if (tg_size (tg->term)) {
      fprintf (file, "term");
      for (i = 0; i < tg_size (tg->term); ++i)
	  fprintf (file, "%s%s[%s]",
		   i == 0 ? " " : ", ",
		   tg_get_name (tg, tg->term[i]->id),
		   tg_get_name (tg, tg->alph[tg->term[i]->tok]->id));
      fprintf (file, ";\n");
  }

  /* nonterm */
  if (tg_size(tg->nonterm)) {
    fprintf (file, "nonterm ");
    for (i = 0; i < tg_size(tg->nonterm); ++i) {
      sym = (tg_sym*) tg->nonterm[i];
      if (i > 0)
	fprintf (file, ", ");
      fprintf (file, tg_get_name (tg, sym->id));
    }
    fprintf (file, ";\n");
  }

  /* Parameters */
  fprintf (file, "/* Parameters */\n");

  /* prob */
  if (tg_size(tg->prob)) {
    fprintf (file, "prob ");
    for (i = 0; i < tg_size(tg->prob); ++i) {
      prob = (tg_prob*) tg->prob[i];
      if (i > 0)
	fprintf (file, ", ");
      fprintf (file, tg_get_name (tg, prob->id));
      for (j = 0; j < tg_size(prob->dim); ++j) {
	subscript_dim = (unsigned int) prob->dim[j];
	subscript_label = (tg_alph*) prob->label[j];
	if (subscript_label) {
	  fprintf (file, "[%s]", tg_get_name (tg, subscript_label->id));
	} else
	  fprintf (file, "[%d]", subscript_dim);
      }
    }
    fprintf (file, ";\n");
  }

  /* rules */
  if (tg_size(tg->rule)) {
      fprintf (file, "/* Rules */\n");

    lastlhs = 0;
    for (i = 0; i < tg_size(tg->rule); ++i) {
      rule = (tg_rule*) tg->rule[i];
      indent = 0;
      /* are we printing a new rule, "lhs -> rhs"
	 or a continuation?          "     | rhs"
      */

      if (rule->lhs == lastlhs) {   /* continuation */
	  /* terminate last line */
	  if (lastlhs)
	      fprintf (file, "\n");

	/* print rule number (as a comment) */
	  fprintf (file, "/*%d*/ ", i + 1);
	  indent += 6;
	  if (i >= 9) ++indent;
	  if (i >= 99) ++indent;

	/* indent */
	for ( ; indent < lastindent - 1; ++indent)
	    fprintf (file, " ");
	fprintf (file, "|");
	++indent;
	if (indent > lastindent)
	    lastindent = indent;
	
      } else {                         /* new rule */
	  /* terminate last line */
	  if (lastlhs)
	      fprintf (file, ";\n");

	/* print rule number (as a comment) */
	  fprintf (file, "/*%d*/ ", i + 1);
	  indent += 6;
	  if (i >= 9) ++indent;
	  if (i >= 99) ++indent;

	/* print LHS */
	fprintf (file, "%s ->", tg_get_name (tg, rule->lhs->id));
	indent += strlen (tg_get_name (tg, rule->lhs->id)) + 3;

	/* retain info about LHS */
	lastlhs = rule->lhs;
	lastindent = indent;
      }

      /* print RHS */
      if (tg_size(rule->rhs))
	for (j = 0; j < tg_size(rule->rhs); ++j) {
	  rhs_entry = (tg_idtab_entry*) rule->rhs[j];
	  fprintf (file, " %s", rhs_entry->name);	
	}
      else
	fprintf (file, " end");

      /* print rule probability, unless it's equal to 1 */
      if (rule->expr->op != tg_Constant
	  || (tg_score) rule->expr->operand1 != 0) {
	fprintf (file, " { ");
	tg_write_expr (tg, rule->expr, file);
	fprintf (file, " }");
      }
    }
    if (lastlhs)
	fprintf (file, ";\n");
  }

  /* Assignments & results */
  if (tg_size(tg->assign)) {
      fprintf (file, "/* Assignments & computations */\n");
      
      for (i = 0; i < tg_size(tg->assign); ++i) {
	tg_write_expr (tg, tg->assign[i], file);
	fprintf (file, ";\n");
      }
  }
}

void tg_write_expr (tg_grammar* tg, tg_expr* expr, FILE* file) {
  int i;
  tg_expr* expr1;
  tg_expr* expr2;
  tg_prob* prob;
  tg_func* func;
  tg_vec vec;
  tg_subscript* sub;
  tg_score score;
  tg_idtab_index id;
  tg_pexpr* pexpr;

  switch (expr->op) {
  case tg_Multiplication:
  case tg_Division:
    expr1 = (tg_expr*) expr->operand1;
    expr2 = (tg_expr*) expr->operand2;
    if (expr1->op == tg_Addition) fprintf (file, "(");
    tg_write_expr (tg, expr1, file);
    if (expr1->op == tg_Addition) fprintf (file, ")");
    fprintf (file, " %c ", expr->op);
    if (expr2->op == tg_Addition) fprintf (file, "(");
    tg_write_expr (tg, expr2, file);
    if (expr2->op == tg_Addition) fprintf (file, ")");
    break;

  case tg_Addition:
    expr1 = (tg_expr*) expr->operand1;
    expr2 = (tg_expr*) expr->operand2;
    tg_write_expr (tg, expr1, file);
    fprintf (file, " + ");
    tg_write_expr (tg, expr2, file);
    break;

  case tg_ProbAssign:
    expr1 = (tg_expr*) expr->operand1;
    expr2 = (tg_expr*) expr->operand2;
    if (expr1->op == tg_Func)
	fprintf (file, "func ");  /* hack */
    tg_write_expr (tg, expr1, file);
    if (expr2->op == tg_Constant
	&& (score = -(tg_score) expr2->operand1) > 0) {
	fprintf (file, " => %d.%3d",
		 score / (int) tg_score2bit_ratio,
		 score % (int) tg_score2bit_ratio);
    } else {
	fprintf (file, " = ");
	tg_write_expr (tg, expr2, file);
    }
    break;

  case tg_ParseAssign:
      id = (tg_idtab_index) expr->operand1;
      pexpr = (tg_pexpr*) expr->operand2;
      fprintf (file, "pset %s = ", tg_get_name (tg, id));
      tg_write_pexpr (tg, pexpr, file);
      break;
    
  case tg_ProbBinding:
    prob = (tg_prob*) expr->operand1;
    vec = (tg_vec) expr->operand2;
    fprintf (file, tg_get_name (tg, prob->id));
    for (i = 0; i < tg_size(vec); ++i) {
      sub = (tg_subscript*) vec[i];
      if (sub->dollar)
	fprintf (file, "[$%d]", sub->index.dollar_subscript + 1);
      else if (prob->label[i])
	fprintf (file, "[%s]", tg_get_name
		 (tg, prob->label[i]->token[sub->index.int_subscript]->id));
      else
	fprintf (file, "[%d]", sub->index.int_subscript);
    }
    break;

  case tg_Func:
    func = (tg_func*) expr->operand1;
    fprintf (file, tg_get_name (tg, func->id));
    break;

  case tg_Constant:
    score = (tg_score) expr->operand1;
    fprintf (file, "%g", tg_score2prob (score));
    break;

  case tg_PSum:
    pexpr = (tg_pexpr*) expr->operand1;
    fprintf (file, "psum (");
    tg_write_pexpr (tg, pexpr, file);
    fprintf (file, ")");
    break;

  case tg_Derivative:
    expr1 = (tg_expr*) expr->operand1;
    expr2 = (tg_expr*) expr->operand2;
    fprintf (file, "dlog (");
    tg_write_expr (tg, expr1, file);
    fprintf (file, ") / dlog (");
    tg_write_expr (tg, expr2, file);
    fprintf (file, ")");
    break;

  default:
    fprintf (file, "???");
    break;
  }
}

void tg_write_parse (tg_grammar* tg, tg_parse_node* node, FILE* file) {
    int i;

    if (node->child) {
	fprintf (file, tg_get_name (tg, node->label.rule->lhs->id));
	if (tg_size (node->child)) {
	    fprintf (file, "(");
	    for (i = 0; i < tg_size (node->child); ++i) {
		if (i > 0)
		    fprintf (file, " ");
		tg_write_parse (tg, node->child[i], file);
	    }
	    fprintf (file, ")");
	}
    } else
	fprintf (file, tg_get_name (tg, node->label.term->id));
}

void tg_write_seq (tg_grammar* tg, tg_parse_node* node, char** seq) {
	int i, j, k, s_i, max; 
	int term_count[tg_size(tg->term)];
	for (k = 0; k < tg_size(tg->term); k++) {
		term_count[k] = 0;
	}
	if (node->child) {
		if (tg_size(node->child)) {
			for (i = 0; i < tg_size(node->child); i++) {	
				if (node->child[i]->child) {
					if (tg_size(node->child[i]->child)) {
						tg_write_seq (tg, node->child[i], seq);
					} else {
						return; /*case for end terminal*/
				   }
				} else {
					for ( ; i < tg_size(node->child); i++) {
						if (node->child[i]->child) {
							if (tg_size(node->child[i]->child)) { /*case for end of terminal block*/
								max = term_count[0];
								/*fill in gaps*/
								for (j = 0; j < tg_size(tg->term); j++) {
									if (term_count[j] > max) 
										max = (term_count[j]);
								}
								for (j = 0; j < tg_size(tg->term); j++) {
									for (k = term_count[j]; k < max; k++) {
										seq[j] = strcat(seq[j], "-");
									}
								}
								--i;
								break;
							} else { /*case for end terminal*/
								return;
							}
						} else { /*case for processing terminal block*/
						s_i = node->child[i]->label.term->alph;
						seq[s_i] = strcat(seq[s_i],
					   tg_get_name (tg, tg->alph[tg->term[s_i]->tok]->token[node->child[i]->label.term->sym]->id));
						term_count[s_i]++;
						}
					}
				}
			}
		}
	}
}	

void tg_write_stockholm (tg_grammar* tg, tg_parse_node* node, FILE* file) {
	/* to hold sequences */
	char** seq = malloc (tg_size(tg->term)*sizeof(char *));
	int i;
	if (node->child)
		fprintf(file, "# STOCKHOLM 1.0 \n");
		for (i = 0; i < tg_size(tg->term); i++) {
			seq[i] = malloc(sizeof(char));
			*seq[i] = '\0';
			/*fills seq w/ empty strings*/
		}
	/* print to file */
	if (tg_size (node->child)) {
		tg_write_seq (tg, node, seq);	
		for (i = 0; i < tg_size(tg->term); i++) {
		fprintf(file, "%s\t%s\n", tg_get_name (tg, tg->term[i]->id), seq[i]);
		}
	}
	fprintf(file, "// \n");
}

void tg_write_pexpr (tg_grammar* tg, tg_pexpr* pexpr, FILE* file) {
    int i, j, flag;
    tg_dynamic_pexpr* dynamic;
    tg_alph* alph;
    tg_subseq_vec seqvec;
    tg_parse_vec pvec;

    switch (pexpr->op) {
    case PSAMPLE:
	fprintf (file, "psample (");
	tg_write_pexpr (tg, (tg_pexpr*) pexpr->operand, file);
	fprintf (file, ")");
	break;

    case PMAX:
	fprintf (file, "pmax (");
	tg_write_pexpr (tg, (tg_pexpr*) pexpr->operand, file);
	fprintf (file, ")");
	break;

    case PSET:
	pvec = (tg_parse_vec) pexpr->operand;
	fprintf (file, "{");
	for (i = 0; i < tg_size (pvec); ++i) {
	    fprintf (file, i > 0 ? ", " : " ");
	    tg_write_parse (tg, pvec[i], file);
	}
	fprintf (file, " }");
	break;

    case PARSE:
	dynamic = (tg_dynamic_pexpr*) pexpr->operand;
	seqvec = dynamic->seq;
	flag = 0;
	fprintf (file, "parse (%s ->",
		 tg_get_name (tg, tg->nonterm[dynamic->start]->id));
	for (i = 0; i < tg_size(dynamic->seq); ++i)
	    if (seqvec[i]) {
		fprintf (file, "%s%s=\"", flag++ ? ", " : " ",
			 tg_get_name (tg, tg->term[i]->id));
		alph = tg->alph[tg->term[i]->tok];
		for (j = 0; j < seqvec[i]->len; ++j)
		    fprintf (file, tg_get_name (tg, alph->token[seqvec[i]->begin[j]]->id));
		fprintf (file, "\"");
	    }
	fprintf (file, ")");
    break;

    case PSETID:
	fprintf (file, tg_get_name (tg, (tg_idtab_index) pexpr->operand));
	break;

    default:
	fprintf (file, "???");
	break;
    }
}

tg_parse_node* tg_new_parse_node (tg_grammar* tg, tg_sym* label,
				  tg_parse_vec child)
{
    tg_parse_node *node, *kid;
    tg_rule* rule;
    int i, j;

    node = malloc (sizeof (tg_parse_node));
    node->child = child;
    if (child) {
	/* We have a (possibly null) child list, so this is a nonterminal node.
	   Let's try to match the rule.
	*/
	rule = 0; /* zero indicates that rule is still unmatched */
	for (j = 0; !rule && j < tg_size(tg->rule); ++j) {
	    rule = (tg_rule*) tg->rule[j];

	    /* check LHS */
	    if (rule->lhs != label) {
		rule = 0;
		continue;
	    }

	    /* check size of RHS */
	    if (tg_size(child) != tg_size(rule->rhs)) {
		rule = 0;
		continue;
	    }

	    /* try to match each kid with corresponding RHS rule symbol */
	    for (i = 0; rule && i < tg_size(child); ++i) {
		/* Check that child node type matches production rule */
		kid = child[i];
		if (kid->child == 0) {
		    /* kid is a terminal symbol */
		    if (rule->rhs[i]->type == TERMID) {
			/* rule expects a specific terminal */
			if ((tg_sym*) rule->rhs[i]->value != kid->label.term)
			    /* kid symbol didn't match */
			    rule = 0;
			
		    } else if (rule->rhs[i]->type == TERMALPHID) {
			/* rule expects a terminal from a specific alphabet */
			if (((tg_alph*) rule->rhs[i]->value)->alph
			    != kid->label.term->alph)
			    /* kid's alphabet didn't match */
			    rule = 0;
			
		    } else
			/* rule expects a nonterminal */
			rule = 0;

		} else {
		    /* child is a nonterminal symbol */
		    if ((tg_sym*) rule->rhs[i]->value
			!= kid->label.rule->lhs)
			/* rule expects a different symbol */
			rule = 0;
		}
	    }
	}

	if (!rule) {
	    yyerror ("Couldn't match rule");
	    for (j = 0; j < 2; ++j) {
		fprintf (stderr, "%s%s%s",
			 j == 0
			 ? "When parsing \""
			 : ")\":\nCouldn't find rule ",
			 tg_get_name (tg, label->id),
			 j == 0 ? " (" : " ->");
		for (i = 0; i < tg_size(child); ++i) {
		    if (i > 0 || j > 0)
			fprintf (stderr, " ");
		    if (child[i]->child)
			fprintf (stderr, tg_get_name
				 (tg, child[i]->label.rule->lhs->id));
		    else
			fprintf (stderr,
				 tg_get_name (tg, child[i]->label.term->id));
		}
	    }
	    if (tg_size(child) == 0)
		fprintf (stderr, " end");
	    fprintf (stderr, "\n");
	    exit(1);
	}

	/* found the rule */
	node->label.rule = rule;

	/* Double-link tree */
	for (i = 0; rule && i < tg_size(child); ++i)
	    child[i]->parent = node;

    } else  /* no child list, so this is a terminal node */
	node->label.term = label;

    return node;
}

void tg_free_parse_tree (tg_parse_node* root)
{
    int i;
    if (root) {
	if (root->child) {
	    for (i = 0; i < tg_size(root->child); ++i)
		tg_free_parse_tree ((tg_parse_node*) root->child[i]);
	    tg_vec_delete ((tg_vec*) &root->child);
	}
	free (root);
    }
}

tg_pexpr* tg_new_pexpr (tg_grammar* tg, int op, void* operand) {
    tg_pexpr* pexpr;
    pexpr = malloc (sizeof (tg_pexpr));
    pexpr->op = op;
    pexpr->operand = operand;
    tg_vec_append ((tg_vec*) &tg->pexpr, pexpr);
    return pexpr;
}

void tg_declare_pset (tg_grammar* grammar, tg_pexpr* pexpr,
		       const char* name)
{
    tg_expr* passign;
    passign = tg_make_binary_expr (grammar, tg_ParseAssign,
				   (void*) 0, (void*) pexpr);
    passign->operand1 = (void*) tg_set_type (grammar, name, PSETID, passign);
    tg_vec_append ((tg_vec*) &grammar->assign, passign);
}

tg_sequence tg_lex_string (tg_grammar* tg, const char* input_string,
			   tg_alph* alphabet) {
  tg_sequence seq;
  int i;
  tg_sym_index sym;
  const char* tokstr;

  seq = (tg_sequence) tg_vec_new(0);
  if (input_string)
      while (*input_string != 0) {
	  for (sym = 0; sym < tg_size (alphabet->token); ++sym) {
	      tokstr = tg_get_name (tg, alphabet->token[sym]->id);
	      for (i = 0; tokstr[i] != 0; ++i)
		  if (tokstr[i] != input_string[i])
		      break;
	      if (tokstr[i] == 0)
		  break;
	  }

	  if (sym == tg_size (alphabet->token)) {
	      yyerror ("Can't lex string");
	      fprintf (stderr, "Unlexable suffix: \"%s\"\n", input_string);
	      exit(1);
	  }
	  tg_vec_append ((tg_vec*) &seq, (void*) sym);
	  input_string += i;
      }
 
 return seq;
}

void tg_assert_is_psum (tg_grammar* tg, tg_expr* result_expr) {
    tg_expr *expr, *expr2;
    tg_func *fpresult, *eresult;
    int i;
    fpresult = (tg_func*) result_expr->operand1;
    for (i = 0; i < tg_size (tg->assign); ++i) {
	/* if not a prob assignment, move on */
	if (tg->assign[i]->op != tg_ProbAssign)
	    continue;
	expr = (tg_expr*) tg->assign[i]->operand1;
	eresult = (tg_func*) expr->operand1;
	/* if results don't match, move on */
	if (eresult != fpresult)
	    continue;
	/* if not a psum, throw an error */
	expr2 = (tg_expr*) tg->assign[i]->operand2;
	if (expr2->op != tg_PSum) {
	    yyerror ("Not a psum");
	    exit(1);
	}
	/* return */
	return;
    }
    yyerror ("psum undefined");
    exit(1);
    /* unreachable */
}

tg_sequence_vec tg_new_binding (tg_grammar* tg) {
  return (tg_sequence_vec) tg_vec_init (tg_size (tg->term), 0);
}

void tg_bind_sequence (tg_grammar* tg, tg_sequence_vec binding,
		       tg_term* term, const char* s)
{
  if (binding[term->term]) {
    yyerror ("Duplicate sequence binding");
    exit(1);
  }
  binding[term->term] = tg_lex_string (tg, s, tg->alph[term->tok]);
}

void tg_check_binding (tg_grammar* tg, tg_sequence_vec binding)
{
  int i;
  if (tg_size (binding) != tg_size (tg->term)) {
    yyerror ("Bad binding");
    exit(1);
  }
  for (i = 0; i < tg_size (tg->term); ++i)
    if (!binding[i]) {
      yyerror ("Sequence unbound");
      fprintf (stderr, "Terminal sequence '%s' is unbound\n",
	       tg_get_name (tg, tg->term[i]->id));
      exit(1);
    }
}

tg_pexpr* tg_make_dynamic_pexpr (tg_grammar* tg, tg_sequence_vec seqs,
				tg_sym_index start)
{
    tg_pexpr* pexpr;
    tg_dynamic_pexpr* dynamic;
    tg_subseq* subseq;
    int i;

    tg_check_binding (tg, seqs);

    dynamic = malloc (sizeof (tg_dynamic_pexpr));
    dynamic->start = start;
    dynamic->seq = (tg_subseq_vec) tg_vec_new (0);
    for (i = 0; i < tg_size (seqs); ++i) {
	subseq = malloc (sizeof (tg_subseq));
	subseq->begin = seqs[i];
	subseq->len = tg_size (seqs[i]);
	tg_vec_append ((tg_vec*) &dynamic->seq, (void*) subseq);
    }

    pexpr = malloc (sizeof (tg_pexpr));
    pexpr->op = PARSE;
    pexpr->operand = (void*) dynamic;
    
    return pexpr;
}
