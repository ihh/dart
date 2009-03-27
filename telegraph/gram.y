%{
#include <stdlib.h>
#include <stdio.h>
#include "tgcompiler.h"
tg_grammar* tg;

/* flex */
FILE* yyin;
int yylex();

/* functions */
void yyerror(char*);

/* The variable TG holds the address of the Telegraph grammar object */
/* #define tg ((tg_grammar*)TG) */
%}

%union {                int intval;
			double realval;
			char* strval;
			char charval;
			tg_sym* sym;
			tg_expr* expr;
			tg_alph* alph;
			tg_term* term;
			tg_prob* prob;
			tg_dim* dim;
			tg_subscript* subscript;
			tg_vec vec;
			tg_idtab_entry* id;
			tg_sequence_vec binding;
			tg_parse_node* parse;
			tg_parse_vec parsevec;
			tg_func* func;
			tg_pexpr* pexpr;
			double* realvalptr;
}



%type <realval> nonnegative_constant 
%type <realvalptr> prob_list_element
%type <intval>  positive_integer_constant
%type <expr>  function sum_function atomic_function
%type <expr>  product_function function_list
%type <expr>  prob_expr fixed_prob_expr dlog_denominator optional_calc
%type <subscript> subscript fixed_subscript
%type <dim> dim
%type <vec>  subscript_list fixed_subscript_list dim_list prob_list 
%type <alph> new_token_list
%type <sym>  lhs_symbol production_rule terminal_symbol
%type <id>   rhs_symbol
%type <vec>  rhs_symbols rhs_symbols_or_end

%type <parse> const_parse
%type <parsevec> const_parse_list child_list child_desc

%type <binding> sequence_binding
%type <expr>  psum_expr deriv_expr func_expr psum_func_expr
%type <pexpr> dynamic_pexpr const_pexpr id_pexpr psample_pexpr pmax_pexpr pexpr

%token <intval> ZERO ONE PLURAL
%token <strval> STRVAL
%token <realval> REAL INF

%token <strval> IDENTIFIER
%token <sym> NONTERMID TOKENID TERMID
%token <alph> TOKALPHID
%token <term> TERMALPHID
%token <prob> PROBID
%token <func> FUNCID
%token <expr> PSETID

%token TOKEN TERM NONTERM PROB
%token PSET PARSE PMAX PSAMPLE
%token FUNC PSUM DLOG
%token PRODUCES EQUALS_BITS END
%token TGBADCHAR


%start file

%%

nonnegative_constant:
		ZERO                       { $$ = 0.; }
	|	positive_integer_constant  { $$ = (double) $1; }
	|	REAL
	|	INF  { $$ = (double) tg_log2_infinity; }
		;

positive_integer_constant:
		ONE    { $$ = 1; }
	|	PLURAL { $$ = $1; }
	;

atomic_function:
		prob_expr

	|	nonnegative_constant
		{ $$ = tg_make_const_expr (tg, $1); }

	|       '(' function_list ')'
		{ $$ = $2; }
		;

prob_expr:
		PROBID subscript_list
		{ $$ = tg_make_prob_expr (tg, $1, $2); }
		;

subscript_list:	
		{ $$ = tg_vec_new (0); }

	|	subscript_list '[' subscript ']'
		{ $$ = $1; tg_vec_append (&$$, (void*) $3); }
	;

subscript:
		fixed_subscript

	|	'$' positive_integer_constant
		{ $$ = tg_new_subscript (1, $2 - 1); }
	;


fixed_prob_expr:
		PROBID fixed_subscript_list
		{$$ = tg_make_prob_expr (tg, $1, $2); }
		;

fixed_subscript_list:
		{$$ = tg_vec_new (0); }

	|	fixed_subscript_list '[' fixed_subscript ']'
		{ $$ = $1; tg_vec_append (&$$, (void*) $3); }
	;

fixed_subscript:
		TOKENID  	      { $$ = tg_new_subscript (0, $1->sym); }
	|	ZERO                      { $$ = tg_new_subscript (0, 0); }
	|	positive_integer_constant { $$ = tg_new_subscript (0, $1); }
		;

product_function:
		atomic_function		{ $$ = $1; }
	| 	product_function '*' atomic_function
		{ $$ = tg_make_binary_expr (tg, '*', (void*) $1, (void*) $3); }
	| 	product_function '/' atomic_function
		{ $$ = tg_make_binary_expr (tg, '/', (void*) $1, (void*) $3); }
		;

function_list:
		function
		{ $$ = $1; }
	|      	function_list function
		{ $$ = tg_make_binary_expr (tg, '*', (void*) $1, (void*) $2); }
	;

sum_function:
		product_function
	| 	sum_function '+' product_function
		{ $$ = tg_make_binary_expr (tg, '+', (void*) $1, (void*) $3); }
	;

function:	sum_function;

alphabet_declaration_statement:
		TOKEN alphabet_declaration_list ';'
	;

alphabet_declaration_list:
		alphabet_declaration
	|	alphabet_declaration_list ',' alphabet_declaration
	;

alphabet_declaration:
		IDENTIFIER '{' new_token_list '}'
		{ tg_declare_alphabet (tg, $3, $1); }
	;

new_token_list:
		IDENTIFIER
		{
		  $$ = tg_new_alphabet (tg);
		  tg_add_token (tg, $$, $1);
		}

	|	new_token_list IDENTIFIER
		{
		  $$ = $1;
		  tg_add_token (tg, $$, $2);
		}
	;

terminal_declaration_statement:
		TERM terminal_declaration_list ';'
	;

terminal_declaration_list:
		terminal_declaration
	|	terminal_declaration_list ',' terminal_declaration
	;

terminal_declaration:
		IDENTIFIER '[' TOKALPHID ']'
		{ tg_declare_terminal (tg, $3, $1); }
	;

nonterminal_declaration_statement:
		NONTERM nonterminal_declaration_list ';'
	;

nonterminal_declaration_list:
		IDENTIFIER
		{ tg_declare_nonterminal (tg, $1); }

	|	nonterminal_declaration_list ',' IDENTIFIER
		{ tg_declare_nonterminal (tg, $3); }
	;

prob_declaration_statement:	PROB prob_declaration_list ';'
		;

prob_declaration_list:
		prob_declaration
	|	prob_declaration_list ',' prob_declaration
		;

prob_declaration:
		IDENTIFIER dim_list
		{tg_declare_prob (tg, $2, $1);	}
	;

dim_list:	
		dim_list '[' dim ']'
		{ $$ = $1; tg_vec_append (&$$, (void*) $3); }
	|	{ $$ = tg_vec_new (0); }
	;

dim:		TOKALPHID
		{ $$ = tg_new_dim ($1, tg_size($1->token)); }

	|	positive_integer_constant
		{ $$ = tg_new_dim ((tg_alph*) 0, $1); }
	;


prob_assignment:
		fixed_prob_expr EQUALS_BITS nonnegative_constant ';'
		{ tg_assign_log2prob (tg, $1, -$3); } 

	|	fixed_prob_expr '=' function
		{ tg_assign_expr (tg, $1, $3); }
		
	|	list_assignment
	;

list_assignment: 
		PROBID '=' prob_list ';'
		{ 
		tg_list_assign_prob(tg, $1, $3);
		 }
		;		

prob_list:
		prob_list ',' prob_list_element  {  $$ = $1; tg_vec_append(&$$, (void*) $3); }
		| prob_list_element {$$ = tg_vec_init(1, (void*) $1);}
		;

prob_list_element:
		nonnegative_constant { 
		double* temp = malloc(sizeof(double));
		*temp = $1;
		$$ = temp; }
		;
		

lhs_symbol:	NONTERMID
		;

terminal_symbol:
		TERMALPHID '[' TOKENID ']'
		{ $$ = tg_get_terminal (tg, $1, $3); }
	;

rhs_symbol:
		NONTERMID		{ $$ = tg->idtab[$1->id]; }
	|	TERMALPHID              { $$ = tg->idtab[$1->id]; }
	|	terminal_symbol		{ $$ = tg->idtab[$1->id]; }
		;

rhs_symbols:
		rhs_symbol
		{ $$ = tg_vec_new (1); tg_vec_append (&$$, (void*) $1); }

	|	rhs_symbols rhs_symbol
		{ $$ = $1; tg_vec_append (&$$, (void*) $2); }
	;

rhs_symbols_or_end:
		rhs_symbols
		{ $$ = $1; }

	|	END
		{ $$ = tg_vec_new (0); }
	;

production_block:
		production_rule ';'
		{ }
	;

production_rule:
		lhs_symbol PRODUCES rhs_symbols_or_end optional_calc
		{ $$ = $1; tg_declare_rule (tg, $1, $3, $4); }

	|	production_rule '|' rhs_symbols_or_end optional_calc
		{ $$ = $1; tg_declare_rule (tg, $1, $3, $4); }
		;

optional_calc:
		{ $$ = tg_make_const_expr (tg, 1.); }
	|	'{' function_list '}'
		{ $$ = $2; }
	|	function_list
		{ $$ = $1; }
		;

func_declaration_statement:	FUNC func_declaration_list ';'
		;

func_declaration_list:
		func_declaration
	|	func_declaration_list ',' func_declaration
		;

func_declaration:
		IDENTIFIER '=' func_expr
		{ tg_declare_func (tg, $1, $3); }
	;

pset_declaration_statement:
		PSET pset_assignment_list ';'
	;

pset_assignment_list:
		pset_assignment
	|	pset_assignment_list ',' pset_assignment
		;

pset_assignment:
		IDENTIFIER '=' pexpr
		{ tg_declare_pset (tg, $3, $1); }
	;

const_parse:	
		NONTERMID child_desc
		{ $$ = tg_new_parse_node (tg, $1, (tg_parse_vec) $2) }

	|	terminal_symbol
		{ $$ = tg_new_parse_node (tg, $1, (tg_parse_vec) 0) }
	;

child_desc:
		'(' child_list ')'  { $$ = $2; }
	|	empty_child_list    { $$ = (tg_parse_vec) tg_vec_new (0); }
	;

child_list:
		const_parse
		{
		  $$ = (tg_parse_vec) tg_vec_new(1);
		  tg_vec_append ((tg_vec*) &$$, $1);
		}

	|	child_list const_parse
		{
		  $$ = $1;
		  tg_vec_append ((tg_vec*) &$$, $2);
		}
		;

empty_child_list:
		'(' ')'
	|	'('END ')'
	|	
	;

sequence_binding:
		TERMALPHID '=' STRVAL
	{ $$ = tg_new_binding (tg); tg_bind_sequence (tg, $$, $1, $3); }

	|	sequence_binding ',' TERMALPHID '=' STRVAL
	{ $$ = $1; tg_bind_sequence (tg, $$, $3, $5); }
		;

deriv_expr:
		DLOG psum_func_expr '/' DLOG dlog_denominator
		{ $$ = tg_make_binary_expr (tg, tg_Derivative, $2, $5); }
		;

psum_func_expr:
		FUNCID
		{
		  $$ = tg_make_func_expr (tg, $1);
		  tg_assert_is_psum (tg, $$);
		}

	|	'(' psum_func_expr ')' { $$ = $2; }
	;

func_expr:	function_list
	|	psum_expr
	|	deriv_expr
	;

dlog_denominator:
		fixed_prob_expr
		{ $$ = $1; }

	|	'(' dlog_denominator ')'
		{ $$ = $2; }
	;

psample_pexpr:	PSAMPLE pexpr
		{ $$ = tg_new_pexpr (tg, PSAMPLE, (void*) $2); }
	;

pmax_pexpr: 	PMAX pexpr
		{ $$ = tg_new_pexpr (tg, PMAX, (void*) $2); }
	;

id_pexpr:	PSETID
		{ $$ = tg_new_pexpr (tg, PSETID, $1->operand1); }
		;

const_pexpr:	'{' const_parse_list '}'
		{ $$ = tg_new_pexpr (tg, PSET, (void*) $2); }
	;

const_parse_list:
		const_parse
		{
		  $$ = (tg_parse_vec) tg_vec_new (0);
		  tg_vec_append ((tg_vec*) &$$, $1);
		}

	|	const_parse_list ',' const_parse
		{
		  $$ = $1;
		  tg_vec_append ((tg_vec*) &$$, $3);
		}
	;

dynamic_pexpr:
		PARSE '(' NONTERMID PRODUCES sequence_binding ')'
		{ $$ = tg_make_dynamic_pexpr (tg, $5, $3->sym); }
		;

psum_expr:
		PSUM '(' id_pexpr ')'
		{ $$ = tg_make_psum_expr (tg, $3); }
	;

pexpr:		const_pexpr
	|	dynamic_pexpr
	|	id_pexpr
	|	psample_pexpr
	|	pmax_pexpr
	|	'(' pexpr ')' { $$ = $2; }
	;

empty_statement:
		';'
		;

declaration:
		empty_statement
	|	alphabet_declaration_statement
	|	terminal_declaration_statement
	|	nonterminal_declaration_statement
	|	prob_declaration_statement
		;

declaration_list:
		declaration
	|	declaration_list declaration
	;

production_list:
		production_block
	|	production_list production_block
	;

action:
		empty_statement
	|	func_declaration_statement
	|	prob_assignment
	|	pset_declaration_statement
		;

action_list:
		action
	|	action_list action
	;

grammar_file:	declaration_list production_list;

script_file:	grammar_file action_list;

file:		grammar_file          { tg_compile (tg); }
        |	script_file           { tg_compile (tg); }
		;

%%

#include <stdio.h>

extern char* yytext;
extern int column;

void yyerror(s)
char *s;
{
  int i;
  fflush(stdout);
  for (i = 0; i <= column; ++i)
    printf ("%c", yytext[i-column]);
  printf("\n%*s\n%*s\n", column, "^", column, s);
}

int main( int argc, char **argv )
{
    ++argv, --argc;  /* skip over program name */
    if ( argc > 0 )
	{
	    yyin = fopen( argv[0], "r" );
	    if (yyin == 0)
		{
		    fprintf (stderr, "Couldn't find file: %s\n", argv[0]);
		    return -1;
		}
	}
    else
	yyin = stdin;

    tg = tg_new_grammar();
    yyparse();
    tg_delete_grammar (tg);

    return 0;
}
