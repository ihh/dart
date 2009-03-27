#ifndef TELEGRAPH_ALPHABET_INCLUDED
#define TELEGRAPH_ALPHABET_INCLUDED

/* Token, terminal and nonterminal symbols */
typedef struct tg_sym_str {
  tg_idtab_index id;    /* index of identifier */
  int            alph;  /* terminal/token alphabet, or -1 for nonterminals */
  tg_sym_index   sym;   /* symbol within alphabet, or nonterminal index */
} tg_sym;

typedef tg_vector(tg_sym*) tg_symbol_vec;  /* elements have type (tg_sym*) */

/* Token alphabets */
typedef struct tg_alph_str {
  tg_idtab_index id;    /* index of identifier */
  tg_alph_index  alph;  /* index of this token alphabet */
  tg_symbol_vec  token; /* the alphabet tokens */
} tg_alph;

typedef tg_vector(tg_alph*) tg_alph_vec;

/* Terminal alphabets */
typedef struct tg_term_str {
  tg_idtab_index id;      /* index of identifier */
  tg_term_index  term;    /* index of this terminal alphabet */
  tg_alph_index  tok;     /* index of token alphabet */
  tg_symbol_vec  termsym; /* the terminal symbols */
} tg_term;

typedef tg_vector(tg_term*) tg_term_vec;

#endif /* TELEGRAPH_ALPHABET_INCLUDED */
