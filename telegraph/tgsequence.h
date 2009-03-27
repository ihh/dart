#ifndef TELEGRAPH_SEQUENCE_INCLUDED
#define TELEGRAPH_SEQUENCE_INCLUDED

/* Sequence */
typedef tg_vector(tg_sym_index) tg_sequence;

/* Vector of sequences */
typedef tg_vector(tg_sequence) tg_sequence_vec;

/* Subsequence */
typedef struct tg_subseq_str {
    tg_sym_index* begin;
    unsigned int len;
} tg_subseq;

/* Vector of subsequences */
typedef tg_vector(tg_subseq*) tg_subseq_vec;

/* Dynamic parse expression */
typedef struct tg_dynamic_pexpr_str {
  tg_sym_index  start;  /* the root nonterminal, or "start state" */
  tg_subseq_vec seq;    /* the observed terminal sequences */
} tg_dynamic_pexpr;

/* Function to lex an individual input string */
tg_sequence tg_lex_string (tg_grammar* tg, const char* input_string,
			   tg_alph* alphabet);

/* Function to create a new sequence binding */
tg_sequence_vec tg_new_binding (tg_grammar* tg);

/* Function to update a sequence binding, checking for duplication */
void tg_bind_sequence (tg_grammar* tg, tg_sequence_vec binding,
		       tg_term* term, const char* s);

/* Function to check a sequence binding for completeness */
void tg_check_binding (tg_grammar* tg, tg_sequence_vec binding);

/* Function to create a subseq */
tg_subseq_vec tg_new_subseq_vec (tg_sequence_vec fullseq_vec,
				 tg_unsigned_vec subseq_offset,
				 tg_unsigned_vec subseq_len);

#endif /* TELEGRAPH_SEQUENCE_INCLUDED */
