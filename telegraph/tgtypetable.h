#ifndef TELEGRAPH_IDENTIFIER_TABLE_INCLUDED
#define TELEGRAPH_IDENTIFIER_TABLE_INCLUDED

/* The Telegraph identifier table.
   A type lookup table for all symbols in a Telegraph file.
*/

  /* ID table entry types, defined in "gram.tab.h", generated from "gram.y" */
  typedef int tg_type;

/* ID table entries */
  /* type       | value type
     -----------+------------
     TERMID     | (tg_sym*)
     NONTERMID  | (tg_sym*)
     TERMALPHID | (tg_alph*)
     PROBID     | (tg_prob*)
     RULEID     | (tg_rule*)
  */
  typedef struct tg_idtab_entry_str
  {
    const char* name;
    tg_type type;
    void* value;
  } tg_idtab_entry;

typedef tg_vector(tg_idtab_entry*) tg_identifier_vec;  /* for RHS of rules */
typedef tg_vector(tg_idtab_entry*) tg_type_table;  /* for grammar */

#endif /* TELEGRAPH_IDENTIFIER_TABLE_INCLUDED */
