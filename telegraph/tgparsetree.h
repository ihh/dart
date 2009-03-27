#ifndef TELEGRAPH_PARSE_TREE_INCLUDED
#define TELEGRAPH_PARSE_TREE_INCLUDED

/* Node in a parse tree */
typedef struct tg_parse_node_str {
  /* Label of this node (rule or terminal) */
    union tg_parse_union {
	tg_rule* rule;  /* internal nodes */
	tg_sym*  term;  /* leaf nodes */
    } label;

    /* Tree structure */
    tg_vector(struct tg_parse_node_str*) child;  /* if 0, then is leaf node */
    struct tg_parse_node_str* parent;  /* tree is doubly linked */

} tg_parse_node;

/* Array of child nodes for a single node in a parse tree */
typedef tg_vector(tg_parse_node*) tg_parse_vec;

/* Parse expressions.
   --------+---------------------+--------------------------
   op      | type of operand     | Interpretation
   --------+---------------------+--------------------------
   PSET    | (tg_parse_vec)      | constant parse expression
   PARSE   | (tg_dynamic_pexpr*) | dynamic parse expression
   PMAX    | (tg_pexpr*)         | pmax(...)
   PSAMPLE | (tg_pexpr*)         | psample(...)
   PSETID  | (tg_idtab_index)    | identifier
   --------+---------------------+--------------------------
*/
typedef struct tg_pexpr_str {
  int op;
  void* operand;
} tg_pexpr;

typedef tg_vector(tg_pexpr*) tg_pexpr_vec;

/* Method to recursively free parse tree nodes */
void tg_free_parse_tree (tg_parse_node* root);

#endif /* TELEGRAPH_PARSE_TREE_INCLUDED */
