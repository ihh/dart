/*Newick parser*/
	
typedef struct nwk_node_str{
	char* label; /*name in Stockholm alignment*/
	double b_len; /*distance to parent*/
	int index; /*index in parent's children*/
	int height; /*distance from root*/
	struct nwk_node_str* parent; /*pointer to parent*/
	tg_vector(struct nwk_node_str*) children; /* children array */
} nwk_node;

typedef tg_vector(nwk_node*) nwk_node_vec;

/*Functions to manipulate Newick trees*/
nwk_node* new_nwk_node();
void add_child(nwk_node* child, nwk_node* parent);
void free_nwk_tree(nwk_node* root);

/*Entrance points to Newick format parser*/
nwk_node* nwk_parse(char* file_name);
nwk_node* nwk_parse_main(FILE *nwk_file);

/*Functions for utility output*/
   /*Input: Stockholm database of pairwise alignments */
   /*Output: A single multiple alignment with inferred guide tree */
void stk_m_align(FILE* db_file, FILE* s_out_file); 
   /* Input: A parsed Newick node and associated Stockholm file*/  
   /* Output: Sequences of subtree rooted at Newick node */
void sub_tree(nwk_node* n_node, FILE* stk_file, FILE* out_file);   
   /*Input: Stockholm multiple alignment with guide tree */
   /*Output: Pairwise Stockholm database */
void guide_tree(FILE* stk_file, FILE* out_file);

/*Misc. functions used in utilities*/
char* tg_char_vec_conv(tg_char_vec* vec); /*converts a char vec to string */
int tree_height(nwk_node* n_node); /* calculates height of deepest leaf from node */
char* nwk_node_to_string(nwk_node* root); /*converts a tree rooted at node to Newick format*/
/*prints tree rooted at node with Newick_node info*/
/*tabs paremeter used for structure: 0 for initial call*/
void print_nwk(nwk_node *n_node, int tabs); 

nwk_node* seq_search(nwk_node_vec* n_array, char* label);

