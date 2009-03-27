/*Newick parser*/
	
typedef struct nwk_seq_str{
	int size; /* total size of sequence */
	int type; /* 0 - DNA; 1 - RNA; 2 - AMINO */
	tg_string_vec seqs; /* array of sequences */
	tg_char_vec RNA_ss; /*to hold RNA secondary structure*/ 

} nwk_seq;

typedef struct nwk_node_str{
	char* label; /*name in Stockholm alignment*/
	double b_len; /*distance to parent*/
	int index; /*index in parent's children*/
	int height; /*distance from root*/
	struct nwk_node_str* parent; /*pointer to parent*/
	int level_index; /*index in level*/ 
	struct nwk_node_str* next; /* pointer to next sibling on this node's level */ 
    nwk_seq* seq; /*pointer to structure for holding sequences*/
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
void stk_m_align(FILE* db_file, FILE* s_out_file);
void sub_tree(nwk_node* node, FILE* stk_file, FILE* out_file);
void guide_tree(FILE* nwk_file, FILE* stk_file, FILE* out_file);

/*Misc. functions used in utilities*/
char* tg_char_vec_conv(tg_char_vec* vec);
int tree_height(nwk_node* n_node);
char* nwk_node_to_string(nwk_node* root); /* writes tree rooted at node in Newick format*/
void print_nwk(nwk_node *n_node, int tabs); 
/*tabs paremeter used for structure: 0 for initial call*/
nwk_node* seq_search(nwk_node_vec* n_array, char* label);

