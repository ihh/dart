#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "tgcompiler.h"

#define max_line_size 1000
#define seq_name_size 5

/*********** TOP CODE FOR PAIRWISE ALIGNMENT TO NEWICK TREE AND MULT. ALIGNMENT ***********/

/*Given as input a Stockholm database of pairwise alignments in breadth first order 
from the root, each corresponding to a single branch of an assumed "guide tree", print 
on the output a single Stockholm multiple alignment (with the guide tree)
*/

/*checks to see if string is composed of only whitespace*/
int is_blank_line(char* line) {
	int c;
	while ((c = *line++)) {
		if (!isspace(c))
			return 0;
		/* if (c == '\n') 		
			 break; */
	}
	return 1;
}

/*returns sequence name from a Stockholm alignment line*/
char* seq_extract(char* line) {
	char* temp = line;
	char* to, *ret;
	int j = 0, i = 0;
	while (!isspace(*line++)) {i++;}
	to = (char*) malloc((i+1)*sizeof(char));
	ret = to;
	while (j < i) {
		*to++ = *temp++;
		j++;
	}
	*to = '\0';
	return ret;
}

/*recompiles a pairwise database of Stockholm alignments into a single multiple align. */
void stk_m_align(FILE* db_file, FILE* s_out_file) {
	fpos_t *holder = malloc(sizeof(fpos_t));
	char *line = malloc(max_line_size*sizeof(char));
	char *label, *label_end;
	int flag = 0, num_of_lines = 0, i, j;
	nwk_node_vec n_array = (nwk_node_vec) tg_vec_new(0);
	nwk_node* parent, *child, *root;			
	/* First pass to extract tree structure */	
	while (fgets(line, max_line_size, db_file)) {
		if (!strncmp("#", line, 1) || is_blank_line(line)){
			continue;
		} else if (!strncmp("//", line, 2)) {
			printf("Error: Empty alignment in database found\n");
			fclose(db_file);
	      fclose(s_out_file);
			return;
		} else { /*New sequence found*/
			label = seq_extract(line);
			if (!(parent = seq_search(&n_array, label))) { /*Seq not in processed nodes*/
				parent = new_nwk_node();
				parent->label = label;
				if (!flag) {
					root = parent;
					flag = 1;
				}
				tg_vec_append((tg_vec*) &n_array, (void*) parent);
			}
			while (fgets(line, max_line_size, db_file)) {
			   label_end = line + strlen(label);
				if (!strncmp("#", line, 1) || 
				   (!strncmp(label, line, strlen(label)) && isspace(*label_end)) ||
				   is_blank_line(line)) {
					continue;
				} else if (!strncmp("//", line, 2)) {
					printf("Error: Unpaired sequence in database found: %s\n", label);
				   fclose(db_file);
	            fclose(s_out_file);
					return;
				} else {
					label = seq_extract(line);				
					if (!(child = seq_search(&n_array, label))) { /*Child seq not in processed nodes*/
						child = new_nwk_node();
						child->label = label;
						tg_vec_append((tg_vec*) &n_array, (void*) child);
					}
					add_child(child, parent);
					while (fgets(line, max_line_size, db_file) && strncmp("//", line, 2)) {}
					break;
				}
			}
		}
		/*End of processing a pairwise alignment*/
	}	
	fseek(db_file, 0, SEEK_SET);
	/* Second pass to write Stockholm file */
	while (fgets(line, max_line_size, db_file)) {	
		if (!strncmp("#", line, 1) || is_blank_line(line)){
		   if (strncmp("#=GF NH", line, 7)) {
  			   fputs(line, s_out_file);
			}
			if (!strncmp ("# STOCKHOLM", line, 11)) {
				fputs("\n", s_out_file);				
				fputs("#=GF NH ", s_out_file);
				fputs(nwk_node_to_string(root), s_out_file);
				fputs("\n", s_out_file);	
			}
		} else if (!strncmp("//", line, 2)) {
			fputs(line, s_out_file);
	      break;
		} else if (!strncmp(root->label, line, strlen(root->label))) {
			fputs(line, s_out_file); /*put the root's seq line first*/
			fgetpos(db_file, holder);
			fseek(db_file, 0, SEEK_SET);				
			for (i = 1; i < tg_size(n_array); i++) { /*root's node at 0 in array*/
				label = n_array[i]->label;
				j = 0;
				while (j <= num_of_lines && fgets(line, max_line_size, db_file) ) {
					if (!strncmp(label, line, strlen(label))) 
						j++;
				}
				fputs(line, s_out_file);	
				fseek(db_file, 0, SEEK_SET);
			}
			fsetpos(db_file, holder);
			num_of_lines++;
		}
	}
	/*Memory management*/
	free(holder);
	free(line);
	tg_vec_delete((tg_vec*) &n_array);
	free_nwk_tree(root);
	fclose(db_file);
	fclose(s_out_file);
}

/***************** BOTTOM CODE FOR NEWICK TREE TO ALIGNMENT ******************/
 
/*	The main for this program takes input in the form:
	# stockholm [path] [seq1] [seq2...seqN]
	where path is the path to a Stockholm formatted file,
	and seq1...seqN are sequence names in the that file.
	It writes the new Stockholm file to the executable's folder.
*/

/* checks for seq in newick tree of stockholm file */
int quick_seq_scan(char* seq, FILE* stk_file) {
	char* line = malloc(max_line_size*sizeof(char));
	if (seq) {
		while (fgets(line, max_line_size, stk_file)) {
			if (!strncmp(seq, line, strlen(seq))) {
				fseek(stk_file, 0, SEEK_SET);			
				free(line);
				return 1;
			}
		}
		fseek(stk_file, 0, SEEK_SET);
		printf("Error Writing to file -- Sequence: %s not named in Stockholm Multiple Alignment\n",seq);
	} else {
		printf("Error Writing to file -- Sequence not named in Newick Tree\n");
	}
	free(line);
	return 0;
}

int put_seq(FILE* stk_file, FILE* out_file, char** seqs, int num_seqs, int flag)  {
	int i;
	char* line = malloc(max_line_size*sizeof(char));
	char* label_end;
	/*Error checking */
	for (i = 0; i < num_seqs; i++) {
		if (!quick_seq_scan(seqs[i], stk_file)) {
			return 0;
		}
	}
	
	while (fgets(line, max_line_size, stk_file)) {
		if (!flag) { /*flag set for only one printing of comments */
			if (!strncmp("#", line, 1)) { 
			   fputs(line, out_file);
			   continue;
			}
		}
		if (!strncmp("//", line, 2) || is_blank_line(line)) {
			fputs(line, out_file);
			continue;
			}
		if (!strncmp("#", line, 1)) {continue;}
		for (i = 0; i < num_seqs; i++) {
		   label_end = line + strlen(seqs[i]);
			if (!strncmp(seqs[i], line, strlen(seqs[i])) && isspace(*label_end)) {
				fputs(line, out_file);
				break;
			}
		}
	}
	free(line);
	return 1;
}

/*
	nwk_file -> Newick tree representation file
	stk_file -> Stockholm alignment file
*/

int guide_tree_help(nwk_node* n_node, FILE* stk_file, FILE* out_file, int flag) {
	char** seq;
	int i;
	seq = malloc(2*sizeof(char *)); /*two sequences per run*/
	if (n_node) {
		for (i = 0; i < tg_size(n_node->children); i++) {
			seq[0] = n_node->label;
			seq[1] = n_node->children[i]->label;
			if (!put_seq(stk_file, out_file, seq, 2, flag))
			   return 0;
			fseek(stk_file, 0, SEEK_SET);
			flag = 1;
		}
		for (i = 0; i < tg_size(n_node->children); i++) { /*seperate from above loop to preserve breadth first order*/
			guide_tree_help(n_node->children[i], stk_file, out_file, flag);
		}
	}
	return 1;
}

/*prints a pairwise database to out_file given a nwk_tree and stk mult. alignment
  returns 0 on error*/
void guide_tree(FILE* stk_file, FILE* out_file) {
	nwk_node* n_node;
	if ((n_node = nwk_parse_main(stk_file))){ /*set n_node to root*/
		if (!guide_tree_help(n_node, stk_file, out_file, 0)) {
			printf("Error: Guide tree not written\n");
		}
		fclose(stk_file);
		fclose(out_file);
	}	
}

/* Helper for sub_tree
	Counts number of nodes rooted at n_node
*/
int count_t_nodes(nwk_node* n_node) {
	int i, count = 0;
	if (n_node) {
		count++;
		for (i = 0; i < tg_size(n_node->children); i++) {
			count += count_t_nodes(n_node->children[i]);
		}
	}
	return count;
}

/* Helper for sub_tree
	Fills seq with labels of nodes rooted at n_node
*/
int fill_seq(nwk_node* n_node, char** seq, int start){
	int i;
	if (n_node) {
		if (!start) { /*initialization*/
			seq[start] = n_node->label;
			start++;
		} 
		for (i = 0; i < tg_size(n_node->children); i++) {
			seq[start] = n_node->children[i]->label;
			start++;
		}
		for (i = 0; i < tg_size(n_node->children); i++) {	
			start = fill_seq(n_node->children[i], seq, start);
		}
	}
	return start;	
}

/* prints subtree */
void sub_tree(nwk_node* n_node, FILE* stk_file, FILE* out_file) {
	char** seqs;
	int size;
	if (n_node) {
		size = count_t_nodes(n_node); /*counts nodes rooted at n_node*/
		seqs = malloc(size*sizeof(char*));
		fill_seq(n_node, seqs, 0); /*fills seq with correct labels in bf order*/
		put_seq(stk_file, out_file, seqs, size, 0); /*prints node's children seq to out_file*/
	}
	fclose(stk_file);
	fclose(out_file);
}
