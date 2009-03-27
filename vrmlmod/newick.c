#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "datastruct.h"
#include "util.h"

#define max_line_size 512
#define pairs_per_run 20

#define DNA 0
#define RNA 1
#define AMINO 2

/*Newick parser*/
/*NOTE: Features of Newick Format not supported :
		Brackets (comments)
		Underscore -> space replacement (because Stockholm format alignments use underscore, not space)
   	Quoted Labels (quotes still show up)
   	Empty nodes (no label or branch length) -> not guaranteed to work	
*/
nwk_node* nwk_fr_file(FILE* nwk_file);

/*Node initializer*/
nwk_node* new_nwk_node() {
	nwk_node* new_node = malloc (sizeof(nwk_node));
	new_node->label = NULL;
	new_node->b_len = 0;
	new_node->parent = NULL;
	new_node->index = 0;
	new_node->height = 0;
	new_node->children = (nwk_node_vec) tg_vec_new(0);
	new_node->seq = malloc(sizeof(nwk_seq));
	new_node->seq->size = 0;
	new_node->seq->type = -1; /* none assigned yet*/
	new_node->seq->seqs = (tg_string_vec) tg_vec_new(0);
	new_node->seq->RNA_ss = (tg_char_vec) tg_vec_new(0);
	return new_node;
}

/*adds a child to a parent node*/
void add_child(nwk_node* child, nwk_node* parent) {
	child->parent = parent;
	child->index = tg_size(parent->children);
	tg_vec_append((tg_vec*) &parent->children, (void*) child);
}

/*Tree destructor, must initially be called from the root*/
void free_nwk_tree(nwk_node* root) {
	int i;
	if (root) {
		free(root->label);
		if (root->children) {
			for (i = 0; i < tg_size(root->children); i++) {
				free_nwk_tree(root->children[i]);
			}
		}
		tg_vec_delete((tg_vec*) &root->children);
	}
}			
			
/*returns nwk_node* with label "label" if it is in n_array, and null otherwise*/
nwk_node* seq_search(nwk_node_vec* n_array, char* label) {
	int i;	
	for (i = 0; i < tg_size(*n_array); i++) {
		if (!strcmp(label, (*n_array)[i]->label)) {
			return (*n_array)[i];
		}
	}
	return NULL;
}
	
/*returns char* of converted tg_char_vec*/
char* tg_char_vec_conv(tg_char_vec* vec) {
	char* ret = malloc(sizeof(char)*tg_size(*vec));  
	int i;
	for (i = 0; i < tg_size(*vec); i++) {
		 ret[i] = (*vec)[i];
	}
	return ret;
}

int tree_height(nwk_node* n_node) {
   int i, max = 0;
   if (n_node) {
      if (tg_size(n_node->children)) {
         for (i = 0; i < (tg_size(n_node->children)); i++) {
		    if (tree_height(n_node->children[i]) > max) {
			   max = tree_height(n_node->children[i]);
			}
		 }
	  } else {
	     return n_node->height;
	  }
	  return max;
   }
   return 0;
}	

void nwk_node_to_string_help(nwk_node* node, tg_char_vec* str) {
	int i;
	if (node) {
		if (tg_size(node->children)) {
			tg_vec_append((tg_vec*) str, (void*) '(');
			for (i = 0; i < tg_size(node->children); i++) {
				nwk_node_to_string_help(node->children[i], str);
				if (i+1 < tg_size(node->children)) {	
					tg_vec_append((tg_vec*) str, (void*) ',');
					tg_vec_append((tg_vec*) str, (void*) ' ');
				} else {
					tg_vec_append((tg_vec*) str, (void*) ')');
				}
			}
		}
		for (i = 0; i < strlen(node->label); i++) {
			tg_vec_append((tg_vec*) str, (void*) ((int) (node->label)[i]));
		}
	}
}

/*returns string representing the Newick tree rooted at "root"*/
char* nwk_node_to_string(nwk_node* root) {
	tg_char_vec str = (tg_char_vec) tg_vec_new(0);
	char* ret;
	nwk_node_to_string_help(root, &str);
	tg_vec_append((tg_vec*) &str, (void*) ';');
	tg_vec_append((tg_vec*) &str, (void*) '\0');
	ret = tg_char_vec_conv(&str);
	tg_vec_delete((tg_vec*) &str);
	return ret;
}		

/*returns 0 if the newick nwk_file has all matching parentheses*/
int parens_check(FILE* nwk_file) {
	int count = 0, next;
	while ((next = fgetc(nwk_file)) !=EOF && next != ';'){
		if (next == '(') {
			count++;
		} else if (next == ')') {
			count--;
		}
	}
	fseek(nwk_file, 0, SEEK_SET); /*reset FILE */
	return count;
}

/*returns 0 if file has a semicolon in it */
int end_check(FILE* nwk_file) {
	int flag = 1, next;
	while ((next = fgetc(nwk_file)) !=EOF){
		if (next == ';')
			flag = 0;
	}
	fseek(nwk_file, 0, SEEK_SET); /*reset FILE */
	return flag;
}

/* returns 0 if empty node found.
	this parser doesn't support empty nodes, this checks for them */
int empty_leaf_check(FILE* nwk_file) {
	int next;
	while ((next = fgetc(nwk_file)) != EOF && next!= ';'){
		if (next == '(' || next == ',') {
			while (isspace(next = fgetc(nwk_file))) {}
			if (next == ')' || next == ';' || next == ',') { /*terminators*/
				return 1;
			} else {
				return empty_leaf_check(nwk_file);
			}
		}
	}
	fseek(nwk_file, 0, SEEK_SET);
	return 0;		
}

/* returns 0 if branch lengths are in correct format */
int number_check(FILE* nwk_file) {
	int next, flag = 0;
	while ((next = fgetc(nwk_file)) != EOF && next!= ';') {
		if (next == ':') { /*entering a branch length seq */
			while (isspace(next = fgetc(nwk_file))) {}
			ungetc(next, nwk_file);
			while ((next = fgetc(nwk_file)) != EOF && next != ';') {
				if (isdigit(next)) {
					continue;
				} else if ((next == '.') && !flag) {
					flag = 1;
					continue;
				} else if (next != ')' && !isspace(next) && next != ',') {
					return 1;
				} else {
					flag = 0;
					break;
				}
			}
		}
	}
	return 0;
}

/* returns 0 if next is a character for a label or branch length */
int info_char(int next) {
	return (!isspace(next) && next != '(' && next != ')' && next != ',' && next != ';');
}

/* checks for correct use of commas in seperating nodes */
int illegal_space_check(FILE* nwk_file) {
	int next;
	while ((next = fgetc(nwk_file)) != EOF && next != ';') {
		if (info_char(next)) {
			while (info_char(next = fgetc(nwk_file))) {}
			ungetc(next, nwk_file);
			while (isspace(next = fgetc(nwk_file))) {}
			if (next != ')' && next != ',' && next != ';') {
				return 1;
			}
		}
	}
	fseek(nwk_file, 0, SEEK_SET);
	return 0;
}

/* returns a char* at the beginning of the elements of the next matching run [name] */
char* get_next_run(char* name, FILE* stk_file) {
   char* line = malloc(max_line_size*sizeof(char));
   while (fgets(line, max_line_size, stk_file)) {
      if (!strncmp(line, "//", 2)) {
	     return NULL;
      } else if(!strncmp(line, name, strlen(name))) {
	     if (isspace(line[strlen(name)])) {
  	        line += strlen(name)-1;
		    while (isspace(*(++line))) {} 
		    return line;
	     }
	  }
   }
   return NULL;
}

/*Can only be one tree represented in stk_file to work !!!!!!!DEAL WITH SEMICOLON*/	
int get_nwk(FILE* stk_file, FILE* nwk_file) {
    char* line = malloc(max_line_size*sizeof(char));
	char* copy;
	int flag1 = 1;
	while (fgets(line, max_line_size, stk_file)) {
	   if (!strncmp("#=GF NH", line, 7)) {
	      flag1 = 0;
		  copy = line;
		  copy += 7; /*skips past mark-up header*/
		  while (*copy != '\n') {
			fputc(*copy, nwk_file);
		    if (*copy == ';') {
		       break;
			}
			copy++;
		  }
		  if (*copy == ';') {
		     break;
	      }
		  fputc(' ', nwk_file);
	   }
    }
	fseek(stk_file, 0, SEEK_SET);
	fseek(nwk_file, 0, SEEK_SET);			
    return flag1;
}

void fill_nwk_seq(nwk_node* n_node,  FILE* stk_file) {
   char* line;
   char* run = malloc((pairs_per_run+1)*sizeof(char));
   char* name; 
   nwk_seq* n_seq;
   int RNA_flag = 0, AMINO_flag = 0; /*for detecting type*/
   int i = 0, size = 0;
   if (n_node) {
      name = n_node->label;
	  n_seq = n_node->seq;
	  while ((line = get_next_run(name, stk_file))) {
	     while (i < pairs_per_run) {
		    if (!isspace(*line)) {
			   switch (*line) {
			      case 'A': case 'C' : case 'G': case 'T': case '-': case '.':
				     break;
				  case 'U': 
				     RNA_flag = 1;
					 break;
				  default:
				     AMINO_flag = 1;
			   }
			   /*Is there any easier way to determine sequence type*/     
			   run[i] = *line;
			   i++;
			   line++;
			   if (i == pairs_per_run) {
			      run[i] = '\0';
			      i = 0;
		          tg_vec_append((tg_vec*) &(n_seq->seqs), (void*) run);
				  run = malloc((pairs_per_run+1)*sizeof(char));
				  size += pairs_per_run;
			   }
		    } else break;
		 }
	  }
	  run[i] = '\0';
	  tg_vec_append((tg_vec*) &(n_seq->seqs), (void*) run);
      fseek(stk_file, 0, SEEK_SET);
	  /*Setting size*/
	  size += i;
	  n_seq->size = size;
	  /*Setting type*/
	  if (AMINO_flag) {
	     n_seq->type = AMINO;
	  } else if (RNA_flag) {
	     n_seq->type = RNA;
	  } else {
	     n_seq->type = DNA;
	  }
	  /*Iterate for children*/
	  for (i = 0; i < tg_size(n_node->children); i++) {
	     fill_nwk_seq(n_node->children[i], stk_file);
	  }
   }
   if (n_node->parent) { 
      if (n_node->parent->seq->type != n_seq->type) {
	     printf("Inconsistent basepair type in seq (%d != %d): %s -- Erroneuous output\n", n_node->parent->seq->type, n_seq->type, n_node->label);
	  }
   }
}		
/* fills in the heights of the nodes in the nwk tree*/
void fill_nwk_height(nwk_node* n_node, int height) {
   int i;
   if (n_node) {
      n_node->height = height;   
      for (i = 0;  i < tg_size(n_node->children); i++) {
         fill_nwk_height(n_node->children[i], height+1);
	  }
   }
}

void fill_RNA_help(nwk_node* n_node) {
   int i;
   if (n_node) {
      for (i = 0; i < tg_size(n_node->children); i++) {
	     n_node->children[i]->seq->RNA_ss = n_node->seq->RNA_ss;
		 fill_RNA_help(n_node->children[i]);
	  }
   }
}

void fill_RNA_ss_cons(nwk_node* n_node, FILE* stk_file) {
   char* line = malloc(max_line_size*sizeof(char));
   char* move;
   if (n_node && n_node->seq->type == RNA) {
      while (fgets(line, max_line_size, stk_file)) {
	     move = line;
		 if (!strncmp("#=GC SS_cons", move, 12)) {
		    move += 12;
			while (isspace(*move)) {move++;}
			while (!isspace(*move)) {
			   switch (*move) {
			   case '.':
			   case '>':
			   case '<':
			      tg_vec_append((tg_vec*)&n_node->seq->RNA_ss, (void*)((int)(*move)));
				  move++;
			      break;
			   default:
			      printf("Illegal character in SS_cons format - Incomplete RNA_ss field\n");
				  return;
			   }
			}
		 }
	  }
	  fseek(stk_file, 0, SEEK_SET);	
	  tg_vec_append((tg_vec*)&n_node->seq->RNA_ss, (void*)((int)'\0'));	  
	  fill_RNA_help(n_node);
   }
   free(line);
}
					  
							

/*returns string representing the label of the 
  node parsed from current position in nwk_file*/
char* nwk_label(FILE* nwk_file) {
	int flag = 0, next;
	tg_char_vec label = (tg_char_vec) tg_vec_new(0);
	char* ret_label;
	
	while (isspace(next = fgetc(nwk_file))); /*skip white space*/
	ungetc(next, nwk_file);	
	while ((next = fgetc(nwk_file)) != EOF && next != ';' 
			&& next != ':' && next != ',' && next != ')' && !isspace(next)) {
		tg_vec_append((tg_vec*) &label, (void*) next);
		flag = 1;
	}
	
	if (next == ';' || next == ':' || next == ')' || next == ',')
		ungetc(next, nwk_file);
	if (flag) {
		tg_vec_append((tg_vec*) &label, (void*) '\0');
		ret_label = tg_char_vec_conv(&label);
		tg_vec_delete((tg_vec*) &label);
		return ret_label;
	} else {
		tg_vec_delete((tg_vec*) &label);
		return NULL; /*no label found*/
	}
}

/*returns string representing the branch length of the 
  node parsed from current position in nwk_file*/
double nwk_b_len(FILE* nwk_file) {
	int i, next;
	tg_char_vec b_str = (tg_char_vec) tg_vec_new(0);
	while (isspace(next = fgetc(nwk_file))) {}
	if (next != ':') { /*b_len has to start with a colon*/
		if (next != ',')
			ungetc(next, nwk_file);
		return 0;	
	}
	for (i = 0; (next = fgetc(nwk_file)) != EOF && (isdigit(next) || next == '.') ; i++) {
		tg_vec_append((tg_vec*) &b_str, (void*) next);
	}
	if (next == ':') {
		printf("Error parsing file: Invalid use of ':'\n");
		return 0;
	} else if (isspace(next)) { /*case of white space after b_length, before comma */
		while (isspace(next = fgetc(nwk_file))) {}
		ungetc(next, nwk_file);
	} else if (next == ';' || next == ')') {
		ungetc(next, nwk_file); 
	}
	tg_vec_append((tg_vec*) &b_str, (void*) '\0');
	return atof(tg_char_vec_conv(&b_str));
}

/*returns nwk_node (subtree) parsed from given position in nwk_file*/
nwk_node* nwk_branch_parse(FILE* nwk_file) {
	int next;
	nwk_node* ret_node = new_nwk_node();
	while ((next = fgetc(nwk_file)) != EOF && next!= ';') {
		if (isspace(next)) {
			continue;
		} else if (next == ')') {
			return NULL;
		} else if (next == '(') { /*next child is a subtree*/
	 		ungetc(next, nwk_file);
			return nwk_fr_file(nwk_file);
		} else { /*next child is a leaf */
			if (next == ':') {
				ungetc(next, nwk_file);
				ret_node->b_len = nwk_b_len(nwk_file);
				if ((next = fgetc(nwk_file)) != ',') {
					ungetc(next, nwk_file);
				}
				return ret_node;				
			} else {
				ungetc(next, nwk_file);
				ret_node->label = nwk_label(nwk_file);
				ret_node->b_len = nwk_b_len(nwk_file);
				if ((next = fgetc(nwk_file)) != ',') {
					ungetc(next, nwk_file);
				}
				return ret_node;
			}
		}
	}
	return NULL;
}

/* Main parser */
nwk_node* nwk_fr_file(FILE* nwk_file) {
	nwk_node* ret_node, *child;
	int next;
	
	ret_node = new_nwk_node(); /*space for a new node*/
	while ((next = fgetc(nwk_file)) != EOF && (next != ';') && (next != ',')) {
		if (isspace(next)) {
			continue;
		} else if (next == '(') { /*Start of a descendant list*/
			while ((child = nwk_branch_parse(nwk_file))) { /*returns NULL when hit right parens*/
			 	add_child(child, ret_node);
			}
		} else { /* at this point, finished dealing with right parens of encountered desc. list */
		ungetc(next, nwk_file);
		ret_node->index = 0;
		ret_node->label = nwk_label(nwk_file);
		ret_node->b_len = nwk_b_len(nwk_file);	
		return ret_node;
		}
	}
	return ret_node;
}
	
/* primary function for parsing a nwk_file
	runs checks on file to prevent crashing*/
nwk_node* nwk_parse_main(FILE* stk_file) {	
    FILE* nwk_file = fopen("NWK_TEMP0000011111.nwk", "w"); /*Open for writing*/
	nwk_node* root;
	if (get_nwk(stk_file, nwk_file)) {
	    printf("Error parsing Newick tree: Illegal or nonexistent representation in file\n");
        fclose(nwk_file);
        remove("NWK_TEMP0000011111.nwk");
		return NULL;
	} else {
	   fclose(nwk_file);
	   nwk_file = fopen("NWK_TEMP0000011111.nwk", "r"); /*Open for reading*/	   
	   if (end_check(nwk_file)) {
		printf("Error parsing Newick tree: Illegal tree terminator (no ';')\n");
		return NULL; }
		if (parens_check(nwk_file)) {
			printf("Error parsing Newick tree: Mismatched parentheses\n");
			return NULL; }
		if (empty_leaf_check(nwk_file)) {
			printf("Error parsing Newick tree: Empty node found\n");		
			return NULL; }
		if (number_check(nwk_file)) {
			printf("Error parsing Newick tree: Malformed branch length\n");
			return NULL; }
		if (illegal_space_check(nwk_file)) {
			printf("Error parsing Newick tree: Nodes not seperated by ','\n");
			return NULL; }
		root = nwk_fr_file(nwk_file);
		fill_nwk_height(root, 0);
		fill_nwk_seq(root, stk_file);
		fill_RNA_ss_cons(root, stk_file);
		fclose(nwk_file);
        remove("NWK_TEMP0000011111.nwk");
		return root;
	}
}

/* Does the same job as nwk_parse_main
	Provides entry point without having to open a FILE* */
nwk_node* nwk_parse(char* file_path){
	FILE* stk_file;
	nwk_node* n_node;
	if ((stk_file = fopen(file_path, "r"))) {
		n_node = nwk_parse_main(stk_file);
		fclose(stk_file);
		return n_node;
	} else {
		printf("Error opening Stockholm format file at: %s\n", file_path);
		return NULL;
	}
}	

/* prints Newick tree to standard out
	tabs paremeter used for structure
	set tabs to 0 in initial call*/
void print_nwk(nwk_node* n_node, int tabs) {
	int i, j;
	char* label; /*name in Stockholm alignment*/
	nwk_seq* n_seq = n_node->seq;
	nwk_node* parent; /*pointer to parent*/
/*	nwk_node* level_sibs = n_node;*/
	
	if (n_node) {
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("*****NODE*****\n");
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("Node Name: %s\n", (label = n_node->label) ? label : "");
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("Child %d of %s\n", n_node->index, (parent = n_node->parent) ? ((label = parent->label) ? label : "") : "no parent");
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("Height: %d\n", n_node->height);
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("Branch Length: %f\n", n_node->b_len);
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("%s[%d]: ", ((n_seq->type == 2) ? "AMINO" : ((n_seq->type == 1) ? "RNA" : "DNA")) , n_seq->size);
		if (!tg_size(n_seq->seqs))
		   printf("none");
		for (i = 0; i < tg_size(n_seq->seqs); i++) {
		   for (j = 0; j < strlen(n_seq->seqs[i]); j++) {
		      printf("%c", n_seq->seqs[i][j]);
		   }
		   printf("   ");
		}
		printf("\n");
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("RNA SS_Cons: %s\n", ((n_seq->type == RNA) ? tg_char_vec_conv(&n_seq->RNA_ss) : "n/a - Not RNA")); 
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
        printf("Number of Children: %d\n", tg_size(n_node->children));
		for (i = 0; i < tabs; i++) {
			printf("\t");
		}
		printf("Children: ");
		if (!tg_size(n_node->children))
			printf("no children");
		for (i = 0; i < tg_size(n_node->children); i++) {
			printf("%s  ", (label = ((nwk_node*) n_node->children[i])->label) ? label: "|no name|");
		}
		printf("\n");
		for (i = 0; i < tg_size(n_node->children); i++) {
			print_nwk(n_node->children[i], tabs+1);
		} 
	}
}



