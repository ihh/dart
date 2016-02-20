#include <stdlib.h>
#include <libguile.h>

#include "guile/guile-keywords.h"
#include "guile/newick-type.h"
#include "tree/tree_alignment.h"

// Smob tag
scm_t_bits newick_tag;

SCM make_newick_smob (const PHYLIP_tree& tree)
{
  SCM smob;
  PHYLIP_tree *newick_smob = new PHYLIP_tree (tree);
  SCM_NEWSMOB (smob, newick_tag, newick_smob);
  return smob;
}

static SCM newick_from_file (SCM s_filename)
{
  SCM smob = SCM_BOOL_F;

  try {
    // read tree from file
    PHYLIP_tree *tree = new PHYLIP_tree();
    char* filename = scm_to_locale_string (s_filename);
    ifstream in(filename);
    tree->read(in);
    free(filename);

    // Create the smob.
    SCM_NEWSMOB (smob, newick_tag, tree);

  } catch (Dart_exception& e) {
    CLOGERR << e.what();
  }

  // Return
  return smob;
}

static SCM newick_from_string (SCM s_string)
{
  SCM smob = SCM_BOOL_F;

  try {
    // read tree from string
    PHYLIP_tree *tree = new PHYLIP_tree();
    char* s = scm_to_locale_string (s_string);
    stringstream ss (s, stringstream::in);
    tree->read(ss);
    free(s);

    // Create the smob.
    SCM_NEWSMOB (smob, newick_tag, tree);

  } catch (Dart_exception& e) {
    CLOGERR << e.what();
  }

  // Return
  return smob;
}

static SCM newick_to_file (SCM tree_smob, SCM s_filename)
{
  PHYLIP_tree *tree = newick_cast_from_scm (tree_smob);

  // write alignment to file
  char* filename = scm_to_locale_string (s_filename);
  ofstream out(filename);
  tree->write(out,0);
  free(filename);

  // Return
  return SCM_UNSPECIFIED;
}

static SCM newick_node_list (SCM tree_smob, vector<Phylogeny::Node> (PHYLIP_tree::*getNodeMethod)())
{
  PHYLIP_tree& tree = *newick_cast_from_scm (tree_smob);
  SCM tree_node_list = SCM_EOL;
  vector<Phylogeny::Node> nodes = (tree.*getNodeMethod)();
  for_const_contents (vector<Phylogeny::Node>, nodes, iter) {
      Phylogeny::Node n = *iter;
      const char* n_name = tree.node_specifier(n).c_str();
      SCM single_element_list = scm_list_1(scm_from_locale_string(n_name));
      tree_node_list = scm_append(scm_list_2(tree_node_list,single_element_list));
  }

  scm_remember_upto_here_1 (tree_smob);

  // Return
  return tree_node_list;
}

static SCM newick_leaf_list (SCM tree_smob)
{
  return newick_node_list (tree_smob, &PHYLIP_tree::leaf_vector);
}

static SCM newick_ancestor_list (SCM tree_smob)
{
  return newick_node_list (tree_smob, &PHYLIP_tree::ancestor_vector);
}

static SCM newick_branch_list (SCM tree_smob)
{
  PHYLIP_tree& tree = *newick_cast_from_scm (tree_smob);
  SCM tree_branch_list = SCM_EOL;
  for_rooted_branches_pre (tree, b) {
    Phylogeny::Node parent = (*b).first, child = (*b).second;
    double length = (*b).length;
    const char* p_name = tree.node_specifier(parent).c_str();
    const char* c_name = tree.node_specifier(child).c_str();

    SCM branch_tuple = scm_list_3 (scm_from_locale_string(p_name),
				   scm_from_locale_string(c_name),
				   length < 0 || parent == tree.root ? SCM_BOOL_F : scm_from_double(length));
    SCM single_element_list = scm_list_1(branch_tuple);
    tree_branch_list = scm_append(scm_list_2(tree_branch_list,single_element_list));
  }

  scm_remember_upto_here_1 (tree_smob);

  // Return
  return tree_branch_list;
}


static SCM newick_unpack (SCM tree_smob)
{
  PHYLIP_tree& tree = *newick_cast_from_scm (tree_smob);
  vector<SCM> node_scm (tree.nodes(), SCM_EOL);
  for_rooted_nodes_post (tree, b) {
    Phylogeny::Node parent = (*b).first, node = (*b).second;
    double length = (*b).length;
    const char* n_name = tree.node_specifier(node).c_str();
    for_children (tree, node, parent, c) {
      node_scm[node] = scm_append (scm_list_2 (node_scm[node],
					       scm_list_1(node_scm[*c])));
    }
    node_scm[node] = scm_append (scm_list_2 (node_scm[node],
					     scm_list_2 (scm_from_locale_string(n_name),
							 parent < 0 || length < 0 ? SCM_BOOL_F : scm_from_double(length))));
  }

  scm_remember_upto_here_1 (tree_smob);

  // Return
  return node_scm[tree.root];
}

static size_t free_newick (SCM tree_smob)
{
  class PHYLIP_tree *tree = (class PHYLIP_tree *) SCM_SMOB_DATA (tree_smob);
  delete tree;
  return 0;
}

static int print_newick (SCM tree_smob, SCM port, scm_print_state *pstate)
{
  class PHYLIP_tree *tree = (class PHYLIP_tree *) SCM_SMOB_DATA (tree_smob);

  SExpr_atom tree_string;
  tree->write(tree_string,0);
  sstring tree_quoted;
  tree_quoted << tree_string;
  scm_puts (tree_quoted.c_str(), port);

  scm_remember_upto_here_1 (tree_smob);

  /* non-zero means success */
  return 1;
}

PHYLIP_tree* newick_cast_from_scm (SCM tree_smob)
{
  scm_assert_smob_type (newick_tag, tree_smob);
  return (PHYLIP_tree *) SCM_SMOB_DATA (tree_smob);
}

// main guile initialization routine
void init_newick_type (void)
{
  newick_tag = scm_make_smob_type ("newick", sizeof (class PHYLIP_tree));
  scm_set_smob_free (newick_tag, free_newick);
  scm_set_smob_print (newick_tag, print_newick);

  // read/write primitives
  scm_c_define_gsubr (GUILE_NEWICK_FROM_STRING, 1, 0, 0, (scm_t_subr) newick_from_string);
  scm_c_define_gsubr (GUILE_NEWICK_FROM_FILE, 1, 0, 0, (scm_t_subr) newick_from_file);
  scm_c_define_gsubr (GUILE_NEWICK_TO_FILE, 2, 0, 0, (scm_t_subr) newick_to_file);
  // primitives to ease migration from xrate macro format
  scm_c_define_gsubr (GUILE_NEWICK_ANCESTOR_LIST, 1, 0, 0, (scm_t_subr) newick_ancestor_list);  // returns list of internal node names (including the root, even if it is a tip node)
  scm_c_define_gsubr (GUILE_NEWICK_LEAF_LIST, 1, 0, 0, (scm_t_subr) newick_leaf_list);  // returns list of leaf node names (excluding the root)
  scm_c_define_gsubr (GUILE_NEWICK_BRANCH_LIST, 1, 0, 0, (scm_t_subr) newick_branch_list);  // returns list of (parent,child,length) tuples representing branches, sorted in preorder
  // convert a Newick tree into a Scheme data structure
  scm_c_define_gsubr (GUILE_NEWICK_UNPACK, 1, 0, 0, (scm_t_subr) newick_unpack);  // returns a tree structure very similar to the Newick file format: ((A:1,B:2)C:3,D:4)E;  --> ((("A" 1.0) ("B" 2.0) "C" 3.0) ("D" 4.0) "E" #f)
}
