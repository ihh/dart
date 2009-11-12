#include <stdlib.h>
#include <libguile.h>

#include "tree/tree_alignment.h"
#include "guile/stockholm-type.h"

scm_t_bits stockholm_tag;

SCM make_stockholm_smob (const Stockholm& stock)
{
  SCM smob;
  Stockholm_smob *stock_smob = new Stockholm_smob (stock);
  SCM_NEWSMOB (smob, stockholm_tag, stock_smob);
  return smob;
}

static SCM stockholm_from_file (SCM s_filename)
{
  SCM smob;

  // read alignment from file
  Stockholm_smob *stock = new Stockholm_smob();
  char* filename = scm_to_locale_string (s_filename);
  stock->read_from_file(filename);
  free(filename);

  // Create the smob.
  SCM_NEWSMOB (smob, stockholm_tag, stock);

  // Return
  return smob;
}

static SCM stockholm_from_string (SCM s_string)
{
  SCM smob;

  // read alignment from file
  Stockholm_smob *stock = new Stockholm_smob();
  char* s = scm_to_locale_string (s_string);
  stock->read_from_string(s);

  // Create the smob.
  SCM_NEWSMOB (smob, stockholm_tag, stock);

  // Return
  return smob;
}

static SCM stockholm_to_file (SCM stock_smob, SCM s_filename)
{
  Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);

  // write alignment to file
  char* filename = scm_to_locale_string (s_filename);
  stock->write_to_file(filename);
  free(filename);

  // Return
  return SCM_UNSPECIFIED;
}

static SCM stockholm_node_list (PHYLIP_tree& tree, vector<Phylogeny::Node> (PHYLIP_tree::*getNodeMethod)())
{
  SCM tree_node_list = SCM_BOOL_F;
  bool list_empty = true;

  vector<Phylogeny::Node> nodes = (tree.*getNodeMethod)();
  for_const_contents (vector<Phylogeny::Node>, nodes, iter) {
      Phylogeny::Node n = *iter;
      const char* n_name = tree.node_specifier(n).c_str();
      SCM single_element_list = scm_list_1(scm_from_locale_string(n_name));
      if (list_empty) {
	tree_node_list = single_element_list;
	list_empty = false;
      } else
	scm_append_x(scm_list_2(tree_node_list,single_element_list));
  }

  // Return
  return tree_node_list;
}

static SCM stockholm_leaf_list (SCM stock_smob)
{
  Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
  Stockholm_tree tree (stock->stock);
  return stockholm_node_list (tree, &PHYLIP_tree::leaf_vector);
}

static SCM stockholm_ancestor_list (SCM stock_smob)
{
  Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
  Stockholm_tree tree (stock->stock);
  return stockholm_node_list (tree, &PHYLIP_tree::ancestor_vector);
}

static SCM stockholm_column_count (SCM stock_smob)
{
  Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
  return scm_from_int (stock->stock.columns());
}

static size_t free_stockholm (SCM stock_smob)
{
  struct Stockholm_smob *stock = (struct Stockholm_smob *) SCM_SMOB_DATA (stock_smob);
  delete stock;
  return 0;
}

static int print_stockholm (SCM stock_smob, SCM port, scm_print_state *pstate)
{
  struct Stockholm_smob *stock = (struct Stockholm_smob *) SCM_SMOB_DATA (stock_smob);

  sstring stock_string;
  stock->write_to_string(stock_string);
  scm_puts (stock_string.c_str(), port);

  /* non-zero means success */
  return 1;
}

Stockholm_smob* Stockholm_smob::cast_from_scm (SCM stock_smob)
{
  scm_assert_smob_type (stockholm_tag, stock_smob);
  return (Stockholm_smob *) SCM_SMOB_DATA (stock_smob);
}

// main guile initialization routine
void init_stockholm_type (void)
{
  stockholm_tag = scm_make_smob_type ("stockholm", sizeof (struct Stockholm_smob));
  scm_set_smob_free (stockholm_tag, free_stockholm);
  scm_set_smob_print (stockholm_tag, print_stockholm);

  scm_c_define_gsubr ("stockholm-from-file", 1, 0, 0, (SCM (*)()) stockholm_from_file);
  scm_c_define_gsubr ("stockholm-to-file", 2, 0, 0, (SCM (*)()) stockholm_to_file);
  scm_c_define_gsubr ("stockholm-from-string", 1, 0, 0, (SCM (*)()) stockholm_from_string);
  scm_c_define_gsubr ("stockholm-column-count", 1, 0, 0, (SCM (*)()) stockholm_column_count);
  scm_c_define_gsubr ("stockholm-ancestor-list", 1, 0, 0, (SCM (*)()) stockholm_ancestor_list);
  scm_c_define_gsubr ("stockholm-leaf-list", 1, 0, 0, (SCM (*)()) stockholm_leaf_list);
}
