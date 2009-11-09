#include <stdlib.h>
#include <libguile.h>

#include "seq/stockholm.h"
#include "tree/tree_alignment.h"
#include "util/sexpr.h"
#include "guile/stockholm-type.h"

static scm_t_bits stockholm_tag;

struct Stockholm_smob {
  // data
  Stockholm stock;
  Sequence_database seqdb;

  // methods (just wrappers for Stockholm methods)
  void read (const char* filename) {
    ifstream infile(filename);
    stock.read_Stockholm (infile, seqdb);
  }

  void write (const char* filename) {
    ofstream outfile(filename);
    stock.write_Stockholm (outfile);
  }
};

static SCM
read_stockholm (SCM s_filename)
{
  SCM smob;

  // read alignment from file
  Stockholm_smob *stock = new Stockholm_smob();
  char* filename = scm_to_locale_string (s_filename);
  stock->read(filename);
  free(filename);

  // Create the smob.
  SCM_NEWSMOB (smob, stockholm_tag, stock);

  // Return
  return smob;
}

static SCM
write_stockholm (SCM stock_smob, SCM s_filename)
{
  scm_assert_smob_type (stockholm_tag, stock_smob);
  Stockholm_smob *stock = (Stockholm_smob *) SCM_SMOB_DATA (stock_smob);

  // write alignment to file
  char* filename = scm_to_locale_string (s_filename);
  stock->write(filename);
  free(filename);

  // Return
  return SCM_UNSPECIFIED;
}

static SCM
stockholm_node_list (PHYLIP_tree& tree, vector<Phylogeny::Node> (PHYLIP_tree::*getNodeMethod)())
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

static SCM
stockholm_leaf_list (SCM stock_smob)
{
  scm_assert_smob_type (stockholm_tag, stock_smob);
  Stockholm_smob *stock = (Stockholm_smob *) SCM_SMOB_DATA (stock_smob);

  Stockholm_tree tree (stock->stock);
  return stockholm_node_list (tree, &PHYLIP_tree::leaf_vector);
}

static SCM
stockholm_ancestor_list (SCM stock_smob)
{
  scm_assert_smob_type (stockholm_tag, stock_smob);
  Stockholm_smob *stock = (Stockholm_smob *) SCM_SMOB_DATA (stock_smob);

  Stockholm_tree tree (stock->stock);
  return stockholm_node_list (tree, &PHYLIP_tree::ancestor_vector);
}

static size_t
free_stockholm (SCM stock_smob)
{
  struct Stockholm_smob *stock = (struct Stockholm_smob *) SCM_SMOB_DATA (stock_smob);
  delete stock;
  return 0;
}

static int
print_stockholm (SCM stock_smob, SCM port, scm_print_state *pstate)
{
  struct Stockholm_smob *stock = (struct Stockholm_smob *) SCM_SMOB_DATA (stock_smob);

  sstring stock_string;
  stock->stock.write_Stockholm(stock_string);
  scm_puts (stock_string.c_str(), port);

  /* non-zero means success */
  return 1;
}

static SExpr*
scm_to_sexpr (SCM scm)
{
  // four guile API calls to get an SCM as a char* string? feel like I'm doing something the hard way here
  const char *s = scm_to_locale_string (scm_object_to_string (scm, scm_variable_ref (scm_c_lookup ("write"))));
  sstring str (s);
  SExpr* sexpr = new SExpr (str.begin(), str.end());
  free((void*) s);
  return sexpr;
}

static SCM
string_to_scm (const char* s)
{
  sstring str;
  str << "(quote " << s << ")";
  SCM scm = scm_c_eval_string(str.c_str());
  return scm;
}

static SCM
sexpr_to_scm (SExpr* sexpr)
{
  sstring str;
  str << *sexpr;
  SCM scm = string_to_scm(str.c_str());
  return scm;
}

// test function that converts an SCM to an SExpr, and back again
static SCM
test_convert_scm (SCM scm) {
  SExpr* sexpr = scm_to_sexpr(scm);
  SCM scm2 = sexpr_to_scm(sexpr);
  delete sexpr;
  return scm2;
}

void
init_stockholm_type (void)
{
  stockholm_tag = scm_make_smob_type ("stockholm", sizeof (struct Stockholm_smob));
  scm_set_smob_free (stockholm_tag, free_stockholm);
  scm_set_smob_print (stockholm_tag, print_stockholm);

  scm_c_define_gsubr ("read-stockholm", 1, 0, 0, (SCM (*)()) read_stockholm);
  scm_c_define_gsubr ("write-stockholm", 2, 0, 0, (SCM (*)()) write_stockholm);
  scm_c_define_gsubr ("stockholm-leaf-list", 1, 0, 0, (SCM (*)()) stockholm_leaf_list);
  scm_c_define_gsubr ("stockholm-ancestor-list", 1, 0, 0, (SCM (*)()) stockholm_ancestor_list);

  scm_c_define_gsubr ("test-convert-scm", 1, 0, 0, (SCM (*)()) test_convert_scm);
}
