#include <stdlib.h>
#include <libguile.h>

#include "tree/tree_alignment.h"
#include "guile/stockholm-type.h"
#include "guile/newick-type.h"

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
  free(s);

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

static SCM stockholm_tree (SCM stock_smob)
{
  SCM scm = SCM_BOOL_F;
  Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
  try {
    Stockholm_tree tree (stock->stock);
    scm = make_newick_smob (tree);
  } catch (Dart_exception& e) {
    CLOGERR << e.what();
  }
  return scm;
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

  // read/write primitives
  scm_c_define_gsubr ("stockholm-from-string", 1, 0, 0, (SCM (*)()) stockholm_from_string);
  scm_c_define_gsubr ("stockholm-from-file", 1, 0, 0, (SCM (*)()) stockholm_from_file);
  scm_c_define_gsubr ("stockholm-to-file", 2, 0, 0, (SCM (*)()) stockholm_to_file);
  // primitives to ease migration from xrate macro format
  scm_c_define_gsubr ("stockholm-column-count", 1, 0, 0, (SCM (*)()) stockholm_column_count);  // returns the number of columns as an integer
  scm_c_define_gsubr ("stockholm-tree", 1, 0, 0, (SCM (*)()) stockholm_tree);  // returns a newick-type smob constructed from the "#=GF NH" tag of the Stockholm alignment, or FALSE if no tree present
  // currently there is no method to convert a Stockholm smob into a Scheme data structure, and the set of accessors is incomplete.
  // for example, you cannot currently read the "#=GC SS_cons" line.
  // if you want anything else you are going to have to either (a) parse the string form or (b) code up a new subroutine here.
  // the recommended subroutine is an analog of newick_branch_list (in newick-type.cc) that returns a Stockholm smob as a flattish Scheme data structure,
  // e.g.
  //    TOP => (GF GC BODY)
  // TAGVAL => (tag value) | TAGVAL TAGVAL | end
  //   BODY => (seqname GS rowdata GR) | BODY BODY | end
  //     GF => (TAGVAL)
  //     GC => (TAGVAL)
  //     GS => (TAGVAL)
  //     GR => (TAGVAL)
}
