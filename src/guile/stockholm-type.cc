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

static SCM stockholm_column_count (SCM stock_smob)
{
  Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
  return scm_from_int (stock->stock->columns());
}

// General function for creating (TAGVAL) SCM from vector<Stockholm::Tag_value> (for #=GF) or Stockholm::Annotation (other #=G's),
// where
// TAGVAL => (tag value) | TAGVAL TAGVAL | end
template<class C>
SCM tagval_list (C& container) {
  SCM tagval_list_scm = SCM_EOL;
  typedef typename C::const_iterator C_iter;
  for (C_iter iter = container.begin(); iter != container.end(); ++iter) {
    SCM tagval_item = scm_list_1 (scm_list_2 (scm_from_locale_string(iter->first.c_str()), scm_from_locale_string(iter->second.c_str())));
    tagval_list_scm = scm_append(scm_list_2(tagval_list_scm,tagval_item));
  }
  // Return
  return tagval_list_scm;
}

// stockholm-unpack creates the following structure
//    TOP => (GF GC BODY)
// TAGVAL => (tag value) | TAGVAL TAGVAL | end
//   BODY => (seqname GS rowdata GR) | BODY BODY | end
//     GF => (TAGVAL)
//     GC => (TAGVAL)
//     GS => (TAGVAL)
//     GR => (TAGVAL)
static SCM stockholm_unpack (SCM stock_smob)
{
  SCM row_list = SCM_EOL;
  Stockholm& stock (*Stockholm_smob::cast_from_scm (stock_smob)->stock);
  try {
    row_list = scm_list_2 (tagval_list(stock.gf_annot), tagval_list(stock.gc_annot));
    for (int r = 0; r < stock.rows(); ++r) {
      const sstring& row_name = stock.row_name[r];
      Stockholm::Row_annotation::const_iterator gs_iter = stock.gs_annot.find(row_name);
      Stockholm::Row_annotation::const_iterator gr_iter = stock.gr_annot.find(row_name);
      SCM row_scm = scm_list_1 (scm_list_4 (scm_from_locale_string(row_name.c_str()),
					    gs_iter == stock.gs_annot.end() ? SCM_EOL : tagval_list(gs_iter->second),
					    scm_from_locale_string(stock.get_row_as_string(r).c_str()),
					    gr_iter == stock.gr_annot.end() ? SCM_EOL : tagval_list(gs_iter->second)));
      row_list = scm_append(scm_list_2(row_list,row_scm));
    }

  } catch (Dart_exception& e) {
    CLOGERR << e.what();
  }

  return row_list;
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

  // convert a Stockholm alignment into a Scheme data structure with the following grammar
  //    TOP => (GF GC BODY)
  // TAGVAL => (tag value) | TAGVAL TAGVAL | end
  //   BODY => (seqname GS rowdata GR) | BODY BODY | end
  //     GF => (TAGVAL)
  //     GC => (TAGVAL)
  //     GS => (TAGVAL)
  //     GR => (TAGVAL)
  scm_c_define_gsubr ("stockholm-unpack", 1, 0, 0, (SCM (*)()) stockholm_unpack);  // returns a Scheme data structure very similar to the Stockholm file format: (GF GC (seq GS row GR) (seq GS row GR) ...)

}
