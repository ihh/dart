#include <stdlib.h>
#include <libguile.h>

#include "seq/stockholm.h"
#include "tree/tree_alignment.h"
#include "util/sexpr.h"
#include "ecfg/ecfgmain.h"
#include "ecfg/ecfgsexpr.h"
#include "guile/guile-inits.h"

static SCM xrate_estimate_tree (SCM stock_smob, SCM alphabet_and_grammar)
{
  SCM scm = SCM_BOOL_F;
  SExpr* sexpr = 0;
  try {
    Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
    sexpr = scm_to_new_sexpr (alphabet_and_grammar);

    ECFG_main xrate;
    Stockholm stock_with_tree = xrate.run_tree_estimation (stock->stock, stock->seqdb, (*sexpr)[0]);

    // make return expression
    scm = make_stockholm_smob (stock_with_tree);

    // delete the temporarily created SExpr
    delete sexpr;
    sexpr = 0;

  } catch (Dart_exception& e) {
    CLOGERR << e.what();
    if (sexpr)
      delete sexpr;
  }

  return scm;
}

static SCM xrate_annotate_alignment (SCM stock_smob, SCM alphabet_and_grammar)
{
  SCM scm = SCM_BOOL_F;
  SExpr* sexpr = 0;
  try {
    Stockholm_smob *stock = Stockholm_smob::cast_from_scm (stock_smob);
    sexpr = scm_to_new_sexpr (alphabet_and_grammar);

    ECFG_main xrate;
    Stockholm stock_with_annotation = xrate.run_alignment_annotation (stock->stock, (*sexpr)[0]);

    // make return expression
    scm = make_stockholm_smob (stock_with_annotation);

    // delete the temporarily created SExpr
    delete sexpr;
    sexpr = 0;

  } catch (Dart_exception& e) {
    CLOGERR << e.what();
    if (sexpr)
      delete sexpr;
  }

  return scm;
}

static SCM xrate_train_grammar (SCM list_of_stock_smobs, SCM alphabet_and_grammar)
{
  SCM scm = SCM_BOOL_F;
  SExpr* sexpr = 0;
  try {
    Stockholm_database stock_db;
    sexpr = scm_to_new_sexpr (alphabet_and_grammar);

    /* Check that first argument is a list */ 
    if (SCM_NFALSEP(scm_list_p(list_of_stock_smobs)))
      {
	/* Iterate through the list */ 
	while (SCM_FALSEP (scm_null_p (list_of_stock_smobs)))
	  {
	    /* Get the head of the list */ 
	    SCM stock_smob = SCM_CAR(list_of_stock_smobs);

	    /* Add it to the Stockholm_database */
	    stock_db.add (Stockholm_smob::cast_from_scm(stock_smob)->stock);

	    /* Discard the head of the list */
	    list_of_stock_smobs = SCM_CDR(list_of_stock_smobs);
	  }
      }
    else   // first arg is not a list, so attempt to treat it as a single alignment
      stock_db.add (Stockholm_smob::cast_from_scm(list_of_stock_smobs)->stock);

    // do the training
    ECFG_main xrate;
    ECFG_scores* trained_grammar = 0;
    ECFG_counts* trained_counts = 0;
    xrate.run_grammar_training (stock_db, (*sexpr)[0], &trained_grammar, &trained_counts);

    // make return expression
    scm = ecfg_to_scm (*trained_grammar, trained_counts);

    // delete the temporarily created SExpr
    delete sexpr;
    sexpr = 0;

  } catch (Dart_exception& e) {
    CLOGERR << e.what();
    if (sexpr)
      delete sexpr;
  }

  return scm;
}

SExpr* scm_to_new_sexpr (SCM scm)
{
  // four guile API calls to get an SCM as a char* string? feel like I'm doing something the hard way here
  const char *s = scm_to_locale_string (scm_object_to_string (scm, scm_variable_ref (scm_c_lookup ("write"))));
  sstring str (s);
  SExpr* sexpr = new SExpr (str.begin(), str.end());
  free((void*) s);
  return sexpr;
}

SCM string_to_scm (const char* s)
{
  sstring str;
  str << "(quote " << s << ")";
  SCM scm = scm_c_eval_string(str.c_str());
  return scm;
}

SCM unparenthesized_string_to_scm (const char* s)
{
  sstring str;
  str << '(' << s << ')';
  return string_to_scm (str.c_str());
}

SCM sexpr_to_scm (SExpr* sexpr)
{
  sstring str;
  str << *sexpr;
  SCM scm = string_to_scm(str.c_str());
  return scm;
}

SCM ecfg_to_scm (const ECFG_scores& ecfg, const ECFG_counts* counts)
{
  sstring grammar_alphabet_string;
  ECFG_builder::ecfg2stream (grammar_alphabet_string, ecfg.alphabet, ecfg, counts);
  ECFG_builder::alphabet2stream (grammar_alphabet_string, ecfg.alphabet);
  return unparenthesized_string_to_scm (grammar_alphabet_string.c_str());
}

// test function that converts a Guile SCM to a Dart SExpr, and back again
static SCM
test_convert_scm (SCM scm)
{
  SExpr* sexpr = scm_to_new_sexpr(scm);
  SCM scm2 = sexpr_to_scm(sexpr);
  delete sexpr;
  return scm2;
}

void init_xrate_primitives (void)
{
  scm_c_define_gsubr ("xrate-estimate-tree", 2, 0, 0, (SCM (*)()) xrate_estimate_tree);
  scm_c_define_gsubr ("xrate-annotate-alignment", 2, 0, 0, (SCM (*)()) xrate_annotate_alignment);
  scm_c_define_gsubr ("xrate-train-grammar", 2, 0, 0, (SCM (*)()) xrate_train_grammar);

  scm_c_define_gsubr ("test-convert-scm", 1, 0, 0, (SCM (*)()) test_convert_scm);
}
