#ifndef SEXPR_VISITOR_INCLUDED
#define SEXPR_VISITOR_INCLUDED

#include "util/sexpr.h"
#include <map>

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include <libguile.h>
#endif /* GUILE_INCLUDED */

#include "util/guile-defs.h"
#include "util/sexpr-keywords.h"

// singleton shorthand map
struct SExpr_macro_aliases
{
  map<sstring,sstring> short2long;
  map<sstring,bool> warned;
  SExpr_macro_aliases();
  sstring expand(const sstring& op);
};
extern SExpr_macro_aliases sexpr_macro_aliases;

// abstract SExpr visitor class
struct SExpr_visitor
{
  // typedefs
  typedef SExpr::SExprIter SExprIter;
  // virtual destructor
  virtual ~SExpr_visitor() { }
  // virtual visitor method
  virtual void visit (SExpr& sexpr) { }
  // preorder & postorder visit wrappers
  void log_visit (SExpr& sexpr);
  void preorder_visit (SExpr& sexpr);
  void postorder_visit (SExpr& sexpr);
};

// SExpr file operations
struct SExpr_file_operations : SExpr_visitor
{
  // methods
  void visit (SExpr& sexpr);
};

// SExpr substitutions, iterators & replacements
struct SExpr_macros : SExpr_visitor
{
  // data
  map<sstring,sstring> replace;  // substitutions
  map<sstring,vector<sstring> > foreach;  // predefined "foreach" lists
  // methods
  void visit (SExpr& sexpr);
  void visit_child (SExpr& sexpr, SExprIter& child_iter, list<SExprIter>& erase_list);
  void visit_and_reap (SExpr& sexpr);  // cleans up erased children
  void handle_replace (SExpr& sexpr);
  void expand_foreach (SExpr& parent_sexpr, SExprIter& parent_pos, unsigned int element_offset, const vector<sstring>& foreach_list, list<SExprIter>& erase_list);
};

// SExpr list & logic operations
struct SExpr_list_operations : SExpr_visitor
{
  // methods
  void visit (SExpr& sexpr);
};

// SExpr Scheme macros using Guile
class SExpr_Scheme_evaluator
{
protected:
  // data
  void* (*register_functions)(void *);  // pointer to function that registers functions
  void* data;  // data that will be passed to *register_functions
  static bool initialized;
public:
  // constructor - sets register_functions and data to dummy values (override in subclasses)
  SExpr_Scheme_evaluator();
  // initialize() - initializes Guile
  // you must call this method before expand_Scheme_expressions.
  // calling it more than once will have no effect.
  void initialize();
  // mark_guile_initialized() - call this if you are intializing guile elsewhere and don't want to do it here
  static void mark_guile_initialized();
  // method to expand all eval/exec blocks in an SExpr tree
  // will throw an exception if an eval/exec block is encountered & program was compiled without Guile
  void expand_Scheme_expressions (SExpr& sexpr, const char* eval_keyword = SEXPR_EVAL, const char* eval_discard_keyword = SEXPR_EXEC) const;
  // helper methods to evaluate an SExpr as a Scheme expression using guile
#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
  static SCM evaluate_SCM (SExpr& sexpr);  // evaluates sexpr; returns a guile SCM object
#endif /* GUILE_INCLUDED */
  static SExpr evaluate_sexpr (SExpr& sexpr);  // evaluates sexpr; returns result as an SExpr
  static SExpr evaluate_values (SExpr& sexpr);  // evaluates all children of sexpr, except for the first child; returns the list of results as an SExpr
};

#endif /* SEXPR_VISITOR_INCLUDED */
