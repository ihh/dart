#ifndef SEXPR_VISITOR_INCLUDED
#define SEXPR_VISITOR_INCLUDED

#include "util/sexpr.h"
#include <map>

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include <libguile.h>
#endif /* GUILE_INCLUDED */

// keywords
#define SEXPR_WARN        "&warn"
#define SEXPR_DEFINE      "&define"

#define SEXPR_FOREACH     "&foreach"
#define SEXPR_FOREACH_INT "&foreach-integer"

#define SEXPR_CONCATENATE "&cat"
#define SEXPR_SUM         "&sum"
#define SEXPR_MULTIPLY    "&mul"
#define SEXPR_DIVIDE      "&div"
#define SEXPR_SUBTRACT    "&sub"
#define SEXPR_MODULUS     "&mod"
#define SEXPR_CONDITIONAL "&if"

#define SEXPR_SHORTHAND_CONCATENATE "&."
#define SEXPR_SHORTHAND_SUM         "&+"
#define SEXPR_SHORTHAND_MULTIPLY    "&*"
#define SEXPR_SHORTHAND_DIVIDE      "&/"
#define SEXPR_SHORTHAND_SUBTRACT    "&-"
#define SEXPR_SHORTHAND_MODULUS     "&%"
#define SEXPR_SHORTHAND_CONDITIONAL "&?"

#define SEXPR_TRUE        "1"
#define SEXPR_FALSE       "0"

#define SEXPR_EQUALS      "&eq"
#define SEXPR_NOT_EQUALS  "&neq"
#define SEXPR_GREATER     "&gt"
#define SEXPR_LESS        "&lt"
#define SEXPR_GEQ         "&geq"
#define SEXPR_LEQ         "&leq"

#define SEXPR_SHORTHAND_EQUALS      "&="
#define SEXPR_SHORTHAND_NOT_EQUALS  "&!="
#define SEXPR_SHORTHAND_GREATER     "&>"
#define SEXPR_SHORTHAND_LESS        "&<"
#define SEXPR_SHORTHAND_GEQ         "&>="
#define SEXPR_SHORTHAND_LEQ         "&<="

#define SEXPR_INCLUDE     "&include"

#define SEXPR_AND         "&and"
#define SEXPR_OR          "&or"
#define SEXPR_NOT         "&not"

#define SEXPR_INT         "&int"
#define SEXPR_CHR         "&chr"
#define SEXPR_ORD         "&ord"

#define SEXPR_EVAL        "&scheme"
#define SEXPR_EXEC        "&scheme-discard"


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
#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
  SCM write_proc;
#endif /* GUILE_INCLUDED */
public:
  // constructor - sets register_functions and data to dummy values (override in subclasses)
  SExpr_Scheme_evaluator();
  // initialize() - initializes Guile, sets write_proc
  // you must call this method before expand_Scheme_expressions
  void initialize();
  // mark_guile_initialized() - call this if you are intializing guile elsewhere and don't want to do it here
  static void mark_guile_initialized();
  // method to expand all eval/exec blocks in an SExpr tree
  // will throw an exception if an eval/exec block is encountered & program was compiled without Guile
  void expand_Scheme_expressions (SExpr& sexpr) const;
  // helper method to evaluate an SExpr as a Scheme expression using guile, then return the result encoded as another SExpr
  SExpr evaluate (SExpr& sexpr) const;
};

#endif /* SEXPR_VISITOR_INCLUDED */
