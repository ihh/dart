#ifndef SEXPR_VISITOR_INCLUDED
#define SEXPR_VISITOR_INCLUDED

#include "util/sexpr.h"
#include <map>

// keywords
#define SEXPR_WARN        "&warn"
#define SEXPR_DEFINE      "&define"

#define SEXPR_FOREACH     "&foreach"
#define SEXPR_FOREACH_INT "&foreach-integer"

#define SEXPR_CONCATENATE "."
#define SEXPR_SUM         "+"
#define SEXPR_MULTIPLY    "*"
#define SEXPR_DIVIDE      "/"
#define SEXPR_SUBTRACT    "-"
#define SEXPR_MODULUS     "%"

#define SEXPR_CONDITIONAL "?"

#define SEXPR_TRUE        "1"
#define SEXPR_FALSE       "0"

#define SEXPR_EQUALS      "="
#define SEXPR_NOT_EQUALS  "!="
#define SEXPR_GREATER     ">"
#define SEXPR_LESS        "<"
#define SEXPR_GEQ         ">="
#define SEXPR_LEQ         "<="

#define SEXPR_INCLUDE     "&include"

#define SEXPR_AND         "&and"
#define SEXPR_OR          "&or"
#define SEXPR_NOT         "&not"

#define SEXPR_INT         "&int"
#define SEXPR_CHR         "&chr"
#define SEXPR_ORD         "&ord"

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

#endif /* SEXPR_VISITOR_INCLUDED */
