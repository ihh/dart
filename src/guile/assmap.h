#ifndef ASS_MAP_INCLUDED
#define ASS_MAP_INCLUDED

#include <libguile.h>
#include <map>

#include "util/sstring.h"
#include "util/sexpr.h"

struct Ass_map : map<sstring,SCM>
{
  // constructors
  Ass_map (SCM scm);
  Ass_map (SExpr& sexpr, int offset = 1);
  // find wrapper
  const_iterator find (const char* tag) const;
  // accessors
  SCM scm_value_or_false (const char* tag) const;  // returns value for a given tag, or SCM_BOOL_F if tag was not foundy
  SCM scm_value_or_die (const char* tag) const;  // returns value for a given tag, or throws an exception if tag was not found
  SExpr sexpr_value_or_empty_list (const char* tag) const;  // returns a newly created SExpr (on the stack), or an empty-list SExpr if tag was not found
  SExpr sexpr_value_or_die (const char* tag) const;  // returns a pointer to a newly created SExpr, or throws an exception if tag was not found
  SExpr* new_sexpr_value_or_null (const char* tag) const;  // returns a pointer to a newly created SExpr, or NULL if tag was not found
  SExpr sexpr_values_list (const char* tag) const;  // returns a newly created SExpr (on the stack), containing a list of all the tag's values
  SExpr sexpr_parent_or_empty_list (const char* tag) const;  // returns a newly created SExpr (on the stack), containing the parent S-expression, or the empty list
};


#endif /* ASS_MAP_INCLUDED */
