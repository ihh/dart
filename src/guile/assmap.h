#ifndef ASS_MAP_INCLUDED
#define ASS_MAP_INCLUDED

#include <libguile.h>
#include <map>

#include "util/sstring.h"
#include "util/sexpr.h"

struct Ass_map : map<sstring,SCM>
{
  Ass_map (SCM scm);
  Ass_map (SExpr& sexpr, int offset = 1);
  SCM value_or_false (const char* tag) const;
  SCM value_or_die (const char* tag) const;
};


#endif /* ASS_MAP_INCLUDED */
