#ifndef SCHEME_KEYWORDS_INCLUDED
#define SCHEME_KEYWORDS_INCLUDED

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

#endif /* SCHEME_KEYWORDS_INCLUDED */
