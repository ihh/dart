D			[0-9]
D1			[1-9]
D2			[2-9]
L			[a-zA-Z_]
E			[Ee][+-]?{D}+
Q                       '

%x comment

%{
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "tgcompiler.h"
#include "gram.tab.h"

#define YY_NO_UNPUT    /* suppresses warning (8/28/2004) */

int line_num = 1;

/* lex functions */
void count();
int yywrap();
char* stripquotes (char* s);

/* Telegraph_compiler interface */
#include "tgcompiler.h"
extern tg_grammar* tg;
%}

%%

"/*"                  { count(); BEGIN(comment); }
<comment>[^*\n]*      { count(); /* eat anything that's not a '*' */ }
<comment>"*"+[^*/\n]* { count(); /* eat up '*'s not followed by '/'s */ }
<comment>\n           { count(); ++line_num; }
<comment>"*"+"/"      { count(); BEGIN(INITIAL); }

{Q}.+{Q}              { count(); yylval.strval = stripquotes(yytext);
                        return(STRVAL); }

\"[^\"]*\"            { count(); yylval.strval = stripquotes(yytext);
                          return(STRVAL); }

"token"  	      { count(); return(TOKEN); }
"term"  	      { count(); return(TERM); }
"nonterm"	      { count(); return(NONTERM); }
"prob"  	      { count(); return(PROB); }

"end"		      { count(); return(END); }

"pset"		      { count(); return(PSET); }
"parse"               { count(); return(PARSE); }
"pmax"                { count(); return(PMAX); }
"psample"             { count(); return(PSAMPLE); }

"func"  	      { count(); return(FUNC); }
"psum"                { count(); return(PSUM); }
"dlog"                { count(); return(DLOG); }

{L}({L}|{D})*	      { count(); return tg_lookup_type
                        (tg, yytext, (void**) &yylval.strval); }

0                    { count(); yylval.intval = 0; return(ZERO); }
1                    { count(); yylval.intval = 1; return(ONE);  }
0+1{D}+		     { count(); yylval.intval = atoi(yytext); return(PLURAL); }
{D2}{D}*	     { count(); yylval.intval = atoi(yytext); return(PLURAL); }
{D}+{E}		     { count(); yylval.realval = atof(yytext); return(REAL); }
{D}*"."{D}+{E}       { count(); yylval.realval = atof(yytext); return(REAL); }
{D}+"."{D}*{E}       { count(); yylval.realval = atof(yytext); return(REAL); }
{D}+"."{D}+          { count(); yylval.realval = atof(yytext); return(REAL); }
"."{D}+              { count(); yylval.realval = atof(yytext); return(REAL); }

"->"		      { count(); return(PRODUCES); }
"=>"                  { count(); return(EQUALS_BITS); }

"$"		      { count(); return('$'); }
":"		      { count(); return(':'); }
";"		      { count(); return(';'); }
"|"		      { count(); return('|'); }
"["		      { count(); return('['); }
"]"		      { count(); return(']'); }
"{"		      { count(); return('{'); }
"}"		      { count(); return('}'); }
"'"		      { count(); return('\''); }
","		      { count(); return(','); }
"="		      { count(); return('='); }
"("		      { count(); return('('); }
")"		      { count(); return(')'); }
"+"		      { count(); return('+'); }
"*"		      { count(); return('*'); }
"/"		      { count(); return('/'); }

[ \t\v\n\f]	      { count(); /* eat up whitespace */ }
.		      { count(); return(isascii(*yytext)?*yytext:TGBADCHAR); }

%%

int yywrap()
{
    return(1);
}

int column = 0;

void count()
{
	int i;

	for (i = 0; yytext[i] != '\0'; i++)
		if (yytext[i] == '\n')
			column = 0;
		else if (yytext[i] == '\t')
			column += 8 - (column % 8);
		else
			column++;

	/* ECHO; */
}

char* stripquotes (char* s)
{
  char* ret;
  ret = strdup (s+1);
  ret[strlen(ret)-1] = '\0';
  return ret;
}
