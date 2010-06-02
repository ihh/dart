#include "ecfg/ecfgsexpr.h"

int main (int argc, char** argv)
{
  // initialise the options parser
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <S-expr file>",
		  "test ECFG S-expression parser\n");

  // parse the command line & do stuff
  try
    {
      // parse command line
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  // do stuff
  try
    {
      // get filename
      sstring ecfg_sexpr_filename = opts.args[0];

      // open file & read it
      ifstream ecfg_sexpr_file (ecfg_sexpr_filename.c_str());
      if (!ecfg_sexpr_file)
	THROWEXPR ("File not found: '" << ecfg_sexpr_filename << "'");

      sstring ecfg_sexpr_string, s;
      while (ecfg_sexpr_file && !ecfg_sexpr_file.eof())
	{
	  s.getline (ecfg_sexpr_file);
	  ecfg_sexpr_string << s;
	}
      ecfg_sexpr_file.close();

      // create SExpr
      SExpr ecfg_sexpr (ecfg_sexpr_string.begin(), ecfg_sexpr_string.end());

      // initialise ECFGs & alphabet
      Alphabet alph ("", 0);
      vector<ECFG_scores*> ecfg;
      ECFG_builder::init_grammars (alph, ecfg, ecfg_sexpr);

      // output ECFGs
      ECFG_builder::grammars2stream (cout, alph, ecfg);

      // delete ECFGs
      for_const_contents (vector<ECFG_scores*>, ecfg, g)
	delete *g;
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }

  return 0;
}
