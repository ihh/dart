#include "ecfg/ecfgsexpr.h"
#include "ecfg/single_chain_ecfg.h"

// main program
int main (int argc, char* argv[])
{
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <hsm filename>", "convert a handalign-format substitution model into an xrate-format grammar");

  // parse the command line
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  // load rate matrix
  SExpr_file irrev_sexpr_file (opts.args[0].c_str());
  Irrev_EM_matrix rate_matrix;
  Alphabet alphabet;
  ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, irrev_sexpr_file.sexpr);

  // create grammar
  Single_chain_ECFG grammar (rate_matrix);

  // output
  ECFG_builder::ecfg2stream (cout, alphabet, grammar);
  ECFG_builder::alphabet2stream (cout, alphabet);

  // return
  return 0;
}
