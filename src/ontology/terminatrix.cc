#include "ecfg/ecfgsexpr.h"
#include "util/unixenv.h"

// path to default chain file, from dart root
#define DEFAULT_CHAIN_FILE "src/ontology/t/test1.tsm"

// default term identifier
#define DEFAULT_TERM_ID "term1"

// wrapper for Alphabet::tok2int
int tok2int (const Alphabet& alphabet, const sstring& s, const char* desc)
{
  if (!alphabet.contains_tok (s))
    THROWEXPR (desc << " token '" << s << "' not found in alphabet");

  return alphabet.token2int (s);
}


// main program
int main (int argc, char** argv)
{
  // initialise the options parser
  INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
		  "read a substitution matrix and exponentiate it\n");

  sstring default_chain_filename;
  default_chain_filename << Dart_Unix::get_DARTDIR() << '/' << DEFAULT_CHAIN_FILE;

  sstring chain_filename, src_tok, dest_tok;
  double time;

  opts.add ("c -chain-filename", chain_filename = default_chain_filename, "chain file");
  opts.add ("s -source-state", src_tok = DEFAULT_TERM_ID, "source term");
  opts.add ("d -dest-state", dest_tok = DEFAULT_TERM_ID, "destination term");
  opts.add ("t -time", time = 1., "rate matrix multiplier");

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
      // read chain file, get SExpr
      SExpr_file sexpr_file (chain_filename.c_str());
      SExpr& sexpr = sexpr_file.sexpr;

      // init alphabet, rate matrix
      Alphabet alphabet ("uninitialized", 1);  // dummy init of alphabet (from seq/biosequence.h)
      Irrev_EM_matrix rate_matrix (1, 1);  // dummy init of rate_matrix (from irrev/irrev_em_matrix.h)
      ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, sexpr, true);

      // convert alphabet symbols to state indices
      const int src_state = tok2int (alphabet, src_tok, "Source");
      const int dest_state = tok2int (alphabet, dest_tok, "Destination");

      // exponentiate matrix
      const array2d<double> cond_submat = rate_matrix.create_conditional_substitution_matrix (time);

      // get equilibrium distribution
      const vector<double> eqm_prob = rate_matrix.create_prior();

      // extract required element
      const double result = cond_submat (src_state, dest_state);

      // print alphabet tokens, in order
      cout << "Alphabet symbols: " << alphabet.tokens() << '\n';

      // print equilibrium distribution
      cout << "Equilibrium distribution: " << eqm_prob << '\n';

      // print required element
      cout << "exp(R*" << time << ")_{" << src_tok << ',' << dest_tok << "} = " << result << '\n';
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }

  return 0;
}
