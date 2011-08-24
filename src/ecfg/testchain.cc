#include "ecfg/ecfgsexpr.h"
#include "util/unixenv.h"

// path to default chain file, from dart root
//#define DEFAULT_CHAIN_FILE "data/handalign/prot.hsm"
// changed here to test tokenized codon model: (OW 8/24/11)
#define DEFAULT_CHAIN_FILE "data/handalign/tcod.eg"

// wrapper for Alphabet::char2int
int char2int (const Alphabet& alphabet, const sstring& s, const char* desc)
{
  if (s.size() < 1)
    THROWEXPR (desc << " state string '" << s << "' has too few characters (must be exactly one char)");

  if (s.size() > 1)
    THROWEXPR (desc << " state string '" << s << "' has too many characters (must be exactly one char)");

  const char c = s[0];
  if (!alphabet.contains_strict (c))
    THROWEXPR (desc << " character '" << c << "' not found in alphabet");

  return alphabet.char2int_strict (c);
}


// main program
int main (int argc, char** argv)
{
  // initialise the options parser
  INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
		  "read a substitution matrix and exponentiate it\n");

  sstring default_chain_filename;
  default_chain_filename << Dart_Unix::get_DARTDIR() << '/' << DEFAULT_CHAIN_FILE;

  sstring chain_filename, src_char, dest_char;
  double time;

  opts.add ("c -chain-filename", chain_filename = default_chain_filename, "chain file");
  opts.add ("s -source-state", src_char = "a", "source state character");
  opts.add ("d -dest-state", dest_char = "a", "destination state character");
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
      SExpr_file ecfg_sexpr_file (chain_filename.c_str());
      SExpr& ecfg_sexpr = ecfg_sexpr_file.sexpr;

      // init alphabet, rate matrix
      Alphabet alphabet;  // dummy init of alphabet (from seq/biosequence.h)
      Irrev_EM_matrix rate_matrix (1, 1);  // dummy init of rate_matrix (from irrev/irrev_em_matrix.h)
      ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, ecfg_sexpr);

//       // print alphabet tokens, in order
      cout << "Alphabet symbols: " << alphabet.tokens() << '\n';


//       // get equilibrium distribution
       const vector<double> eqm_prob = rate_matrix.create_prior();

      // print equilibrium distribution
      cout << "Equilibrium distribution: " << eqm_prob << '\n';

      // convert alphabet symbols to state indices
	  const int src_state = 0;
	  const int dest_state = 1;

	  //       // exponentiate matrix
	  const array2d<double> cond_submat = rate_matrix.create_conditional_substitution_matrix (time);
	  
	  
	  
	  //       // extract required element
	  const double result = cond_submat (src_state, dest_state);
	  
	  cout << result<<endl;
      // print required element
	  //      cout << "exp(R*" << time << ")_{" << src_char << ',' << dest_char << "} = " << result << '\n';
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }

  return 0;
}
