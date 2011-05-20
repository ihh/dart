#include "ontology/onto_sexpr.h"
#include "util/unixenv.h"
#include "ecfg/guile-ecfg.h"

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

  sstring chain_filename, write_filename;
  bool compute_evidence;

  opts.print_title ("Modeling options");

  opts.add ("c -chain-filename", chain_filename = default_chain_filename, "file to load model from");
  opts.add ("w -write-filename", write_filename = "", "file to save model to", false);
  opts.add ("e -evidence", compute_evidence = false, "compute log-evidences");

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

      // create Scheme context
      ECFG_Scheme_evaluator scheme;
      scheme.initialize();

      // init the Terminatrix
      Terminatrix term (scheme);
      Terminatrix_builder::init_terminatrix (term, sexpr);
      term.eval_funcs();

      // calculate evidences, if asked
      if (compute_evidence)
	{
	  Terminatrix_log_evidence log_ev (term);
	  SExpr log_ev_sexpr = log_ev.map_reduce_sexpr();
	  cout << ";; log-evidence\n";
	  cout << log_ev_sexpr.to_parenthesized_string() << '\n';
	}

      // save
      if (write_filename.size())
	{
	  ofstream write_file (write_filename.c_str());
	  Terminatrix_builder::terminatrix2stream (write_file, term);
	}
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }

  return 0;
}
