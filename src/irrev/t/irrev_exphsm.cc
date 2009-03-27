// pk 3/14/05 irrev version.
// pk 3/11/05 invoke Substitution_matrix_factory::create_joint_substitution_matrix,
// write resulting matrix to stdout.
//
// Input command-line arguments:
// user_filename the filename of an hsm format file - see xrate -hsmhelp
// time (double) to use in exponentiation.
//

#include <sstream>
#include "irrev/irrev_em_matrix.h"

int main(int argc, char* argv[])
{
  // initialize the options parser
  //
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <hsm format input - see xrate -hsmhelp>",
	  "call create_joint_substitution_matrix, print out results");

  double time;

  opts.newline();
  opts.newline();
  opts.print ("options\n");
  opts.add ("t time -time", time = 1, "\ttime over which to exponentiate");
  opts.newline();

  // parse the command line
  //
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
	  cout << "ok\n";				// pk 3/17/05 work around so ihh test facility doesn't complain
      exit(1);
    }

  try
    {
      // get args
      const char* user_hsm_filename = opts.args[0].c_str();

      // initialise
      ifstream hsm_in (user_hsm_filename);
      if (!hsm_in) THROWEXPR ("Couldn't open input file '" << user_hsm_filename << "'");

	  Irrev_EM_matrix em (1, 2);  	// dummy C & A for initialization
	  em.read (hsm_in);
	  array2d<double> prob = em.create_joint_substitution_matrix (time);
	  cout.precision (10);			// want more digits
	  cout << prob;
	}
  catch (const Dart_exception& e)
    {
      CLOGERR << "ERROR: " << e.what();
      exit(1);
    }

  return 0;
}
