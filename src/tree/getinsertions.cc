#include "tree/tree_alignment.h"

// maximum length for lines in Stockholm input alignments
#define MAX_LINE_LEN 12345678

// main program
int main (int argc, char** argv)
{
  // create Opts_list object using a macro
  INIT_OPTS_LIST (opts, argc, argv, -1, "[options] [<alignment database(s) in Stockholm format>]",
		  "split an alignment into insertion sub-alignments");

  opts.parse_or_die();  // parse the command-line options
  try
    {
      // initialise Stockholm_database
      Sequence_database seq_db;
      Stockholm_database stock_db;
      if (!opts.args.size())
	{
	  // no alignment filenames specified; read from stdin
	  CLOGERR << "[waiting for alignments on standard input]\n";
	  stock_db.read (cin, seq_db, MAX_LINE_LEN);
	}
      else
	for_const_contents (vector<sstring>, opts.args, align_db_filename)
	{
	  ifstream align_db_in ((const char*) align_db_filename->c_str());
	  if (!align_db_in) THROWEXPR ("Couldn't open alignment file '" << *align_db_filename << "'");
	  stock_db.read (align_db_in, seq_db, MAX_LINE_LEN);
	}

      // loop through alignments, splitting each into insertion subalignments
      for (int n_align = 0; n_align < stock_db.size(); n_align++)
	{
	  Stockholm& stock = *stock_db.align_index[n_align];
	  Insertion_database insertion_database (stock);
	  insertion_database.write (cout);
	}
    }
  catch (const Dart_exception& e)  // exception; bail out gracefully
    {
      CLOGERR << e.what();
      exit(1);
    }
  return 0;
}
