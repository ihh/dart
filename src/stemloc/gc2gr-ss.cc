#include <math.h>
#include <stack>

#include "seq/stockholm.h"
#include "scfg/foldenv.h"

int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <alignment database>",
		  "utility to map consensus secondary structures to by-sequence structures in a Stockholm alignment database");

  bool over;
  opts.add ("over", over = TRUE, "override existing structures for each sequence, if any exist");

  opts.parse_or_die();
  try
    {
      // get alignment database filename
      const sstring stock_db_filename = opts.args[0];

      // read in alignments
      Sequence_database seq_db;
      Stockholm_database stock_db;
      ifstream stock_db_file (stock_db_filename.c_str());
      stock_db.read (stock_db_file, seq_db);

      // do it
      stock_db.propagate_consensus_folds (over);

      // write alignments back out
      stock_db.write (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
