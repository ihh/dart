#include <sstream>
#include "seq/stockholm.h"

int main (int argc, char** argv)
{
  try
    {
      const sstring stockfile ("# STOCKHOLM 1.0\n\nseq1 AAA\nseq2 AAG\nseq3 GTC\n");
      istringstream stockfile_inp (stockfile);
      Stockholm_database stock_db;
      Sequence_database seq_db;
      stock_db.read (stockfile_inp, seq_db);
      if (stock_db.size() != 1) THROWEXPR ("Number of alignments = " << stock_db.size() << ", should be 1");
      const Stockholm& stock = *stock_db.align.begin();
      if (stock.rows() != 3) THROWEXPR ("Number of rows in first alignment = " << stock.rows() << ", should be 3");

      // if an argument was specified, try reading in a specified file
      if (argc > 1)
	{
	  cerr << "(reading alignment database from '" << argv[1] << "')\n";
	  Stockholm_database sdb;
	  ifstream inp (argv[1]);
	  sdb.read (inp, seq_db);
	  cerr << "(database has " << sdb.size() << " alignments)\n";
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cerr << "(Stockholm databases tested)\n";
  cout << "ok\n";
  return 0;
}
