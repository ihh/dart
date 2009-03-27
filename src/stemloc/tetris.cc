#include <math.h>
#include "stemloc/tetrastem.h"
#include "scfg/paircfgdp.h"
#include "util/logfile.h"
#include "util/vector_output.h"
#include "scfg/wehits.h"

int main (int argc, char** argv)
{
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <sequence file>",
		  "tetraloop finder");

  int min_stem_len;
  int max_subseq_len;
  int max_hits;
  bool allow_GNRA;
  bool allow_UNCG;
  bool allow_XUUY;
  bool allow_UYU;
  bool allow_ARU;

  const int min_loop_len = 5;  // minimum length of "generic" loops (i.e. non-sequence-specific)
  int max_loop_len;

  opts.newline();
  opts.add ("l -minstemlen", min_stem_len = 4, "minimum stem length");
  opts.add ("m -maxlooplen", max_loop_len = 0, "maximum generic (i.e. >4 bases) loop length (by default such loops are not allowed; set this to 5+ to allow them)", FALSE);
  opts.add ("gnra", allow_GNRA = TRUE, "\t\tallow GNRA tetraloops");
  opts.add ("uncg", allow_UNCG = TRUE, "\t\tallow UNCG tetraloops");
  opts.add ("xuuy", allow_XUUY = TRUE, "\t\tallow XUUY tetraloops");
  opts.add ("uyu", allow_UYU = FALSE, "\t\tallow UYU triloops");
  opts.add ("aru", allow_ARU = FALSE, "\t\tallow ARU triloops");

  opts.newline();
  opts.add ("s -maxsubseq", max_subseq_len = 30, "maximum subsequence length");
  opts.add ("x -maxhits", max_hits = 100, "maximum number of hits to report per sequence");

  opts.parse_or_die();
  try
    {
      // get sequence filename
      const sstring seq_filename = opts.args[0];

      // read in sequences
      FASTA_sequence_database seq_db (seq_filename.c_str(), 0, Profile_flags_enum::NONE);
      const Alphabet& alphabet = CFG_alphabet;
      seq_db.seqs2dsqs (alphabet);

      // create Tetra_stem
      Tetra_stem tetra_stem (+1, -InfinityScore,
			     min_loop_len, max_loop_len,
			     allow_GNRA, allow_UNCG, allow_XUUY,
			     allow_UYU, allow_ARU);

      // loop through database, printing out GFF hit-sets
      for_const_contents (Sequence_database, seq_db, np)
	{
	  WE_hits hits;
	  hits.get_hits (tetra_stem, tetra_stem.paired_states, *np, TRUE, max_subseq_len, min_stem_len, max_hits);
	  cout << hits;
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
