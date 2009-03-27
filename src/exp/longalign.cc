#include "ecfg/swexmid.h"
#include "util/opts_list.h"

// main program
int main(int argc, char* argv[])
{
  // initialise the options parser
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <amino acid sequence database>",
		  "optimal-accuracy neighbor-joining progressive alignment using a two-event long indel approximation\n");

  // parse the command line
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }

      // get args
      const char* seqdb_filename = opts.args[0].c_str();

      // set alphabet
      const Alphabet& alphabet = Protein_alphabet;

      // read sequences
      FASTA_sequence_database seqdb (seqdb_filename, &alphabet);

      // make SWEXMID_params
      Protein_default_SWEXMID_params params;

      // make SWEXMID_AXY_PHMM
      SWEXMID_AXY_PHMM hmm (params.indelTypes(), alphabet);
 
      // make alignment
      Tree_alignment tree_align = hmm.progressive_alignment (seqdb.index, params);

      // show alignment
      tree_align.align.write_MUL (cout, alphabet);
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  return 0;
}
