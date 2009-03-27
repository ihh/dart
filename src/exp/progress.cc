#include "handel/proghmm.h"
#include "handel/progalign.h"

int main (int argc, char** argv)
{
  // initialise the options parser
  INIT_OPTS_LIST (opts, argc, argv, -1, "[options] <amino acid sequence database>",
		  "max-expected-SPS multiple alignment using the algorithm of Do, Brudno and Batzoglou\n");

  int indelTypes;
  int zones;

  bool train_only;
  sstring load_filename;
  sstring train_filename;
  sstring save_filename;

  opts.newline();
  opts.print ("Alignment options\n");
  opts.print ("-----------------\n");
  opts.newline();

  opts.add ("i -indeltypes", indelTypes = 2, "number of indel states");
  opts.add ("z -zones", zones = 1, "number of zones");

  opts.newline();
  opts.add ("lp -loadparams", load_filename = "", "load parameters from file (Telegraph format)", FALSE);
  opts.add ("tp -trainparams", train_filename = "", "train parameters from alignment file (Stockholm format)", FALSE);
  opts.add ("sp -saveparams", save_filename = "", "save parameters to file (Telegraph format)", FALSE);
  opts.add ("to -trainonly", train_only = FALSE, "don't do alignment, just do any specified training");

  // parse the command line & do stuff
  try
    {
      // parse command line
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
      if (opts.args.size() != 1 && !train_only) { cerr << opts.short_help(); exit(1); }

      // set alphabet
      const Alphabet& alphabet = Protein_alphabet;

      // make Prog_HMM
      Prog_HMM hmm (indelTypes, zones, alphabet);

      // make Optimal_accuracy_progressive_aligner
      Optimal_accuracy_progressive_aligner aligner (hmm, hmm.prior);

      // load params
      if (load_filename.size())
	aligner.load_params (load_filename.c_str());
      else
	aligner.init_default_params();

      // train params
      if (train_filename.size())
	aligner.train_params (train_filename.c_str());

      // save params
      if (save_filename.size())
	aligner.save_params (save_filename.c_str());

      // stop here if we're just training
      if (!train_only)
	{
	  // get sequence filename
	  const char* seqdb_filename = opts.args[0].c_str();
	  
	  // read sequences
	  FASTA_sequence_database seqdb (seqdb_filename, &alphabet);

	  // make alignment
	  Tree_alignment tree_align = aligner.make_progressive_alignment (seqdb.index);
	  
	  // show alignment
	  tree_align.align.write_MUL (cout, alphabet);
	}
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  return 0;
}
