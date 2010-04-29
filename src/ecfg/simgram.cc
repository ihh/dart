#include "ecfg/ecfgsexpr.h"
#include "util/vector_output.h"
#include "irrev/irrev_em_matrix.h"

// main program
int main (int argc, char* argv[])
{
  INIT_OPTS_LIST (opts, argc, argv, 0, "<options>", "sample alignments from an xrate grammar file");

  opts.print_title ("Simulation");

  sstring grammar_filename;
  sstring tree_filename;
  int n_align;
  opts.add ("g -grammar", grammar_filename, "xrate format grammar file (required)", false);
  opts.add ("t -tree", tree_filename, "Newick format tree file (required)", false);
  opts.add ("n -nalign", n_align = 1, "number of alignments to generate");

  sstring gap_chars, gap_chars_help;
  gap_chars_help << "set character(s) that are used to denote gaps in alignment (default is \"" << DEFAULT_GAP_CHARS << "\")";
  opts.add ("gc -gapchar", gap_chars = "", gap_chars_help.c_str(), false);

  // parse the command line
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  try
    {
      // get args
      if (!grammar_filename.size() || !tree_filename.size())
	THROWEXPR ("Missing grammar or tree file\n"
		   << opts.help()
		   << "\n\nPlease specify both a grammar file (-g) and a tree file (-t)\n");

      // set up the gap characters
      if (gap_chars.size())
	Alignment::set_gap_chars (gap_chars);

      // load grammar
      Empty_alphabet alph;
      vector<ECFG_scores*> grammar_vec;
      ECFG_builder::load_xgram_alphabet_and_grammars (grammar_filename, alph, grammar_vec);
      if (grammar_vec.size() != 1)
	THROWEXPR (grammar_vec.size() << " grammars in file -- expected one");

      // load tree
      ifstream tree_stream (tree_filename.c_str());
      if (!tree_stream)
	THROWEXPR ("Can't open tree file " << tree_filename);
      PHYLIP_tree tree;
      tree.read (tree_stream);

      // loop over N alignments
      for (int n = 0; n < n_align; ++n)
	{
	  // log
	  CTAG(6,SIMGRAM) << "Generating alignment #" << n+1 << '\n';

	  // generate alignment
	  ECFG_scores& gram = *grammar_vec.front();
	  ECFG_simulation ecfg_sim (gram, tree);
	  ecfg_sim.add_counts_to_stockade();

	  // output it
	  ecfg_sim.align.write_Stockholm (cout);
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << "ERROR: " << e.what();
      exit(1);
    }
}
