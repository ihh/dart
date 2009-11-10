#include "ecfg/ecfgmain.h"
#include "ecfg/pfold.h"
#include "ecfg/protgrammar.h"

typedef ECFG_main::ECFG_map ECFG_map;

// main program
int main (int argc, char* argv[])
{
  // create ECFG_main
  Sequence_database seq_db;
  ECFG_main ecfg_main (argc, argv, seq_db);

  // initialise grammars
  ecfg_main.add_grammar ("rev", new Null_ECFG (DNA_alphabet, true));
  ecfg_main.add_grammar ("irrev", new Null_ECFG (DNA_alphabet, false));
  ecfg_main.add_grammar ("dinuc", new Nearest_neighbor_ECFG());
  ecfg_main.add_grammar ("codon", new Codon_ECFG (true));
  ecfg_main.add_grammar ("aa", new Protein_grammar(1,1));
  ecfg_main.add_grammar ("aa2", new Protein_grammar(1,2));
  ecfg_main.add_grammar ("aa3", new Protein_grammar(1,3));
  ecfg_main.add_grammar ("aa4", new Protein_grammar(1,4));

  // set default grammars
  ecfg_main.default_grammars = "rev";

  // initialise the options parser
  ecfg_main.init_opts ("annotate alignments and estimate evolutionary rates\n");

  // run
  try {
    ecfg_main.run_xrate (cout);
  } catch (const Dart_exception& e)
    {
      CLOGERR << "ERROR: " << e.what();
      exit(1);
    }

  // delete grammars
  for_contents (ECFG_map, ecfg_main.ecfg_map, sg)
    delete (*sg).second;

  // return
  return 0;
}
