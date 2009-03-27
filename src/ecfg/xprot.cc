#include "ecfg/ecfgmain.h"
#include "ecfg/protgrammar.h"
#include "tree/pam.h"

typedef ECFG_main::ECFG_map ECFG_map;

// main program
int main (int argc, char* argv[])
{
  // create ECFG_main
  ECFG_main ecfg_main (argc, argv);

  // initialise grammars
  ecfg_main.add_grammar ("prot4", new Protein_grammar (4));
  ecfg_main.add_grammar ("nullprot", new Protein_grammar (1));

  // set default grammars
  ecfg_main.default_grammars = "prot4";

  // initialise the options parser
  ecfg_main.init_opts ("annotate protein alignments using evolutionary grammars\n");

  // run
  ecfg_main.run();

  // delete grammars
  for_contents (ECFG_map, ecfg_main.ecfg_map, sg)
    delete (*sg).second;

  // return
  return 0;
}
