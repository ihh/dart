#include "ecfg/ecfgmain.h"

typedef ECFG_main::ECFG_map ECFG_map;

// main program
int main (int argc, char* argv[])
{
  // create ECFG_main
  ECFG_main ecfg_main (argc, argv);

  // initialise grammars
  ecfg_main.add_grammar ("pfold", new PFOLD_ECFG());
  ecfg_main.add_grammar ("ifold", new IFOLD_ECFG());
  ecfg_main.add_grammar ("nullrna", new Null_RNA_ECFG());
  ecfg_main.add_grammar ("inullrna", new Null_indel_RNA_ECFG());
  ecfg_main.add_grammar ("dinuc", new Nearest_neighbor_ECFG());
  ecfg_main.add_grammar ("codon", new Codon_ECFG (true));

  // set default grammars
  ecfg_main.default_grammars = "pfold";

  // initialise the options parser
  ecfg_main.init_opts ("annotate RNA alignments using evolutionary grammars\n");

  // run
  ecfg_main.run();

  // delete grammars
  for_contents (ECFG_map, ecfg_main.ecfg_map, sg)
    delete (*sg).second;

  // return
  return 0;
}
