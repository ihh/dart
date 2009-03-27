#include "ecfg/ecfgmain.h"
#include "ecfg/protgrammar.h"
#include "tree/pam.h"

typedef ECFG_main::ECFG_map ECFG_map;

// main program
int main (int argc, char* argv[])
{
  // create ECFG_main
  ECFG_main ecfg_main (argc, argv);

  // add a null grammar
  ecfg_main.add_grammar ("null_grammar", new Null_ECFG (Roman_alphabet));

  // initialise the options parser
  ecfg_main.init_opts ("annotate multiple alignments using evolutionary grammars\n");

  // run
  ecfg_main.run();

  // delete grammars
  for_contents (ECFG_map, ecfg_main.ecfg_map, sg)
    delete (*sg).second;

  // return
  return 0;
}
