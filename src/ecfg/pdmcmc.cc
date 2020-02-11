#include "ecfg/pdmcmc_main.h"
#include "ecfg/pfold.h"
#include "ecfg/protgrammar.h"

typedef PDMCMC_main::ECFG_map ECFG_map;

// main program
int main (int argc, char* argv[])
{
  // create PDMCMC_main
  PDMCMC_main ecfg_main (argc, argv);

  // initialise grammars
  ecfg_main.add_standard_grammars("rev");

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
  ecfg_main.delete_grammars();

  // return
  return 0;
}
