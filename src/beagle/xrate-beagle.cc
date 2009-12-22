#include "ecfg/ecfgmain.h"
#include "ecfg/pfold.h"
#include "ecfg/protgrammar.h"

// TODO: uncomment the following line, and #define BEAGLE_INCLUDE_PATH somewhere
// #include BEAGLE_INCLUDE_PATH "/beagle.h"

typedef ECFG_main::ECFG_map ECFG_map;

// main program
int main (int argc, char* argv[])
{
  // create main xrate program object
  // TODO: define ECFG_main-derived class that uses Beagle, and use that class here instead of ECFG_main
  ECFG_main ecfg_main (argc, argv);

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
