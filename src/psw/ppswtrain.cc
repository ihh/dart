#include "psw/pswmain.h"

int main (int argc, char** argv)
{
  // create Opts_list object using a macro
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <list-of-alignments file>",
		  "Baum-Welch parameter estimation for simple Protein Probabilistic Smith-Waterman alignment");

  PPSW_params params;  // protein PSW params
  PSW_trainer trainer (opts, params);
  trainer.run();
  return 0;
}
