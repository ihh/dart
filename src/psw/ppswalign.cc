#include "psw/pswmain.h"

int main (int argc, char** argv)
{
  // create Opts_list object using a macro
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <sequence file>",
		  "simple Protein Probabilistic Smith-Waterman alignment using a Pair HMM");
  
  PPSW_params params;  // protein PSW params
  PSW_aligner aligner (opts, params);
  aligner.run();
  return 0;
}
