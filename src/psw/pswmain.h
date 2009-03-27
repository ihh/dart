#ifndef PSWMAIN_INCLUDED
#define PSWMAIN_INCLUDED

#include <math.h>

#include "psw/simplepsw.h"
#include "util/logfile.h"
#include "util/vector_output.h"
#include "hmm/postpairhmm.h"

struct PSW_aligner
{
  // references
  Opts_list& opts;
  Simple_PSW_params& params;
  const Alphabet& alphabet;

  // command-line options
  bool odds_ratio;  // flag indicating whether to use odds-ratios rather than likelihoods
  bool opt_acc;  // flag indicating whether to use "optimal accuracy" algorithm, rather than Viterbi
  sstring subtab_prefix;  // filename prefix for saving substitution tables
  sstring posttab_prefix;  // filename prefix for saving posterior tables
  sstring submat_filename;  // filename for loading substitution matrix
  sstring params_filename;  // filename for loading all parameters

  // constructor
  PSW_aligner (Opts_list& opts, Simple_PSW_params& params);

  // run method
  void run();  // parses command-line options, etc

  // method to output an alignment
  void write_alignment (const Alignment& align, const Alphabet& alphabet, double score, const char* score_units);
};

struct PSW_trainer
{
  // references
  Opts_list& opts;
  Simple_PSW_params& params;
  const Alphabet& alphabet;
  
  // command-line parameters
  int max_iter;  // max EM iterations
  bool symmetrise;  // flag indicating whether to forcibly "symmetrise" model by flipping alignments
  bool ignore_align;   // flag indicating whether to ignore alignments provided

  // constructor
  PSW_trainer (Opts_list& opts, Simple_PSW_params& params);

  // run method
  void run();
};

#endif /* PSWMAIN_INCLUDED */
