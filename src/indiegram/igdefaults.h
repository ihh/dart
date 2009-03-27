#ifndef IGDEFAULTS_INCLUDED
#define IGDEFAULTS_INCLUDED

#include "telegraph/tgvar.h"
#include "indiegram/tripletscfg.h"

struct Indiegram_defaults
{
  /// Default grammar for building the fold envelopes.
  /*
   * Relies on stemloc.
   */
  static void init_superstem_params (Telegraph_PScores_adaptor& fold_tgio);

  /// default grammar for constrained dp to triplet scfg
  /*
   * 
   */
  static Triplet_SCFG init_tkfst_scfg (double t, double u, double v);

  /// The scfg used in testtripletscfgdp.cc.
  /*
   * Included for testing and debugging purposes.
   */
  static Triplet_SCFG init_test_scfg();

};

#endif /* IGDEFAULTS_INCLUDED */
