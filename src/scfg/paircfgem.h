#ifndef PAIRCFG_TRAINER_INCLUDED
#define PAIRCFG_TRAINER_INCLUDED

#include "scfg/pairpcfg.h"
#include "scfg/paircfgdp.h"
#include "seq/stockholm.h"
#include "seq/dirichlet.h"

// class to do EM for a pair CFG on a training set of structurally annotated pairwise alignments
struct Pair_CFG_trainer
{
  // references
  const Stockholm_database& stock;
  Pair_PCFG& cfg;
  PScores& pscore;
  Dirichlet_prior prior;
  bool local;
  // constructors
  Pair_CFG_trainer (const Stockholm_database& stock, Trainable_PCFG& cfg, bool local = FALSE);
  Pair_CFG_trainer (const Stockholm_database& stock, Pair_PCFG& cfg, Dirichlet_prior& prior, PScores& pscore, bool local = FALSE);
  // EM methods
  void train_pairwise (double min_inc = .001, int max_rounds = -1, int forgive = 0);  // requires alignments to be pairwise
  void train_single (double min_inc = .001, int max_rounds = -1, int forgive = 0);  // requires alignments to be single-sequence
};

#endif /* PAIRCFG_TRAINER_INCLUDED */
