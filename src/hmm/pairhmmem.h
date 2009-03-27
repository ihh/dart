#ifndef PAIRHMM_TRAINER_INCLUDED
#define PAIRHMM_TRAINER_INCLUDED

#include "hmm/pairphmm.h"
#include "seq/stockholm.h"
#include "seq/dirichlet.h"

// class to do EM for a pair HMM on a training set of structurally annotated pairwise alignments
struct Pair_HMM_trainer
{
  // references
  const Stockholm_database& stock;
  Pair_PHMM& hmm;
  Dirichlet_prior& prior;
  PScores& pscore;
  // flags that can be changed after construction
  bool symmetrise;  // FALSE by default
  // constructors
  Pair_HMM_trainer (const Stockholm_database& stock, Pair_PHMM& hmm, Dirichlet_prior& prior, PScores& pscore);
  // EM methods
  void train (double min_inc = .001, int max_rounds = -1, int forgive = 0);  // requires alignments to be pairwise
};

#endif /* PAIRHMM_TRAINER_INCLUDED */
