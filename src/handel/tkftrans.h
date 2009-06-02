#ifndef TKF_TRANSDUCER_INCLUDED
#define TKF_TRANSDUCER_INCLUDED

#include "handel/alitrans.h"

// TKF91_transducer_factory
// Implements a transducer for the TKF91 model
struct TKF91_transducer_factory : Transducer_alignment_with_subst_model
{
  // data
  double lambda, mu;  // insertion & deletion rates

  // constructor
  TKF91_transducer_factory (double lambda = 0., double mu = 1.);


  // Handel_base virtual methods

  // gamma method
  double gamma() const { return lambda / mu; }

  // clone method
  TKF91_transducer_factory* clone();


  // Transducer_alignment virtual methods

  // virtual factory methods
  Pair_transducer_scores prior_pair_trans_sc();
  Pair_transducer_scores branch_pair_trans_sc (double time);

  // virtual banding methods
  double gap_rate();
  double mean_gap_size();
};

#endif /* TKF_TRANSDUCER_INCLUDED */
