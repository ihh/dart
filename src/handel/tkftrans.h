#ifndef TKF_TRANSDUCER_INCLUDED
#define TKF_TRANSDUCER_INCLUDED

#include "handel/alitrans.h"

struct TKF91_transducer_factory_param_container
{
  double lambda, mu;  // insertion & deletion rates
};

// TKF91_transducer_factory
// Implements a transducer for the TKF91 model
struct TKF91_transducer_factory : TKF91_transducer_factory_param_container, Transducer_alignment_with_subst_model
{
  // prior
  // global prior on delete_rate
  static double alpha, beta;  // hyperparameters of beta prior on delete_rate

  // global prior on 'gamma', the equilibrium sequence extension probability (probably better renamed 'eqm_extend')
  static double gamma_alpha, gamma_beta;  // hyperparameters of beta prior on 'gamma'

  // constructor
  TKF91_transducer_factory (double lambda = 0., double mu = 1.);

  // Handel_base virtual methods
  // gamma method
  double gamma() const { return lambda / mu; }

  void sample_indel_params();
  TKF91_transducer_factory_param_container propose_indel_params();

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
