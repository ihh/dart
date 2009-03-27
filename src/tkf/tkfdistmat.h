#ifndef TKF_DISTMAT_INCLUDED
#define TKF_DISTMAT_INCLUDED

#include "seq/distmat.h"
#include "tkf/tkfparams.h"

// concrete Dist_func_factory for TKF model
struct TKF_dist_func_factory : Dist_func_factory
{
  // Dist_func class
  struct TKF_dist_func : Dist_func
  {
    const TKF_dist_func_factory& factory;
    const Alignment& align;
    // constructor
    TKF_dist_func (const TKF_dist_func_factory& factory, const Alignment& align) : factory (factory), align (align) { }
    // method to evaluate pairwise distance
    double operator() (int i, int j);
  };
  // data
  const TKF_params& params;
  bool use_indels, use_subst;
  double tres, tmax;
  // Time-likelihood function f(t) depends on n_realign as follows:
  //              { < 0: f = TKF_unaligned_counts_function
  //    n_realign { = 0: f = TKF_aligned_counts_function, using supplied Alignment
  //              { > 0: f = TKF_aligned_counts_function, using ML pairwise alignment at default_time
  int n_realign;
  double default_time;
  // constructor
  TKF_dist_func_factory (const TKF_params& params,
			 bool use_indels = FALSE, bool use_subst = TRUE,
			 double tres = .01, double tmax = 10.,
			 int n_realign = 0, double default_time = 1.);
  // factory method
  TKF_dist_func* create_dist_func (const Alignment& align);
};

#endif /* TKF_DISTMAT_INCLUDED */
