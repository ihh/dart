#ifndef DISCRETE_GAMMA_INCLUDED
#define DISCRETE_GAMMA_INCLUDED

#include <vector>

using namespace std;

// A vector of equal-probability bins for the gamma distribution
// alpha: shape parameter
// beta: scale parameter
// K: number of rate categories
// median: a flag indicating how to choose the rates
//   if true, use median rate for each bin, rescaled so mean rate taken over all bins is 1.0
//   if false, use mean rate for each bin
struct Discrete_gamma : vector<double> {
  // data
  const double alpha, beta;
  const int K;
  const bool median;
  // constructor
  Discrete_gamma (double alpha, double beta, int K, bool median = true);
};

// mathematical functions, taken from MrBayes
#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta))

double  PointChi2 (double prob, double v);
double  LnGamma (double alp);
double  PointNormal (double prob);
double  IncompleteGamma (double x, double alpha, double LnGamma_alpha);
double  GammaDensity (double x, double alpha, double beta);

#endif /* DISCRETE_GAMMA_INCLUDED */
