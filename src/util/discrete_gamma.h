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
  const double alpha, beta;
  const int K;
  const bool median;
  Discrete_gamma (double alpha, double beta, int K, bool median = true);
};

#endif /* DISCRETE_GAMMA_INCLUDED */
