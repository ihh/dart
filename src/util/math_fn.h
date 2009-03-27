#ifndef MATH_FN_INCLUDED
#define MATH_FN_INCLUDED

#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <complex>

#include "util/newmat_adaptors.h"

using namespace std;

// complex numbers
typedef std::complex<double> Complex;

// Much of this code is shamelessly ripped off from HMMER

struct Math_fn
{
  // gamma function
  static double log_gamma (double x);

  // combinatorial functions
  static unsigned int factorial (unsigned int x);
  static double log_factorial (unsigned int x);
  static unsigned int nCk (unsigned int n, unsigned int k);

  // Dirichlet distribution
  static double log_dirichlet (const vector<double>& p, const vector<double>& alpha);
  static double log_dirichlet_normaliser (const vector<double>& alpha);

  // log, exp, pow wrappers
  static inline double math_log (double a) { return log(a); }
  static inline double math_exp (double a) { return exp(a); }
  static inline double math_pow (double a, double b) { return pow(a,b); }

  // helpful stuff
  // p_extend() calculates a loop probability from a waiting time
  // NB total length = mean_extend_len + 1
  static double p_extend (double mean_extend_len);
};

// Numerical functions
//

template<class N> inline N sgn (N n) { return n > 0 ? 1 : (n < 0 ? -1 : 0); }
template<class N> inline N abs (N n) { return n > 0 ? n : -n; }
template<class N> inline N minmax (N n, N min_bound, N max_bound) { return n < min_bound ? min_bound : (n > max_bound ? max_bound : n); }
template<class N> inline N acos(N x) { return atan2(sqrt(1.-x*x),x); }
template<class N> inline N asin(N x) { return atan2(x,sqrt(1.-x*x)); }

#define Pi 3.141592654
#define PI Pi

// Schwartzian transform sort functor
// Sorts by increasing cost
template<class N>
struct Schwartzian : binary_function <int, int, bool>
{
  const vector<N>& cost;
  Schwartzian (const vector<N>& cost) : cost (cost) { }
  inline bool operator() (int a, int b) { return cost[a] < cost[b]; }
  vector<int> sorted()
  {
    const int sz = cost.size();
    vector<int> n (sz);
    for (int i = 0; i < sz; ++i) n[i] = i;
    sort (n.begin(), n.end(), *this);
    return n;
  }
};

// Reverse Schwartzian transform sort functor
// Sorts by decreasing value
template<class N>
struct Reverse_Schwartzian : binary_function <int, int, bool>
{
  const vector<N>& value;
  Reverse_Schwartzian (const vector<N>& value) : value (value) { }
  inline bool operator() (int a, int b) { return value[a] > value[b]; }
  vector<int> sorted()
  {
    const int sz = value.size();
    vector<int> n (sz);
    for (int i = 0; i < sz; ++i) n[i] = i;
    sort (n.begin(), n.end(), *this);
    return n;
  }
};

// The value of TINY is used by various numerical routines.
// It should be larger than the rounding error of the machine.
//

#define TINY 1.0e-20

#endif
