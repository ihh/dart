#ifndef RND_INCLUDED
#define RND_INCLUDED

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <numeric>
#include "util/opts_list.h"
#include "util/logfile.h"

class Rnd
{
 private:
  static int    seed_val;
  static bool   seeded;

  static void   seed();

 public:
  // initialization & seeding methods
  static void   set_seed (int val);

  static bool   force_seed (Opts_list* ol);
  static bool   seed_on_time (Opts_list* ol);  // actually seeds on (microseconds since 1/1/1970) XOR (current PID)
  static bool   log_seed (Opts_list* ol);
  
  static void   add_opts (Opts_list& ol);

  // the following two methods are the only ones that actually generate random numbers
  static double prob() { Rnd::seed(); return drand48(); }
  static long   long_int() { Rnd::seed(); return lrand48(); }

  // the remaining methods are just helpers, that manipulate the numbers
  static bool   decide (double prob) { return Rnd::prob() < prob; }
  static int    choose (const vector<double>& weight_vec)
    {
      const double norm = accumulate (weight_vec.begin(), weight_vec.end(), 0.0);
      double prob = Rnd::prob() * norm;
      for (int i = 0; i < ((int) weight_vec.size()) - 1; i++) if ((prob -= weight_vec[i]) <= 0) return i;
      return weight_vec.size() - 1;
    }

  static int rnd_int (int n) { return min ((int) (prob() * ((double) n)), n-1); }      // returns integer in range 0 ... n-1
  static double rnd_double (double scale_factor);

  static double sample_gamma (const double alpha, const double beta);
  static double sample_beta (const double p0_count, const double p1_count);   // returns sampled p0

  static vector<double> sample_dirichlet (const vector<double>& counts);
  static vector<double> dirichlet_max (const vector<double>& counts);                  // this just returns the normalised counts vector
};

#endif
