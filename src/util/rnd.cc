#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#include "util/rnd.h"
#include "randlib/randlib.h"

#define DEFAULT_RND_SEED 1234567890

int  Rnd::seed_val = DEFAULT_RND_SEED;
bool Rnd::seeded   = 0;

void Rnd::seed()
{
  if (!seeded)
    {
      srand48 (seed_val);
      seeded = TRUE;
      log_seed (0);
    }
}

void Rnd::set_seed (int val)
{
  seed_val = val;
  seeded = FALSE;
}

void Rnd::add_opts (Opts_list& ol)
{
  ol.add ("rndseed", Rnd::seed_val, "seed random number generator");
  ol.add ("rndtime", Rnd::seed_on_time, "\tseed random number generator on current clock time");
  ol.post_parse_callback.push_back (&force_seed);
}

bool Rnd::seed_on_time (Opts_list* ol)
{
  struct timeval tv;
  struct timezone tz;
  struct tm *tm;
  gettimeofday(&tv, &tz);
  set_seed (((long) getpid()) ^ ((long) tv.tv_usec));
  return TRUE;
}

bool Rnd::force_seed (Opts_list* ol)
{
  if (seeded)
    CTAG(8,RND) << "Warning -- used random number generator during initialisation\n";
  seeded = FALSE;
  return 1;
}

bool Rnd::log_seed (Opts_list* ol)
{
  CTAG(8,RND) << "Random number seed = " << Rnd::seed_val << "\n";
  return 1;
}

vector<double> Rnd::sample_dirichlet (const vector<double>& counts)    // not totally convinced this is working correctly
{
  vector<double> result (counts.size());
  double norm = 0;
  for (int i = 0; i < (int) counts.size(); ++i)
    norm += (result[i] = sample_gamma (counts[i], 1.0));
  if (norm > 0)
    for (int i = 0; i < (int) counts.size(); ++i)
      result[i] = result[i] / norm;
  return result;
}

vector<double> Rnd::dirichlet_max (const vector<double>& counts)
{
  vector<double> result = counts;
  double norm = accumulate (counts.begin(), counts.end(), 0.0);
  if (norm > 0)
    for (int i = 0; i < (int) counts.size(); ++i)
      result[i] = result[i] / norm;
  return result;
}

double Rnd::sample_beta (const double p0_count, const double p1_count)
{
  vector<double> counts (2);
  counts[0] = p0_count;
  counts[1] = p1_count;

  const vector<double> p = sample_dirichlet (counts);
  return p[0];
}

double Rnd::sample_gamma (const double alpha, const double beta)
{
  return gengam (alpha, beta);
}
