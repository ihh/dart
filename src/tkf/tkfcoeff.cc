#include <math.h>
#include "tkf/tkfcoeff.h"

/* Function definitions */

TKF_coeff::TKF_coeff (double dt, double T) : dt(dt), T(T)
{
  N = 1 + (int) (T / dt);
  if (dt * (double) (N-1) < T) ++N;

  T = dt * (double) (N-1);

  alpha = vector<double> (N);
  beta = vector<double> (N);
  gamma = vector<double> (N);
  lambdaD = vector<double> (N);
  lambdaI = vector<double> (N);
  lambdaE = vector<double> (N);
  I0 = vector<double> (N);
  I1 = vector<double> (N);
  I01 = vector<double> (N);
  I11 = vector<double> (N);
  Z0 = vector<double> (N);
  Z1 = vector<double> (N);
  Z01 = vector<double> (N);
  Z11 = vector<double> (N);
}

void TKF_coeff::initialise (const TKF_params& params)
{
  CLOG(6) << "Initialising TKF EM coefficients: dt= " << dt << " T= " << T << " lambda= " << params.lambda << " mu= " << params.mu << "\n";

  /* tmp params */
  double lambda_dt, mu_dt;

  /* tmp variables for outer loop */
  int n;
  double tmp_T;
  double x;  /* exp((lambda-mu) * tmp_T) */
  double y;  /* exp(-mu * tmp_T) */
  double z;  /* (1-x) / (mu - lambda * x) */
  double tmp_alpha;
  double noalpha_gamma;  /* (1-alpha) * gamma */
  double noalpha_nogamma;  /* (1-alpha) * (1-gamma) */
  
  /* tmp variables for inner loop */
  int m;
  double t;
  double v;
  double a, b, c;
  double tmp_lambdaD, tmp_lambdaI, tmp_lambdaE;
  double tmp_I0, tmp_I1, tmp_I01;
  double tmp_Z0, tmp_Z1, tmp_Z01;

  /* main function body */

  lambda_dt = params.lambda * dt;
  mu_dt = params.mu * dt;

  g = params.lambda / params.mu;

  /* outer loop: fills arrays */
  for (n = 0; n < N; ++n) {
    /* calculate tmp_T */
    tmp_T = dt * (double) n;

    /* calculate alpha, beta, gamma */
    x = exp ((params.lambda - params.mu) * tmp_T);
    y = exp (-params.mu * tmp_T);
    z = (1-x) / (params.mu - params.lambda * x);
    alpha[n] = y;
    beta[n] = params.lambda * z;
    if (n == 0)
      gamma[n] = 0;
    else
      gamma[n] = 1 - params.mu * z / (1-y);

    /* set up temp variables */
    tmp_alpha = alpha[n];
    noalpha_gamma = (1 - alpha[n]) * gamma[n];
    noalpha_nogamma = (1 - alpha[n]) * (1 - gamma[n]);

    /* calculate effective insertion rates */
    lambdaD[n] = lambda_dt * tmp_alpha;
    lambdaI[n] = lambda_dt * noalpha_gamma;
    lambdaE[n] = lambda_dt * noalpha_nogamma;

    /* calculate expected insertions and site time by zone type */
    I0[n] = I1[n] = I01[n] = I11[n] = 0;
    Z0[n] = Z1[n] = Z01[n] = Z11[n] = 0;
    /* inner loop: does integrals */
    v = 1.;
    for (m = 0; m < n; ++m) {

      /* calculate t */
      t = dt * (double) m;

      /* make temporary copies of effective insertion rates */
      tmp_lambdaD = lambdaD[n-m];
      tmp_lambdaI = lambdaI[n-m];
      tmp_lambdaE = lambdaE[n-m];

      /* update v */
      v *= 1 - (mu_dt + tmp_lambdaD + tmp_lambdaI);

      /* calculate a, b, c */
      if (n < 2) {  /* protect against division by zero at n=1 */
	a = 1;
	b = 0;
      } else {
	a = v * (1 - alpha[n-m-1]) * (1 - gamma[n-m-1]) / noalpha_nogamma;
	b = alpha[m] * (1 - alpha[n-m-1]) * (1 - gamma[n-m-1]) / noalpha_gamma;
      }
      c = (1 - beta[n-m-1]) / beta[n];

      /* update expected insertions */
      tmp_I0 = I0[n-m-1];
      tmp_I1 = I1[n-m-1];
      tmp_I01 = I01[n-m-1];
      I0[n] += a * tmp_lambdaE * (1 + tmp_I0);
      I1[n] += tmp_lambdaE * (1 + tmp_I0);
      I01[n] += b * (tmp_lambdaD * (1 + tmp_I0) + tmp_lambdaI * (1 + tmp_I0 + tmp_I01));
      I11[n] += c * (tmp_lambdaD * (1 + tmp_I1) + tmp_lambdaI * (1 + tmp_I1 + tmp_I01));

      /* update expected site time */
      tmp_Z0 = Z0[n-m-1];
      tmp_Z1 = Z1[n-m-1];
      tmp_Z01 = Z01[n-m-1];
      Z0[n] += a * (dt + tmp_lambdaE * tmp_Z0);
      Z1[n] += dt + tmp_lambdaE * tmp_Z0;
      Z01[n] += b * (tmp_lambdaD * tmp_Z0 + tmp_lambdaI * (tmp_Z0 + tmp_Z01));
      Z11[n] += c * (tmp_lambdaD * tmp_Z1 + tmp_lambdaI * (tmp_Z1 + tmp_Z01));
    }
  }

  CLOG(6) << "Finished initialising TKF EM coefficients\n";
}

TKF_counts::TKF_counts()
{
  reset_counts();
}

void TKF_counts::reset_counts()
{
  T = Z = I = D = 0;
}

void TKF_coeff::accumulate_counts (const TKF_transition_counts& trans, double T, TKF_counts& counts)
{
  int n;
  double deltaZ, deltaD, deltaI;  /* increments for Z, D and I */
  double obsI, obsD;  /* observed insertions & deletions */

  /* get time index for TKF_coeff arrays */
  n = (int) (T / dt);
  if (n < 0 || n >= N)
    {
      fprintf (stderr, "Bad time value in accumulate_TKF_counts: %g is outside range [0,%g] given by TKF_coeff struct\n", T, T);
    }
  /* count observed indels */
  obsD = trans.n0 + trans.n01;
  obsI = trans.n01 + trans.n11;
  /* zero deltas */
  deltaI = deltaZ = deltaD = 0;
  /* get expected insertions delta */
  deltaI += trans.n0 * I0[n];
  deltaI += trans.n1 * I1[n];
  deltaI += trans.n01 * I01[n];
  deltaI += trans.n11 * I11[n];
  /* get expected deletions delta */
  deltaD = deltaI + obsD - obsI;
  /* get expected site-time delta */
  deltaZ += trans.n0 * Z0[n];
  deltaZ += trans.n1 * Z1[n];
  deltaZ += trans.n01 * Z01[n];
  deltaZ += trans.n11 * Z11[n];
  deltaZ -= T;  /* subtract site-time contribution of immortal link */
  /* add deltas */
  counts.T += T;
  counts.Z += deltaZ;
  counts.D += deltaD;
  counts.I += deltaI;
}

void TKF_coeff::accumulate_naive_counts (const TKF_transition_counts& trans, double T, TKF_counts& counts)
{
  double obsI, obsD;  /* observed insertions & deletions */
  double initLen, finalLen;  /* initial & final length */

  /* count observed indels & lengths */
  obsD = trans.n0 + trans.n01;
  obsI = trans.n01 + trans.n11;
  initLen = trans.n0 + trans.n1 - 1;
  finalLen = trans.n1 + trans.n11 - 1;
  /* add deltas */
  counts.T += T;
  counts.Z += (initLen + finalLen) * T / 2.;
  counts.D += obsD;
  counts.I += obsI;
}

void TKF_counts::set_params (TKF_params& params)
{
  params.lambda = I / (Z + T);
  params.mu = D / Z;
}

TKF_transition_counts::TKF_transition_counts() : n0(0), n1(0), n01(0), n11(0) { }

int TKF_transition_counts::TKF_pair_state (const Pairwise_path& path, int pos)
{
  if (path.parent()[pos])
    return path.child()[pos] ? MAT : DEL;
  if (!path.child()[pos])
    THROWEXPR ("All-gap column in TKF pairwise alignment -- illegal");
  return INS;
}

void TKF_transition_counts::add_path_counts (const Pairwise_path& path)
{
  if (path.columns() == 0)
    ++n1;
  else
    {
      if (TKF_pair_state (path, 0) == INS)
	++n11;
      else
	++n1;
      /* intra-model transitions */
      for (int n = 0; n < path.columns() - 1; ++n) {
	if (TKF_pair_state (path, n) == DEL) {
	  if (TKF_pair_state (path, n+1) == INS)
	    ++n01;
	  else
	    ++n0;
	} else if (TKF_pair_state (path, n+1) == INS)
	  ++n11;
	else
	  ++n1;
      }
      /* transition to end */
      if (TKF_pair_state (path, path.columns() - 1) == DEL)
	++n0;
      else
	++n1;
    }
}

void TKF_transition_counts::show (ostream& o) const
{
  o << "TKF transition counts: n0=" << n0 << " n1=" << n1 << " n01=" << n01 << " n11=" << n11 << "\n";
}

void TKF_counts::show (ostream& o) const
{
  o << "TKF EM counts: T=" << T << " Z=" << Z << " I=" << I << " D=" << D << "\n";
}
