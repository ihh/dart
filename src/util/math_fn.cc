#include "util/math_fn.h"

// Much of this code is shamelessly ripped off from HMMER

/* Function: log_gamma()
 *
 * Returns the natural log of the gamma function of x.
 * x is > 0.0.  
 *
 * Adapted from a public domain implementation in the
 * NCBI core math library. Thanks to John Spouge and
 * the NCBI. (According to the NCBI, that's Dr. John
 * "Gammas Galore" Spouge to you, pal.)
 */
double Math_fn::log_gamma(double x)
{
  int i;
  double xx, tx;
  double tmp, value;
  static double cof[11] = {
    4.694580336184385e+04,
    -1.560605207784446e+05,
    2.065049568014106e+05,
    -1.388934775095388e+05,
    5.031796415085709e+04,
    -9.601592329182778e+03,
    8.785855930895250e+02,
    -3.155153906098611e+01,
    2.908143421162229e-01,
    -2.319827630494973e-04,
    1.251639670050933e-10
  };
  
  /* Protect against x=0. We see this in Dirichlet code,
   * for terms alpha = 0. This is a severe hack but it is effective
   * and safe. (due to GJM)
   */ 
  if (x <= 0.0) return 999999.; 

  xx       = x - 1.0;
  tx = tmp = xx + 11.0;
  value    = 1.0;
  for (i = 10; i >= 0; --i)	/* sum least significant terms first */
    {
      value += cof[i] / tmp;
      tmp   -= 1.0;
    }
  value  = log(value);
  tx    += 0.5;
  value += 0.918938533 + (xx+0.5)*log(tx) - tx;
  return (double) value;
}


/* Function: factorial()
 * 
 * Purpose:  Calculate the factorial of an integer
 *           
 * Return:   x!
 */          
unsigned int Math_fn::factorial (unsigned int x)
{
  unsigned int f = 1;
  while (x > 0) f *= x--;
  return f;
}

/* Function: log_factorial()
 * 
 * Purpose:  Calculate the log-factorial of an integer
 *           
 * Return:   log (x!)
 */          
double Math_fn::log_factorial (unsigned int x)
{
  double logf = 0;
  while (x > 1) logf += log ((double) (x--));
  return logf;
}

/* Function: nCk()
 * 
 * Purpose:  Calculate the combinatorial coefficient for two ints
 *           
 * Return:   n!/(k!(n-k)!)
 */          
unsigned int Math_fn::nCk (unsigned int n, unsigned int k)
{
  if (k > n) return 0;
  return factorial(n) / (factorial(k) * factorial(n-k));
}

/* Function: log_dirichlet()
 * 
 * Purpose:  Calculate the log probability of a probability
 *           vector given a single Dirichlet component, alpha.
 *           Follows Sjolander (1996) appendix, lemma 2.
 *           
 * Return:   log P(p | alpha)
 */          
double Math_fn::log_dirichlet (const vector<double>& p, const vector<double>& alpha)
{
  double sum;		        /* for log_gamma(|alpha|) in Z  */
  double logp;			/* RETURN: log P(p|alpha)       */
  int x;

  sum = logp = 0.0;
  for (x = 0; x < (int) p.size(); x++)
    if (p[x] > 0.0)		/* any param that is == 0.0 doesn't exist */
      {
	logp += (alpha[x]-1.0) * log(p[x]);
	logp -= log_gamma(alpha[x]);
	sum  += alpha[x];
      }
  logp += log_gamma (sum);
  return logp;
}


/* Function: log_dirichlet_normaliser()
 * 
 * Purpose:  Calculate the normaliser for a single Dirichlet component, alpha.
 *           Follows Durbin et al (1998), section 11.1
 *           
 * Return:   log Z(alpha)
 */          
double Math_fn::log_dirichlet_normaliser (const vector<double>& alpha)
{
  double sum;		        /* for log_gamma(|alpha|) in denominator  */
  double logz;			/* RETURN: log Z(alpha)       */
  int x;

  sum = logz = 0.0;
  for (x = 0; x < (int) alpha.size(); x++)
    {
      logz += log_gamma(alpha[x]);
      sum  += alpha[x];
    }
  logz -= log_gamma (sum);
  return logz;
}

// p_extend(len) calculates a loop probability from a waiting time

double Math_fn::p_extend (double mean_extend_len)
{
  return 1.0 - 1.0 / (1.0 + mean_extend_len);
}
