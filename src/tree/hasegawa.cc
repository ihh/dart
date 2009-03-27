#include <math.h>
#include <numeric>
#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "util/newmat_adaptors.h"
#include "tree/hasegawa.h"
#include "util/logfile.h"

Hasegawa_likelihood_matrix::Hasegawa_likelihood_matrix (double transition_rate,
							double transversion_rate,
							double time,
							const vector<double>& base_frequency,
							bool differentiate,
							bool exponentiate,
							bool do_ivratio_correction)
  : array2d<double> (4, 4, 0.)
{
  double norm = accumulate (base_frequency.begin(), base_frequency.end(), 0.0);

  // 'transition_rate' parameter corresponds to the expected total rate of i's
  // 'transversion_rate' parameter corresponds to half the expected total rate of v's
  // ... because this is intuitive for the caller
  // but we need the actual rate for each base divided by the base frequency
  // so we have to correct for this
  //
  // Oct 26, 2000:
  // This correction has caused problems in situations where the base frequencies preclude
  // a particular kind of mutation.
  // e.g. when f[a]=f[g]=0.5 and f[c]=f[t]=0 then transversions never happen and so vmul is infinite.
  //
  // To protect against this, and because the correction is itself non-biological in cases
  // where base frequency is an additional level of selection (such as synonymous substitution with codon bias),
  // I've made this correction optional.
  //
  if (do_ivratio_correction)
    {
      double imul = 0, vmul = 0;
      for (int i=0; i<4; i++)
	{
	  imul += base_frequency[i] * base_frequency[i^2] / (norm * norm);
	  vmul += base_frequency[i] * (base_frequency[i^1] + base_frequency[i^3]) / (norm * norm);
	}
      transition_rate = transition_rate / imul;
      transversion_rate = (transversion_rate * 2) / vmul;
    }

  const double transitions = transition_rate * time;
  const double transversions = transversion_rate * time;

  vector<int> base_type (4);          // base_type[] == 0 for purine, 1 for pyrimidine
  base_type[0] = base_type[2] = 0;
  base_type[1] = base_type[3] = 1;

  vector<double> type_frequency (2, (double) 0);  // purine & pyrimidine frequencies
  for (int i = 0; i < 4; i++)
    type_frequency [base_type[i]] += base_frequency[i] / norm;

  const double transversion_prob = exp (-transversions);
  const double transversion_deriv = -transversion_rate * transversion_prob;

  for (int i = 0; i < 4; i++)
    {
      const int bti = base_type[i];
      const double tfi = type_frequency[bti];
      const double tfo = type_frequency[1 - bti];

      const double transition_prob = exp (- tfo*transversions - tfi*transitions);
      const double transition_deriv = (- tfo*transversion_rate - tfi*transition_rate) * transition_prob;

      const double bfo = base_frequency [i^2] / norm;

      for (int j = 0; j < 4; j++)
	{
	  const double bfj = base_frequency[j] / norm;

	  if (differentiate)
	    {
	      if (i == j)
		(*this)(i,j) = (bfj*tfo/tfi)*transversion_deriv + (bfo/tfi)*transition_deriv;
	      
	      else if (bti == base_type[j])
		(*this)(i,j) = (bfj*tfo/tfi)*transversion_deriv - (bfo/tfi)*transition_deriv;
	      
	      else
		(*this)(i,j) = - bfj*transversion_deriv;
	    }
	  else if (exponentiate)
	    {
	      if (i == j)
		(*this)(i,j) = bfj + (bfj*tfo/tfi)*transversion_prob + (bfo/tfi)*transition_prob;
	      
	      else if (bti == base_type[j])
		(*this)(i,j) = bfj + (bfj*tfo/tfi)*transversion_prob - (bfo/tfi)*transition_prob;
	      
	      else
		(*this)(i,j) = bfj - bfj*transversion_prob;
	    }
	  else
	    {
	      if (i != j)
		{
		  if (bti == base_type[j])
		    (*this)(i,j) = transversion_rate * bfj;
		  else
		    (*this)(i,j) = transition_rate * bfj;
		  (*this) (i,i) -= (*this) (i,j);
		}
	    }
	}
    }
}

array2d<double> Hasegawa_matrix_factory::rate_matrix() const
{
  return Hasegawa_likelihood_matrix (transition_rate, transversion_rate, 0., base_frequency, false, false, do_ivratio_correction);
}

array2d<double> Hasegawa_matrix_factory::create_conditional_substitution_matrix (double time)
{
  array2d<double> sub = Hasegawa_likelihood_matrix (transition_rate, transversion_rate, time, base_frequency, false, true, do_ivratio_correction);
  if (CLOGGING(-1)) clog_stream << "Hasegawa substitution matrix at t=" << time << ":\n" << sub;
  if (CLOGGING(-2))
    {
      array2d_to_Matrix_adaptor submat (sub);
      vector_to_RowVector_adaptor p (base_frequency);
      clog_stream << "Prior: " << p;
      clog_stream << "Hasegawa substitution matrix * prior: " << p * submat.t();
    }
  return sub;
}


array2d<double> Hasegawa_matrix_factory::differentiate_conditional_substitution_matrix (double time)
{
  array2d<double> sub = Hasegawa_likelihood_matrix (transition_rate, transversion_rate, time, base_frequency, true, true, do_ivratio_correction);
  if (CLOGGING(-1)) clog_stream << "Time derivative of Hasegawa substitution matrix at t=" << time << ":\n" << sub;
  return sub;
}


