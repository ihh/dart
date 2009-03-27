#ifndef HASEGAWA_INCLUDED
#define HASEGAWA_INCLUDED

#include "tree/substitution_matrix_factory.h"

struct Hasegawa_likelihood_matrix : array2d<double>
{
  Hasegawa_likelihood_matrix (double transition_rate,
			      double transversion_rate,
			      double time,
			      const vector<double>& base_frequency,
			      bool differentiate = 0,
			      bool exponentiate = 1,
			      bool do_ivratio_correction = 1);
};

struct Hasegawa_matrix_factory : Substitution_matrix_factory
{
  double          transition_rate;
  double          transversion_rate;
  vector<double>  base_frequency;
  bool            do_ivratio_correction;

  array2d<double> create_conditional_substitution_matrix (double time);
  array2d<double> differentiate_conditional_substitution_matrix (double time);
  vector<double>  create_prior() { return base_frequency; }

  array2d<double> rate_matrix() const;

  const Alphabet& alphabet() { return DNA_alphabet; }

  Hasegawa_matrix_factory() : transition_rate(0), transversion_rate(0), base_frequency (4, 0.0), do_ivratio_correction(1) { }

  Hasegawa_matrix_factory (double transition_rate,
			   double transversion_rate,
			   const vector<double>& base_frequency,
			   bool do_ivratio_correction = 1)
    : transition_rate(transition_rate),
      transversion_rate(transversion_rate),
      base_frequency(base_frequency),
      do_ivratio_correction(do_ivratio_correction) {}
};

#endif
