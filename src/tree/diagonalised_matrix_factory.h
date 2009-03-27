#ifndef DIAGONALISED_MATRIX_FACTORY_INCLUDED
#define DIAGONALISED_MATRIX_FACTORY_INCLUDED

#include "tree/substitution_matrix_factory.h"

struct Diagonalised_matrix_factory : Substitution_matrix_factory_with_expected
{
  virtual const Alphabet& alphabet() = 0;          // inherited from Substitution_matrix_factory

  // min_prob is the minimum allowable probability for times > 0
  double min_prob;

  // NB the following are the eigenvectors/values for the underlying symmetric matrix, S_ij = sqrt(pi_i/pi_j) R_ij

  vector<double> eigenvalue;
  vector<vector<double> > eigenvector;  // eigenvector[i] == i'th eigenvector

  void truncate_eigenvalues();  // rounds down positive eigenvalues, with a warning
  int zeroth_eigenvalue_index() const;

  array2d<double> rate_matrix() const;  // R(x,y) = rate from x to y
  vector<double> prior() const;
  double expected_substitution_rate() const;  // expected no. of substitutions per unit of time
  void normalise_substitution_rate();  // divide by expected subst.rate

  // concrete methods
  
  array2d<Prob> create_conditional_substitution_matrix (double time);         // entry (x,y) is prob. FROM x TO y
  array2d<Prob> differentiate_conditional_substitution_matrix (double time);  // entry (x,y) is derivative w.r.t. time of prob. from x to y
  vector<Prob>  create_prior();

  Diagonalised_matrix_factory();
};

#endif
