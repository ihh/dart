#ifndef IR_DIAGONALISED_MATRIX_FACTORY_INCLUDED
#define IR_DIAGONALISED_MATRIX_FACTORY_INCLUDED

#include "tree/substitution_matrix_factory.h"
#include "newmat/newmat.h"

struct Irrev_diagonalised_matrix_factory : Substitution_matrix_factory_with_expected
{
  virtual const Alphabet& alphabet() = 0;          // inherited from Substitution_matrix_factory

  // min_prob is the minimum allowable probability for times > 0
  double min_prob;

  // the following are the eigenvalues & right/left eigenvectors for the rate matrix
  vector<Complex> mu;
  array2d<Complex> U;
  array2d<Complex> Uinv;

  // The irrev_prior vector describes the initial probability distribution over states.
  // It is an exact copy of EM_matrix_params::pi (which is wasteful...)
  vector<double> irrev_prior;

  void truncate_eigenvalues();  // rounds down positive eigenvalues, with a warning

  // method to construct rate matrix from eigenvalues
  virtual array2d<double> rate_matrix() const;  // R(x,y) = rate from x to y
  vector<double> prior() const;
  double expected_substitution_rate() const;  // expected no. of substitutions per unit of time

  // concrete methods
  
  array2d<Prob> create_conditional_substitution_matrix (double time);         // entry (x,y) is prob. FROM x TO y
  array2d<Prob> differentiate_conditional_substitution_matrix (double time);  // entry (x,y) is derivative w.r.t. time of prob. from x to y
  vector<Prob>  create_prior();

  // constructor, initialiser
  Irrev_diagonalised_matrix_factory();
  void init_matrix_factory (int states);
};

#endif
