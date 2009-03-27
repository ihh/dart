#ifndef SUBSTITUTION_MATRIX_FACTORY_INCLUDED
#define SUBSTITUTION_MATRIX_FACTORY_INCLUDED

#include "util/array2d.h"
#include "seq/biosequence.h"

struct Substitution_matrix_factory
{
  virtual array2d<Prob> create_conditional_substitution_matrix (double time) = 0;            // entry (x,y) is P(y|x) = prob. FROM x TO y (conditioned on x, i.e. column sum (over y) is 1)
  virtual array2d<Prob> differentiate_conditional_substitution_matrix (double time) = 0;     // entry (x,y) is d/dt P(y|x)
  virtual vector<Prob>  create_prior() = 0;
  virtual const Alphabet& alphabet() = 0;

  virtual array2d<double> rate_matrix() const = 0;  // R(x,y) = rate from x to y

  array2d<Prob> create_joint_substitution_matrix (double time);                  // entry (x,y) is P(x) * P(y|x)
  array2d<Prob> differentiate_joint_substitution_matrix (double time);           // entry (x,y) is P(x) * d/dt P(y|x)
  array2d<Prob> create_substitution_odds_ratio_matrix (double time);             // entry (x,y) is P(y|x) / P(y)

  vector<Prob> create_state_vector (double time);  // entry [y] is P(x_t=y) = \sum_x P(y|x) P(x)

  virtual ~Substitution_matrix_factory() {}
};

struct Substitution_matrix_factory_with_expected : Substitution_matrix_factory
{
  virtual double expected_substitution_rate() const = 0;  // expected no. of substitutions per unit of time. Implemented in Diagonalised_matrix_factory and Irrev_diagonalised_matrix_factory
};
#endif
