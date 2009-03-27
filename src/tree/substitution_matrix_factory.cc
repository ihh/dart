#include "tree/substitution_matrix_factory.h"

array2d<double> Substitution_matrix_factory::create_joint_substitution_matrix (double time)
{
  array2d<double> a = create_conditional_substitution_matrix (time);
  vector<double> v = create_prior();

  for (int i = 0; i < a.xsize(); i++)
    for (int j = 0; j < a.ysize(); j++)
      a(i,j) *= v[i];

  if (CTAGGING(-2,JOINTSUBMAT))
    CL << "Joint substitution matrix at t=" << time << ":\n" << a;

  return a;
}

array2d<double> Substitution_matrix_factory::differentiate_joint_substitution_matrix (double time)
{
  array2d<double> a = differentiate_conditional_substitution_matrix (time);
  vector<double> v = create_prior();
  
  for (int i = 0; i < a.xsize(); i++)
    for (int j = 0; j < a.ysize(); j++)
      a(i,j) *= v[i];

  return a;
}

array2d<double> Substitution_matrix_factory::create_substitution_odds_ratio_matrix (double time)
{
  array2d<double> a = create_conditional_substitution_matrix (time);
  vector<double> v = create_prior();

  for (int i = 0; i < a.xsize(); i++)
    for (int j = 0; j < a.ysize(); j++)
      a(i,j) /= v[j];

  return a;
}

vector<Prob> Substitution_matrix_factory::create_state_vector (double time)
{
  const array2d<double> a = create_conditional_substitution_matrix (time);
  const vector<double> v = create_prior();

  vector<Prob> p (a.xsize(), 0.);
  for (int i = 0; i < a.xsize(); ++i)
    for (int j = 0; j < a.ysize(); ++j)
      p[j] += v[i] * a(i,j);

  return p;
}
