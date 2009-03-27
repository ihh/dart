#include "tree/irrev_diag_matrix_factory.h"
#include "util/dexception.h"
#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "util/logfile.h"
#include "util/newmat_adaptors.h"

Irrev_diagonalised_matrix_factory::Irrev_diagonalised_matrix_factory() : min_prob (1e-10) { }

void Irrev_diagonalised_matrix_factory::init_matrix_factory (int states)
{
  mu = vector<Complex> (states);
  U = array2d<Complex> (states, states);
  Uinv = array2d<Complex> (states, states);
  irrev_prior = vector<double> (states);
}

array2d<double> Irrev_diagonalised_matrix_factory::create_conditional_substitution_matrix (double time)
{
  truncate_eigenvalues();
  const int s = mu.size();
  vector<Complex> exp_eval (s);
  for (int i = 0; i < s; ++i)
    exp_eval[i] = exp (mu[i] * time);

  array2d<double> M (s, s, 0.);
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      for (int k = 0; k < s; ++k)
	M(i,j) += std::real (U(i,k) * exp_eval[k] * Uinv(k,j));

  int cutoffs = 0;
  double t_min_prob = time == 0 ? 0.0 : min_prob;
  double lowest_neg_prob = 0.;
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      if (M(i,j) < min_prob)
	{
	  if (M(i,j) < lowest_neg_prob)
	    lowest_neg_prob = M(i,j);
	  M(i,j) = min_prob;
	  ++cutoffs;
	}
  if (cutoffs > 0)
    CLOG(3) << "Warning: " << cutoffs << " probabilities rounded up to " << t_min_prob << " at time " << time
	    << " (lowest was " << lowest_neg_prob << ")\n";

  if (CTAGGING(-1,COND_SUBMAT))
    CL << "Conditional substitution matrix at t=" << time << ":\n" << M;
  if (CLOGGING(-2))
    {
      array2d_to_Matrix_adaptor submat = M;
      RowVector pi(s);
      for (int i = 0; i < s; ++i)
        pi(i+1) = irrev_prior[i];
      CL << "Pi: " << pi << "Pi * substitution matrix: " << pi * submat.t();
    }

  return M;
}

array2d<double> Irrev_diagonalised_matrix_factory::differentiate_conditional_substitution_matrix (double time)
{
  truncate_eigenvalues();

  const int s = mu.size();
  vector<Complex> exp_eval_deriv (s);
  for (int i = 0; i < s; ++i)
    exp_eval_deriv[i] = mu[i] * exp (mu[i] * time);

  array2d<double> deriv (s, s, 0.0);
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      for (int k = 0; k < s; ++k)
	deriv(i,j) += std::real (U(i,k) * exp_eval_deriv[k] * Uinv(k,j));

  if (CLOGGING(-1))
    CL << "Time derivative of substitution matrix at t=" << time << ":\n" << deriv;

  return deriv;
}

vector<double> Irrev_diagonalised_matrix_factory::create_prior() { return prior(); }

vector<double> Irrev_diagonalised_matrix_factory::prior() const
{
  return irrev_prior;
}

void Irrev_diagonalised_matrix_factory::truncate_eigenvalues()
{
  const int s = mu.size();
  for (int i = 0; i < s; i++)
	if (std::real(mu[i]) > 0)
	{
	  CLOGERR << "Warning: zeroing rate eigenvalue with positive real component (was " << mu[i] << ")\n";
	  mu[i] = 0;
	}
}

array2d<double> Irrev_diagonalised_matrix_factory::rate_matrix() const
{
  const int s = mu.size();
  array2d<double> R (s, s, 0.0);
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      if (j != i)
	{
	  for (int k = 0; k < s; ++k)
	    R(i,j) += std::real (U(i,k) * mu[k] * Uinv(k,j));
	  if (R(i,j) < 0)
	    R(i,j) = 0;
	  R(i,i) -= R(i,j);
	}
  return R;
}

double Irrev_diagonalised_matrix_factory::expected_substitution_rate() const
{
  const int s = mu.size();
  const array2d<double> R = rate_matrix();
  double srate = 0;
  for (int i = 0; i < s; ++i)
    srate -= irrev_prior[i] * R(i,i);
  return srate;
}
