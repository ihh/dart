#include "tree/diagonalised_matrix_factory.h"
#include "util/dexception.h"
#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "util/logfile.h"
#include "util/newmat_adaptors.h"

Diagonalised_matrix_factory::Diagonalised_matrix_factory() : min_prob (1e-10) { }

array2d<double> Diagonalised_matrix_factory::create_conditional_substitution_matrix (double time)
{
  truncate_eigenvalues();

  const int s = eigenvalue.size();
  vector<double> exp_eval (s);
  for (int i = 0; i < s; ++i) exp_eval[i] = exp (eigenvalue[i] * time);

  const vector<double>& sqrt_prior = eigenvector[zeroth_eigenvalue_index()];

  array2d<double> a (s, s, 0.0);
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      for (int k = 0; k < s; ++k)
	a(i,k) += eigenvector[j][i] * exp_eval[j] * eigenvector[j][k] * sqrt_prior[k] / sqrt_prior[i];

  int cutoffs = 0;
  double t_min_prob = time == 0 ? 0.0 : min_prob;
  for (int i = 0; i < s; ++i)
    for (int k = 0; k < s; ++k)
      if (a(i,k) < min_prob)
	{
	  a(i,k) = min_prob;
	  ++cutoffs;
	}
  if (cutoffs > 0)
    CLOG(4) << "Warning: " << cutoffs << " probabilities rounded up to " << t_min_prob << " at time " << time << "\n";

  if (CTAGGING(-1,CONDSUBMAT))
    CL << "Conditional substitution matrix at t=" << time << ":\n" << a;
  if (CLOGGING(-2))
    {
      array2d_to_Matrix_adaptor submat = a;
      vector_to_RowVector_adaptor p = create_prior();
      CL << "Prior: " << p << "Substitution matrix * prior: " << p * submat.t();
    }

  return a;
}

array2d<double> Diagonalised_matrix_factory::differentiate_conditional_substitution_matrix (double time)
{
  truncate_eigenvalues();

  const int s = eigenvalue.size();
  vector<double> exp_eval_deriv (s);
  for (int i = 0; i < s; ++i) exp_eval_deriv[i] = eigenvalue[i] * exp (eigenvalue[i] * time);

  const vector<double>& sqrt_prior = eigenvector[zeroth_eigenvalue_index()];

  array2d<double> a (s, s, 0.0);
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      for (int k = 0; k < s; ++k)
	a(i,k) += eigenvector[j][i] * exp_eval_deriv[j] * eigenvector[j][k] * sqrt_prior[k] / sqrt_prior[i];

  if (CLOGGING(-1)) clog_stream << "Time derivative of substitution matrix at t=" << time << ":\n" << a;

  return a;
}

vector<double> Diagonalised_matrix_factory::create_prior() { return prior(); }

vector<double> Diagonalised_matrix_factory::prior() const
{
  const int z = zeroth_eigenvalue_index();
  vector<double> p (eigenvalue.size());
  for (int j = 0; j < (int) eigenvalue.size(); ++j) p[j] = eigenvector[z][j] * eigenvector[z][j];
  return p;
}

void Diagonalised_matrix_factory::truncate_eigenvalues()
{
  for_contents (vector<double>, eigenvalue, eval)
    if (*eval > 0)
      {
	CLOGERR << "Warning: zeroing rate eigenvalue (was " << *eval << ")\n";
	*eval = 0;
      }
}

int Diagonalised_matrix_factory::zeroth_eigenvalue_index() const
{
  return max_element (eigenvalue.begin(), eigenvalue.end()) - eigenvalue.begin();
}

array2d<double> Diagonalised_matrix_factory::rate_matrix() const
{
  const int s = eigenvalue.size();
  const vector<double>& sqrt_prior = eigenvector[zeroth_eigenvalue_index()];
  array2d<double> R (s, s, 0.0);
  for (int i = 0; i < s; ++i)
    for (int j = 0; j < s; ++j)
      for (int k = 0; k < s; ++k)
	R(i,k) += eigenvector[j][i] * eigenvalue[j] * eigenvector[j][k] * sqrt_prior[k] / sqrt_prior[i];
  return R;
}

double Diagonalised_matrix_factory::expected_substitution_rate() const
{
  const int s = eigenvalue.size();
  const array2d<double> R = rate_matrix();
  const vector<double> pi = prior();
  double srate = 0;
  for (int i = 0; i < s; ++i)
    srate -= pi[i] * R(i,i);
  return srate;
}
