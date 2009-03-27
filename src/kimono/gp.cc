#include <math.h>
#include <iomanip>

#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "kimono/gp.h"
#include "util/macros.h"
#include "util/newmat_adaptors.h"
#include "util/vector_output.h"
#include "randlib/randlib.h"

void Gaussian_cluster_params::read_ColumnVector (istream& in, ColumnVector& v, const Phonebook& index, int size)
{
  v.ReSize (size);
  
  sstring s;
  s.getline (in);
  vector<sstring> experiment_labels = s.split ("\t", 0);
  if ((int) experiment_labels.size() != size + 1) THROW Format_exception (in, "First row of vector has the wrong number of elements\n");

  vector<int> col_index (size, -1);
  for (int i = 0; i < size; ++i)
    {
      int table_index = index.lookup (experiment_labels[i+1]);
      if (table_index < 0)
	{
	  sstring e;
	  e << "Can't find a data point with label " << experiment_labels[i+1];
	  THROW Standard_exception (e);
	}
      else if (col_index[table_index] >= 0)
	{
	  sstring e;
	  e << "Multiple data points with label " << experiment_labels[i+1];
	  THROW Standard_exception (e);
	}
      col_index[i] = table_index + 1;
    }
  
  s.getline (in);
  vector<sstring> row_data = s.split ("\t", 0);
  if ((int) row_data.size() != size + 1)
    THROW Format_exception (in, "Second row of vector has the wrong number of elements\n");
  for (int col = 0; col < size; ++col)
    v (col_index[col]) = atof (row_data[col+1].c_str());
}

void Gaussian_cluster_params::read_SymmetricMatrix (istream& in, SymmetricMatrix& m, const Phonebook& index, int size)
{
  m.ReSize (size);
  
  sstring s;
  s.getline (in);
  s.chomp();
  vector<sstring> experiment_labels = s.split ("\t", 0);
  if ((int) experiment_labels.size() != size + 1)
    THROW Format_exception (in, "Covariance matrix is the wrong size\n");

  vector<int> col_index (size, -1);
  for (int i = 0; i < size; ++i)
    {
      int table_index = index.lookup (experiment_labels[i+1]);
      if (table_index < 0)
	{
	  sstring e;
	  e << "Can't find a data point with label " << experiment_labels[i+1];
	  THROW Standard_exception (e);
	}
      else if (col_index[table_index] >= 0)
	{
	  sstring e;
	  e << "Multiple data points with label " << experiment_labels[i+1];
	  THROW Standard_exception (e);
	}
      col_index[i] = table_index + 1;
    }
  
  for (int row = 0; row < size; ++row)
    {
      s.getline (in);
      s.chomp();
      vector<sstring> row_data = s.split ("\t", 0);
      if ((int) row_data.size() != row + 2)
	{
	  sstring e;
	  e << "Row #" << row+1 << " of covariance matrix is the wrong size\n";
	  THROW Format_exception (in, e);
	}
      if (row_data[0] != experiment_labels[row+1]) { sstring e; e << "Row #" << row+1 << " of covariance matrix has the wrong label\n"; THROW Format_exception (in, e); }
      for (int col = 0; col <= row; ++col)
	m (col_index[row], col_index[col]) = atof (row_data[col+1].c_str());
    }
}

void Gaussian_cluster_params::calculate_inverses()
{
  CLOG(6) << "Calculating inverses of covariance matrices\n";

  double dlog2Pi = ((double) table.experiments()) * (log(2.) + log(Pi));
  
  invC = C.i();
  LogAndSign det_invC = invC.LogDeterminant();
  log_invC_norm = 0.5 * (det_invC.LogValue() - dlog2Pi);

  invN = N.i();
  LogAndSign det_invN = invN.LogDeterminant();
  log_invN_norm  = 0.5 * (det_invN.LogValue() - dlog2Pi);

  if (log_invC_norm > log_invN_norm) CLOGERR << "Warning: intra-cluster variance-per-point is " << exp ((det_invC.LogValue() - det_invN.LogValue()) / table.experiments()) << " * inter-cluster variance-per-point\n";
}


void Gaussian_cluster_params::read (const sstring& Emu_file, const sstring& C_file, const sstring& N_file, double N_scale)
{
  CLOG(7) << "Setting expression profile clustering parameters\n";

  if (Emu_file != "")
    {
      CLOG(6) << "Reading profile expectation vector from file '" << Emu_file << "'\n";
      ifstream Emu_stream (Emu_file.c_str());
      read_ColumnVector (Emu_stream, Emu, table.expt_idx, table.experiments());
    }
  else
    {
      CLOG(6) << "Using default profile expectation vector (zero)\n";
      Emu.ReSize (table.experiments());
      for (int i = 0; i < table.experiments(); ++i) Emu(i+1) = 0;
    }

  if (C_file != "")
    {
      CLOG(6) << "Reading profile covariance matrix from file '" << C_file << "'\n";
      ifstream C_stream (C_file.c_str());
      read_SymmetricMatrix (C_stream, C, table.expt_idx, table.experiments());
    }
  else
    {
      CLOG(6) << "Using default profile covariance matrix (identity matrix)\n";
      C.ReSize (table.experiments());
      for (int i = 0; i < table.experiments(); ++i)
	{
	  for (int j = 0; j < i; ++j) C(i+1,j+1) = 0;
	  C(i+1,i+1) = 1;
	}
    }
  
  if (N_file != "")
    {
      CLOG(6) << "Reading within-cluster covariance matrix from file '" << N_file << "'\n";
      ifstream N_stream (N_file.c_str());
      read_SymmetricMatrix (N_stream, N, table.expt_idx, table.experiments());
      N = N * N_scale;
    }
  else
    {
      CLOG(6) << "Using default within-cluster covariance matrix (identity matrix)\n";
      N.ReSize (table.experiments());
      for (int i = 0; i < table.experiments(); ++i)
	{
	  for (int j = 0; j < i; ++j) N(i+1,j+1) = 0;
	  N(i+1,i+1) = N_scale;
	}
    }
  debug_display();
  calculate_inverses();
  debug_display_inv();
}

void Gaussian_cluster_params::debug_display() const
{
  if (CLOGGING(-1)) {
    Matrix_to_array2d_adaptor C_a (C);
    Matrix_to_array2d_adaptor N_a (N);
    ColumnVector_to_vector_adaptor Emu_a (Emu);
    CLOG(-1) << "C:\n" << C_a;
    CLOG(-1) << "N:\n" << N_a;
    CLOG(-1) << "Emu:\n" << Emu_a;
  }
}

void Gaussian_cluster_params::debug_display_inv() const
{
  if (CLOGGING(-1)) {
    Matrix_to_array2d_adaptor invC_a (invC);
    Matrix_to_array2d_adaptor invN_a (invN);
    CL << "invC:\n" << invC_a;
    CL << "invN:\n" << invN_a;
    CL << "log_invC_norm = " << Nats2Bits (log_invC_norm) << " bits\n";
    CL << "log_invN_norm = " << Nats2Bits (log_invN_norm) << " bits\n";
  }
}


Gaussian_cluster_null_model::Gaussian_cluster_null_model (const Scaled_expr_table& table, const Gaussian_cluster_params& params)
  : table (table),
    params (params),
    profile (table.probes()),
    profile_null_exponent (table.probes())
{
  CLOG(6) << "Making null model for expression vectors\n";

  for (int n = 0; n < table.probes(); ++n)
    {
      profile[n] = new vector_to_ColumnVector_adaptor (table.profile[n]);
      ColumnVector p_dev = *profile[n];
      for (int j = 0; j < table.experiments(); ++j) {
	p_dev(j+1) += table.probe_offset[n] + table.experiment_offset[j];
	p_dev(j+1) *= table.probe_scaling[n] * table.experiment_scaling[j];
	p_dev(j+1) -= params.Emu(j+1);
      }
      profile_null_exponent[n] = - 0.5 * (p_dev.t() * params.invC * p_dev).AsScalar();
    }
}

void Gaussian_cluster::display_all_scores (ostream& o, const Gaussian_cluster_null_model& null, const vector<Gaussian_cluster*>& model, double kappa)
{
  Stream_saver ss;  // temporary variable

  int old_prec = o.precision(3);
  ss.save_flags (o);
  ss.left_align (o);

  o << "Gaussian cluster log-likelihoods (in bits):\n";
  o << "[gauss] Profile    Null      ";
  for (int m = 0; m < (int) model.size(); ++m) o << " Model #" << setw(3) << m+1;
  o << "\n";
  
  for (int g = 0; g < null.table.probes(); ++g)
    {
      o << "[gauss] ";
      sstring n = null.table.probe_name[g];
      if (n.size() > 10) n.erase (n.begin() + 10, n.end());
      o << setw(11) << n;
      ss.right_align (o);
      o << setw(10) << Score2Bits (null.profile_null_score (g, kappa));
      for (int m = 0; m < (int) model.size(); ++m)
	o << " " << setw(10) << Score2Bits (model[m]->log_likelihood_sc (g, kappa));
      ss.left_align (o);
      o << "\n";
    }

  ss.restore_flags (o);
  o.precision (old_prec);
}

Gaussian_cluster::Gaussian_cluster (const Gaussian_cluster_null_model& null)
  : params (null.params),
    null_model (null),
    x (null.profile),
    mu (params.Emu.Nrows()),
    weight (x.size(), (double) 0),
    invK_valid (0),
    sqrtK_valid (0)
{
}


void Gaussian_cluster::calculate_invK()
{
  if (!invK_valid)
    {
      CLOG(4) << "Gaussian_cluster: finding inverse of " << experiments() << "*" << experiments() << " transformation matrix\n";
      K = params.invC + effective_vectors() * params.invN;
      invK = K.i();
      if (CLOGGING(-1))
	{
	  Matrix_to_array2d_adaptor k_a (K);
	  Matrix_to_array2d_adaptor invk_a (invK);
	  CL << "K:\n" << k_a << "inv(K):\n" << invk_a;
	}
      invK_valid = 1;
      sqrtK_valid = 0;
    }
}

void Gaussian_cluster::calculate_sqrtK()
{
  calculate_invK();
  if (!sqrtK_valid)
    {
      CLOG(4) << "Gaussian_cluster: finding eigenvectors of " << experiments() << "*" << experiments() << " transformation matrix\n";
      DiagonalMatrix D (K.Nrows());
      Matrix V;
      Jacobi_catch (K, D, V);
      for (int i = 1; i <= D.Nrows(); ++i) D(i,i) = sqrt(D(i,i));
      sqrtK << V * D * V.t();    // have to use <<, not =, to prevent a newmat exception
      if (CLOGGING(-1))
	{
	  Matrix_to_array2d_adaptor sqrtk_a (sqrtK);
	  CL << "sqrt(K):\n" << sqrtk_a;
	}
      sqrtK_valid = 1;
    }
}


void Gaussian_cluster::calculate_mu_mean (ColumnVector& mu_mean)
{
  ColumnVector weighted_centroid (experiments());
  weighted_centroid = 0;

  for (int i = 0; i < vectors(); ++i) {
    for (int j = 0; j < experiments(); ++j) {
      double d = (*x[i]) (j+1);
      d += params.table.probe_offset[i] + params.table.experiment_offset[j];
      d *= params.table.probe_scaling[i] * params.table.experiment_scaling[j];
      weighted_centroid(j+1) += weight[i] * d;
    }
  }

  if (CLOGGING(0)) { ColumnVector_to_vector_adaptor wc (weighted_centroid); CL << "Weighted centroid = {" << wc << "}\n"; }

  calculate_invK();
  mu_mean = invK * (params.invC * params.Emu + params.invN * weighted_centroid);           // doing the invC * Emu here is slow - should probably move it elsewhere
}

void Gaussian_cluster::set_membership_probability (int point, int score)
{
  weight[point] = Score2Prob (score);
  dirty_K();
}

int Gaussian_cluster::log_likelihood_sc (int point, double kappa) const
{
  ColumnVector x_dev = *x[point];
  for (int j = 0; j < experiments(); ++j) {
    x_dev(j+1) += params.table.probe_offset[point] + params.table.experiment_offset[j];
    x_dev(j+1) *= params.table.probe_scaling[point] * params.table.experiment_scaling[j];
    x_dev(j+1) -= mu(j+1);
  }
  return Nats2Score (params.log_invN_norm - 0.5*log(kappa) - (1/kappa) * 0.5 * (x_dev.t() * params.invN * x_dev).AsScalar());
}

int Gaussian_cluster::log_parameter_prior_sc (double kappa) const
{
  const ColumnVector mu_dev = mu - params.Emu;
  return Nats2Score (params.log_invC_norm - 0.5*log(kappa) - (1/kappa) * 0.5 * (mu_dev.t() * params.invC * mu_dev).AsScalar());
}

void Gaussian_cluster::sample_parameters (double kT)
{
  ColumnVector mu_mean;
  calculate_mu_mean (mu_mean);
  for (int d = 0; d < experiments(); ++d) mu(d) = gennor (0.0, 1.0);
  calculate_sqrtK();
  mu = sqrtK * mu + mu_mean;
}

void Gaussian_cluster::optimise_parameters()
{
  calculate_mu_mean (mu);
}

void Gaussian_cluster::display (ostream& o)
{
  ColumnVector_to_vector_adaptor v (mu);
  o << "Mean: " << v << "\n";
  o << "SD:";
  for (int i = 0; i < experiments(); ++i)
    {
      double m = mu (i+1);
      double s = 0, n = 0;
      for (int p = 0; p < vectors(); ++p)
	{
	  double d = (*x[p]) (i+1) - m;
	  d = (d + params.table.probe_offset[p] + params.table.experiment_offset[i]) * params.table.probe_scaling[p] * params.table.experiment_scaling[p];
	  double w = weight[p];
	  s += w*d*d;
	  n += w;
	}
      o << " " << (n>0 ? sqrt(s/n) : 0);
    }
  o << "\n";
  o << "\n";
}

void Gaussian_cluster::display_scores (ostream& o, double kappa, double threshold)
{
  for (int i = 0; i < vectors(); ++i)
    if (weight[i] >= threshold)
      {
	o << "Probe '" << params.table.probe_name[i] << "' log-likelihood = " << Score2Bits (log_likelihood_sc (i, kappa)) << " bits";
	o << "  (null = " << Score2Bits (null_model.profile_null_score (i, kappa)) << " bits)\n";
      }
  o << "Log-parameter prior = " << Score2Bits(log_parameter_prior_sc()) << "\n";
}

Gaussian_cluster_factory::Gaussian_cluster_factory (const Gaussian_cluster_params& params, const Gaussian_cluster_null_model& null)
  : params (params), null (null)
{
  CLOG(7) << "Gaussian_cluster_factory constructed\n";
}

Gaussian_cluster* Gaussian_cluster_factory::new_model() const
{
  return new Gaussian_cluster (null);
}

