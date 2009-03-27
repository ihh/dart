#include "hsm/em_matrix.h"
#include "hsm/phylo_em.h"
#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "util/newmat_adaptors.h"
#include "util/vector_output.h"

EM_matrix::EM_matrix (int C,
		      int A,
		      int max_fork,
		      const Tree_alignment_database* align_db,
		      double timepoint_res)
  : EM_matrix_base (C, A, max_fork, align_db, timepoint_res)
{
  // hardwire the tolerance & max iterations for NR & EM. this is HACKY
  nr_tol = .1 * (double) n_params();
  nr_max_iter = n_params() * 5000;

  // hardwire against doing equilibration
  do_equilibrate = false;
}

void EM_matrix::init_alphabet (const Alphabet& base_alphabet)
{
  EM_matrix_base::init_alphabet (base_alphabet);	// does all but init_param_labels()
  init_param_labels (base_alphabet);
}

void EM_matrix::update_pi(const Update_statistics& stats)
{
  bool found_zero_up_triang = false;	// generate error message if protection against zero upper(i+1,i+1) invoked
  bool found_zero_pi = false;		// also protect against zero pi[i]
  // decompose R.t() by QR
  CTAG(3,RATE_EM RATE_EM_NEWMAT) << "Solving for equilibrium probabilities by QR decomposition\n";
  Matrix Q = R.t();
  UpperTriangularMatrix upper;
  QRZ (Q, upper);
  // solve upper triangular decomposition of R to find pi
  // this is the solution to the matrix equation pi * R = 0
  pi[m()-1] = 1;
  double pi_norm = pi[m()-1];
  for (int i = m() - 2; i >= 0; --i)
    {
      double x = 0;
      for (int j = i + 1; j < m(); ++j)
	x -= upper(i+1,j+1) * pi[j];
      if (abs(upper(i+1,i+1)) > TINY)
	  pi[i] = x / upper(i+1,i+1);
      else
	{
	  found_zero_up_triang = true;
	  break;	// recover using start counts
	}
      if (pi[i] < TINY)
	{
	  found_zero_pi = true;
	  break;
	}
      pi_norm += pi[i];
    }
  if (found_zero_up_triang || found_zero_pi)	// recovery
  {
    pi_norm = 0;
    for (int ci = 0; ci < C; ++ci)
      for (int ai = 0; ai < A; ++ai)
      {
        const int i = ca (ci, ai);
        pi[i] = stats.s[i];
        pi_norm += pi[i];
      }
  }
  // normalise, compute log_pi
  for (int i = 0; i < m(); ++i)
    {
      pi[i] /= pi_norm;
      log_pi[i] = Prob2Nats (pi[i]);
    }
  if (found_zero_up_triang)
    CLOGERR << "Warning: attempt to find equilibrium distribution directly from rate matrix\n"
	    << "failed (possibly due to noncommunicating states or other sparse matrix problem).\n"
	    << "Setting initial distribution from observed counts instead.\n"
	    << "(Consider using pseudocounts to avoid problems like this)\n";
  if (found_zero_pi)
    CLOGERR << "Warning: attempt to find equilibrium distribution directly from rate matrix\n"
	    << "yielded zero or negative probabilities. Setting initial distribution from\n"
	    << "observed counts instead. (This may have been caused by a sparse matrix;\n"
	    << "consider using pseudocounts to avoid problems like this)\n";
  // print out debugging message
  if (CTAGGING(-1,RATE_EM_CHECK_EQM))
    {
      vector_to_RowVector_adaptor pi_row (pi);
      CL << "QR upper diagonal of R.t():\n" << upper;
      CL << "pi:\n" << pi_row;
      CL << "pi * R:\n" << pi_row * R;
    }
}

void EM_matrix::diagonalize()
{
  // calculate sqrt_pi_transform
  vector<double> sqrt_pi_transform (m());
  for (int i = 0; i < m(); ++i)
    if (pi[i] > TINY)
      sqrt_pi_transform[i] = sqrt(pi[i]);
    else
      sqrt_pi_transform[i] = 1.0;

  // calculate symmetric matrix S
  SymmetricMatrix S(m());
  double sym_sq = 0;
  double dev_sq = 0;
  double n = 0;
  for (int i = 0; i < m(); ++i)
    for (int j = i; j < m(); ++j)
    {
	  const double S_ij = R(i+1,j+1) * sqrt_pi_transform[i] / sqrt_pi_transform[j];
	  const double S_ji = R(j+1,i+1) * sqrt_pi_transform[j] / sqrt_pi_transform[i];
	  const double sym = (S_ij + S_ji) / 2;
	  S(i+1,j+1) = sym;
	  if (j > i)
	  {
	    const double dev = S_ij - S_ji;
	    dev_sq += dev * dev;
	    sym_sq += sym * sym;
	    n += 1;
	  }
    }
  double rms_sym = sqrt(sym_sq/n);
  double rms_dev = sqrt(dev_sq/n);
  if (CTAGGING(1,RATE_EM_SYM))
    CL << "Rate matrix R:\n" << R << "Matrix S:\n" << S << "Equilibrium pi: (" << pi << ")\n";
  if (CTAGGING(2,RATE_EM RATE_EM_SYM RATE_EM_CHECK_EQM))
  {
    CL << "RMS asymmetry = " << rms_dev;
    if (rms_sym != 0) CL << " (" << 100*rms_dev/rms_sym << "% of element RMS)";
    CL << "\n";
  }

  // call Newmat to do eigenvector decomposition
  Matrix U_mx (m(), m());
  DiagonalMatrix mu_diag (m());
  CTAG(3,RATE_EM RATE_EM_NEWMAT) << "Finding eigenvectors using Jacobi algorithm\n";
  Jacobi_catch (S, mu_diag, U_mx);

  // copy into our eigensystem data structures
  for (int i = 0; i < m(); ++i)
    {
      mu[i] = mu_diag(i+1,i+1);
      for (int j = 0; j < m(); ++j)
	Uinv(j,i) = U(i,j) = U_mx(i+1,j+1);
    }

  // do inverse sqrt_pi_transform
  for (int i = 0; i < m(); ++i)
    for (int j = 0; j < m(); ++j)
      {
	U(i,j) /= sqrt_pi_transform[i];
	Uinv(i,j) *= sqrt_pi_transform[j];
      }
}

void EM_matrix::transform_symmetrised_waits_transitions (Update_statistics& stats, bool symmetrise) const
{
  transform_waits_transitions (stats, symmetrise);
}

Update_statistics EM_matrix::get_stats (bool infer_class_labels) const
{
  Update_statistics stats = get_stats_unequilibrated (TRUE, infer_class_labels);  // symmetrise
  if (do_equilibrate)  // IH, 4/20/2005
    equilibrate (stats);
  if (CTAGGING(3,RATE_EM_STATS)) CL << "Equilibrated update statistics:\n" << stats;
  return stats;
}

void EM_matrix::equilibrate (Update_statistics& stats) const
{
  // incorporate s[] counts into w[] and u()
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca (ci, ai);
	stats.w[i] -= stats.s[i] / R(i+1,i+1);
	for (int aj = 0; aj < A; ++aj)
	  if (aj != ai)
	    {
	      const int j = ca(ci,aj);
	      const double p = stats.s[j] * R(j+1,i+1) / R(j+1,j+1);
	      stats.u(i,j) -= p / 2;
	      stats.u(j,i) -= p / 2;
	    }
	for (int cj = 0; cj < C; ++cj)
	  if (cj != ci)
	    {
	      const int j = ca(cj,ai);
	      const double p = stats.s[j] * R(j+1,i+1) / R(j+1,j+1);
	      stats.u(i,j) -= p / 2;
	      stats.u(j,i) -= p / 2;
	    }
      }
}

map<int,double> EM_matrix::param_minima() const
{
  const double pmin = .1 / (double) m();
  map<int,double> xmin;
  for (int c = 0; c < C; ++c)
    for (int a = 0; a < A; ++a)
      xmin[pi_idx(c,a)] = pmin;
  for (int c = 0; c < C; ++c)
    for (int ai = 0; ai < A; ++ai)
      for (int aj = ai + 1; aj < A; ++aj)
	xmin[X_idx(c,ai,aj)] = pmin;
  for (int a = 0; a < A; ++a)
    for (int ci = 0; ci < C; ++ci)
      for (int cj = ci + 1; cj < C; ++cj)
	xmin[Y_idx(ci,cj,a)] = pmin;
  return xmin;
}

Diff EM_matrix::pi_diff (int c, int a) const
{
  Diff_param D (pi_idx (c, a));
  return D;
}

Diff EM_matrix::X_diff (int c, int ai, int aj) const
{
  Diff D;
  if (aj > ai)
    D = Diff_param (X_idx (c, ai, aj));
  else if (aj < ai)
    {
      D = Diff_param (X_idx (c, aj, ai));
      D *= Diff_param (pi_idx (c, aj));
      D /= Diff_param (pi_idx (c, ai));
    }
  else  // ai == aj
    for (int ak = 0; ak < A; ++ak)
      if (ak != ai)
	D -= X_diff (c, ai, ak);
  return D;
}

Diff EM_matrix::Y_diff (int ci, int cj, int a) const
{
  Diff D;
  if (cj > ci)
    D = Diff_param (Y_idx (ci, cj, a));
  else if (cj < ci)
    {
      D = Diff_param (Y_idx (cj, ci, a));
      D *= Diff_param (pi_idx (cj, a));
      D /= Diff_param (pi_idx (ci, a));
    }
  else  // ci == cj
    for (int ck = 0; ck < C; ++ck)
      if (ck != ci)
	D -= Y_diff (ci, ck, a);
  return D;
}

Diff EM_matrix::lagrangian (const Update_statistics& stats)
{
  Diff L;
  Diff pi_sum;
  double s_sum = 0;
  vector<Diff> pi_R (m(), Diff());
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca(ci,ai);
	
	Diff_const s_i (stats.s[i]);
	Diff_const w_i (stats.w[i]);

	Diff pi_i = pi_diff (ci, ai);
	Diff X_ii = X_diff (ci, ai, ai);
	Diff Y_ii = Y_diff (ci, ci, ai);

	Diff log_pi_i = pi_i.log();
	Diff R_ii = X_ii + Y_ii;
	Diff pi_R_ii = pi_i * R_ii;

	L += s_i * log_pi_i;
	L += w_i * R_ii;

	s_sum += stats.s[i];
	pi_sum += pi_i;

	pi_R[i] += pi_R_ii;

	for (int aj = 0; aj < A; ++aj)
	  if (ai != aj)
	    {
	      const int j = ca(ci,aj);

	      Diff_const u_ij (stats.u(i,j));

	      Diff X_ij = X_diff (ci, ai, aj);
	      Diff pi_j = pi_diff (ci, aj);

	      Diff log_X_ij = X_ij.log();
	      Diff pi_X_ij = pi_i * X_ij;

	      L += u_ij * log_X_ij;

	      pi_R[j] += pi_X_ij;
	    }

	for (int cj = 0; cj < C; ++cj)
	  if (ci != cj)
	    {
	      const int j = ca(cj,ai);

	      Diff_const u_ij (stats.u(i,j));

	      Diff_param kappa_j (kappa_idx (cj, ai));

	      Diff Y_ij = Y_diff (ci, cj, ai);
	      Diff pi_j = pi_diff (cj, ai);

	      Diff log_Y_ij = Y_ij.log();
	      Diff pi_Y_ij = pi_i * Y_ij;

	      L += u_ij * log_Y_ij;

	      pi_R[j] += pi_Y_ij;
	    }
      }
  L += Diff_const (-s_sum) * pi_sum.log();

  Diff kappa_sum;
  for (int c = 0; c < C; ++c)
    for (int a = 0; a < A; ++a)
      kappa_sum += Diff_param (kappa_idx (c, a)) * pi_R [ca (c, a)];
  // dividing kappa_sum by pi_sum gives many 2nd derivatives that otherwise would vanish...
  //  L += kappa_sum / pi_sum;
  // ...so we don't
  L += kappa_sum;

  L += Diff_param (beta_idx()) * (pi_sum + Diff_const (-1));
  
  L.param_name = &param_label;
  return L;
}

vector<double> EM_matrix::unconstrained_optimum (const Update_statistics& stats) const
{
  vector<double> x (n_params(), (double) 0);
  double beta = 0;
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca(ci,ai);
	const double s_i = stats.s[i];
	const double w_i = stats.w[i];
	beta += s_i;
	for (int aj = ai + 1; aj < A; ++aj)
	  {
	    const int j = ca(ci,aj);
	    const double Xu_ij = stats.u(i,j);
	    //	    const double Xu_ji = stats.u(j,i);
	    double& X_ij = x[X_idx(ci,ai,aj)];
	    //	    X_ij = (Xu_ij + Xu_ji) / (2 * w_i);
	    X_ij = Xu_ij / w_i;
	  }
	for (int cj = ci + 1; cj < C; ++cj)
	  {
	    const int j = ca(cj,ai);
	    const double Yu_ij = stats.u(i,j);
	    //	    const double Yu_ji = stats.u(j,i);
	    double& Y_ij = x[Y_idx(ci,cj,ai)];
	    //	    Y_ij = (Yu_ij + Yu_ji) / (2 * w_i);
	    Y_ij = Yu_ij / w_i;
	  }
      }
  for (int c = 0; c < C; ++c)
    for (int a = 0; a < A; ++a)
      {
	const int i = ca(c,a);
	const double s_i = stats.s[i];
	double& pi_i = x[pi_idx(c,a)];
	pi_i = s_i / beta;
      }
  return x;
}

void EM_matrix::set_params_from_vector (const vector<double>& x)
{
  double beta = 0;
  for (int c = 0; c < C; ++c)
    for (int a = 0; a < A; ++a)
      beta += x[pi_idx(c,a)];
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca(ci,ai);
	const double pi_i = x[pi_idx(ci,ai)] / beta;
	pi[i] = pi_i;
	for (int aj = ai + 1; aj < A; ++aj)
	  {
	    const double X_ij = x[X_idx(ci,ai,aj)];
	    const double pi_j = x[pi_idx(ci,aj)] / beta;
	    X[ci](ai,aj) = X_ij;
	    X[ci](aj,ai) = X_ij * pi_i / pi_j;
	  }
	double X_ii = 0;
	for (int aj = 0; aj < A; ++aj)
	  if (ai != aj)
	    X_ii -= X[ci](ai,aj);
	X[ci](ai,ai) = X_ii;

	for (int cj = ci + 1; cj < C; ++cj)
	  {
	    const double Y_ij = x[Y_idx(ci,cj,ai)];
	    const double pi_j = x[pi_idx(cj,ai)] / beta;
	    Y[ai](ci,cj) = Y_ij;
	    Y[ai](cj,ci) = Y_ij * pi_i / pi_j;
	  }
	double Y_ii = 0;
	for (int cj = 0; cj < C; ++cj)
	  if (ci != cj)
	    Y_ii -= Y[ai](ci,cj);
	Y[ai](ci,ci) = Y_ii;
      }
  update();
}

void EM_matrix::init_param_labels (const Alphabet& base_alphabet)
{
  param_label = vector<sstring> (n_params(), sstring());
  for (int ci = 0; ci < C; ++ci)
    {
      const char ci_char = ci < 10 ? ('0' + ci) : ('a' + ci - 10);
      for (int ai = 0; ai < A; ++ai)
	{
	  const char ai_char = toupper (base_alphabet.int2char(ai));
	  param_label[pi_idx(ci,ai)] << "Pi" << ai_char << ci_char;
	  param_label[kappa_idx(ci,ai)] << "kappa" << ai_char << ci_char;
	  for (int aj = ai + 1; aj < A; ++aj)
	    {
	      const char aj_char = toupper (base_alphabet.int2char(aj));
	      param_label[X_idx(ci,ai,aj)] << "X" << ai_char << aj_char << ci_char;
	    }
	  for (int cj = ci + 1; cj < C; ++cj)
	    {
	      const char cj_char = cj < 10 ? ('0' + cj) : ('a' + cj - 10);
	      param_label[Y_idx(ci,cj,ai)] << "Y" << ai_char << ci_char << cj_char;
	    }
	}
    }
  param_label[beta_idx()] = "beta";
}

Update_statistics EM_matrix::single_EM (bool infer_class_labels)
{
  CTAG(5,RATE_EM RATE_EM_PROGRESS) << "Computing update statistics (E-phase)\n";
  const Update_statistics stats = get_stats (infer_class_labels);
  const Diff L = lagrangian (stats);
  Laplacian laplace (L);
  const vector<double> nr_seed = unconstrained_optimum (stats);
  const map<int,double> nr_min = param_minima ();
  CTAG(5,RATE_EM RATE_EM_PROGRESS) << "Starting Newton-Raphson (M-phase)\n";
  const vector<double> nr_max = Newton_Raphson::iterate (laplace, nr_seed, nr_min, nr_tol, nr_max_iter);
  set_params_from_vector (nr_max);
  return stats;
}

Loge EM_matrix::iterate_EM (bool infer_class_labels)
{
  Loge best_log_likelihood = 0;
  EM_matrix_params best_params;
  for (int iter = 0; ; ++iter)
    {
      if (iter >= em_max_iter)
	{
	  CTAG(7,RATE_EM RATE_EM_PROGRESS) << "EM hit " << em_max_iter << " iterations; stopping\n";
	  break;
	}
      const EM_matrix_params old_params = *this;
      Update_statistics stats = single_EM (infer_class_labels);
      const Loge prev_best = best_log_likelihood;
      if (iter == 0 || stats.log_likelihood > best_log_likelihood)
	{
	  best_log_likelihood = stats.log_likelihood;
	  best_params = old_params;
	}
      CTAG(6,RATE_EM RATE_EM_PROGRESS) << "EM iteration #" << iter+1 << ": log-likelihood = " << Nats2Bits(stats.log_likelihood) << " bits\n";
      if (iter > 0)
	{
	  const double inc = -(stats.log_likelihood - prev_best) / prev_best;
	  if (inc < em_min_inc)
	    {
	      if (stats.log_likelihood < prev_best)
		CTAG(7,RATE_EM RATE_EM_PROGRESS) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_best) << " to " << Nats2Bits(stats.log_likelihood) << " bits during EM\n";
	      CTAG(7,RATE_EM RATE_EM_PROGRESS) << "Fractional log-likelihood increase below EM improvement threshold; stopping\n";
	      break;
	    }
	}
    }
  // since we only update the best log-likelihood the round *after* the new params are calculated, we need to check for an increase one last time
  CTAG(6,RATE_EM RATE_EM_PROGRESS) << "Checking post-iteration log-likelihood\n";
  const Loge final_log_likelihood = log_likelihood();
  CTAG(6,RATE_EM RATE_EM_PROGRESS) << "Post-iteration log-likelihood = " << Nats2Bits(final_log_likelihood) << " bits\n";
  if (em_max_iter != 0 && final_log_likelihood < best_log_likelihood)
    {
      ((EM_matrix_params&) *this) = best_params;
      update();
    }
  else
    best_log_likelihood = final_log_likelihood;
  // and return
  return best_log_likelihood;
}

Update_statistics EM_matrix::single_quick_EM (bool intra, bool inter, bool infer_class_labels)
{
  CTAG(5,RATE_EM RATE_EM_PROGRESS) << "Computing update statistics (E-phase)\n";
  const Update_statistics stats = get_stats (infer_class_labels);
  CTAG(5,RATE_EM RATE_EM_PROGRESS) << "Finding nearest reversible rate matrix (M-phase)\n";
  quick_M (stats, intra, inter);
  return stats;
}

void EM_matrix::quick_M (const Update_statistics& stats, bool intra, bool inter)
{
  // get X & Y
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
    {
      const int i = ca (ci, ai);

      if (intra)
      {
        X[ci](ai,ai) = 0;
        for (int aj = 0; aj < A; ++aj)
          if (aj != ai)
	    {
	      if (X_update_flag[ci](ai, aj))
		X[ci](ai, aj) = stats.u(i,ca(ci,aj)) / stats.w[i];
	      X[ci](ai,ai) -= X[ci](ai, aj);
	    }
      }

      if (inter && C > 1)
      {
        Y[ai](ci,ci) = 0;
        for (int cj = 0; cj < C; ++cj)
          if (cj != ci)
	    {
	      if (Y_update_flag[ai](ci,cj))
		Y[ai](ci,cj) = stats.u(i,ca(cj,ai)) / stats.w[i];
	      Y[ai](ci,ci) -= Y[ai](ci,cj);
	    }
      }
    }

  // need update_R() before update_pi()
  update_R();
  update_pi(stats);
  update();
}
