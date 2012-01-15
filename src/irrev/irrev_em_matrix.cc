#include "irrev/irrev_em_matrix.h"
#include "newmat/newmatap.h"
#include "newmat/nm_exp.h"
#include "newmat/linalg.h"
#include "newmat/newmatio.h"
#include "util/newmat_adaptors.h"
#include "util/vector_output.h"
#include "hsm/phylo_em.h"

Irrev_EM_matrix::Irrev_EM_matrix (int C,
				  int A,
				  int max_fork,
				  const Tree_alignment_database* align_db,
				  double timepoint_res)
  : EM_matrix_base (C, A, max_fork, align_db, timepoint_res)
{ }

Irrev_EM_matrix::Irrev_EM_matrix(void) // added by OW 7-15-2011
 :EM_matrix_base (1, 1, 1, 0, DEFAULT_TIMEPOINT_RES)
{ }


// matrix diagonalization - S is the symmetrized version of Sb
void Irrev_EM_matrix::diagonalize()
{
  // calculate sqrt_pi_transform
  vector<double> sqrt_pi_transform (m());
  for (int i = 0; i < m(); ++i)
    if (pi[i] > TINY)
      sqrt_pi_transform[i] = sqrt(pi[i]);
    else
      sqrt_pi_transform[i] = 1.0;

  // apply sqrt_pi_transform to calculate symmetric matrix S and asymmetric Sb
  SymmetricMatrix S(m());
  Matrix Sb (m(), m());
  for (int i = 0; i < m(); ++i)
    for (int j = i; j < m(); ++j)
    {
      const double S_ij = R(i+1,j+1) * sqrt_pi_transform[i] / sqrt_pi_transform[j];
      const double S_ji = R(j+1,i+1) * sqrt_pi_transform[j] / sqrt_pi_transform[i];
      const double sym = (S_ij + S_ji) / 2;
      S(i+1,j+1) = sym;
      Sb(i+1,j+1) = S_ij;
      Sb(j+1,i+1) = S_ji;
    }
  if (CTAGGING(1,RATE_EM_SYM))
    CL << "Rate matrix R:\n" << R << "Matrix Sb:\n" << Sb << "Equilibrium pi: (" << pi << ")\n";

  // call Newmat for eigenvector decomposition
  bool used_MatrixExpEigen = true;
  int MatrixExp_retries = 0;
  Matrix Sb_tmp = Sb;		// MatrixExpEigenPrepare destroys the input matrix, so use Sb_tmp as a temporary matrix

  SimpleIntArray cplx (m());
  Matrix U_mx (m(), m());
  Matrix Uinv_mx (m(), m());
  DiagonalMatrix mu_diag (m());
  ColumnVector mu_column;

  CTAG(3,RATE_EM RATE_EM_NEWMAT) << "Finding eigenvectors using MatrixExpEigenPrepare\n";

  int MatrixExp_status;
  while (true)
    {
      CTAG(2,RATE_EM_NEWMAT NEWMAT_BUG) << "Calling MatrixExpEigenPrepare on the following matrix\n" << Sb_tmp;
      MatrixExp_status = MatrixExpEigenPrepare (Sb_tmp, mu_column, U_mx, Uinv_mx, cplx);

      if (MatrixExp_status == NO_ERROR)
	break;
      if (++MatrixExp_retries > MATRIXEXP_RETRIES)
	break;

      // perturb if matrix is singular
      for (int i = 0; i < m(); i++)
	for (int j = 0; j < m(); j++)
	  Sb(i+1,j+1) += MATRIXEXP_PERTURB * (S(i+1,j+1) - Sb(i+1,j+1));	// small perturbation toward symmetric matrix

      Sb_tmp = Sb;
    }

  if (MatrixExp_status != NO_ERROR)	// reached retry count
    {
      CTAG(3,RATE_EM RATE_EM_NEWMAT) << "MatrixExpEigenPrepare perturb failed - use Jacobi\n";
      Jacobi_catch (S, mu_diag, U_mx);
      Uinv_mx = U_mx.t();
      Sb = S;  // precaution, just in case we later refer to Sb as the matrix that was actually diagonalized
      used_MatrixExpEigen = false;
    }

  if (MatrixExp_retries > 0)
    {
      double perturb = MATRIXEXP_PERTURB;
      CTAG(3,RATE_EM RATE_EM_NEWMAT) << "Retried MatrixExpEigenPrepare " << MatrixExp_retries << " time(s); perturbation ~ " << perturb << "\n";
    }

  // copy into our eigensystem data structures
  vector<int> cplx_pairs;
  for (int k = 0; k < m(); k++)	// did we wander off onto the complex plane?
    if (used_MatrixExpEigen)
      {
	if (cplx[k])
	  {
	    CTAG(3,RATE_EM) << "MatrixExpEigenPrepare eigenvalues #"
			    << k+1 << ',' << k+2 << ": ("
			    << mu_column(k+1) << " +/- i*" << mu_column(k+2)
			    << " [complex]\n";

	    mu[k] = Complex (mu_column(k+1), mu_column(k+2));
	    mu[k+1] = std::conj (mu[k]);

	    for (int i = 0; i < m(); ++i)
	      {
		U(i,k) = Complex (U_mx(i+1,k+1), U_mx(i+1,k+2));
		U(i,k+1) = std::conj (U(i,k));

		Uinv(k,i) = Complex (Uinv_mx(k+1,i+1) / 2, -Uinv_mx(k+2,i+1) / 2);   // weird EISPACK format for Uinv...
		Uinv(k+1,i) = std::conj (Uinv(k,i));
	      }

	    ++k;  // skip conjugate
	  }
	else  // used_MatrixExpEigen && cplx[k]==0
	  {
	    CTAG(3,RATE_EM) << "MatrixExpEigenPrepare eigenvalue #"
			    << k+1 << ": " << mu_column(k+1)
			    << " [real]\n";

	    mu[k] = mu_column(k+1);

	    for (int i = 0; i < m(); ++i)
	      {
		U(i,k) = U_mx(i+1,k+1);
		Uinv(k,i) = Uinv_mx(k+1,i+1);
	      }

	    if (CTAGGING(3,RATE_EM))
	      {
		CL << "MatrixExpEigenPrepare right eigenvector #" << k+1 << ": [";
		for (int i = 0; i < m(); ++i)
		  CL << ' ' << U(i,k);
		CL << " ]\n";

		CL << "MatrixExpEigenPrepare left eigenvector #" << k+1 << ": [";
		for (int i = 0; i < m(); ++i)
		  CL << ' ' << Uinv(k,i);
		CL << " ]\n";
	      }
	  }
      }
    else  // !used_MatrixExpEigen
      {
	mu[k] = mu_diag(k+1,k+1);
	for (int i = 0; i < m(); ++i)
	  {
	    U(i,k) = U_mx(i+1,k+1);
	    Uinv(k,i) = Uinv_mx(k+1,i+1);
	  }
      }

  // apply inverse sqrt_pi_transform
  for (int i = 0; i < m(); ++i)
    for (int j = 0; j < m(); ++j)
      {
	U(i,j) /= sqrt_pi_transform[i];
	Uinv(i,j) *= sqrt_pi_transform[j];
      }

}

void Irrev_EM_matrix::transform_symmetrised_waits_transitions (Update_statistics& stats, bool symmetrise) const
{
  if (symmetrise)
    CTAG(-1,RATE_EM) << "Warning: ignoring 'symmetrise' flag in Irrev_EM_matrix::transform_symmetrised_waits_transitions\n";
  transform_waits_transitions (stats, false);  // ignore symmetrise flag
}

Update_statistics Irrev_EM_matrix::single_quick_EM (bool intra, bool inter, bool infer_class_labels)
{
  CTAG(5,RATE_EM RATE_EM_PROGRESS) << "Computing update statistics (E-phase)\n";
  const Update_statistics stats = get_stats_unequilibrated (FALSE, infer_class_labels);	// don't equilibrate
  CTAG(5,RATE_EM RATE_EM_PROGRESS) << "Finding nearest reversible rate matrix (M-phase)\n";
  quick_M (stats, intra, inter);
  return stats;
}

void Irrev_EM_matrix::quick_M (const Update_statistics& stats, bool intra, bool inter)
{
  double pi_sum = 0;
  bool found_w_zero = false;
  // get X & Y
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca (ci, ai);
	pi[i] = stats.s[i];
	pi_sum += pi[i];

	if (intra)
	  {
	    X[ci](ai,ai) = 0;
	    for (int aj = 0; aj < A; ++aj)
	      if (aj != ai)
		{
		  if (X_update_flag[ci](ai, aj))
		    {
		      double newX;
		      if (stats.w[i] != 0)
			newX = stats.u(i,ca(ci,aj)) / stats.w[i];
		      else
			{
			  newX = 0;
			  found_w_zero = true;
			}
		      X[ci](ai, aj) = newX;
		    }
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
		    {
		      double newY;
		      if (stats.w[i] != 0)
			newY = stats.u(i,ca(cj,ai)) / stats.w[i];
		      else
			{
			  newY = 0;
			  found_w_zero = true;
			}
		      Y[ai](ci,cj) = newY;
		    }
		  Y[ai](ci,ci) -= Y[ai](ci,cj);
		}
	  }
      }
  if (found_w_zero)
    CLOGERR << "Warning: zero wait time found\n";

  // normalize pi[i]
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca (ci, ai);
	pi[i] /= pi_sum;
      }

  // update eigenvalues & eigenvectors
  update();
}
