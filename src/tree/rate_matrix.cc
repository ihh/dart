#include <math.h>
#include <algorithm>
#include "tree/rate_matrix.h"
#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "util/input_string.h"
#include "util/logfile.h"
#include "util/newmat_adaptors.h"


Rate_matrix::Rate_matrix (const Alphabet& a, const vector<double>& freq, const Matrix& evec_cols, const vector<double>& eval)
  : a (a.name.c_str(), a.nondegenerate_chars()),
    frequency(freq),
    sqrt_freq(freq.size()),
    D(U.Ncols()),
    U(evec_cols)
{
  for (int i = 0; i < (int) eval.size(); ++i) D(i+1,i+1) = eval[i];
  for (int i = 0; i < (int) freq.size(); ++i) sqrt_freq[i] = sqrt (freq[i]);
}

void Rate_matrix::read (istream& in)
{
  CLOG(6) << "Reading in alphabet, residue frequencies, unbiased rate matrix\n";

  // read in the alphabet
  //
  Input_string first_line;
  in >> first_line;
  vector<sstring> chars = first_line.split();
  for_contents (vector<sstring>, chars, s)
    if ((*s).size() != 1)
      THROW Format_exception ("Alphabet symbols may only be one character wide");
  sstring alphabet_chars = sstring::join (chars, "");
  a = Simple_alphabet ("user",alphabet_chars.c_str());

  CLOG(0) << "Alphabet: " << alphabet_chars << "\n";

  // read in residue frequencies
  //
  sqrt_freq = frequency = vector<double> (a.size(), (double) 0);
  for (int i = 0; i < a.size(); i++)
    {
      if (in.eof()) THROW Format_exception ("Too few elements in residue frequency vector");
      in >> frequency[i];
      sqrt_freq[i] = sqrt (frequency[i]);
    }
  
  // read in lower triangle of unbiased rate matrix S (excluding diagonal elements) row by row
  // actual rate matrix R is obtained by multiplying columns of S by frequency vector (R = S * F)
  // the convention used here is that dP/dt = R*P where P is the residue probability vector
  // so that R_ij is the rate of changing from j to i, and columns add up to zero
  //
  SymmetricMatrix S (a.size());
  for (int row = 1; row <= a.size(); row++)
    for (int col = 1; col < row; col++) 
      {
	if (in.eof()) THROW Format_exception ("Too few elements in unbiased rate matrix");
	in >> S(row,col);
      }

  // force rows of R to add up to zero
  //
  for (int row = 1; row <= a.size(); row++)
    {
      double row_total = 0;
      for (int col = 1; col <= a.size(); col++)
	if (col != row) row_total += S(row,col) * frequency[col-1];
      S(row,row) = -row_total / frequency[row-1];
    }

  CLOG(6) << "Diagonalising symmetrised rate matrix\n";
  
  // find the symmetric matrix Z = sqrt_freq * R * (sqrt_freq^-1) = sqrt_freq * S * sqrt_freq
  //
  DiagonalMatrix sqrt_freq_diag (a.size());
  for (int i = 1; i <= a.size(); i++) sqrt_freq_diag(i,i) = sqrt_freq[i-1];
  SymmetricMatrix Z (a.size());
  Z << sqrt_freq_diag * S * sqrt_freq_diag;    // have to use <<, not =, to prevent a newmat exception

  CLOG(-1) << "S matrix:\n" << S;
  CLOG(-1) << "Z matrix:\n" << Z;

  // diagonalise Z
  //
  Jacobi_catch (Z, D, U);

  CLOG(6) << "Rate matrix initialisation completed\n";

  for (int j = 1; j <= a.size(); j++)
    {
      CLOG(0) << "Eigenvalue #" << j << ": " << D(j,j) << "\tEigenvector: (";
      for (int i = 1; i <= a.size(); i++)
	{
	  clog_stream << U(i,j);
	  if (i < a.size()) clog_stream << ", ";
	}
      clog_stream << ")\n";
    }

  CLOG(-1) << "Z * U:\n" << Z * U << "U * D:\n" << U * D << "U * Ut:\n" << U * U.t();
}

array2d<double> Rate_matrix::create_conditional_substitution_matrix (double time)
{
  vector<double> exp_ev (a.size(), (double) 0);
  for (int k = 0; k < a.size(); k++) exp_ev[k] = exp(D(k+1,k+1) * time);

  array2d<double> sub (a.size(), a.size(), 0);
  for (int i = 0; i < a.size(); i++)
    for (int j = 0; j < a.size(); j++)
      {
	for (int k = 0; k < a.size(); k++)
	  sub(j,i) += U(i+1,k+1) * U(j+1,k+1) * exp_ev[k];
	sub(j,i) *= sqrt_freq[i] / sqrt_freq[j];
	sub(j,i) = max (sub(j,i), 0.0);
      }

  if (CLOGGING(-1)) clog_stream << "Substitution matrix at t=" << time << ":\n" << sub;
  if (CLOGGING(-2))
    {
      array2d_to_Matrix_adaptor submat = sub;
      vector_to_RowVector_adaptor p = create_prior();
      clog_stream << "Prior: " << p << "Substitution matrix * prior: " << p * submat.t();
    }
  
  return sub;
}

array2d<double> Rate_matrix::differentiate_conditional_substitution_matrix (double time)
{
  vector<double> dexp_ev (a.size(), (double) 0);
  for (int k = 0; k < a.size(); k++) dexp_ev[k] = D(k+1,k+1) * exp(D(k+1,k+1) * time);

  array2d<double> sub (a.size(), a.size(), 0);
  for (int i = 0; i < a.size(); i++)
    for (int j = 0; j < a.size(); j++)
      {
	for (int k = 0; k < a.size(); k++)
	  sub(j,i) += U(i+1,k+1) * U(j+1,k+1) * dexp_ev[k];
	sub(j,i) *= sqrt_freq[i] / sqrt_freq[j];
	sub(j,i) = max (sub(j,i), 0.0);
      }

  if (CLOGGING(-1)) clog_stream << "Time derivative of substitution matrix at t=" << time << ":\n" << sub;
  
  return sub;
}

bool Rate_matrix::rate_matrix_help (Opts_list* ol)
{
  CLOGERR << "\n";
  CLOGERR << "Format for rate matrix file\n";
  CLOGERR << "===========================\n";
  CLOGERR << "First line specifies the symbol alphabet, currently one-character symbols separated by whitespace, e.g. \"a c g t\".\n";
  CLOGERR << "Second line specifies equilibrium frequencies for the alphabet symbols. These need not be normalised\n";
  CLOGERR << "(although the scale factor will modulate the overall substitution rate).\n";
  CLOGERR << "Following this are the elements of the off-diagonal lower triangle of the symmetric part of the rate matrix.\n";
  CLOGERR << "The reversible rate matrix R (from the equation dP/dt = R*P, where P is the column vector representing the\n";
  CLOGERR << "state occupancy probabilities) is decomposable into a symmetric matrix S and a diagonal matrix F,\n";
  CLOGERR << "i.e. R = S * F, where the diagonal elements of F are the equilibrium symbol frequencies.\n";
  CLOGERR << "Entry R_ij (row i, column j) is the rate of change of symbol j into symbol i.\n";
  CLOGERR << "Since the format is specified in terms the symmetric matrix S, the distinction between S_ij and S_ji vanishes,\n";
  CLOGERR << "but to get the order right it's important to realise that the lower triangle is being entered.\n";
  CLOGERR << "Since diagonal elements are omitted, the number of matrix entries that must be specified is 0.5 * N * (N-1)\n";
  CLOGERR << "where N is the alphabet size.\n";
  CLOGERR << "For example, a DNA substitution rate matrix file might look like this:\n";
  CLOGERR << "\n";
  CLOGERR << "A    C    G    T\n";
  CLOGERR << "F_A  F_C  F_G  F_T\n";
  CLOGERR << "S_CA\n";
  CLOGERR << "S_GA S_GC\n";
  CLOGERR << "S_TA S_TC S_TG\n";
  CLOGERR << "\n";
  CLOGERR << "... where the F_i and S_ij are floating-point values.\n";
  CLOGERR << "For Kimura's two-parameter model, F_A=F_C=F_G=F_T = 1, S_CA=S_GC=S_TA=S_TG = 1/(r+2), S_GA=S_TC = r/(r+2)\n";
  CLOGERR << "where r is the transition/transversion ratio.\n";
  CLOGERR << "\n";
  exit(0);
  return 0;
}

