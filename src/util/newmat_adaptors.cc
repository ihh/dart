#include "util/newmat_adaptors.h"
#include "util/logfile.h"
#include "newmat/newmatio.h"

void Jacobi_catch(const SymmetricMatrix& X, DiagonalMatrix& D, Matrix& A)
{
  try
    {
      Jacobi (X, D, A);
    }
  catch (ConvergenceException& e)
    {
      CLOGERR << "Newmat exception: Jacobi(X,D,A) fails to converge\n";
      CLOGERR << "X:\n" << X;
      CLOGERR << "D:\n" << D;
      CLOGERR << "A:\n" << A;
    }
  catch (...)
    {
      CLOGERR << "Newmat exception: Jacobi(X,D,A)\n";
      CLOGERR << "X:\n" << X;
      CLOGERR << "D:\n" << D;
      CLOGERR << "A:\n" << A;
    }
}
