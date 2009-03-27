#define WANT_STREAM

#include "newmatap.h"
#include "newmatio.h"
#include "nm_exp.h"
#include "linalg.h"

// pk 3/3/05 modify MatrixExpEigenPrepare for use with dart
// pk 4/05 pass return of Gerton's MatrixExpEigenPrepare to caller.
class C_sub2
{
   Real** a;
public:
   C_sub2(Matrix& A);
   ~C_sub2();
   operator Real**() { return a; }
};

C_sub2::C_sub2(Matrix& A)
{
   int n = A.ncols();
   int m = A.nrows();
   a = new Real*[m];
   if (!a) Throw(Bad_alloc("No space"));
   Real* d = A.data();
   for (int i = 0; i < n; ++i) a[i] = d + i * m;
} 

C_sub2::~C_sub2()
{
   delete[] a;
}

void MatrixExp(Matrix& A, Real t, Matrix& Exp)
{
   Tracer et("MatrixExp ");
   int n = A.nrows();
   if (n != A.ncols()) Throw(Logic_error("Not square"));
   Exp.ReSize(n,n);
   Matrix W1(n,n);
   Matrix W2(n,n);
   Matrix W3(n,n);
   SimpleIntArray W4(n);
   int i = MatrixExp(
      C_sub2(A),
      n,
      t,
      C_sub2(Exp),
      C_sub2(W1),
      C_sub2(W2),
      C_sub2(W3),
      W4.Data() );
   if (i) Throw(Runtime_error("Fails"));
}

void Eigenvector(Matrix& A, Real ev, ColumnVector& C, Matrix& A_inv, bool right)
{
   Tracer et("Eigenvector ");
   int n = A.nrows();
   if (n != A.ncols()) Throw(Logic_error("Not square"));
   C.ReSize(n);
   SimpleIntArray Indx(n);
   A_inv.ReSize(n, n);
   int i = Eigenvector(
      C_sub2(A),
      n,
      ev,
      C.data(),
      Indx.Data(),
      C_sub2(A_inv),
      (char)right );  
   if (i) Throw(Runtime_error("Fails"));
}

int MatrixExpEigenPrepare(Matrix& A, ColumnVector& V, Matrix& U, Matrix& Uinv,	// pk
   SimpleIntArray& Cplx)
{
   Tracer et("MatrixExpEigenPrepare ");
   int n = A.nrows();
   if (n != A.ncols()) Throw(Logic_error("Not square"));
   V.ReSize(n);
   U.ReSize(n, n);
   Uinv.ReSize(n, n);
   Cplx.ReSize(n);
   ColumnVector W1(n), W2(n);
   int i = MatrixExpEigenPrepare(
      C_sub2(A),
      n,
      V.data(),
      C_sub2(U),
      C_sub2(Uinv),
      Cplx.Data(),
      W1.data(),
      W2.data() );  
     /*     arguments of MatrixExpEigenPrepare */
     /*     **a = matrix to exponentiate (destroyed)  */  
     /*       n = order of matrix */  
     /*      *v = output vector (eigenvalues) */  
     /*     **u = output matrix (eigenvectors) */  
     /*  **uinv = output matrix (inverse of eigenvector matrix) */  
     /*    *cplx = output vector (booleans, true for real/imag eigenvalue pairs) */  
     /*  *work1 = work vector */  
     /*  *work2 = work vector */
   return i;	// pk
}

void MatrixExpEigen(Matrix& Mexp, ColumnVector& V, Matrix& U, Matrix& Uinv,
   SimpleIntArray& Cplx, Real t)
{
   Tracer et("MatrixExpEigen ");
   int n = Mexp.nrows();
   if (n != Mexp.ncols()) Throw(Logic_error("Not square"));
   if (n != V.nrows()) Throw(Logic_error("Incompatible dimension"));
   if (n != U.nrows()) Throw(Logic_error("Incompatible dimension"));
   if (n != U.ncols()) Throw(Logic_error("Incompatible dimension"));
   if (n != Uinv.nrows()) Throw(Logic_error("Incompatible dimension"));
   if (n != Uinv.ncols()) Throw(Logic_error("Incompatible dimension"));
   if (n != Cplx.Size()) Throw(Logic_error("Incompatible dimension"));
   ColumnVector W1(n);
     
   MatrixExpEigen(
      C_sub2(Mexp),
      n,
      V.data(),
      C_sub2(U),
      C_sub2(Uinv),
      Cplx.Data(),
      t,
      W1.data() );  
}





void LUDecompose(Matrix& A, SimpleIntArray& Indx, double& pd)
{
   Tracer et("LUDecompose ");
   int n = A.nrows();
   if (n != A.ncols()) Throw(Logic_error("Not square"));
   ColumnVector V(n);
   Indx.ReSize(n);
   int i = LUDecompose(C_sub2(A), n, V.data(), Indx.Data(), &pd);
   if (i) Throw(Runtime_error("Singular matrix"));
}
   
void InvertMatrix(Matrix& A, Matrix& A_inv)
{
   int n = A.nrows();
   if (n != A.ncols()) Throw(Logic_error("Not square"));
   SimpleIntArray Indx(n);
   ColumnVector V(n);
   A_inv.ReSize(n,n);
   int i = InvertMatrix(
      C_sub2(A),
      n,
      V.data(),
      Indx.Data(),
      C_sub2(A_inv) );  
   if (i) Throw(Runtime_error("Singular matrix"));
}

void LUBackSubst(Matrix& A, SimpleIntArray& Indx, ColumnVector& B)
{
   Tracer et("LUBackSubst ");
   int n = A.nrows();
   if (n != A.ncols()) Throw(Logic_error("Not square"));
   if (n != B.nrows()) Throw(Logic_error("Incompatible dimension"));
   if (n != Indx.Size()) Throw(Logic_error("Incompatible dimension"));
   LUBackSubst(C_sub2(A), n, Indx.Data(), B.data());
}   

