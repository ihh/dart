// pk 4/05 int MatrixExpEigenPrepare

void MatrixExp(Matrix& A, Real t, Matrix& Exp);

void Eigenvector(Matrix& A, Real ev, ColumnVector& C, Matrix& A_inv, bool right);

int MatrixExpEigenPrepare(Matrix& A, ColumnVector& V, Matrix& U, Matrix& Uinv,
   SimpleIntArray& Cplx);

void MatrixExpEigen(Matrix& Mexp, ColumnVector& V, Matrix& U, Matrix& Uinv,
   SimpleIntArray& Cplx, Real t);

void LUDecompose(Matrix& A, SimpleIntArray& Indx, double& pd);
   
void InvertMatrix(Matrix& A, Matrix& A_inv);

void LUBackSubst(Matrix& A, SimpleIntArray& Indx, ColumnVector& B);

// body file: nm_exp.cpp

