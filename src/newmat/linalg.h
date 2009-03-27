/*  
|  
|  linalg.h  
|  
|       Gerton Lunter, lunter@stats.ox.ac.uk.  
|  
|	Prototypes for matrix exponentiation and inversion, and eigensystem functions.  
|  
|  
|  
|  1983:  EISPACK linear algebra package  (Burton S. Garbow, Argonne National Laboratory)  
|  
|  1998:  EISPACK Fortran package translated into C by David L. Swofford (Smithsonian Institution).  
|         This version is taken from MrBayes 3.0 by John Huelsenbeck (UCSD) and   
|         Frederik Ronquist (Uppsala University)  
|  
|  2003:  Added MatrixExp, Eigenvector, MatrixExpEigenPrepare, MatrixExpEigen  (Gerton Lunter, Oxford University)  
|         MatrixExp is a translation of a Matlab-implementation in EXPOKIT by Roger B. Sidje, rbs@maths.uq.edu.au.  
|  
|         Refs:-  
|           ACM - Transactions On Mathematical Software, 24(1):130-156, 1998  
|           Moler & Van Loan: Nineteen Dubious ways to compute the Exponential of a Matrix, Twenty-Five Years Later,  
|                             SIAM Review 45 (2003) No. 1, pp. 3--49.  
|  
|  
*/  
  
  
#ifndef __linalg_h__  
#define __linalg_h__  
  
  
int  MatrixExp (double **aaA, int iN, double iT, double **aaExp, double** aaW1, double** aaW2, double** aaW3, int* aW4);  
int  InvertMatrix (double **a, int n, double *col, int *indx, double **a_inv);  
int  Eigenvector (double **a, int n, double ev, double *col, int *indx, double **a_inv, char right);  
int  MatrixExpEigenPrepare (double **a, int n, double *v, double **u, double **uinv, int *cplx, double *work1, double *work2);  
void MatrixExpEigen(double **mexp, int n, double *v, double **u, double **uinv, int *cplx, double t, double *work1);  
int  LUDecompose (double **a, int n, double *vv, int *indx, double *pd);  
void LUBackSubst (double **a, int n, int *indx, double *b);  
  
  
#define NO_ERROR 0  
#define ERROR 1  
  
  
#define TINY		1.0e-20          /* Used to handle singular matrices in LUDecompose / MatrixInverse                        */  
#define EVECTORTOL      1.0e-10          /* Used in stopping condition for Eigenvector iteration                                   */  
#define PERTURB         1.0e-7           /* Used to perturb matrix away from generalized eigenvector case (MatrixExpEigenPrepare)  */  
  
  
#define VERBOSE 1                        /* Makes MatrixExpEigenPrepare shout if eigensystem is (nearly) degenerate                */  
  
  
#endif  
  
// body file: linalg.cpp
