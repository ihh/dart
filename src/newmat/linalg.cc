/*  
|  
|  linalg.c  
|  
|       Gerton Lunter, lunter@stats.ox.ac.uk.  
|  
|	Routines for matrix exponentiation and inversion, and eigensystem functions.  
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
|  2005:  Pete K float to double in MatrixExp, include Gerton's bug fix in MatrixExpEigenPrepare, throw exceptions
|         on errors detected in MatrixExpEigenPrepare, dart logging replaces printf.  Limit number of times perturb
|         matrix in MatrixExpEigenPrepare
|
|  2008: Ian Holmes added bounds test to BalBak to prevent segfault, 11/19/2008
|  
*/  
  
#include <stdio.h>  
#include <stdlib.h>  
#include <math.h>  
#include <float.h>  
  
#include "linalg.h"  
#include "newmat.h"
#include "util/logfile.h"

#define PERTURB_TRIES 20
  
static int      EigenRG (int n, double **a, double *wr, double *wi, double **z, int *iv1, double *fv1);  
static void     Balanc (int n, double **a, int *pLow, int *pHigh, double *scale);  
static void     Exchange (int j, int k, int l, int m, int n, double **a, double *scale);  
static void     ElmHes (int n, int low, int high, double **a, int *intchg);  
static void     ElTran (int n, int low, int high, double **a, int *intchg, double **z);  
static int      Hqr2 (int n, int low, int high, double **h, double *wr, double *wi, double **z);  
static void     BalBak (int n, int low, int high, double *scale, int m, double **z);  
static void     CDiv (double ar, double ai, double br, double bi, double *cr, double *ci);  
static double   D_sign (double a, double b);  
  
  
  
#if !defined(MAX)  
#	define MAX(a,b)	(((a) > (b)) ? (a) : (b))  
#endif  
  
#if !defined(MIN)  
#	define MIN(a,b)	(((a) < (b)) ? (a) : (b))  
#endif  
  
  
  
#define FALSE 0  
#define TRUE 1  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|       MatrixExp  
|  
|       Commputes  exp(t * A),  with  t  real and  A  a real matrix, using the Pad\'e approximation.  
|       Returns NO_ERROR if successful, ERROR if not  (unlikely to happen for 6th order Pad\'e).  
|       This version is based on code from EXPOKIT - see top of file for references.  
|  
|       Matrices are in row-column format; aaA[i][j] = A_{ij}.  
|  
|       aaA                   the input matrix A.  (Not destroyed)  
|       iN                    matrix dimension  
|       iT                    t  
|       aaExp                 the output matrix.  (Must not be equal to aaA.)  
|       aaW1, aaW2, aaW3      workspace; matrices of dimension iN  
|       aW4                   workspace; integer array of length iN  
|  
|       (Gerton Lunter, 27/9/03)  
|  
*/  
  
int MatrixExp(double **aaA, int iN, double iT, double **aaExp, double** aaW1, double** aaW2, double** aaW3, int* aW4)  
{  
  const int iP = 6;         /* Order of approximation. */  
  static double iC[iP+1];  
  double iNorm,iFac,iD;  
  double** aaP = aaW1;  
  double** aaQ = aaW2;  
  double** aaR = aaW3;  
  double** aaTemp;  
  int iS;  
  int iSign;  
  int i,j,k,p;  
  
  /* Compute Pad\'e coefficients   */  
  
  iC[0] = 1.0;  
  for (i=0; i<iP; i++) {  
    iC[i+1] = iC[i] * (iP-i) / (double)((i+1)*(2*iP-i));  
  }  
  
  /* Scaling */  
  
  iNorm = 0.0;  
  for (i=0; i<iN; i++) {  
    iD = 0.0;  
    for (j=0; j<iN; j++) {  
      iD += fabs(aaA[i][j]);  
    }  
    iNorm = MAX(iNorm,iD);  
  }  
  iS = 0;  
  iFac = iT;  
  iNorm *= fabs(iT);  
  while (iNorm > 0.5) {  
    iS++;  
    iFac /= 2;  
    iNorm /= 2;  
  }  
  
  /* Horner evaluation of the irreducible fraction - see   */  
  /* ACM - Transactions On Mathematical Software, 24(1):130-156, 1998 */  
  
  for (i=0; i<iN; i++) {  
    for (j=0; j<iN; j++) {  
      aaP[i][j] = 0.0;  
      aaQ[i][j] = 0.0;  
      iD = 0.0;  
      for (k=0; k<iN; k++) {  
	iD += aaA[i][k]*aaA[k][j];  
      }  
      aaExp[i][j] = iD * iFac * iFac;  /* A^2 */  
    }  
    aaP[i][i] = iC[iP-1];      /* P */  
    aaQ[i][i] = iC[iP];        /* Q */  
  }  
  
  iSign = -1;  
  for (p=iP-2; p>=0; --p) {  
    if (iSign == -1) {  
      for (i=0; i<iN; i++) {  
	for (j=0; j<iN; j++) {  
	  if (i==j)  
	    iD = iC[p];  
	  else  
	    iD = 0.0;  
	  for (k=0; k<iN; k++) {  
	    iD += aaQ[i][k] * aaExp[k][j];   /* Q=Q*A^2 + c(p)*I */  
	  }  
	  aaR[i][j] = iD;  
	}  
      }  
      aaTemp = aaQ;  
      aaQ = aaR;  
      aaR = aaTemp;  
    } else {  
      for (i=0; i<iN; i++) {  
	for (j=0; j<iN; j++) {  
	  if (i==j)  
	    iD = iC[p];  
	  else  
	    iD = 0.0;  
	  for (k=0; k<iN; k++) {  
	    iD += aaP[i][k] * aaExp[k][j];   /* P=P*A^2 + c(p)*I */  
	  }  
	  aaR[i][j] = iD;  
	}  
      }  
      aaTemp = aaP;  
      aaP = aaR;  
      aaR = aaTemp;  
    }  
    iSign = -iSign;  
  }  
  if (iSign == -1) {  
    for (i=0; i<iN; i++) {  
      for (j=0; j<iN; j++) {  
	iD = 0.0;  
	for (k=0; k<iN; k++) {  
	  iD += aaQ[i][k] * aaA[k][j];   /* Q = Q*A */  
	}  
	aaR[i][j] = iD * iFac;  
      }  
    }  
    aaTemp = aaQ;  
    aaQ = aaR;  
    aaR = aaTemp;  
  } else {  
    for (i=0; i<iN; i++) {  
      for (j=0; j<iN; j++) {  
	iD = 0.0;  
	for (k=0; k<iN; k++) {  
	  iD += aaP[i][k] * aaA[k][j];   /* P = P*A */  
	}  
	aaR[i][j] = iD * iFac;  
      }  
    }  
    aaTemp = aaP;  
    aaP = aaR;  
    aaR = aaTemp;  
  }  
  for (i=0; i<iN; i++) {  
    for (j=0; j<iN; j++) {  
      aaQ[i][j] -= aaP[i][j];              /* Q = Q-P */  
    }  
  }  
  
  /* Compute Q^-1 in R using aaExp[0] and aW4 as workspace */  
  
  if (InvertMatrix(aaQ, iN, aaExp[0], aW4, aaR) != NO_ERROR)  
    return ERROR;  
  
  /* Compute iSign * (I + 2 P Q^-1) in Q */  
  
  for (i=0; i<iN; i++) {  
    for (j=0; j<iN; j++) {  
      if (i==j)  
	iD = 0.5;  
      else   
	iD = 0.0;  
      for (k=0; k<iN; k++) {  
	iD += aaP[i][k] * aaR[k][j];  
      }  
      aaQ[i][j] = iD * 2.0 * iSign;  
    }  
  }  
  
  /* Squaring */  
  
  for (p=0; p<iS; p++) {  
    for (i=0; i<iN; i++) {  
      for (j=0; j<iN; j++) {  
	iD = 0.0;  
	for (k=0; k<iN; k++) {  
	  iD += aaQ[i][k] * aaQ[k][j];  
	}  
	aaP[i][j] = iD;  
      }  
    }  
    aaTemp = aaP;  
    aaP = aaQ;  
    aaQ = aaTemp;  
  }  
  
  /* Copy result */  
  
  for (i=0; i<iN; i++) {  
    for (j=0; j<iN; j++) {  
      aaExp[i][j] = aaQ[i][j];  
    }  
  }  
  return NO_ERROR;  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|       MatrixExpEigenPrepare  
|  
|       From a general real matrix, computes eigenvalues and eigenvectors and the inverse of the  
|       eigenvector matrix, in preparation of matrix exponentiation by MatrixExpEigen.  
|  
|       Deals with non-semisimple matrices (which would require generalized eigenvectors) by  
|       perturbing the matrix and trying again.   
|  
|       This algorithm is unstable!  Nearly confluent and confluent eigenvalues lead to a loss of   
|       precision.  Both cases are detected and generate warnings (if VERBOSE is #define'd to 1).  
|       If you need stability, use the Pad\'e approximation; however, the method is fast if  exp(A t)    
|       is needed for fixed  A  and many  t  values.  
|  
|       (Gerton Lunter, 22/8/03)  
|  
*/  
  
int MatrixExpEigenPrepare(double **a, int n, double *v, double **u, double **uinv, int *cplx, double *work1, double *work2)  
     /*     **a = matrix to exponentiate (destroyed)  */  
     /*       n = order of matrix */  
     /*      *v = output vector (eigenvalues) */  
     /*     **u = output matrix (eigenvectors) */  
     /*  **uinv = output matrix (inverse of eigenvector matrix) */  
     /*    *cplx = output vector (booleans, true for real/imag eigenvalue pairs) */  
     /*  *work1 = work vector */  
     /*  *work2 = work vector */  
  
{   
  int i,j,tries;  
  char bad;  
  double iM;  
  
  /* copy initial matrix into a safe place */  
  for (i=0; i<n; i++)  
    for (j=0; j<n; j++)  
      uinv[i][j] = a[i][j];  
  
  /* compute eigenvectors and eigenvalues, perturbing matrix if necessary  */  
  tries = 0;
  do {  
  
    /* work1 receives imaginary parts; cplx and work2 are temporary */  
    if (EigenRG (n, a, v, work1, u, cplx, work2) != NO_ERROR)  
	  Throw(Runtime_error("Error return from EigenRG"));	// pk
  
    /* check whether all eigenvalues are simple (i.e. all eigenvectors nonzero) */  
    for (i=0, bad=FALSE; i<n && !bad; i++) {  
      bad = TRUE;  
      if (work1[i] != 0.0) {  
		/* complex eigenvalue */  
		for (j=0; j<n && bad; j++) {  
		  if (fabs(u[i][j]) > PERTURB || fabs(u[i+1][j]) > PERTURB) {  
			bad = FALSE;  
		  }
		}
		i+=1;   /* skip complex conjugate (G.A.L., 25/3/05) */
      } else {  
		/* real eigenvalue */  
		for (j=0; j<n && bad; j++) {  
		  if (fabs(u[i][j]) > PERTURB) {  
			bad = FALSE;  
		  }
		}
      }  
    }  
  
    if (bad) {  
      /* retrieve original matrix, and perturb it to get rid of non-simple eigenvalues */  
//      if (VERBOSE)  
//	printf("MatrixExpPrepare: Perturbing matrix.  Losing precision!\n");  
      CLOGERR << "MatrixExpPrepare: Perturbing matrix.  Losing precision!\n";

      for (i=0; i<n; i++) {  
	for (j=0; j<n; j++) {  
	  if (uinv[i][j] == 0.0)  
	    a[i][j] = (rand()%5) * PERTURB;  
	  else   
	    a[i][j] = uinv[i][j] * (1.0+((rand()%5)-2) * PERTURB);  
	}  
      }  
    }  
  
  } while (bad && ++tries < PERTURB_TRIES);  
  if (bad)
  {
	CTAG(3,RATE_EM_NEWMAT) << "MatrixExpPrepare:Reached maximum number of matrix perturbations\n";	// pk 
    return ERROR;  
  }
  
  /* copy eigenvector matrix u into a */  
  for (i=0; i<n; i++)  
    for (j=0; j<n; j++)  
      a[i][j] = u[i][j];  
  
  /* invert I,R matrix. work2 and cplx are temporary */  
  if (InvertMatrix(a, n, work2, cplx, uinv) != NO_ERROR)  
	Throw(Runtime_error("Error return from InvertMatrix"));		// pk
  
  /* check that matrix is sufficiently non-degenerate... */  
  if (VERBOSE) {  
    iM = 0.0;  
    for (i=0; i<n; i++) {  
      for (j=0; j<n; j++) {  
	iM = MAX(iM,fabs(uinv[i][j]));  
      }  
    }  
    if (iM > 1e6) {  
//      printf("MatrixExpPrepare: Warning: eigensystem is nearly degenerate, leading to a loss of precision.\n");  
	  CTAG(7,RATE_EM RATE_EM_PROGRESS) << "MatrixExpPrepare: Warning: eigensystem is nearly degenerate, leading to a loss of precision.\n";  
    }  
  }  
  
  /* sort real and imaginary parts */  
  for (i=0; i<n; i++) {  
    if (work1[i] != 0.0) {  
      v[i+1] = work1[i];  
      cplx[i++] = TRUE;  
    }  
    cplx[i] = FALSE;  
  }  
  
  return NO_ERROR;  
  
}  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|       MatrixExp  
|  
|       Returns exp(t*a), where a is the matrix given to MatrixExpPrepare  
|  
|       (Gerton Lunter, 22/8/03)  
|  
*/  
  
void MatrixExpEigen(double **mexp, int n, double *v, double **u, double **uinv, int *cplx, double t, double *work1)  
{  
  int i,j,k;  
  long double lam,s,c;  				// changed from float to double (pk) and then "long double" (Peter Schattner, IH)
  
  /* Multiply exp(t*diag(v)) with uinv, taking complex eigenvalues into account. */  
  for (i=0; i<n; i++) {  
    lam = exp( v[i]*t );  
    if (cplx[i]) {  
      /* Pair of complex conjugate eigenvalues.  Deal with two rows at once. */  
      s = sin( v[i+1]*t );  
      c = cos( v[i+1]*t );  
      for (j=0; j<n; j++) {  
	mexp[i][j] = lam * (c * uinv[i][j] + s * uinv[i+1][j]);  
	mexp[i+1][j] = lam * (-s * uinv[i][j] + c * uinv[i+1][j]);  
      }  
      i++;  
    } else {  
      /* Simple real eigenvalue */  
      for (j=0; j<n; j++) {  
	mexp[i][j] = lam * uinv[i][j];  
      }  
    }  
  }  
  
  /* multiply result with u */  
  for (j=0; j<n; j++) {  
    for (i=0; i<n; i++) {  
      lam = 0.0;  
      for (k=0; k<n; k++) {  
	lam += u[i][k]*mexp[k][j];  
      }  
      work1[i] = lam;  
    }  
    for (i=0; i<n; i++) {  
      mexp[i][j] = work1[i];  
    }  
  }  
}  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	Eigenvector  
|  
|	Finds eigenvector for eigenvalue (close to) 'ev', using inverse iteration (Numerical Recipes, 11.7)  
|	Input matrix 'a' is destroyed.  Returns eigenvector in 'col' (left eigenvector if 'right' == 0),  
|       normalized so that sum of absolute values is 1.    
|       Returns inverse of matrix (a - ev I) in workspace 'a_inv'.  
|       Returns ERROR if matrix is zero.  
|  
|       (Gerton Lunter, 29/8/03)  
|  
*/  
  
int Eigenvector (double **a, int n, double ev, double *col, int *indx, double **a_inv, char right)  
{  
  int i,j;  
  double d,norm;  
  for (i=0; i<n; i++)  
    a[i][i] -= ev;  
  if (InvertMatrix (a, n, col, indx, a_inv) == ERROR)  
    return ERROR;  
  for (i=0; i<n; i++)  
    col[i] = 1.0/(i+1);  
  do {  
    norm = 0;  
    for (i=0; i<n; i++) {  
      d = 0;  
      if (right) {  
	for (j=0; j<n; j++) {  
	  d += a_inv[i][j] * col[j];  
	}  
      } else {  
	for (j=0; j<n; j++) {  
	  d += a_inv[j][i] * col[j];  
	}  
      }  
      norm += fabs(d);  
      a[0][i] = d;  
    }  
    d = 0;  
    for (i=0; i<n; i++) {  
      a[0][i] /= norm;  
      if (fabs(fabs(col[i]) - fabs(a[0][i])) > d)  
	d = fabs(fabs(col[i]) - fabs(a[0][i]));  
      col[i] = a[0][i];  
    }  
  } while (d > EVECTORTOL);  
  return NO_ERROR;  
}  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	InvertMatrix  
|  
|	Invert matrix 'a' using LU-decomposition technique, storing inverse in 'a_inv'.  Matrix 'a'  
|	is destroyed.  Returns ERROR if matrix is singular, NO_ERROR otherwise.  
|  
*/  
  
int InvertMatrix (double **a, int n, double *col, int *indx, double **a_inv)  
	/*     **a = matrix represented as vector of row pointers      */  
	/*       n = order of matrix                                   */  
	/*    *col = work vector of size n                             */  
	/*   *indx = work vector of size n                             */  
	/* **a_inv = inverse of input matrix a (matrix a is destroyed) */  
{  
  int			rc, i, j;  
    
  rc = LUDecompose(a, n, col, indx, (double *)NULL);  
  if (rc == NO_ERROR) {  
    for (j = 0; j < n; j++) {  
      for (i = 0; i < n; i++) {  
	col[i] = 0.0;  
      }  
      col[j] = 1.0;  
      LUBackSubst(a, n, indx, col);  
      for (i = 0; i < n; i++) {  
	a_inv[i][j] = col[i];  
      }  
    }  
  }  
  return rc;  
}  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	LUDecompose  
|  
|	Replace matrix 'a' with its LU-decomposition.  Returns ERROR if matrix is singular, NO_ERROR  
|	otherwise.  
*/  
  
int LUDecompose (double **a, int n, double *vv, int *indx, double *pd)  
	/*   **a = the matrix whose LU-decomposition is wanted                    */  
	/*     n = order of a                                                     */  
	/*   *vv = work vector of size n (stores implicit scaling of each row)    */  
	/* *indx => row permutation according to partial pivoting sequence        */  
	/*   *pd => 1 if number of row interchanges was even, -1 if odd (NULL OK) */  
{  
	int			i, imax=0, j, k;  
	double		big, dum, sum, temp, d;  
  
	d = 1.0;  
	for (i = 0; i < n; i++)  
		{  
		big = 0.0;  
		for (j = 0; j < n; j++)  
			{  
			if ((temp = fabs(a[i][j])) > big)  
				big = temp;  
			}  
		if (big == 0.0)  
			{  
			printf("singular matrix in routine LUDecompose");  
			return ERROR;  
			}  
		vv[i] = 1.0 / big;  
		}  
	for (j = 0; j < n; j++)  
		{  
		for (i = 0; i < j; i++)  
			{  
			sum = a[i][j];  
			for (k = 0; k < i; k++)  
				sum -= a[i][k] * a[k][j];  
			a[i][j] = sum;  
			}  
		big = 0.0;  
		for (i = j; i < n; i++)  
			{  
			sum = a[i][j];  
			for (k = 0; k < j; k++)  
				sum -= a[i][k] * a[k][j];  
			a[i][j] = sum;  
			dum = vv[i] * fabs(sum);  
			if (dum >= big)  
				{  
				big = dum;  
				imax = i;  
				}  
			}  
		if (j != imax)  
			{  
			for (k = 0; k < n; k++)  
				{  
				dum = a[imax][k];  
				a[imax][k] = a[j][k];  
				a[j][k] = dum;  
				}	  
			d = -d;  
			vv[imax] = vv[j];  
			}  
		indx[j] = imax;  
		if (a[j][j] == 0.0)  
			a[j][j] = TINY;  
		if (j != n - 1)  
			{  
			dum = 1.0 / (a[j][j]);  
			for (i = j + 1; i < n; i++)  
				a[i][j] *= dum;  
			}  
		}  
  
	if (pd != NULL)  
		*pd = d;  
	return NO_ERROR;  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	LUBackSubst  
|  
|	Perform back-substition into LU-decomposed matrix in order to obtain inverse.  
*/  
  
void LUBackSubst (double **a, int n, int *indx, double *b)  
  
{  
	int			i, ip, j,  
				ii = -1;  
	double		sum;  
  
	for (i = 0; i < n; i++)  
		{  
		ip = indx[i];  
		sum = b[ip];  
		b[ip] = b[i];  
		if (ii >= 0)  
			{  
			for (j = ii; j <= i - 1; j++)  
				sum -= a[i][j] * b[j];  
			}  
		else if (sum != 0.0)  
			ii = i;  
		b[i] = sum;  
		}  
	for (i = n - 1; i >= 0; i--)  
		{  
		sum = b[i];  
		for (j = i + 1; j < n; j++)  
			sum -= a[i][j] * b[j];  
		b[i] = sum / a[i][i];  
		}  
}  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	EigenRG  
|  
|	This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine  
|	package (EISPACK) to find the eigenvalues of a real general matrix.  It was converted from  
|	Fortran to C by David Swofford.  
|  
|	ON INPUT:  
|  
|		n  is the order of the matrix 'a'  
|  
|		a  contains the real general matrix  
|  
|	ON OUTPUT:  
|  
|		wr  and  wi  contain the real and imaginary parts, respectively, of the eigenvalues.  
|		Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the  
|		positive imaginary part first.  
|  
|		z  contains the real and imaginary parts of the eigenvectors.  If the j-th eigenvalue is  
|		real, the j-th column of  z  contains its eigenvector.  If the j-th eigenvalue is complex  
|		with positive imaginary part, the j-th and (j+1)-th columns of  z  contain the real and  
|		imaginary parts of its eigenvector.  The conjugate of this vector is the eigenvector for  
|		the conjugate eigenvalue.  
|  
|		ierr  is an integer output variable set equal to an error completion code described in the  
|		documentation for Hqr and Hqr2.  The normal completion code is zero.  
|  
|		iv1  and  fv1  are temporary storage vectors of size n  
*/  
  
int EigenRG (int n, double **a, double *wr, double *wi, double **z, int *iv1, double *fv1)  
  
{  
	static int	is1, is2;  
	int			ierr;  
  
	Balanc (n, a, &is1, &is2, fv1);  
	ElmHes (n, is1, is2, a, iv1);  
	ElTran (n, is1, is2, a, iv1, z);  
	ierr = Hqr2 (n, is1, is2, a, wr, wi, z);  
	if (ierr == 0)  
		BalBak (n, is1, is2, fv1, n, z);  
  
	return ierr;  
}  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	Balanc  
|  
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments  
|	follow.  
|  
|	This subroutine is a translation of the algol procedure BALANCE, Num. Math. 13, 293-304(1969)  
|	by Parlett and Reinsch. Handbook for Auto. Comp., Vol. II-Linear Algebra, 315-326( 1971).  
|  
|	This subroutine balances a real matrix and isolates eigenvalues whenever possible.  
|  
|	ON INPUT:  
|  
|	   n is the order of the matrix.  
|  
|	   a contains the input matrix to be balanced.  
|  
|	ON OUTPUT:  
|  
|	   a contains the balanced matrix.  
|  
|	   low and high are two integers such that a(i,j) is equal to zero if  
|	      (1) i is greater than j and  
|	      (2) j=1,...,low-1 or i=high+1,...,n.  
|  
|	   scale contains information determining the permutations and scaling factors used.  
|  
|	Suppose that the principal submatrix in rows low through high has been balanced, that p(j)  
|	denotes the index interchanged with j during the permutation step, and that the elements of the  
|	diagonal matrix used are denoted by d(i,j).  Then  
|	   scale(j) = p(j),    for j = 1,...,low-1  
|	            = d(j,j),      j = low,...,high  
|	            = p(j)         j = high+1,...,n.  
|	The order in which the interchanges are made is n to high+1,  then 1 to low-1.  
|  
|	Note that 1 is returned for high if high is zero formally.  
*/  
  
void Balanc (int n, double **a, int *pLow, int *pHigh, double *scale)  
  
{  
	double		c, f, g, r, s, b2;  
	int			i, j, k, l, m, noconv;  
  
	b2 = FLT_RADIX * FLT_RADIX;  
	k = 0;  
	l = n - 1;  
   	  
	/* search for rows isolating an eigenvalue and push them down */  
  
	for (j = l; j >= 0; j--)  
		{  
		for (i = 0; i <= l; i++)  
			{  
			if (i != j)  
				{  
				if (a[j][i] != 0.0)  
					goto next_j1;  
				}  
			}  
#		if 0  /* bug that dave caught */  
		m = l;  
		Exchange(j, k, l, m, n, a, scale);  
		if (l < 0)  
			goto leave;  
		else  
			j = --l;  
#		else  
		m = l;  
		Exchange(j, k, l, m, n, a, scale);  
		if (--l < 0)  
			goto leave;  
#		endif  
		  
		next_j1:  
			;  
		}  
  
	/* search for columns isolating an eigenvalue and push them left */  
  
	for (j = k; j <= l; j++)  
		{  
		for (i = k; i <= l; i++)  
			{  
			if (i != j)  
				{  
				if (a[i][j] != 0.0)  
					goto next_j;  
				}  
			}  
  
		m = k;  
		Exchange(j, k, l, m, n, a, scale);  
		k++;  
  
		next_j:  
			;  
		}  
  
	/* now balance the submatrix in rows k to l */  
	for (i = k; i <= l; i++)  
		scale[i] = 1.0;  
  
	/* iterative loop for norm reduction */  
	  
	do	{  
		noconv = FALSE;  
	  
		for (i = k; i <= l; i++)  
			{  
			c = 0.0;  
			r = 0.0;  
		  
			for (j = k; j <= l; j++)  
				{  
				if (j != i)  
					{  
					c += fabs(a[j][i]);  
					r += fabs(a[i][j]);  
					}  
				}  
			/* guard against zero c or r due to underflow */  
			if ((c != 0.0) && (r != 0.0))  
				{  
				g = r / FLT_RADIX;  
				f = 1.0;  
				s = c + r;  
	  
				while (c < g)  
					{  
					f *= FLT_RADIX;  
					c *= b2;  
					}  
	  
				g = r * FLT_RADIX;  
	  
				while (c >= g)  
					{  
					f /= FLT_RADIX;  
					c /= b2;  
					}  
		  
				/* now balance */  
	  
				if ((c + r) / f < s * .95)  
					{  
					g = 1. / f;  
					scale[i] *= f;  
					noconv = TRUE;				  
					for (j = k; j < n; j++)  
						a[i][j] *= g;  
					for (j = 0; j <= l; j++)  
						a[j][i] *= f;  
					}  
				}  
			}	  
		}  
		while (noconv);  
  
	leave:  
		*pLow = k;  
		*pHigh = l;  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	Exchange  
|  
|	Support function for EISPACK routine Balanc.  
*/  
  
void Exchange (int j, int k, int l, int m, int n, double **a, double *scale)  
  
{  
	int			i;  
	double		f;  
  
	scale[m] = (double)j;  
	if (j != m)  
		{  
		for (i = 0; i <= l; i++)  
			{  
			f = a[i][j];  
			a[i][j] = a[i][m];  
			a[i][m] = f;  
			}	  
		for (i = k; i < n; i++)  
			{  
			f = a[j][i];  
			a[j][i] = a[m][i];  
			a[m][i] = f;  
			}  
		}  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	ElmHes  
|  
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments  
|	follow.  
|  
|	This subroutine is a translation of the algol procedure ELMHES, Num. Math. 12, 349-368(1968) by  
|	Martin and Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 339-358 (1971).  
|  
|	Given a real general matrix, this subroutine reduces a submatrix situated in rows and columns  
|	low through high to upper Hessenberg form by stabilized elementary similarity transformations.  
|  
|	ON INPUT:  
|  
|		n is the order of the matrix.  
|  
|		low and high are integers determined by the balancing subroutine BALANC.  If BALANC has not  
|		been used, set low=1, high=n.  
|  
|		a contains the input matrix.  
|  
|	ON OUTPUT:  
|  
|		a contains the Hessenberg matrix.  The multipliers which were used in the reduction are  
|		stored in the remaining triangle under the Hessenberg matrix.  
|  
|		int contains information on the rows and columns interchanged in the reduction.  Only  
|		elements low through high are used.  
*/  
  
void ElmHes (int n, int low, int high, double **a, int *intchg)  
  
{  
	int			i, j, m;  
	double		x, y;  
	int			la, mm1, kp1, mp1;  
	  
	la = high - 1;  
	kp1 = low + 1;  
	if (la < kp1)  
		return;  
  
	for (m = kp1; m <= la; m++)  
		{  
		mm1 = m - 1;  
		x = 0.0;  
		i = m;  
	  
		for (j = m; j <= high; j++)  
			{  
			if (fabs(a[j][mm1]) > fabs(x))  
				{  
				x = a[j][mm1];  
				i = j;  
				}  
			}  
	  
		intchg[m] = i;  
		if (i != m)  
			{  
			/* interchange rows and columns of a */  
			for (j = mm1; j < n; j++)  
				{  
				y = a[i][j];  
				a[i][j] = a[m][j];  
				a[m][j] = y;  
				}  
			for (j = 0; j <= high; j++)  
				{  
				y = a[j][i];  
				a[j][i] = a[j][m];  
				a[j][m] = y;  
				}  
			}  
  
		if (x != 0.0)  
			{  
			mp1 = m + 1;  
		  
			for (i = mp1; i <= high; i++)  
				{  
				y = a[i][mm1];  
				if (y != 0.0)  
					{  
					y /= x;  
					a[i][mm1] = y;  
					for (j = m; j < n; j++)  
						a[i][j] -= y * a[m][j];  
					for (j = 0; j <= high; j++)  
						a[j][m] += y * a[j][i];  
					}  
				}  
			}  
		}  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	ElTran  
|  
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments  
|	follow.  
|  
|	This subroutine is a translation of the algol procedure ELMTRANS,  Num. Math. 16, 181-204 (1970)  
|	by Peters and Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 372-395 (1971).  
|  
|	This subroutine accumulates the stabilized elementary similarity transformations used in the  
|	reduction of a  real general matrix to upper Hessenberg form by  ElmHes.  
|  
|	ON INPUT:  
|  
|		n is the order of the matrix.  
|  
|		low and high are integers determined by the balancing subroutine Balanc.  if  Balanc has  
|		not been used, set low=1, high=n.  
|  
|		a contains the multipliers which were used in the reduction by ElmHes in its lower triangle  
|		below the subdiagonal.  
|  
|		intchg contains information on the rows and columns interchanged in the reduction by ElmHes.  
|		Only elements low through high are used.  
|  
|	ON OUTPUT:  
|  
|	   z contains the transformation matrix produced in the reduction by ElmHes.  
*/  
  
void ElTran (int n, int low, int high, double **a, int *intchg, double **z)  
  
{  
	int			i, j, mp;  
  
	/* initialize z to identity matrix */  
	for (j = 0; j < n; j++)  
		{  
		for (i = 0; i < n; i++)  
			z[i][j] = 0.0;  
		z[j][j] = 1.0;  
		}  
  
	for (mp = high - 1; mp >= low + 1; mp--)  
		{  
		for (i = mp + 1; i <= high; i++)  
			z[i][mp] = a[i][mp-1];  
	  
		i = intchg[mp];  
		if (i != mp)   
			{  
			for (j = mp; j <= high; j++)  
				{  
				z[mp][j] = z[i][j];  
				z[i][j] = 0.0;  
				}  
			z[i][mp] = 1.0;  
			}  
		}  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	Hqr2  
|  
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments  
|	follow.  
|  
|	This subroutine is a translation of the algol procedure HQR2, Num. Math. 16, 181-204 (1970) by  
|	Peters and Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 372-395 (1971).  
|  
|	This subroutine finds the eigenvalues and eigenvectors of a real upper Hessenberg matrix by  
|	the QR method.  The eigenvectors of a real general matrix can also be found if ElmHes and  
|	ElTran or OrtHes  and  OrTran  have been used to reduce this general matrix to Hessenberg form  
|	and to accumulate the similarity transformations.  
|  
|	ON INPUT:  
|  
|		n is the order of the matrix  
|  
|		low and high are integers determined by the balancing subroutine Balanc.  If Balanc has not  
|		been used, set low=0, high=n-1.  
|  
|		h contains the upper Hessenberg matrix  
|  
|		z contains the transformation matrix produced by ElTran after the reduction by ElmHes, or  
|		by OrTran after the reduction by OrtHes, if performed.  If the eigenvectors of the  
|		Hessenberg matrix are desired, z must contain the identity matrix.  
|  
|	ON OUTPUT:  
|  
|		h has been destroyed  
|  
|		wr and wi contain the real and imaginary parts, respectively, of the eigenvalues.  The  
|		eigenvalues are unordered except that complex conjugate pairs of values appear consecutively  
|		with the eigenvalue having the positive imaginary part first.  If an error exit is made, the  
|		eigenvalues should be correct for indices ierr,...,n-1.  
|  
|		z contains the real and imaginary parts of the eigenvectors.   If the i-th eigenvalue is  
|		real, the i-th column of z contains its eigenvector.  If the i-th eigenvalue is complex with  
|		positive imaginary part, the i-th and (i+1)-th columns of z contain the real and imaginary  
|		parts of its eigenvector.  The eigenvectors are unnormalized.  If an error exit is made,  
|		none of the eigenvectors has been found.   
|  
|		Return value is set to:  
|			zero	for normal return,  
|			j		if the limit of 30*n iterations is exhausted while the j-th eigenvalue is  
|					being sought.  
|  
|	Calls CDiv for complex division.  
*/  
  
int Hqr2 (int n, int low, int high, double **h, double *wr, double *wi, double **z)  
  
{  
	int			i, j, k, l, m, na, en, notlas, mp2, itn, its, enm2, twoRoots;  
	double		norm, p=0, q=0, r=0, s=0, t, w=0, x, y=0, ra, sa, vi, vr, zz=0, tst1, tst2;  
  
	/* store roots isolated by Balanc and compute matrix norm */  
	norm = 0.0;  
	k = 0;  
	for (i = 0; i < n; i++)  
		{  
		for (j = k; j < n; j++)  
			norm += fabs(h[i][j]);  
  
		k = i;  
		if ((i < low) || (i > high))  
			{  
			wr[i] = h[i][i];  
			wi[i] = 0.0;  
			}  
		}  
  
	en = high;  
	t = 0.0;  
	itn = n * 30;  
  
	/* search for next eigenvalues */  
  
	while (en >= low)  
		{  
		its = 0;  
		na = en - 1;  
		enm2 = na - 1;  
		twoRoots = FALSE;  
  
		/* look for single small sub-diagonal element */  
		for (;;)  
			{  
			for (l = en; l > low; l--)  
				{  
				s = fabs(h[l-1][l-1]) + fabs(h[l][l]);  
				if (s == 0.0)  
					s = norm;  
				tst1 = s;  
				tst2 = tst1 + fabs(h[l][l-1]);  
				if (tst2 == tst1)  
					break;  
				}  
	  
			/* form shift */  
		  
			x = h[en][en];  
			if (l == en)  
				break;  
			y = h[na][na];  
			w = h[en][na] * h[na][en];  
			if (l == na)  
				{  
				twoRoots = TRUE;  
				break;  
				}  
  
			if (itn == 0)  
				{  
				/* set error -- all eigenvalues have not converged after 30*n iterations */  
				return en;  
				}  
			if ((its == 10) || (its == 20))  
				{  
				/* form exceptional shift */  
				t += x;  
			  
				for (i = low; i <= en; i++)  
					h[i][i] -= x;  
			  
				s = fabs(h[en][na]) + fabs(h[na][enm2]);  
				x = s * 0.75;  
				y = x;  
				w = s * -0.4375 * s;  
				}  
	  
			its++;  
			--itn;  
	  
			/* look for two consecutive small sub-diagonal elements */  
			for (m = enm2; m >= l; m--)  
				{  
				zz = h[m][m];  
				r = x - zz;  
				s = y - zz;  
				p = (r * s - w) / h[m+1][m] + h[m][m+1];  
				q = h[m+1][m+1] - zz - r - s;  
				r = h[m+2][m+1];  
				s = fabs(p) + fabs(q) + fabs(r);  
				p /= s;  
				q /= s;  
				r /= s;  
				if (m == l)  
					break;  
				tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));  
				tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));  
				if (tst2 == tst1)  
					break;  
				}  
		  
			mp2 = m + 2;  
			for (i = mp2; i <= en; i++)  
				{  
				h[i][i-2] = 0.0;  
				if (i != mp2)  
					h[i][i-3] = 0.0;  
				}  
	  
			/* double qr step involving rows l to en and columns m to en */  
			for (k = m; k <= na; k++)  
				{  
				notlas = (k != na);  
				if (k != m)  
					{  
					p = h[k][k-1];  
					q = h[k+1][k-1];  
					r = 0.0;  
					if (notlas)  
						r = h[k+2][k-1];  
					x = fabs(p) + fabs(q) + fabs(r);  
					if (x == 0.0)  
						continue;  
					p /= x;  
					q /= x;  
					r /= x;  
					}  
	  
				s = D_sign(sqrt(p*p + q*q + r*r), p);  
				if (k != m)  
					h[k][k-1] = -s * x;  
				else if (l != m)  
					h[k][k-1] = -h[k][k-1];  
				p += s;  
				x = p / s;  
				y = q / s;  
				zz = r / s;  
				q /= p;  
				r /= p;  
				if (!notlas)  
					{  
					/* row modification */  
					for (j = k; j < n; j++)  
						{  
						p = h[k][j] + q * h[k+1][j];  
						h[k][j] -= p * x;  
						h[k+1][j] -= p * y;  
						}   
				  
					j = MIN(en, k + 3);  
					/* column modification */  
					for (i = 0; i <= j; i++)  
						{  
						p = x * h[i][k] + y * h[i][k+1];  
						h[i][k] -= p;  
						h[i][k+1] -= p * q;  
						}  
					/* accumulate transformations */  
					for (i = low; i <= high; i++)  
						{  
						p = x * z[i][k] + y * z[i][k+1];  
						z[i][k] -= p;  
						z[i][k+1] -= p * q;  
						}  
					}  
				else  
					{  
					/* row modification */  
					for (j = k; j < n; j++)  
						{  
						p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];  
						h[k][j] -= p * x;  
						h[k+1][j] -= p * y;  
						h[k+2][j] -= p * zz;  
						}  
				  
					j = MIN(en, k + 3);  
					/* column modification */  
					for (i = 0; i <= j; i++)  
						{  
						p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];  
						h[i][k] -= p;  
						h[i][k+1] -= p * q;  
						h[i][k+2] -= p * r;  
						}  
					/* accumulate transformations */  
					for (i = low; i <= high; i++)  
						{  
						p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];  
						z[i][k] -= p;  
						z[i][k+1] -= p * q;  
						z[i][k+2] -= p * r;  
						}  
					}  
				}  
			}  
  
		if (twoRoots)  
			{  
			/* two roots found */  
			p = (y - x) / 2.0;  
			q = p * p + w;  
			zz = sqrt(fabs(q));  
			h[en][en] = x + t;  
			x = h[en][en];  
			h[na][na] = y + t;  
			/* DLS 28aug96: Changed "0.0" to "-1e-12" below.  Roundoff errors can cause this value  
			                to dip ever-so-slightly below zero even when eigenvalue is not complex.  
			*/  
			if (q >= -1e-12)  
				{  
				/* real pair */  
				zz = p + D_sign(zz, p);  
				wr[na] = x + zz;  
				wr[en] = wr[na];  
				if (zz != 0.0)  
					wr[en] = x - w/zz;  
				wi[na] = 0.0;  
				wi[en] = 0.0;  
				x = h[en][na];  
				s = fabs(x) + fabs(zz);  
				p = x / s;  
				q = zz / s;  
				r = sqrt(p*p + q*q);  
				p /= r;  
				q /= r;  
				/* row modification */  
				for (j = na; j < n; j++)  
					{  
					zz = h[na][j];  
					h[na][j] = q * zz + p * h[en][j];  
					h[en][j] = q * h[en][j] - p * zz;  
					}  
				/* column modification */  
				for (i = 0; i <= en; i++)  
					{  
					zz = h[i][na];  
					h[i][na] = q * zz + p * h[i][en];  
					h[i][en] = q * h[i][en] - p * zz;  
					}  
				/* accumulate transformations */  
				for (i = low; i <= high; i++)  
					{  
					zz = z[i][na];  
					z[i][na] = q * zz + p * z[i][en];  
					z[i][en] = q * z[i][en] - p * zz;  
					}  
				}  
			else  
				{  
				/* complex pair */  
				wr[na] = x + p;  
				wr[en] = x + p;  
				wi[na] = zz;  
				wi[en] = -zz;  
				}  
			en = enm2;  
			}  
		else  
			{  
			/* one root found */  
			h[en][en] = x + t;  
			wr[en] = h[en][en];  
			wi[en] = 0.0;  
			en = na;  
			}  
		}  
	  
	/* All roots found.  Backsubstitute to find vectors of upper triangular form */  
  
	if (norm == 0.0)  
		return 0;  
  
	for (en = n - 1; en >= 0; en--)  
		{  
		p = wr[en];  
		q = wi[en];  
		na = en - 1;  
		/* DLS 28aug96: Changed "0.0" to -1e-12 below (see comment above) */  
		if (q < -1e-12)  
			{  
			/* complex vector */  
			m = na;  
			/* last vector component chosen imaginary so that eigenvector matrix is triangular */  
			if (fabs(h[en][na]) > fabs(h[na][en]))  
				{  
				h[na][na] = q / h[en][na];  
				h[na][en] = -(h[en][en] - p) / h[en][na];  
				}  
			else  
				CDiv(0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);  
  
			h[en][na] = 0.0;  
			h[en][en] = 1.0;  
			enm2 = na - 1;  
			if (enm2 >= 0)  
				{  
				for (i = enm2; i >= 0; i--)  
					{  
					w = h[i][i] - p;  
					ra = 0.0;  
					sa = 0.0;  
			  
					for (j = m; j <= en; j++)  
						{  
						ra += h[i][j] * h[j][na];  
						sa += h[i][j] * h[j][en];  
						}  
			  
					if (wi[i] < 0.0)  
						{  
						zz = w;  
						r = ra;  
						s = sa;  
						}  
					else  
						{  
						m = i;  
						if (wi[i] == 0.0)  
							CDiv(-ra, -sa, w, q, &h[i][na], &h[i][en]);  
						else  
							{  
							/* solve complex equations */  
							x = h[i][i+1];  
							y = h[i+1][i];  
							vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;  
							vi = (wr[i] - p) * 2.0 * q;  
							if ((vr == 0.0) && (vi == 0.0))  
								{  
								tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));  
								vr = tst1;  
								do	{  
									vr *= .01;  
									tst2 = tst1 + vr;  
									}  
									while (tst2 > tst1);  
								}  
							CDiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);  
							if (fabs(x) > fabs(zz) + fabs(q))  
								{  
								h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;  
								h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;  
								}  
							else  
								CDiv(-r - y * h[i][na], -s - y * h[i][en], zz, q, &h[i+1][na], &h[i+1][en]);  
							}  
				  
						/* overflow control */  
						tst1 = fabs(h[i][na]);  
						tst2 = fabs(h[i][en]);  
						t = MAX(tst1, tst2);  
						if (t != 0.0)  
							{  
							tst1 = t;  
							tst2 = tst1 + 1.0 / tst1;  
							if (tst2 <= tst1)  
								{  
								for (j = i; j <= en; j++)  
									{  
									h[j][na] /= t;  
									h[j][en] /= t;  
									}  
								}  
							}  
						}  
					}  
				}  
			/* end complex vector */  
			}  
		else if (q == 0.0)  
			{  
			/* real vector */  
			m = en;  
			h[en][en] = 1.0;  
			if (na >= 0)  
				{  
				for (i = na; i >= 0; i--)  
					{  
					w = h[i][i] - p;  
					r = 0.0;  
			  
					for (j = m; j <= en; j++)  
						r += h[i][j] * h[j][en];  
			  
					if (wi[i] < 0.0)  
						{  
						zz = w;  
						s = r;  
						continue;  
						}  
					else  
						{  
						m = i;  
						if (wi[i] == 0.0)  
							{  
							t = w;  
							if (t == 0.0)  
								{  
								tst1 = norm;  
								t = tst1;  
								do	{  
									t *= .01;  
									tst2 = norm + t;  
									}  
									while (tst2 > tst1);  
								}			  
							h[i][en] = -r / t;  
							}  
						else  
							{  
							/* solve real equations */  
							x = h[i][i+1];  
							y = h[i+1][i];  
							q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];  
							t = (x * s - zz * r) / q;  
							h[i][en] = t;  
							if (fabs(x) > fabs(zz))  
								h[i+1][en] = (-r - w * t) / x;  
							else  
								h[i+1][en] = (-s - y * t) / zz;  
							}  
				  
						/* overflow control */  
						t = fabs(h[i][en]);  
						if (t != 0.0)  
							{  
							tst1 = t;  
							tst2 = tst1 + 1. / tst1;  
							if (tst2 <= tst1)  
								{  
								for (j = i; j <= en; j++)  
									h[j][en] /= t;  
								}  
							}  
						}  
					}  
				}  
			/* end real vector */  
			}  
		}  
	/* end back substitution */  
	  
	/* vectors of isolated roots */  
	for (i = 0; i < n; i++)  
		{  
		if ((i < low) || (i > high))  
			{  
			for (j = i; j < n; j++)  
				z[i][j] = h[i][j];  
			}  
		}  
  
	/* multiply by transformation matrix to give vectors of original full matrix */  
	for (j = n - 1; j >= low; j--)  
		{  
		m = MIN(j, high);  
		for (i = low; i <= high; i++)  
			{  
			zz = 0.0;  
			for (k = low; k <= m; k++)  
				zz += z[i][k] * h[k][j];  
			z[i][j] = zz;  
			}  
		}  
  
	return 0;  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	BalBak  
|  
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments  
|	follow.  
|  
|	This subroutine is a translation of the algol procedure BALBAK, Num. Math. 13, 293-304 (1969)  
|	by Parlett and Reinsch.  Handbook for Auto. Comp., vol. II-Linear Algebra, 315-326 (1971).  
|  
|	This subroutine forms the eigenvectors of a real general matrix by back transforming those of  
|	the corresponding balanced matrix determined by  Balanc.  
|  
|	ON INPUT:  
|  
|		n is the order of the matrix.  
|  
|		low and high are integers determined by Balanc.  
|  
|		scale contains information determining the permutations and scaling factors used by Balanc.  
|  
|		m is the number of columns of z to be back transformed.  
|  
|		z contains the real and imaginary parts of the eigenvectors to be back transformed in its  
|		first m columns.  
|  
|	ON OUTPUT:  
|  
|		z contains the real and imaginary parts of the transformed eigenvectors in its first m  
|		columns.  
*/  
  
void BalBak (int n, int low, int high, double *scale, int m, double **z)  
  
{  
	int			i, j, k, ii;  
	double		s;  
  
	if (m != 0)  
		{  
		if (high != low)  
			{  
			for (i = low; i <= high; i++)  
				{  
				s = scale[i];	/* left hand eigenvectors are back transformed if this statement is  
								   replaced by  s = 1.0/scale[i] */  
				for (j = 0; j < m; j++)  
					z[i][j] *= s;  
				}  
			}  
		for (ii = 0; ii < n; ii++)  
			{  
			i = ii;  
			if ((i < low) || (i > high))  
				{  
				if (i < low)  
					i = low - ii;  
				k = (int)scale[i];  
				if (k != i && k >= 0 && k < m)   /* (k>=0 && k<m) test added to prevent bounds-overrun segfault, IH 11/19/2008 */
					{  
					for (j = 0; j < m; j++)  
						{  
						s = z[i][j];  
						z[i][j] = z[k][j];  
						z[k][j] = s;  
						}  
					}  
				}  
			}  
		}  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	CDiv  
|  
|	Complex division, (cr,ci) = (ar,ai)/(br,bi)  
*/  
  
void CDiv (double ar, double ai, double br, double bi, double *cr, double *ci)  
  
{  
	double		s, ais, bis, ars, brs;  
  
	s = fabs(br) + fabs(bi);  
	ars = ar / s;  
	ais = ai / s;  
	brs = br / s;  
	bis = bi / s;  
	s = brs*brs + bis*bis;  
	*cr = (ars*brs + ais*bis) / s;  
	*ci = (ais*brs - ars*bis) / s;  
}  
  
  
  
  
  
/*--------------------------------------------------------------------------------------------------  
|  
|	D_sign  
|  
|	"Sign" function.  
*/  
  
double D_sign (double a, double b)  
  
{  
	double		x;  
  
	x = (a >= 0 ? a : -a);  
	return (b >= 0 ? x : -x);  
}
