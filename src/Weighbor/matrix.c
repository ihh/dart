/*
 *
 * WEIGHBOR 0.2.13
 *
 * Copyright (c) 1997-1998 
 *    Los Alamos National Laboratories
 *
 * This file is part of the WEIGHBOR program that can be used free of
 * charge in academic research and teaching. Any commercial use of
 * this software requires prior permission and possibly a license. For
 * commercial use contact: W. Bruno at billb@t10.lanl.gov Los Alamos
 * National Laboratories makes no representations about the
 * suitability of this software for any purpose. It is provided 
 * "as is" without express or implied warranty
 *	
 * For those using this program for research please use the following 
 * citation for this program:
 *
 *    Bruno, W. J., Socci, N. D. and Halpern, A. L., ``Weighted 
 *        Neighbor Joining: A Fast Approximation to Maximum-Likelihood
 *        Phylogeny Reconstruction,'' Submitted to Molecular Biology
 *        and Evolution, (1998).
 *
 * 
 */

/*
  
  \input cnoweb

  Program: weighbor                             \hfill\break
  File: matrix.c                                \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 4-Feb-97                             \hfill\break
  
  \title{matrix.c}
  \job{Weighbor}
  
  \synopsis{Routines for matrix and vector creation and manipulation}
  

  */

/* $Id: matrix.c,v 1.3 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: matrix.c,v $
  Revision 1.3  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.2  2003/12/15 12:28:06  ihh
  Fixed so gcc can compile with -Wall.

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.13  1998-06-15 19:56:17-04  nds
  First Release Version!

  Revision 1.12  1998-06-09 18:33:39-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.11  1997-10-14 22:16:41-04  nds
  Added \|RECT| function with rectifies its argument; ie, max(0,x)

  Revision 1.10  1997-09-18 16:56:06-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.3  1997-08-19 13:18:07-04  nick
  Bug in printVector routine. Printing the last element
  incorrectly

  Revision 1.2  1997/08/15 03:44:28  nds
  Moved the macros SQR, DMAX, DMIN here and made them fuctions. The
  functions declarations are in weighbor.h

  Revision 1.1  1997-08-08 13:59:06-04  nds
  Initial revision


  */


/*

  \section{Introduction}

  This file contains functions to manipulate the basic vector and
  matrix data types.

  */

#include <stdio.h>
#include <stdlib.h>
#include "Weighbor/matrix.h"


/*
  
  \section{Constructors and Destructors}

  Here are the routines for allocating and deallocating matrices and
  vectors. Note all matrices are square; all matrices and vectors have
  indices from $0\ldots N-1$ and the matrices are stored in the C
  convention (row first). 

  */


/*

  \subsection{matrix}

  A \|MatrixT| is actually an array of pointers into an $N\times N$
  block of memory. The pointers, point to the beginning of each
  row. It waste some space (ie $N$ pointers) but has some
  conveniences. 

  Note the matrix idicies go from $0\ldots N-1$; ie the normal \"C"
  convention.

  */

MatrixT
matrix(int N)
{

  int i,j;
  MatrixT m;

  /* Allocate the row pointers */

  m = (MatrixT)malloc((size_t)((N)*sizeof(double*)));
  if(!m) {
    perror("weighbor::matrix::matrix(1)");
    exit(1);
  }

  /* Allocate the $N\times N$ block of memory and set pointers*/

  m[0] = (double *)malloc((size_t)((N*N)*sizeof(double)));
  if(!m[0]) {
    perror("weighbor::matrix::matrix(2)");
    exit(1);
  }
  
  for(i=1;i<N;++i)
    m[i] = m[i-1]+N;
  
  /* Initial to zero to save any grief */

  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      m[i][j] = 0.0;

  return(m);

}

void
freeMatrix(MatrixT m)
{

  free(m[0]);
  free(m);

}


/*

  \subsection{vector}

  A \|VectorT| is just an pointer to a block of doubles. 

  */

VectorT
vector(int N)
{

  int i;
  VectorT v;

  v=(VectorT)malloc((size_t)((N)*sizeof(double)));
  if (!v) {
    perror("weighbor::matrix::vector");
    exit(1);
  }

  for(i=0;i<N;++i)
    v[i] = 0.0;

  return(v);

}  

void
freeVector(VectorT v)
{
  free(v);
}

/*

  \section{Basic Matrix Manipulations}

  Some basic function to manipulate matrices and vectors

  */

/*

  \subsection{setMM}

  Set ${\bf B} = {\bf A}$

  There are certainly better ways to do this, but this will suffice
  for now. Maybe better to rewrite every thing in C++ and use BLAS.

  */

void
setMM(int N, MatrixT A, MatrixT B)
{

  int i, j;

  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      B[i][j] = A[i][j];

}

/*

  \subsection{Delete column or row}

  Delete the $i^{\rm th}$ row and column of the matrix. We delete them
  by setting the off diagonal elments to zero and setting the diagonal
  element $A_{ii}=1$ this will keep the matrix non-singular. This is
  necessary if you are going to refactor the matrix with a Cholesky
  factorization. 

  */

void
delRowCol(int N, MatrixT A, int i)
{
  
  int j;
  
  for(j=0;j<N;++j)
      A[i][j] = A[j][i] = 0.0;

  A[i][i] = 1.0;

}

/*

  \section{Printing}

  These are some simple functions for printing matrices and vectors;
  mostly for debugging purposes. 

  */

void
printMatrix(int n, MatrixT m)
{
  int i, j;
  printf("\t[");
  for(i=0;i<n;++i) {
    for(j=0;j<n;++j)
      if(j<n-1)
	printf("%10.7f, ", m[i][j]);
      else if(i<n)
	printf("%10.7f; ", m[i][j]);
      else
	printf("%10.7f]", m[i][j]);
    printf("\n\t ");
  }

  printf("\n");

}

void
printVector(int n, VectorT v)
{

  int i;

  printf("\t[");
  for(i=0;i<n-1;++i)
    printf("%10.7f, ", v[i]);
  printf("%10.7f]\n\n", v[i]);

}

/*

  \section{Miscellaneous Functions}

  This are some miscellaneous functions that have nothing to do 
  with vectors or matrices but I was too lazy to start a new file so
  I am sticking them here. Also these functions use to be macros but it
  seems gcc has some sort of bug that causes an error with the \|DMIN| 
  
  ie \|z = DMIN(f(a),f(b))+DMIN(f(c),f(d))| evalutes to \|a+d| even if 
  \|a>b| and \|c>d|

  */

double SQR(double a)
{

  return(a*a);

}

double DMAX(double a, double b)
{

  return( a>b ? a : b );

}

double DMIN(double a, double b)
{

  return( a<b ? a : b);

}

double RECT(double a)
{
  return( a>0.0 ? a : 0.0);
}

/* \endc */ 


