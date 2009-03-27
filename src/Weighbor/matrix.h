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
  File: matrix.h                                \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 4-Feb-97                             \hfill\break
  
  \title{matrix.h}
  \job{Weighbor}
  
  \synopsis{Include file for various linear algebra data types: matrix
  and vector} 
  

  */

/* $Id: matrix.h,v 1.1.1.1 2003/11/10 19:07:59 ihh Exp $ */

/*
  
  $Log: matrix.h,v $
  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.13  1998-06-15 19:56:17-04  nds
  First Release Version!

  Revision 1.12  1998-06-09 18:33:40-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.11  1997-11-15 00:31:32-05  nds
  Changes for weighor version 0.2.7.1

  Revision 1.10  1997-09-18 16:56:22-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.1  1997-08-08 13:58:30-04  nds
  Initial revision

  
  */


#ifndef __MATRIX_H
#define __MATRIX_H

#include <math.h>

/*

  \section{Introduction}

  This file contains a rudimentary vector and square matrix datatype.

  */

/*
  
  \section{Data Types}

  Here are the basic data types: \|MatrixT| which is a 2D dimensional
  array of doubles. 
  It will have the \"C" index convention; ie $i=0\ldots N-1$. 
  The data will be stored in C (row-first)
  ordering. \|VectorT| is a 1D array of doubles again with the index
  going from 0 to $N-1$.

  I also include the ever useful BooleanT here.

  */

typedef double **MatrixT;
typedef double *VectorT;

typedef enum {False=0, True=-1} BooleanT;


/*

  \section{Function Definitions}

  \subsection{Allocation and Deallocation}

  This functions allocate \& deallocate matrices and vectors. 
  Note all the matrices used by weighbor are square so matrix takes
  just one argument.


  */

MatrixT matrix(int);
VectorT vector(int);

void freeVector(VectorT);
void freeMatrix(MatrixT);

/*
  
  \subsection{Linear Algebra Routines}

  Here are the function def's for any linear algebra routines that
  weighbor will used. They will simply be wrappers around other
  standard library routines
  */

void setMM(int, MatrixT, MatrixT);
void delRowCol(int, MatrixT, int);


/* \section{Printing} */

void printMatrix(int, MatrixT);
void printVector(int, VectorT);

/* \section{Miscellaneous Functions} */

/* 

   Complimentary error function: in source file \|calerf.c| Using our
   own code since the clib of gcc only computes to single precision
   even though it returns a double.

   */

double derf(double x); 
double derfc(double x); 
double derfcx(double x); 

#endif /* __MATRIX_H */

/* \endc */ 






