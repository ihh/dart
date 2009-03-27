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
  Version: II
  File: weighbor.h                              \hfill\break
  Author: N. D. Socci                           \hfill\break
  
  \title{weighbor.h}
  \job{Weighbor}
  
  \synopsis{\"Weigh"ted neigh\"bor"-joining.
This program will created a evolutionary tree based on distance
between sequences. The tree built attempts to minimize the distances
in it versus the actually distances inputted to the program. The
output is the tree built along with the branch lengths.}


*/

/* $Id: weighbor.h,v 1.3 2005/04/26 21:00:24 ihh Exp $ */

/*

  $Log: weighbor.h,v $
  Revision 1.3  2005/04/26 21:00:24  ihh
  tkfalign_tree_sampling

  Revision 1.2  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.1.1.1  2003/11/10 19:08:00  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 2.23  1998-06-15 19:56:18-04  nds
  First Release Version!

  Revision 2.22  1998-06-09 18:33:41-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 2.21  1998-06-04 18:35:22-04  nds
  Update for version 2.12
      Changed z score
      Change renormalization if z-score negative
      non-deterministic warning to stderr
      -Xn flag to change R(i,j,j') normalization

  Revision 2.20  1998-05-20 17:14:14-04  nsocci
  Initial release 1.0.0 alpha

  Revision 2.19  1998/05/20 20:00:04  nsocci
  Added -b option to change the value of B and changed the default
  setting of the -S option (now default is off, -S turn on barred
  values).

  Revision 2.18  1998/05/20 01:44:00  nsocci
  Added the extended tournament option (-e)

  Revision 2.17  1998/05/19 17:12:15  nsocci
  Added -S and -T options do control the use of the average (barred)
  quanitites in the calculation of R(i,j,j')

  Revision 2.16  1998/05/18 21:26:33  nsocci
  Added function \|sigma2tinv| which is used to calculate the $c_i$

  Revision 2.15  1998/05/01 20:57:20  nds
  Made numerous changes to log file and adjusted what was printed at
  various verbose levels. To print the new LLR value of R(q(i),i,j) the
  calcR funtion need to be called from buildTree which require the
  moving of the deltaB and associated matrices to global level. This
  should improve performance since now they are allocated only once per
  tree rather than a every iteration

  Revision 2.14  1998-04-29 18:47:09-04  nds
  Added a compiler constant DEFAULT_LENGTH to set the default sequence
  length

  Revision 2.13  1998-04-14 16:55:01-04  nds
  Added new flag to control which z score method to use.

  Revision 2.12  1998-04-08 19:23:44-04  nds
  $\Delta b$ optimization: calculation in $N^2$ time.

  Revision 2.11  1998-03-20 00:02:17-05  nds
  Added -r option to turn of the recalculation of $\Delta b$ and
  associated variables.

  Revision 2.10  1998-03-13 16:52:44-05  nds
  Added -v command line option
  BETA release version

  Revision 2.9  1998-02-06 16:04:52-05  nds
  Added the constant \|MINB| used in the calculation of the b's (eq
  0.13-0.15, 0.22-0.24)

  Revision 2.8  1998-01-15 15:23:17-05  nds
  Version 0.2.8.5.1

  Added a new noise function $\sigma^2_{ij;k}{}'$ equations 0.11 and
  0.12.

  Revision 2.7  1997-12-18 22:13:12-05  nds
  version 0.2.8.4
   + Check to make sure the winner of the tournament is the true winner
     print warning if not or tie.

   + -q option to turn of checking of q[q[i]]=i

   + New R(i,j,j') function (eq 0.11)

  Revision 2.6  1997-12-16 15:20:49-05  nds
  Added the checkQQI boolean variable and the -q command line switch to
  turn of the checking of q[q[i]]=i in selection the best pair.

  Revision 2.5  1997-12-10 18:52:30-05  nds
  Version 0.2.8.2

  Revision 2.4  1997-11-29 21:00:27-05  nds
  Version 0.2.8.2-; still did not implement equation 0.11
  and the q[i] check.

  Revision 2.3  1997-11-15 00:32:06-05  nds
  Changes for weighor version 0.2.7.1

  Revision 2.2  1997-10-17 19:22:30-04  nds
  Version 0.2.6: new renormalization code

  Revision 2.1  1997-10-14 22:16:41-04  nds
  Added \|RECT| function with rectifies its argument; ie, max(0,x)

  Revision 2.0  1997-10-09 19:23:47-04  nds
  Version 0.2.x initial revision

  Revision 1.12  1997-10-04 18:46:21-04  nds
  Added \|calc_b| declaration. Last 0.1.x version.

  Revision 1.11  1997-09-19 18:48:38-04  nds
  Added LOGFILE compile flag to switch the logging of extra info to the
  weighbor.out file on and off.

  Revision 1.10  1997-09-18 16:58:17-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.8  1997-09-11 22:26:51-04  nds
  Changed \|EPSILON| to 1e-9

  Revision 1.7  1997-08-23 17:19:14-04  nsocci
  Added compiler flags
          CHECK_QQOFI
    to control whether the code to check that $q(q(i))==i$, and
          USEBMIN
    to control whether we use $b_min_i$ for branch lengths or if not set
    then use $b_{i;j}$.

  Set or unset these flag in weighbor.h file

  Revision 1.6  1997/08/23 01:16:21  nsocci
  Renormalization code added. Alpha version still needs to be tested
  Note the sigma2_3 function's calling parameteres has been
  changed. Added two booleans to determine whether to renormalize or
  not.

  Revision 1.5  1997/08/20 21:25:02  nick

  */

#ifndef __WEIGHBOR_H
#define __WEIGHBOR_H

#include "Weighbor/matrix.h"
#include "Weighbor/tree.h"


#define VERSION_INFO "Weighbor ver 1.0.0 alpha [20-May-98]\n"

/*
  \section{Macros}
  */

/* 

   \subsection{Distance \|D(a,b)|}

   This macro accesses the distance array given indices into the node
   array.  
   
   Note it assumes that the distance matrix is called \|mD| and that the
   array of nodes is called \|nodes|

  */

#define D(a,b) (mD[nodes[(a)]->ind][nodes[(b)]->ind])

/*
  Same for the \|C| array $c_i$ (eq 0.39) 
  */

#define C(a) (vC[nodes[(a)]->ind])

/*
  This macros are for access to the \|mS|, \|mDelB| and \|mDel2B|
  arrays
  */

#define S(a,b)     (mS[nodes[(a)]->ind][nodes[(b)]->ind])
#define DelB(a,b)  (mDelB[nodes[(a)]->ind][nodes[(b)]->ind])
#define Del2B(a,b) (mDel2B[nodes[(a)]->ind][nodes[(b)]->ind])

/*
  
  \section{Global Varibles}

  The \|nodes| global varibles is the list of currently active nodes (taxa) 
  which still need to be joined.

  The matrix \|mD| is the distance matrix. The indices into this matrix 
  stored in each taxon's node varible. Ie to find the distance between 
  taxon $i$ and $j$ use: 

  \qquad\|mD[nodes[i]->ind][nodes[j]->ind]|

  The macro \|D(a,b)| is used for this purpose.

  */

/*

  We need a small positive number $\epsilon$ in cases where all the
  $\sigma$'s turn out to be zero. However, we can not make epsilon too
  small as the algorithm requires that the following be true: 
  $$
  (x+{1\over \epsilon})-{1\over\epsilon} = x
  $$
  when $x$ is of order 1. This means that $\epsilon$ must be greater
  then the machine precision and probably a lot greater. We set
  $\epsilon=10^{-9}$ which is about the square root of most double
  precision and also is roughly the minimal noise for sequence the
  size of the human genome

  */
  
/*

   These are the unrenormalized vaules for \|EPSILON| and \|MINB|; ie,
   before we renormalize the vaules for sequence length. \|MINB| is
   used in the calculation of the b's (eq 0.13-0.15 and 0.22-0.24)

   */

#define EPS_BARE (1e-9)
#define MINB_BARE (0.5)

#ifdef __WEIGHBOR_C   /* only declare globals once */

NodeT **nodes;    /* array of nodes left to be joined */
MatrixT mD;       /* Distance Matrix */
VectorT vC;       /* Renormalization vector */

/* 

   The following global matricies hold the old vaules of the $\Delta
   b$ and other related matrices. They are used in the new $N^3$
   updateing of this matrices by the function \|calcb|. They are
   stored in the format as the \|mD| matrix and are access via the
   \|S|, \|DelB| and \|Del2B| macros which use the node structures to
   get the right indicies.

   */

MatrixT mS;        /* Old value of $s_{ij}$ matrix (eq 0.9) */
MatrixT mDelB;     /* $\Delta b_{ij}$ (eq 0.8) */
MatrixT mDel2B;    /* $\Delta^2b_{ij}$ */

/*

  Make these varibles globals to make life easier so we do not
  constantly have to pass them around

  Also need them hanging around so we can calculate $R(q(i),i,j')$
  after we have decided which pair to join in \|build.c|

  */

MatrixT s;       /* $s_{ij}$ eq 0.9 */
MatrixT deltaB;  /* $\Delta b_{ij}$ eq 0.8 */
MatrixT delta2B; /* $\Delta^2 b_{ij}$ */
MatrixT oldDeltaB; /* Save orginal value of \|deltaB| */

/*

  These varible keep track of the last joined pair \|(ta,tb)| and the
  newest taxon \|tc|. They are needed in \|calcb| to update the
  $\Delta b$ values in constant time.

  N.B. these are really indicies into the \|nodes|. 

  */
  
int ta, tb, tc;

#define DEFAULT_LENGTH (500.0)
double L=DEFAULT_LENGTH;        /* Length of sequences */

/*
  
  Use a default sequence length of 500

  */

double EPSILON=(EPS_BARE)/500.0; 
double MINB=(MINB_BARE)/500.0; 

#define DEFAULT_B (4.0)  /* Default number of DNA base types used in the noise and
			    dissimiliarty functions */

double B=DEFAULT_B;

FILE *outfile;    /* File for auxillary output */

int printLevel=0; /* Controls the amount of output */

BooleanT checkQQI=True; /* Turns on checking of $q[q[i]]=i$ */
BooleanT recalcB=False;  /* controls whether $\Delta b$ is recalculated with 
			    the improved $\sigma$'s */

BooleanT oldZflag=False; /* Determines which algorithm is used to calculate
			    $z(i,j)$ */

BooleanT useSigmaBar=False; /* Controls whether the barred version of 
			      $\sigma^2_{i,j;j'}$ is used in the calculation 
			      of R(i,j,j') */

BooleanT useBarValues=True; /* Controls whether the barred version of 
			       $d_{ij}$, $\Delta b$ and $s$ are used in
			       the calculation of R(i,j,jp) */

BooleanT extendedTourn=False; /* Activate the extended tournament */

BooleanT warnFlag=False; /* use to limit the printing of the non-deterministic
			    warning to once per data file */

BooleanT n_Flag=False; /* Change the normalization of $R(i,j,j')$ */
BooleanT w_Flag=False; /* Controls the printing of warnings */
BooleanT x_Flag=False; /* Change normalization of $R(i,j)$ */

#ifdef HASH 
double *hash; /* sigma23 hash array */
int hashsize;
int total=0;
int hashcnt=0;
#endif

#else

extern NodeT **nodes;    /* array of nodes left to be joined */
extern MatrixT mD;        /* Distance Matrix */
extern VectorT vC;

extern MatrixT mS;        /* Old value of $s_{ij}$ matrix (eq 0.9) */
extern MatrixT mDelB;     /* $\Delta b_{ij}$ (eq 0.8) */
extern MatrixT mDel2B;    /* $\Delta^2b_{ij}$ */

extern MatrixT s;       /* $s_{ij}$ eq 0.9 */
extern MatrixT deltaB;  /* $\Delta b_{ij}$ eq 0.8 */
extern MatrixT delta2B; /* $\Delta^2 b_{ij}$ */
extern MatrixT oldDeltaB; /* Save orginal value of \|deltaB| */

extern int ta, tb, tc;

extern double L;
extern double EPSILON;
extern double MINB;

extern double B;

extern FILE *outfile;

extern int printLevel;

extern BooleanT checkQQI;
extern BooleanT recalcB;
extern BooleanT oldZflag; 
extern BooleanT useSigmaBar;
extern BooleanT useBarValues;
extern BooleanT extendedTourn;
extern BooleanT n_Flag;
extern BooleanT w_Flag;
extern BooleanT x_Flag;

extern BooleanT warnFlag;

#ifdef HASH
extern double *hash;
extern int hashsize;
extern int total;
extern int hashcnt;
#endif

#endif

/*
  
  \subsection{Compiler Node_pairives}

  Here are the flags to control the compliation of various options

  \|CHECK_QQOFI| controls whether the code to check that $q(q(i))=i$
  is used in \"build.c". Set flag to turn this code on.

  \|LOGFILE| controls whether extra info is dumped to the weighbor.out
  file

  */


#undef CHECK_QQOFI
#define LOGFILE 1

/*

  \section{function defs}

  */

void calc_b(int N, MatrixT b, MatrixT deltaB, MatrixT delta2B, MatrixT s);
void recalc_b(int N, MatrixT b, MatrixT deltaB, MatrixT delta2B, MatrixT s);
void calc_q(int N, int* q, VectorT R, MatrixT b, 
	    int* q2, VectorT LLR, VectorT Zscore);
double calcR(int N, int i, int j, int jp);
double calcR2(int N, int i, int j, int k);
void calc_z1(int N, VectorT z, int* q);
void calc_z2(int N, VectorT z, int* q, int* q2);
double calcZ2(int N, int i, int q, int q2);
double calcPhi(int N, int i, int ip);

double sigma2_3(int, int, int);
double sigma_na(double, double);
double sigma2_3p(int, int, int, double);
double sigma2t(double);
double sigma2tinv(double);
int maxVector(int, VectorT);

double SQR(double);
double DMAX(double, double);
double DMIN(double, double);
double RECT(double);


#endif

/* \endc */ 
