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
  File: noise.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{noise.c}
  \job{Weighbor}
  

*/


/* $Id: noise.c,v 1.2 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: noise.c,v $
  Revision 1.2  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 2.14  1998-06-15 19:56:17-04  nds
  First Release Version!

  Revision 2.13  1998-06-09 18:33:40-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 2.12  1998-06-04 18:35:22-04  nds
  Update for version 2.12
      Changed z score
      Change renormalization if z-score negative
      non-deterministic warning to stderr
      -Xn flag to change R(i,j,j') normalization

  Revision 2.11  1998-05-18 17:26:13-04  nsocci
  Added function \|sigma2tinv| which is used to calculate the $c_i$

  Revision 2.10  1998/04/08 23:23:43  nds
  $\Delta b$ optimization: calculation in $N^2$ time.

  Revision 2.9  1998-03-20 00:00:55-05  nds
  Simple optimization of \|sigma2t| function

  Revision 2.8  1998-03-13 16:52:44-05  nds
  Added -v command line option
  BETA release version

  Revision 2.7  1998-02-06 16:07:05-05  nds
  Change the calcultion of the b's in \|sigma2_3p| to make sure the
  smallest value is >= MINB

  Revision 2.6  1998-02-06 15:41:54-05  nds
  Check for negative sigma in \|sigma_na|

  Check to make sure b_i, b_j, and b_k where non-negative. Set to
  EPSILON if so.

  Revision 2.5  1998-01-15 15:23:16-05  nds
  Version 0.2.8.5.1

  Added a new noise function $\sigma^2_{ij;k}{}'$ equations 0.11 and
  0.12.

  Revision 2.4  1997-12-18 22:13:12-05  nds
  version 0.2.8.4
   + Check to make sure the winner of the tournament is the true winner
     print warning if not or tie.

   + -q option to turn of checking of q[q[i]]=i

   + New R(i,j,j') function (eq 0.11)

  Revision 2.3  1997-11-08 15:56:40-05  nds
  Fixed sigma2t to use the length of the sequences L

  Revision 2.2  1997-10-17 19:22:30-04  nds
  Version 0.2.6: new renormalization code

  Revision 2.1  1997-10-14 23:00:09-04  nds
  Last version of the old renormalization code.

  Revision 2.0  1997-10-09 19:23:22-04  nds
  Version 0.2.x initial revision

  Revision 1.12  1997-10-04 17:05:46-04  nds
  Added some debugging statements to \|sigma2_3|

  Revision 1.11  1997-09-19 18:52:28-04  nds
  Bug fix: extra inversion in the calculation of sig4oo (ie 1/sig4oo
  instead of just sig4oo).

  Revision 1.10  1997-09-18 16:56:41-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.6  1997-08-23 18:29:11-04  nsocci
  Added compile flag RENORM to control whether the renormalization code
  in \|sigma2_3| is used.

  Revision 1.5  1997/08/23 01:16:19  nsocci
  Renormalization code added. Alpha version still needs to be tested
  Note the sigma2_3 function's calling parameteres has been
  changed. Added two booleans to determine whether to renormalize or
  not.

  Revision 1.4  1997/08/22 21:52:56  nsocci
  Added DEBUG flag

  Revision 1.3  1997/08/20 21:23:21  nick
  Added the B def and changed BM1OB (=(b-1)/b) to use the B macro
  instead of hardcoding its value

  Revision 1.2  1997/08/12 00:24:01  nds
  Added forth case to triangle inequality test
  	if \|d_ik > d_ij + d_jk| then \|d_il=d_ij|

  Revision 1.1  1997-08-08 13:59:06-04  nds
  Initial revision


  */

#include <math.h>
#include <assert.h>
#include "Weighbor/weighbor.h"

/*

  \section{Dissimilarity---\|D(d)|}

  Equation 0.1

  */

#define BM1OB ((B-1.0)/B)  /* $(b-1)/b$ */

#define DISS(d) (BM1OB*(1-exp(-(d)/BM1OB)))

/*
  
  \section{Total Variance---\|sigma2t(d)|}

  Equation 0.3

  $$
  \sigma^2_t(d) = n*e^{2d b/(b-1)}D(d)(1-D(d))
  $$

  and its inverse (eq 0.29)
  $$
  \sigma^2_{\rm inv}(x) = {b-1 \over b} \ln\left( 2[x b^2 L + (b-1)^2] \over
  b\sqrt{4 x (b-1) L+(b-1)^2}+ (b-1)(b-2) \right)
  $$
  
  

  */

#define sigB  (B)
#define sigBi (1.0/B)
#define sigBB ((B-1.0)/B)
#define sigBBi (B/(B-1.0))

double
sigma2t(double d)
{

  /*  return( exp(2.0*d/BM1OB)*DISS(d)*(1-DISS(d))/((double)L) ); */

  double sub1 = exp(d*sigBBi)-1.0;

  return( sigBB*(sub1*sub1*sigBi + sub1)/((double)L));
      

}

double
sigma2tinv(double x)
{

  double dL = (double)L;
  double bb = (double)(sigB);

  return(
    sigBB*log(2.0*(x*bb*bb*dL+(bb-1.0)*(bb-1.0))/
	      (bb*sqrt(4.0*x*(bb-1.0)*dL+(bb-1.0)*(bb-1.0))
	       +(bb-1.0)*(bb-2.0))));
  
}
  

/*

  \section{\|sigma2_3|---$\sigma^2_{ij;k}$}

  Non-additive noise for the joining of taxa $i$ and $j$ about taxon $k$.
  See equation 0.4

  $$
  \sigma^2_{ij;k} = \sigma^2(d_{ij}) -  \sigma^2(d_{il}) -  \sigma^2(d_{lj})
  $$

  The two boolean varibles \|bareI| and \|bareJ| are used to control the
  recursion and prevent the calculation of the renormalized sigma when we
  want the bare sigma for a specific index.
  
  */ 

double
sigma2_3(int i, int j, int k)
{

  double d_il, d_lj;

  double sigma;

#ifdef HASH
  int hashidx, a, b, c;
#endif

  if(i==j || i==k || k==j)
    return(0.0);

#ifdef HASH
  ++total;
  a = nodes[i]->ind;
  b = nodes[j]->ind;
  c = nodes[k]->ind;
  hashidx = a+hashsize*b+hashsize*hashsize*c;
  if(hash[hashidx]>0.0)
    {
      ++hashcnt;
      return(hash[hashidx]);
    }
#endif
  
#ifdef DEBUG

  fprintf(stdout, ">>sigma[%d,%d,%d]\n", i, j, k);
  fprintf(stdout, "indices %d=%d %d=%d %d=%d\n",
	  i, nodes[i]->ind,
	  j, nodes[j]->ind,
	  k, nodes[k]->ind);
  
  fprintf(stdout, "D(%d,%d)=%lg\n", i, j, D(i,j));
  fprintf(stdout, "D(%d,%d)=%lg\n", i, k, D(i,k));
  fprintf(stdout, "D(%d,%d)=%lg\n", j, k, D(j,k));
  
#endif


  /*
    Equation 0.5
    */


  if( D(j,k) > (D(i,j)+D(i,k)) ) 
    {
      d_il = 0.0;
    }
  else if( D(i,j) > (D(i,k)+D(j,k))) 
    {
      d_il = D(i,k);
    }
  else if( D(i,k) > (D(i,j)+D(j,k)) )
    {
      d_il = D(i,j);
    }
  else
    {
      d_il = (D(i,j) + D(i,k) - D (j,k))/2;
    }

  d_lj = D(i,j) - d_il;

  sigma = sigma2t(D(i,j)+C(i)+C(j)) - sigma2t(d_il+C(i)) - sigma2t(d_lj+C(j));

  if(sigma<0.0) {
    sigma=0.;
    /*
    if(-sigma<1e-14)
      sigma=0.0;
    else
      {
	printf("ERROR: sigma2_3(%d,%d,%d)=%g < 0.0\n",
	       i,j,k, sigma);
	exit(0);
      }
    */
  }
  
#ifdef DEBUG
  fprintf(stdout, "d_ij=%6lg  %.16lg\n", D(i, j), sigma2t(D(i,j)) );
  fprintf(stdout, "d_il=%6lg  %.16lg\n", d_il, sigma2t(d_il) ); 
  fprintf(stdout, "d_lj=%6lg  %.16lg\n", d_lj, sigma2t(d_lj) );
  fprintf(stdout, "sigma2_3(%d,%d;%d)=%g\n",
	  i, j, k, sigma);

#endif
  
#ifdef HASH
  hash[hashidx] = sigma;
#endif

#ifdef DEBUG
  fprintf(stdout, "sigma2_3(%d,%d;%d)=%g\n",
	  i, j, k, sigma);
#endif

  return(sigma);
  
}


double
sigma_na(double X, double Y)
{

  double sigma = sigma2t(X+Y)-sigma2t(X)-sigma2t(Y);

  if(sigma<0.0) {
    sigma=0.;
    /*
    if(-sigma<1e-14)
      sigma=0.0;
    else
      {
	printf("ERROR: sigma_na(%g,%g)=%g < 0.0\n",
	       X, Y, sigma);
	exit(0);
      }
    */
  }

  return(sigma);

}

double
sigma2_3p(int i, int k, int j, double deltaBij)
{

  double bi, bj, bk;

  bi = DMAX(0.5*(D(i,j)+deltaBij),MINB);
  bj = DMAX(0.5*(D(i,j)-deltaBij),MINB);
  
  bk = (
	DMAX(D(i,k)-bi,MINB)/(sigma2_3(i,k,j)+EPSILON)
	+
	DMAX(D(j,k)-bj,MINB)/(sigma2_3(j,k,i)+EPSILON)
	)
    /( 1.0/(sigma2_3(i,k,j)+EPSILON) + 1.0/(sigma2_3(j,k,i)+EPSILON) );

  return(sigma_na(bi+C(i),bk+C(k)));

}

/* \endc */ 
