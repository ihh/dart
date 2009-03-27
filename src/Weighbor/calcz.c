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
  File: calcz.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{calcz.c}
  \job{Weighbor}
  
  Calculate the $z(i,q(i))$ array equation 0.23


*/


/* $Id: calcz.c,v 1.3 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: calcz.c,v $
  Revision 1.3  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.2  2003/12/15 12:28:06  ihh
  Fixed so gcc can compile with -Wall.

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.19  1998-06-15 19:56:16-04  nds
  First Release Version!

  Revision 1.18  1998-06-09 18:33:39-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.17  1998-06-04 18:35:22-04  nds
  Update for version 2.12
      Changed z score
      Change renormalization if z-score negative
      non-deterministic warning to stderr
      -Xn flag to change R(i,j,j') normalization

  Revision 1.16  1998-05-19 21:43:20-04  nsocci
  Added the extended tournament option (-e)

  Revision 1.15  1998/05/14 21:32:53  nds
  Fixed documentation errors (TeX problems)

  Revision 1.14  1998-05-01 16:55:17-04  nds
  Major Bug, z score was being truncated to int!!

  Revision 1.13  1998-04-14 16:53:43-04  nds
  New way to calculate the z score in order N time. Old method can be
  selected with the -z command line switch.

  Revision 1.12  1997-12-10 18:52:30-05  nds
  Version 0.2.8.2

  Revision 1.11  1997-11-29 21:00:27-05  nds
  Version 0.2.8.2-; still did not implement equation 0.11
  and the q[i] check.

  Revision 1.10  1997-09-18 16:54:14-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.7  1997-09-13 17:52:43-04  nds
  Used the new SIG2_3 macro for the renormalized sigma2_3

  Revision 1.6  1997-09-11 22:34:54-04  nds
  Added EPSILON to the calculation of sigPQ. Very bad trees (ie with
  triangle inequality violations) may have sigPQ=0

  Revision 1.5  1997-08-22 21:16:18-04  nsocci
  Renormalization code added. Alpha version still needs to be tested
  Note the sigma2_3 function's calling parameteres has been
  changed. Added two booleans to determine whether to renormalize or
  not.

  Revision 1.4  1997/08/21 16:54:15  nick
  Changed the equation numbers in the comments to correspond to the
  lastest version (97-Aug-21) of the program.

  Revision 1.3  1997/08/20 20:39:59  nick
  Error in sigma2PQ function when all sigmas are zero. Returning wrong
  value.

  Revision 1.2  1997/08/15 04:04:16  nds
  Bug, was not doing the sum of the second term in eq 0.25
  correctly. The test to not duplicate k,l pairs was skipping terms in
  this sum which should sum over both the k,l and l,k pair. Moved sum to
  before test for duplicate.

  Revision 1.1  1997-08-08 13:59:06-04  nds
  Initial revision


  */


#include <math.h>
#include <assert.h>

#include "Weighbor/tree.h"
#include "Weighbor/matrix.h"
#include "Weighbor/weighbor.h"


/*

  \section{\|calc_z1|}

  Calculate the $Z$-score of each $(i, j)$ pair. This is the older
  $N^3$ method and will be called if the \|oldZflag| is set.

  */

void
calc_z1(int N, VectorT z, int *q)
{

  int i;
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  VectorT klpair=vector(N);

  for(i=0;i<N;++i) 
    {
      
      int k, j=q[i];
      double PQ_ij=0.0;           /* $PQ(ij)$ (eq 0.16)            */
      double s_PQ_ij=0.0;         /* $s_{PQ}(ij)$                  */
      double sigma_PQ_ij=0.0;     /* $\sigma^2_{PQ}(ij)$ (eq 0.25) */
      
      
      /*
	reset pair matrix
	*/
      
      for(k=0;k<N;++k)
	klpair[k] = -1;
      
      /*
	
	\subsection{Equation 0.16}
	
	$$
	PQ(ij) \equiv {1\over s_{PQ}(ij)}
	\sum_{k} {1\over \sigma^2_{PQ}(ij)}
	\min_l[PQ(ij,kl)]
	$$
	
	\'N.B.' $j\equiv q(i)$, ie we want $PQ(i\,q(i))$
	
	Note we calculate 
	$$
	\sum_{k} {1\over \sigma^2_{PQ}(ij)}
	\min_l[PQ(ij,kl)]
	$$
	and
	$$
	\sum_{k} {1\over \sigma^2_{PQ}(ij)} \equiv s_{PQ}(ij)
	$$
	the normalization factor and also
	$$
	\sigma^2_{PQ}(ij)\equiv{1\over s_{PQ}(ij)}
	+{1\over N-2}\sum_{k\ne i,j}\sigma^2_{ij;k}
	$$
	
	
	*/
      
      for(k=0;k<N;++k)
	{
	  int l;
	  double sigPQ;
	  double minPQ;
	  int minl;
	  
	  /*
	    $k \ne i$ or $j$
	    */
	  
	  if(k==i || k==j)
	    continue;
	  
	  /*
	    get an inital $l\ne i$ or $j=q(i)$ or $k$
	    */
	  
	  l=0;
	  while(l==i || l==j || l==k)
	    ++l;
	  
	  assert(l<N); /* there had better be at least one l left */
	  
	  minPQ=PQ(i,j,k,l);
	  minl = l;
	  
	  /*
	    
	    \subsection{$\min_l[PQ(ij,kl)]$}
	    
	    Find the minimum over $l$ of $PQ(ij,kl)$
	    
	    */
	  
	  for(++l;l<N;++l)
	    if(l!=i && l!=j && l!=k && PQ(i,j,k,l)<minPQ)
	      {
		minPQ=PQ(i,j,k,l);
		minl=l;
	      }
	  
	  /*
	    Accumulate the second term of eq. 0.25
	    $$
	    {1\over N-2}\sum_{k\ne i,j}\sigma^2_{ij;k}
	    $$
	    where $j=q(i)$
	    */
	  
	  if(k!=i && k!=j)
	    {
	      sigma_PQ_ij += sigma2_3(i,j,k);
	    }
	  
	  /*
	    Keep track of the (k,l) pairs to avoid duplicates; i.e. pair 
	    (l,k)
	    */
	  
	  klpair[k] = minl;
	  
	  /*
	    Now check for duplicates
	    */
	  
	  if(minl<k && klpair[minl] == k)
	    {
	      continue; /* skip rest of loop; do not add in this term */
	    }
	  
	  /*
	    
	    \|sigPQ| = $\sigma^2_{PQ}(ij,kl)$. Add $\epsilon$ to
	    handle cases when it equals zero
	    
	    */
	  
	  sigPQ = sigma2PQ(i,j,k,minl) + EPSILON;
	  
	  assert(sigPQ>0.0);
	  
	  PQ_ij += minPQ/sigPQ;
	  s_PQ_ij += 1.0/sigPQ;
	  
	} /* \|for(k=0;k<N;++k)| */
      
      
      /*
	
	Finish calculation of 
	$$
	\sigma^2_{PQ}(ij) \equiv{1\over s_{PQ}(ij)}
	+{1\over N-2}\sum_{k\ne i,j}\sigma^2_{ij;k}
	$$
	
	*/
      
      sigma_PQ_ij = 1.0/s_PQ_ij + (0.25)*(1.0/((double)N-2.0))*sigma_PQ_ij;
      
      /*
	
	And now the $Z$-score (eq 0.23)
	$$
	z(i,q(i)) \equiv PQ(i\,q(i))/\sqrt{\sigma^2_{PQ}(ij)}
	$$
	
	\'N.B.' Remember \|PQ_ij|=$s_{PQ}(ij) PQ(ij)$ so need to divide
	by normalization \|s_PQ_ij|
	
	*/

      calcPhi(N, i, j);

      z[i] = PQ_ij/(s_PQ_ij*sqrt(sigma_PQ_ij));
	  
    }
  

  freeVector(klpair);

}

/*

  \section{\|calc_z2|}

  This is the newer method that runs in $N^2$ time and which does a
  heuristic min over both $(k,l)$.

  */

void
calc_z2(int N, VectorT z, int* q, int* q2)
{

  int i, j, k, l;
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);

  for(i=0;i<N;++i) 
    {

      int minl, mink;
      double mArgl = 0, mArgk = 0;

      j = q[i];
      k = q2[i];

      /*

	min over $l$
	$$
	\min_l PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
	$$

	*/
      
      minl=-1;
      for(l=0;l<N;++l)
	if(l!=i && l!=j && l!=k)
	  {
	    double argl = PQ(i,j,k,l)/(sqrt(sigma2PQ(i,j,k,l)+EPSILON));
	    
	    if(minl==-1 || argl<mArgl)
	      {
		minl = l;
		mArgl = argl;
	      }

	  }

      /*

	min over $k$
	$$
	\min_k PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
	$$
      
	*/

      mink=-1;
      for(k=0;k<N;++k)
	if(k!=i && k!=j && k!=minl)
	  {
	    double argk = PQ(i,j,k,minl)/(sqrt(sigma2PQ(i,j,k,minl)+EPSILON));
	    
	    if(mink==-1 || argk<mArgk)
	      {
		mink = k;
		mArgk = argk;
	      }

	  }

      z[i] = mArgk;

    }
}

/*

  \subsection{\|calcZ2|}

  This is the newer method that runs in $N^2$ time and which does a
  heuristic min over both $(k,l)$. This function just calculates Z for
  just one (i,q(i)) pair.

  */

double
calcZ2(int N, int i, int q, int q2)
{

  int j, k, l;
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  int minl, mink;
  double mArgl = 0, mArgk = 0;
  
  j = q;
  k = q2;
  
  /*

    min over $l$
    $$
    \min_l PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
    $$
    
  */
  
  minl=-1;
  for(l=0;l<N;++l)
    if(l!=i && l!=j && l!=k)
      {
	double argl 
	  = PQ(i,j,k,l)
	  /(sqrt(
		 sigma2PQ(i,j,k,l)
		 +(sigma2_3(i,j,k)+sigma2_3(i,j,l))/8.0
		 +EPSILON));
	
	if(minl==-1 || argl<mArgl)
	  {
	    minl = l;
	    mArgl = argl;
	  }
	
      }
  
  /*
    
    min over $k$
    $$
    \min_k PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
    $$
    
  */
  
  mink=-1;
  for(k=0;k<N;++k)
    if(k!=i && k!=j && k!=minl)
      {
	double argk 
	  = PQ(i,j,k,minl)
	  /(sqrt(
		 sigma2PQ(i,j,k,minl)
		 +(sigma2_3(i,j,k)+sigma2_3(i,j,minl))/8.0
		 +EPSILON));
	
	if(mink==-1 || argk<mArgk)
	  {
	    mink = k;
	    mArgk = argk;
	  }
	
      }
  
  return(mArgk);

   
}


/*
  
  \section{\|PQ(ij;kl)|}

  Equation 0.22. Note this equation can be written has the following form
  $$
  {1\over 2}\left[{ {{\alpha\over z_1}+{\beta\over z_2}}
       \over{{1\over z_1}+{1\over z_2}}}-d_{ij}-d_{kl}\right]
  $$
  where
  $\alpha=(d_{ik}+d_{jl})$, $\beta=(d_{il}+d_{jk})$, 
  $z_1=\min(\sigma^2_{ik;j},\sigma^2_{ik;l})
      +\min(\sigma^2_{jl;i},\sigma^2_{jl;k})$ and 
  $z_2=\min(\sigma^2_{il;j},\sigma^2_{il;k})
      +\min(\sigma^2_{jk;i},\sigma^2_{jk;l})$
  
  However if this expression is unstable if $z_1$ or $z_2$ is zero even 
  though the expression is still well defined.

  The expression is evaluated in the following form
  $$
  {1\over 2}\left[{{\alpha z_2 + \beta z_1}\over{z_2+z_1}}-d_{ij}-d_{kl}\right]
  $$
  which is finite as long as either is greater than zero. If both are zero 
  then the expression has the following limit
  $$
  {1\over 2}\left[{{\alpha + \beta}\over{2}}-d_{ij}-d_{kl}\right]
  $$

  */


double
PQ(int i, int j, int k, int l)
{

  double alpha, beta, z1, z2;

  alpha = D(i,k)+D(j,l);
  beta  = D(i,l)+D(j,k);

  z1 = DMIN(sigma2_3(i,k,j),sigma2_3(i,k,l))
     + DMIN(sigma2_3(j,l,i),sigma2_3(j,l,k));

  z2 = DMIN(sigma2_3(i,l,j),sigma2_3(i,l,k))
     + DMIN(sigma2_3(j,k,i),sigma2_3(j,k,l));

  if( (z1+z2) == 0.0 )
    return( (((alpha+beta)/2.0)-D(i,j)-D(k,l))/2.0 );
  else
    return( (((alpha*z2+beta*z1)/(z2+z1)) - D(i,j) - D(k,l))/2.0 );

}


/*

  \section{\|sigma2PQ(ij,kl)|}

  As in the case of $PQ(ij, kl)$, $\sigma^2_{PQ}(ij,kl)$ (equation 0.16)
  can also be rewritten in the following form:
  $$
  {1\over 4}\left[ {{z_1 z_2}\over{z_1+z_2}} 
              + (\sigma^2_{kl;i}+\sigma^2_{kl;j})/2 \right]
  $$
  where $z_1$ and $z_2$ are as defined in the comment for the $PQ(ij, kl)$
  function. Again this will be finite as long as either is non-zero. If
  both are then the limit is just zero.

  */

double 
sigma2PQ(int i, int j, int k, int l)
{

  double z1, z2;

  z1 = DMIN(sigma2_3(i,k,j),sigma2_3(i,k,l))
     + DMIN(sigma2_3(j,l,i),sigma2_3(j,l,k));

  z2 = DMIN(sigma2_3(i,l,j),sigma2_3(i,l,k))
     + DMIN(sigma2_3(j,k,i),sigma2_3(j,k,l));

  if( (z1+z2) == 0.0 )
    return( (sigma2_3(k,l,i)+sigma2_3(k,l,j))/8.0 );
  else
    return( ((z1*z2)/(z1+z2)
	      + (sigma2_3(k,l,i)+sigma2_3(k,l,j))/2.0)/4.0 );

}

/*
  \|calcPhi| is used in the renormalization phase and is called from 
 \|build.c|
 */

double
calcPhi(int N, int i, int ip)
{

  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  VectorT klpair=vector(N);
  int k;
  double PQ_iip=0.0;           /* $PQ(ii')$ (eq 0.16)      */
  double s_PQ_iip=0.0;         /* $s_{PQ}(ii')$            */
  double sigma_iip=0.0;        /* \sigma^2_{ii'} (eq 0.27) */

  /*
    reset pair matrix
    */
  

  for(k=0;k<N;++k)
    klpair[k] = -1;
  
  
  for(k=0;k<N;++k)
    {
      int l;
      double sigPQ;
      double minPQ;
      int minl;
      
      /*
	$k \ne i$ or $i'$
	*/
      
      if(k==i || k==ip)
	continue;
      
      /*
	get an inital $l\ne i$ or $i'=q(i)$ or $k$
	*/
      
      l=0;
      while(l==i || l==ip || l==k)
	++l;
      
      assert(l<N); /* there had better be at least one l left */
      
      minPQ=PQ(i,ip,k,l);
      minl = l;
      
      /*
	
	\subsection{$\min_l[PQ(ii',kl)]$}
	
	Find the minimum over $l$ of $PQ(ii',kl)$
	
	*/
      
      for(++l;l<N;++l)
	if(l!=i && l!=ip && l!=k && PQ(i,ip,k,l)<minPQ)
	  {
	    minPQ=PQ(i,ip,k,l);
	    minl=l;
	  }
      
      
      /*
	Keep track of the (k,l) pairs to avoid duplicates; i.e. pair 
	(l,k)
	*/
      
      klpair[k] = minl;
      
      /*
	Now check for duplicates
	*/
      
      if(minl<k && klpair[minl] == k)
	{
	  continue; /* skip rest of loop; do not add in this term */
	}
      
      /*
	
	\|sigPQ| = $\sigma^2_{PQ}(ii',kl)$. Add $\epsilon$ to
	handle cases when it equals zero
	
	*/
      

      sigPQ = sigma2PQ(i,ip,k,minl) + EPSILON;
      
      assert(sigPQ>0.0);
      
      PQ_iip += minPQ/sigPQ;
      s_PQ_iip += 1.0/sigPQ;

      
    } /* \|for(k=0;k<N;++k)| */
  
  
  /*
    
    Calculate $\sigma^2_{ii'}$ equation (0.27)
    
    */
  
  sigma_iip=0.0;
  for(k=0;k<N;++k)
    if(k!=i && k!=ip)
      sigma_iip += sigma2_3(i,ip,k);

  sigma_iip /= ((double)N)-2.0;

  /*
    
    And now $\phi$ eq (0.26)

    */


  freeVector(klpair);

  return( (-PQ_iip/2.0)
	  /( (1.0/(sigma_iip+EPSILON)) + s_PQ_iip/4.0 + EPSILON ) );
  
}



/* \endc */



