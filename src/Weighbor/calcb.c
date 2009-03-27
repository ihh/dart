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
  File: calcb.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{calcb.c}
  \job{Weighbor}
  
  Calculate the $b_{i;j}$ matrix.

*/


/* $Id: calcb.c,v 1.2 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: calcb.c,v $
  Revision 1.2  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.9  1998-06-15 19:56:15-04  nds
  First Release Version!

  Revision 1.8  1998-06-09 18:33:38-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.7  1998-05-14 17:32:53-04  nds
  Fixed documentation errors (TeX problems)

  Revision 1.6  1998-04-08 19:23:43-04  nds
  $\Delta b$ optimization: calculation in $N^2$ time.

  Revision 1.5  1998-03-13 16:52:43-05  nds
  Added -v command line option
  BETA release version

  Revision 1.4  1998-02-06 16:10:33-05  nds
  Added some more debugging printf's

  Revision 1.3  1998-01-15 15:23:16-05  nds
  Version 0.2.8.5.1

  Added a new noise function $\sigma^{2'}_{ij;k}$ equations 0.11 and
  0.12.

  Revision 1.2  1997-10-17 19:22:30-04  nds
  Version 0.2.6: new renormalization code

  Revision 1.1  1997-10-04 18:42:08-04  nds
  Initial revision


  */


#include <math.h>
#include <stdlib.h>

#include "Weighbor/tree.h"
#include "Weighbor/matrix.h"
#include "Weighbor/weighbor.h"

#include <assert.h>

/*
  
  \section{calc_b}

  Calculate the $b_{i;j}$ matrix eq 0.7 and also calculate and return
  $\Delta b_{ij}$ (eq 0.8), $s_{ij}$ (eq 0.9), and the new terms;
  $\Delta^2b_{ij}$ (eq 0.14), and $\overline{\sigma^2_{ij}}$ (eq 0.15)

  */

void
calc_b(int N, MatrixT b, MatrixT deltaB,
       MatrixT delta2B, MatrixT s)
{

  int i, j, k;

  if(printLevel>3)
    fprintf(outfile, "Previous pair (%d,%d)-->%d\n", 
	    ta, tb, tc);

  for(i=0;i<N-1;++i)
    for(j=i+1;j<N;++j) 
      {

	s[i][j]       = 0.0;
	deltaB[i][j]  = 0.0;
	delta2B[i][j] = 0.0;
	
	if(S(i,j)<0) 
	  {
	    for(k=0;k<N;++k)
	      if(k!=i && k!=j)
		{
		  
		  /*
		    epsilon added to handle the cases when the two 
		    $\sigma^2_{ik;j}$ are zero (ie \|sigfac=0|
		    */
		  
		  double weight=1.0/(sigma2_3(i,k,j)
				     +sigma2_3(j,k,i)+EPSILON);
		  
		  s[i][j]       += weight;
		  deltaB[i][j]  += weight*(D(i,k)-D(j,k));
		  delta2B[i][j] += weight*(D(i,k)-D(j,k))*(D(i,k)-D(j,k));
				
		  if(s[i][j]<0) {
		    printf("i %d j %d k %d weight %lf s_ij %lf\n", 
			   i, j, k, weight, s[i][j]);
		    exit(0);
		    assert(s[i][j]>0.0);
		  }
		  
		}
	  }
	else
	  {

	    /* 

	       We have already calculated $s_{ij}$ for this $(i,j)$
	       pair.  So now just update by subtracting out the
	       previously joined pair (\|ta|,\|tb|) and adding in the
	       new taxon \|tc|.

	       */

		double weighta=1.0/(sigma2_3(i,ta,j)
				   +sigma2_3(j,ta,i)+EPSILON);
		double weightb=1.0/(sigma2_3(i,tb,j)
				   +sigma2_3(j,tb,i)+EPSILON);
		double weightc=1.0/(sigma2_3(i,tc,j)
				   +sigma2_3(j,tc,i)+EPSILON);

		
		s[i][j] = S(i,j) - weighta - weightb + weightc;

		deltaB[i][j] 
		  = DelB(i,j)
		  - weighta*(D(i,ta)-D(j,ta)) 
		  - weightb*(D(i,tb)-D(j,tb))
		  + weightc*(D(i,tc)-D(j,tc));

		delta2B[i][j] 
		  = Del2B(i,j)
		  - weighta*(D(i,ta)-D(j,ta))*(D(i,ta)-D(j,ta)) 
		  - weightb*(D(i,tb)-D(j,tb))*(D(i,tb)-D(j,tb))
		  + weightc*(D(i,tc)-D(j,tc))*(D(i,tc)-D(j,tc));


	  }

      }

  for(i=0;i<N-1;++i)
    for(j=i+1;j<N;++j) 
      {
	S(i,j) = S(j,i) = s[j][i] = s[i][j];

	DelB(i,j) = deltaB[i][j];
	DelB(j,i) = -deltaB[i][j];

	deltaB[i][j]  /= s[i][j];
	deltaB[j][i] = -deltaB[i][j];
	
	Del2B(i,j) = delta2B[i][j];
	Del2B(j,i) = delta2B[i][j];

	delta2B[i][j] /= s[i][j];
	delta2B[j][i] = delta2B[i][j];
      }
  
  if(printLevel>2) {
    fprintf(outfile,"Delta b_{ij}=\n");
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) 
	fprintf(outfile,"%9.4g ", deltaB[i][j]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    
    fprintf(outfile,"Delta^2 b_{ij}=\n");
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) 
	fprintf(outfile,"%9.4g ", delta2B[i][j]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    
    fprintf(outfile,"s_{ij}=\n");
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) 
	fprintf(outfile,"%9.4g ", s[i][j]);
    fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
  }
    
  /*
    Finally, calculate $b_{i;j}$ (eq 0.7)
    */

   for(i=0;i<N;++i)
    for(j=0;j<N;++j) 
      {
	
	if( D(i,j) >= fabs(deltaB[i][j]) )

	  b[i][j] = (deltaB[i][j] + D(i,j))/2.0;

	else if( D(i,j) < -deltaB[i][j] )
	  
	  b[i][j] = 0.0;

	else
	  
	  b[i][j] = D(i,j);

      }

}

/*

  \section{\|recalc_b|}

  Redo the calculation of $b_{i;j}$ (eq 0.7) $\Delta b_{ij}$ (eq 0.8),
  $s_{ij}$ (eq 0.9), $\Delta^2b_{ij}$ (eq 0.14) using the new
  $\sigma^{2'}_{ij;k}$ and the old values of $\Delta b_{ij}$.

  */

void
recalc_b(int N, MatrixT b, MatrixT deltaB,
       MatrixT delta2B, MatrixT s)
{

  int i, j, k;

  MatrixT old_deltaB = matrix(N); /* Tmp matrix to hold old values of 
				     $\delta b_{ij}$ */

  setMM(N, deltaB, old_deltaB);


  for(i=0;i<N;++i)
    for(j=0;j<N;++j) 
      {

	s[i][j]       = 0.0;
	deltaB[i][j]  = 0.0;
	delta2B[i][j] = 0.0;
	
	if(i!=j) 
	  {
	    for(k=0;k<N;++k)
	      if(k!=i && k!=j)
		{

		  /*
		    epsilon added to handle the cases when the two 
		    $\sigma^2_{ik;j}$ are zero (ie \|sigfac=0| */

		  /*
		    N.B. We are using the new \|sigma2_3p| function
		    */


		  double weight=1.0/(sigma2_3p(i,k,j,old_deltaB[i][j])
				     +sigma2_3p(j,k,i,old_deltaB[j][i])
				     +EPSILON);
		  
		  s[i][j]       += weight;
		  deltaB[i][j]  += weight*(D(i,k)-D(j,k));
		  delta2B[i][j] += weight*(D(i,k)-D(j,k))*(D(i,k)-D(j,k));

		  if(s[i][j]<0) {
		    printf("Recalc N %d i %d j %d k %d weight=%g sij=%lg " 
			   "sigma %lg %lg [%g %g] deltaB %lg %lg\n", 
			   N, i, j, k, weight, s[i][j], 
			   sigma2_3p(i,k,j,old_deltaB[i][j]),
			   sigma2_3p(j,k,i,old_deltaB[j][i]),
			   sigma2_3p(k,i,j,old_deltaB[k][j]),
			   sigma2_3p(k,j,i,old_deltaB[k][i]),
			   old_deltaB[i][j], old_deltaB[j][i]);
		    exit(0);
		    assert(s[i][j]>0.0);
		  }
		}
	    
	    deltaB[i][j]  /= s[i][j];
	    delta2B[i][j] /= s[i][j];

	    
	  }
      }

  
  if(printLevel>2)
    {
      fprintf(outfile,"Delta b_{ij}=\n");
      for(i=0;i<N;++i) {
	for(j=0;j<N;++j) 
	  fprintf(outfile,"%9.4g ", deltaB[i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
      
      fprintf(outfile,"Delta^2 b_{ij}=\n");
      for(i=0;i<N;++i) {
	for(j=0;j<N;++j) 
	  fprintf(outfile,"%9.4g ", delta2B[i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
      
      fprintf(outfile,"s_{ij}=\n");
      for(i=0;i<N;++i) {
	for(j=0;j<N;++j) 
	  fprintf(outfile,"%9.4g ", s[i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }

  /*
    Finally, calculate $b_{i;j}$ (eq 0.7)
    */

   for(i=0;i<N;++i)
    for(j=0;j<N;++j) 
      {
	
	if( D(i,j) >= fabs(deltaB[i][j]) )

	  b[i][j] = (deltaB[i][j] + D(i,j))/2.0;

	else if( D(i,j) < -deltaB[i][j] )
	  
	  b[i][j] = 0.0;

	else
	  
	  b[i][j] = D(i,j);

      }


   freeMatrix(old_deltaB);

}

/* \endc */ 
