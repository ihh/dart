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
  File: build.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{build.c}
  \job{Weighbor}
  
  \synopsis{\|buildTree| will be a tree given a distance matrix
  \"D". It uses the weighted neighbor joining algorithm of W. Bruno
  (\'refs').}


*/


/* $Id: build.c,v 1.3 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: build.c,v $
  Revision 1.3  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.2  2003/12/15 12:28:06  ihh
  Fixed so gcc can compile with -Wall.

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 2.25  1998-06-15 19:56:14-04  nds
  First Release Version!

  Revision 2.24  1998-06-09 18:33:38-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 2.23  1998-06-04 18:37:35-04  nds
  Minor bug fix: changed a stderr to outfile

  Revision 2.22  1998-06-04 18:35:22-04  nds
  Update for version 2.12
      Changed z score
      Change renormalization if z-score negative
      non-deterministic warning to stderr
      -Xn flag to change R(i,j,j') normalization

  Revision 2.21  1998-05-19 21:42:19-04  nsocci
  Added the extended tournament option (-e)

  Revision 2.20  1998/05/19 17:10:29  nsocci
  Added tie warning message to weighbor.out file

  Revision 2.19  1998/05/18 21:27:41  nsocci
  Changed how the $c_i$ are calculate. They now use the $\sigma^2_{\rm inv}$
  function

  Revision 2.18  1998/05/14 21:32:53  nds
  Fixed documentation errors (TeX problems)

  Revision 2.17  1998-05-08 15:27:52-04  nds
  Now print just the max of the two p-values
  Moved the swap of i and ip to after the printing of which pair to keep
  the names straight

  Revision 2.16  1998-05-05 15:24:45-04  nds
  Fixed bug when joining nodes (if i=Nleft-1). Make sure i<ip to prevent
  this problem.

  Adjusted what is printed at various verbose levels.

  Now printing p-values instead of LLR (also fixed LLR, factor of two
  off).

  Revision 2.15  1998-05-01 19:31:49-04  nds
  Fixed sign error in LLR and LLRp calculation

  Revision 2.14  1998-05-01 16:57:20-04  nds
  Made numerous changes to log file and adjusted what was printed at
  various verbose levels. To print the new LLR value of R(q(i),i,j) the
  calcR funtion need to be called from buildTree which require the
  moving of the deltaB and associated matrices to global level. This
  should improve performance since now they are allocated only once per
  tree rather than a every iteration

  Revision 2.13  1998-04-08 19:23:42-04  nds
  $\Delta b$ optimization: calculation in $N^2$ time.

  Revision 2.12  1998-03-13 16:52:43-05  nds
  Added -v command line option
  BETA release version

  Revision 2.11  1998-02-25 21:47:34-05  nds
  Fixed equation 0.32: changed the form of the MAX(*,BMIN) terms

  Revision 2.10  1998-02-20 20:43:26-05  nds
  Fixed problem with negative branch lengths for equation 0.31

  Revision 2.9  1998-02-19 23:52:10-05  nds
  Version 0.2.8.5.2

  Using new formula for $\sigma_\infty$. No longer need to store vaules
  of it in the nodes so removed \|sigma2inf| variables from node
  structure.

  Revision 2.8  1998-01-15 18:01:05-05  nds
  Added test to make sure that the new branch lengths from the newly
  joined taxa are not less than zero. If so then set them to zero.

  Revision 2.7  1997-12-16 15:20:49-05  nds
  Added the checkQQI boolean variable and the -q command line switch to
  turn of the checking of q[q[i]]=i in selection the best pair.

  Revision 2.6  1997-12-10 18:52:29-05  nds
  Version 0.2.8.2

  Revision 2.5  1997-11-29 21:00:26-05  nds
  Version 0.2.8.2-; still did not implement equation 0.11
  and the q[i] check.

  Revision 2.4  1997-11-17 23:10:36-05  nds
  Fix bug in an assertion test. Values less than 1e-14 are considered
  identically zero.

  Revision 2.3  1997-11-15 00:31:32-05  nds
  Changes for weighor version 0.2.7.1

  Revision 2.2  1997-10-17 19:22:29-04  nds
  Version 0.2.6: new renormalization code

  Revision 2.1  1997-10-15 16:51:18-04  nds
  Put in test for negative branch lengths. Set to zero and print warning
  to logfile

  Revision 2.0  1997-10-09 19:22:44-04  nds
  Version 0.2.x initial revision

  Revision 1.12  1997-10-04 18:45:03-04  nds
  Fixed some TeX errors in comments. Last 0.1.x version

  Revision 1.11  1997-09-19 18:49:56-04  nds
  Fixed two bugs:
  (1) Added EPSILON to sigma_inf to handle cases when it equals zero
  (2) Using N when should have been using Nleft in a number of places

  Revision 1.10  1997-09-18 16:47:48-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.9  1997-09-11 22:42:14-04  nds
  Check for negative lengths in the final three taxa join and set any
  negative lengths to zero

  Revision 1.8  1997-08-23 18:31:37-04  nsocci
  FIxed bug in the node joining algorithm. When joining nodes i and j
  where i is the last node in the list then there was a problem that
  node I was clobbered by the code the moves the nodes around. Fix is to
  make sure that i<j if not then swap them.

  Revision 1.7  1997/08/23 20:49:58  nsocci
  Added compiler flags
  	CHECK_QQOFI
    to control whether the code to check that $q(q(i))==i$, and
  	USEBMIN
    to control whether we use $b_min_i$ for branch lengths or if not set
    then use $b_{i;j}$.

  Set or unset these flag in weighbor.h file

  Revision 1.6  1997/08/23 01:16:16  nsocci
  Renormalization code added. Alpha version still needs to be tested
  Note the sigma2_3 function's calling parameteres has been
  changed. Added two booleans to determine whether to renormalize or
  not.

  Revision 1.5  1997/08/22 18:21:42  nsocci
  Trival change, added code to print if q[q[i]]==i condition was met

  Revision 1.4  1997/08/21 16:59:41  nsocci
  Changed the algorithm that picks the best pair to first consider only
  those pairs what have the property q(q(i))=i; and pick the best z
  score. If no pair has that property then just pick the best z score of
  all of them.

  Also the branch lengths in eq 26 are now bmin[i] instead of b[i][j]

  Revision 1.3  97/08/21  16:54:14  nick
  Changed the equation numbers in the comments to correspond to the
  lastest version (97-Aug-21) of the program.
  
  Revision 1.2  1997/08/20 22:28:50  nick
  Added code to calcualate $\sigma^2_\infty$ for both of the added
  branches and also fixed the calculation of the distances from the new
  taxon to the remaining ones.

  Revision 1.1  1997/08/08 17:59:06  nds
  Initial revision


  */


#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "Weighbor/tree.h"
#include "Weighbor/matrix.h"
#include "Weighbor/weighbor.h"
#include "Weighbor/io.h"


RootNodeT *
buildTree(int N, MatrixT Dorig, char **taxon)
{

  /* 
     \section{Initialize} 
     Initialize main variables 
     */

  int i;
  RootNodeT *root;
  int Nleft, Nnext; /* number of Nodes left to be joined 
		       and the next index to be used */
  

  MatrixT b=matrix(N);      /* the $b_{i;j}$ matrix (eq 0.7) */ 

  /* $q(i)$ array: value which minimizes $R(i,q(i),j)\,\forall j\ne i,q(i)$ */
  int *q;
  
  int *q2;                  /* Second best value */ 
  VectorT R=vector(N);      /* $R(i,q(i))$ (eq 0.10) */
  VectorT LLR=vector(N);    /* $R(i,q(i),q2(i))$ */
  VectorT Zscore=vector(N); /* $z(i,q(i))$ */

  /*

    This auxilary matrices are globally defined in \|weighbor.h| we do
    this to make it simplier so we do not always have to pass these
    around. Note that the need to be visible here as we will be
    calling \|calcR| later in this function and \|calcR| needs these
    values

    */

  s       = matrix(N);      /* $s_{ij}$ eq 0.9 */
  deltaB  = matrix(N);      /* $\Delta b_{ij}$ eq 0.8 */
  delta2B = matrix(N);      /* $\Delta^2 b_{ij}$ */
  if(recalcB)
    oldDeltaB = matrix(N);

  /* 

     This will hold this orignal $N$ distances plus any distances from
     the $N-3$ internal nodes. Note we do not care about the root node
     so $N-3$ and not $N-2$

     */

  mD=matrix(2*N-3); 

  /*

    This is the renormalization vector $c_i$ (eq 0.39) and matrix
    $c_{i;j}$ (eq 0.43 ver0.2.5); again it must
    be large enough to hold both the original and the new joined taxa

    N.B. \|vector| sets all elements to zero.

    */

  vC=vector(2*N-3);


  /*
    
    This matrices hold the previous iterations values of $s_{ij}$,
    $\Delta b_{ij}$, etc. They are used to speed up the next
    iterations calcultions of these quantities.

    */

  mS     = matrix(2*N-3);
  mDelB  = matrix(2*N-3);
  mDel2B = matrix(2*N-3);

  /*

    Init \|mS| to -1 to keep track of which entries have not yet been
    computed.  */

  for(i=0;i<2*N-3;++i)
    {
      int j;
      for(j=0;j<2*N-3;++j)
	mS[i][j] = -1.0;
    }

  /*

    Make a copy of the original distance matrix; embed it in the
    larger matrix which will hold the new distance from the added
    internal nodes of the tree.

    */

  setMM(N, Dorig, mD);

  /*
    
    Allocate and initialize the \|q|, \|q2| and \|nodes| array. $2N-3$
    nodes to hold both the original and the added nodes.
    
    */

  q = (int *)malloc(N*sizeof(int));
  if(!q) printError("build::buildTree:out of memory-q\n");
  q2 = (int *)malloc(N*sizeof(int));
  if(!q2) printError("build::buildTree:out of memory-q2\n");

  nodes = (NodeT **)malloc( (2*N-3)*sizeof(NodeT *));
  if(!nodes) printError("build::buildTree:out of memory-nodes");

  for(i=0;i<N;++i) {
    nodes[i] = createNode();
    nodes[i]->name = taxon[i];
    nodes[i]->ind  = i;
  }


  Nleft = N;
  Nnext = N;

  /*
    
    \section{Loop until 3 taxa left}

    While we have more than 3 nodes left do the neighbor joining algorithm. 
    Each pass of the algorithm will join 2 nodes replacing them with one.

    */

  while(Nleft>3)
    {

      int j, k, ip, ip2;
      double minR = 0, min2R = 0;
      NodeT *newNode, *tmpNode;
      double sigma_inf_i, sigma_inf_ip;
      double sig_r, sig_l;
      int jj, jjmin;
      double LLRp = 0, tR;

      /* \subsection{Calculate Residual} */

      calc_q(Nleft, q, R, b, q2, LLR, Zscore);

      if(printLevel>2)
	for(k=0;k<Nleft;++k)
	  fprintf(outfile, "q[%d]=%d R(%d,%d)=%g\n",
		  k, q[k], k, q[k], R[k]); 

      /*

	Find $i$ than minimizes $R(i,q(i))$. With the constraint that
	$q(q(i))=i$ first if no pair found then find the best $i$
	without this constraint.

	Note: the \|checkQQI| flag determines if we will use the
	$q(q(i))=i$ constraint.

	Note: j will hold the next best pair

	*/

      i = -1;
      j = -1;

      if(checkQQI) { 
	for(k=0;k<Nleft;++k)
	  if(q[q[k]]==k) {
 	    if(R[k]<minR || i==-1)
	      {
		if(printLevel>3)
		  fprintf(outfile, 
			  "ij=%d%d k=%d q[k]=%d minR = %.16g R[k] = %.16g\n",
			  i,j,k, q[k], minR, R[k]);
		j = i;
		min2R = minR;
		i = k;
		minR = R[k];
	      }
	    else if(R[k]>minR && (R[k]<min2R || j==-1) )
	      {
		j = k;
		min2R = R[k];
	      }
	  }
      }

      if(i==-1) /* No pair had $q(q(i))=i$ */
	{
	  if(R[0]<R[1]) {
	    i = 0;
	    minR = R[0];
	    j = 1;
	    min2R = R[1];
	  } else {
	    i = 1;
	    minR = R[1];
	    j = 0;
	    min2R = R[0];
	  }	    
	  for(k=1;k<Nleft;++k)
	    if(R[k]<minR)
	      {
		j = i;
		min2R = minR;
		i = k;
		minR = R[k];
	      }
	    else if(R[k] < min2R && R[k] > minR)
	      {
		j = k;
		min2R = R[k];
	      }

	  if(checkQQI && printLevel>1)
	      fprintf(outfile, "No pair with q[q[i]]==i ");
	  else
	    if(q[q[i]]!=i && printLevel>1)
	      fprintf(outfile, 
		      "The pair does not satisfy q[q[i]]==i (checking is off)"
		      );
	}

      ip = q[i];
      ip2 = j;

      /*
	
	If the extended tournament option is set (-e) then run two more tournaments for 
	(i,q[i]) to see who really wins. 

	*/

      if(extendedTourn)
	{

	  double minR1 = 0, minR2 = 0, tmpR, oldR=R[i];
	  int jmin=-1, jpmin=-1;

	  /* 
	     First fine the j the minimizes R(i,j)
	     */

	  for(j=0;j<Nleft;++j) 
	    if(j!=i && j!=q[i])
	      {
		if(j!=q2[i])
		  tmpR = calcR2(Nleft, i, j, q2[i]);
		else
		  tmpR = calcR2(Nleft, i, j, q[i]);

		if(tmpR<minR1 || jmin==-1)
		  {
		    minR1=tmpR;
		    jmin = j;
		  }
	      }

	  /* 
	     and now the $j'$ that minimizes $R(j',q[i])$
	     */

	  for(j=0;j<Nleft;++j) 
	    if(j!=i && j!=q[i])
	      {
		if(j!=q2[i])
		  tmpR = calcR2(Nleft, j, q[i], q2[i]);
		else
		  tmpR = calcR2(Nleft, j, q[i], i);

		if(tmpR<minR2 || jpmin==-1)
		  {
		    minR2=tmpR;
		    jpmin = j;
		  }
	      }

	  /*
	    Now fnd which of the three is the smallest
	    */
	  
	  if(minR1<minR2 && minR1<R[i])
	    {
	      ip = jmin;
	      if(printLevel>1)
		fprintf(outfile, 
			"Extended Tournament New Winner(A): (%d, %d) R=%g\n",
			i, ip, minR1);
	    }
	  else if(minR2<minR1 && minR2<R[i])
	    {
	      i = jpmin;
	      if(printLevel>1)
		fprintf(outfile, 
			"Extended Tournament New Winner(B): (%d, %d) R=%g\n",
			i, ip, minR2);
	    }	    

	  if(printLevel>3)
	    fprintf(outfile, "R=%g, R1=%g, R2=%g\n", oldR, minR1, minR2);

	}

      /*
	
	Find the $jj$ that minimizes $R(q(i),i,jj)$ and then print out 
	the LLR and LLR' values.
	
      */
      
      jjmin=-1;
      
      for(jj=0;jj<Nleft;++jj)
	if(jj!=i && jj!=ip 
	   && (((tR=calcR(Nleft, ip, jj, i))<LLRp) || jjmin==-1))
	  {
	    jjmin = jj;
	    LLRp = tR;
	  }
	
      LLRp *= 0.5;

      if( (LLR[i]<1e-6) && (LLRp<1e-6) ) {
	fprintf(stderr, 
		"warning: tie scores encountered; topology may depend on sequence order!\n");
	fprintf(stderr, "taxon %s and taxon %s\n\n",
		nodes[i]->name, nodes[ip]->name);

	if(printLevel>1) {
	  fprintf(outfile, 
		  "warning: tie scores encountered; topology may depend on sequence order!\n");
	  fprintf(outfile, "taxon %s and taxon %s\n\n",
		  nodes[i]->name, nodes[ip]->name);
	}
      }
      
      if(printLevel>0) {
	fprintf(outfile, 
		"\nJoin taxon %s to taxon %s (%s next best choice)\n",
		nodes[i]->name, nodes[ip]->name, nodes[q2[i]]->name);
	

	fprintf(outfile, "     p-value = %g\n", 
		DMAX(1.0/(exp(LLR[i])+1.0), 1.0/(exp(LLRp)+1.0)));
		
	if(printLevel>1) {
	  fprintf(outfile,"\nJoin taxon %s to taxon %s; R=%g\n", 
		  nodes[i]->name, nodes[ip]->name, minR);
	  
	  if(ip2!=-1 && ip2!=i && ip2!=ip)
	    fprintf(outfile, "Second best pair (%s, %s); R=%g\n",
		    nodes[ip2]->name, nodes[q[ip2]]->name, min2R);
	  else
	    fprintf(outfile, "No second best pair\n");
	}
      }
      

      /* 

	 Note due to the way we shuffle around nodes after joining:
	 i->Nnext, New->i, ip<->Nleft-1, if ip is less than i and
	 i=Nleft-1 then the new node will be in position ip not i!!
	 But tc (the global that is suppose to point to the position
	 of the new node for calcb) is set to i so this will screw us
	 up. The simpliest solution is to make sure i<ip; swap if they
	 are not.

	*/


      if(ip<i) {
	int tt;
	tt=i;
	i=ip;
	ip=tt;
      }

      /*
	
	Need to calculate the new branch lengths $\bar b_{i;i'}$ and
	$\bar b_{i';i}$, eq. 0.19.

	Note if the z-score is negative then we calculate $\phi$ eq
	(0.26) and use it to renormalize $d_{i,i'}$ and recompute 
	$b_{i;i'}$ and $b_{i';i}$.

	*/

      if(Zscore[i]<0.0) {

	double phi_iip, dBar_iip;
	
	phi_iip = calcPhi(Nleft, i, ip);
	if(printLevel>2)
	  fprintf(outfile, "Renormalizing z[%d,%d] = %g\n", i, ip, Zscore[i]);
	if(phi_iip>0)
	  {
	    dBar_iip = D(i,ip)-phi_iip;
	    if(printLevel>2)
	      fprintf(outfile, "phi=%g dBar_iip=%g\n", phi_iip, dBar_iip);
	    
	    /* renormalize the b's */
	    
	    if( dBar_iip >= fabs(deltaB[i][ip]) )
	      b[i][ip] = (deltaB[i][ip] + dBar_iip)/2.0;
	    else if( dBar_iip < -deltaB[i][ip] )
	      b[i][ip] = 0.0;
	    else
	      b[i][ip] = dBar_iip;
	    
	    
	    if( dBar_iip >= fabs(deltaB[ip][i]) )
	      b[ip][i] = (deltaB[ip][i] + dBar_iip)/2.0;
	    else if( dBar_iip < -deltaB[ip][i] )
	      b[ip][i] = 0.0;
	    else
	      b[ip][i] = dBar_iip;
	  }
      }

      nodes[i ]->rho = b[i][ip];
      nodes[ip]->rho = b[ip][i];

      if(nodes[i ]->rho < 0.0) 
	{
	  if(printLevel>0)
	    fprintf(outfile, 
		    "WARNING: Negative branch length %lg set to zero\n", 
		    nodes[i ]->rho);
	  nodes[i ]->rho = 0.0;
	  nodes[ip]->rho = D(i,ip);
	}
      else if(nodes[ip]->rho < 0.0) 
	{
	  if(printLevel>0)
	    fprintf(outfile, 
		    "WARNING: Negative branch length %lg set to zero\n", 
		    nodes[ip]->rho);
	  nodes[ip]->rho = 0.0;
	  nodes[i ]->rho = D(i,ip);
	}

      if(printLevel>3)
	{
	  fprintf(outfile, "\\bar b_[%d%d] = %g b_[%d%d]=%g\n",
		  i, ip, nodes[i]->rho, i, ip, b[i][ip]);
	  fprintf(outfile, "\\bar b_[%d%d] = %g b_[%d%d]=%g\n\n",
		  ip, i, nodes[ip]->rho, ip, i, b[ip][i]);
	}
      
      newNode = createNode();
      
      newNode->ind = Nnext;
      newNode->child_r = nodes[i];
      newNode->child_l = nodes[ip];
      newNode->name = nodes[i]->name;
      nodes[Nnext] = newNode;

      /*

	Calculate $\sigma^2_\infty(i\bar\imath)$ (eq. 0.27) for each
	of the joined taxa.

	*/

      sigma_inf_i  = 0.0;
      sigma_inf_ip = 0.0;

      for(j=0;j<Nleft;++j) 
	{
	  if(j!=i && j!=ip)
	    {
	      sigma_inf_i  
		+= sigma_na(DMAX(b[i][ip],MINB)+C(i), 
			    DMAX(D(i,j)-b[i][ip],MINB)+C(j) );
	      sigma_inf_ip 
		+= sigma_na(DMAX(b[ip][i],MINB)+C(ip), 
			    DMAX(D(ip,j)-b[ip][i],MINB)+C(j) );
	    }
	}

      /*
	
	Add \|EPSILON| here to make the following formulae a bit simplier

	*/

      sigma_inf_i  += EPSILON;
      sigma_inf_ip += EPSILON;


      /*

	Calculate the new distances from eq. 0.24
	$$
	d_{\bar\imath k} = {{(d_{ik}-b_{i;i'}+\phi_i)/\sigma^2_\infty(i\bar\imath)+
	               (d_{i'j}-b_{i';i}+\phi_{i'})/\sigma^2_\infty(i'\bar\imath)}
		   \over{
		       {1\over\sigma^2_\infty(i'\bar\imath)} +
		       {1\over\sigma^2_\infty(i'\bar\imath)}}}
	$$
	where\hfill\break
	$i=$ \|newNode->child_r->ind|,\hfill\break
	$i'=$ \|newNode->child_l->ind|,\hfill\break
	$b_{i;i'}=$ \|newNode->child_r->rho|,\hfill\break
	$b_{i';i}=$ \|newNode->child_l->rho|


	Also calcuate the renormalization terms $c_{i;j}$ (eq 0.43 ver0.2.5)
	and $c_i$

	*/
      
      for(j=0;j<Nleft;++j) 
	{
	  if(j!=i && j!=ip)
	    {
	      
	      /* $1/\sigma^2_\infty(i\bar\imath)+1/\sigma^2_\infty(i'\bar\imath)$ */
	      double norm = 
		1.0/( 1.0/sigma_inf_i + 1.0/sigma_inf_ip);

	      /*
		First calcuate the new distances
		*/
	      
	      D(Nnext,j) = D(j,Nnext) = 
		norm *
		(
		 (D(i,j)-RHO(newNode->child_r))/(sigma_inf_i)
		 + 
		 (D(ip,j)-RHO(newNode->child_l))/(sigma_inf_ip)
		 );

	      if(D(Nnext,j)<0.0)
		D(Nnext,j) = D(j,Nnext) = 0.0;

	    }
	}
      
      D(Nnext,Nnext) = 0.0;

      /*
	
	And now the new renormalization quantity $c_{\bar\imath}$
	
	N.B. eq 0.30 has been rewritten from
	$$
	{1\over{{1\over X}+{1\over Y}}}
	$$
	to
	$$
	{XY\over X+Y}
	$$
	which is better behaved numerically when $X$ or $Y$ is
	small (and cheeper since it only has one division).
	
	*/

      sig_r = sigma2t(C(i)+DMAX(RHO(newNode->child_r), MINB));
      sig_l = sigma2t(C(ip)+DMAX(RHO(newNode->child_l), MINB));

      if(sig_r*sig_l>0.0)
	{
	  C(Nnext) = sigma2tinv(sig_r*sig_l/(sig_r+sig_l));
	}
      else
	C(Nnext) = sigma2tinv(0.0);

      if(!
	 (C(Nnext)<=DMIN(DMAX(RHO(newNode->child_r),MINB)+C(i)+1e-14,
			    DMAX(RHO(newNode->child_l),MINB)+C(ip)+1e-14)))
	{
	  printf("C(Nnext=%d)=%g\n"
		 "RHO_R=%g C(i=%d)=%g sig_r=%g\nRHO_L=%g C(ip=%d)=%g sig_l=%g -- %g\n",
		 Nnext, C(Nnext),
		 RHO(newNode->child_r), i, C(i), sig_r,
		 RHO(newNode->child_l), ip, C(ip), sig_l,
		 sig_r*sig_l/(sig_r+sig_l));
	} 

      assert((C(Nnext)<=DMIN(DMAX(RHO(newNode->child_r),MINB)+C(i)+1e-14,
			    DMAX(RHO(newNode->child_l),MINB)+C(ip)+1e-14)));
            
      /*
	Swap $i$ node to the empty node at the end of the list and
	place the new node in position $i$ */
      
      nodes[Nnext] = nodes[i];
      nodes[i] = newNode;

      /*
	Swap the $ip$ node and the last node on the list this moves
	$ip$ to the end. When we decrease \|Nleft| by one there will be
	on less node and the two joined nodes $i$ and $ip$ will now be
	after then end (\|Nleft|) of the list
	*/
      
      tmpNode = nodes[ip];
      nodes[ip] = nodes[Nleft-1];
      nodes[Nleft-1] = tmpNode;

      /*
	In the new node set the child indecies to the
	new indexes of the the joined nodes. This info
	will be used by \|sigma2_3| in the renormalization
	step
	*/
	
      newNode->cind_r=Nnext;
      newNode->cind_l=Nleft-1;

      /*
	Set up the \|ta|, \|tb| and \|tc| node array indices.  \|ta|
	and \|tb| point to the two taxa that where just joined, and
	\|tc| points to the newly created taxon.

	These globals will be used in the next call to \|calcb|.
	*/

      ta = Nnext;
      tb = Nleft - 1;
      tc = i;
      
      --Nleft;
      ++Nnext;

      /* 
	 Print out the values of the various variables
	 */

      if(printLevel>2)
	{
	  int a, b;
	  fprintf(outfile, "\nReduced d_ij=\n");
	  for(a=0;a<Nleft;++a)
	    {
	      for(b=0;b<Nleft;++b)
		fprintf(outfile,"%7.4lg ", D(a,b));
	      fprintf(outfile,"\n");
	    }
	  fprintf(outfile,"\n");
	}


      if(printLevel>3) {
	int a, b;

	for(a=0;a<Nnext;++a)
	    {
	      for(b=0;b<Nnext;++b)
		fprintf(outfile,"%7.4lg ", mD[a][b]);
	      fprintf(outfile,"\n");
	    }
	  fprintf(outfile,"\n");
	  
	  fprintf(outfile, "c_i = ");
	  for(a=0;a<Nleft;++a)
	    {
	      fprintf(outfile,"%7.4lg ", C(a));
	    }
	  fprintf(outfile,"\n");
	  
	  for(a=0;a<Nnext;++a)
	    {
	      fprintf(outfile,"%7.4lg ", vC[a]);
	    }
	  fprintf(outfile,"\n");
	  
	  
	  fprintf(outfile, "\n");
	}
    }	    
  
  /*
    
    \section{Final three taxa}
    
    Now there are just three taxa left. They will join to the root
    node of our tree. Find their branch lengths (which we can do
    exactly) and set up the root node to be passed back on return from
    this functin.
    
    */
  
  root = createRootNode();
  if(!root) printError("build::buildTree:out of memory-root");

  root->child_l = nodes[0];
  root->child_m = nodes[1];
  root->child_r = nodes[2];

  /*
    
    Now get the root branch lengths. We can solve this exactly since
    we have three equations and three unknows. The equations to solve
    are:
    $$
    \rho_0+\rho_1 = d_{01},
    \rho_0+\rho_2 = d_{02},
    \rho_1+\rho_2 = d_{12}
    $$
    And the solution is:
    $$
    \rho_0={1 \over 2}\left(d_{01}+d_{02}-d_{12}\right),
    \rho_1={1 \over 2}\left(d_{01}-d_{02}+d_{12}\right),
    \rho_2={1 \over 2}\left(-d_{01}+d_{02}+d_{12}\right)
    $$
    
    */

  root->child_l->rho = 0.5*( D(0,1)+D(0,2)-D(1,2));
  root->child_m->rho = 0.5*( D(0,1)-D(0,2)+D(1,2));
  root->child_r->rho = 0.5*(-D(0,1)+D(0,2)+D(1,2));

  /* check for negative lengths and set to zero if found and decrease
    the other each by half the the negative length (note + a neg
    number is a decrease) */

  if(root->child_l->rho < 0.0)
    {
      root->child_m->rho += 0.5*root->child_l->rho;
      root->child_r->rho += 0.5*root->child_l->rho;
      root->child_l->rho=0.0;
    }
  
  if(root->child_m->rho < 0.0)
    { 
      root->child_l->rho += 0.5*root->child_m->rho;
      root->child_r->rho += 0.5*root->child_m->rho;
      root->child_m->rho=0.0;
    }
  if(root->child_r->rho < 0.0) 
    {
      root->child_l->rho += 0.5*root->child_r->rho;
      root->child_m->rho += 0.5*root->child_r->rho;
      root->child_r->rho=0.0;
    }

  /*
    Clean up
    */

  freeMatrix(mD);

  freeMatrix(b);
  freeMatrix(delta2B);
  freeMatrix(deltaB);
  if(recalcB)
    freeMatrix(oldDeltaB);
  freeMatrix(s);

  freeVector(R);
  freeVector(LLR);
  freeVector(Zscore);

  freeVector(vC);

  freeMatrix(mS);
  freeMatrix(mDelB);
  freeMatrix(mDel2B);

  free(nodes);
  free(q);
  free(q2);

  return(root);

}




/*

  \section{maxVector}

  Return the index of the maximum element of a \|VectorT|

  */

int
maxVector(int num, VectorT v)
{

  int i, maxind = 0;
  double max = v[maxind];

  for(i=1;i<num;++i)
    if(v[i]>max)
      {
	maxind = i;
	max = v[i];
      }

  return(maxind);

}


/* \endc */ 
