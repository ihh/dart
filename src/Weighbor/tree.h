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
  \input epsf

  Program: weighbor                             \hfill\break
  File: tree.h                                  \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 4-Feb-97                             \hfill\break
  
  \title{tree.h}
  \job{Weighbor}
  
  \synopsis{This file contains the various data structures used by the
  weighbor program. In particular those to build the best distance tree}
  

  */

/* $Id: tree.h,v 1.2 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: tree.h,v $
  Revision 1.2  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.15  1998-06-15 19:56:17-04  nds
  First Release Version!

  Revision 1.14  1998-06-09 18:33:41-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.13  1998-02-19 23:52:10-05  nds
  Version 0.2.8.5.2

  Using new formula for $\sigma_\infty$. No longer need to store vaules
  of it in the nodes so removed \|sigma2inf| variables from node
  structure.

  Revision 1.12  1997-10-17 19:22:30-04  nds
  Version 0.2.6: new renormalization code

  Revision 1.11  1997-09-19 18:53:03-04  nds
  Added RHO(n) macro the get the branch length of a node

  Revision 1.10  1997-09-18 16:57:26-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.5  1997-08-23 18:29:53-04  nsocci
  Added some more tree macros (needs to be completed)

  Revision 1.4  1997/08/23 01:16:20  nsocci
  Renormalization code added. Alpha version still needs to be tested
  Note the sigma2_3 function's calling parameteres has been
  changed. Added two booleans to determine whether to renormalize or
  not.

  Revision 1.3  1997/08/22 22:50:50  nsocci
  Added cind_r and cind_l varibles to the NodeT structure to keep track
  of the indexes into the node array of the childern of a node. This is
  necessary for calculating the recursion of sigma2_3 efficiently.

  Revision 1.2  1997/08/22 21:57:40  nsocci
  Added macro to test if node is a leaf

  Revision 1.1  1997/08/08 17:58:30  nds
  Initial revision

  
  */

#ifndef __TREE_H
#define __TREE_H

#include <stdio.h>
#include "Weighbor/matrix.h"

/* \section{NodeT}

  \|NodeT| is the basic binary tree node. It will hold not only the
  pointer to its childern but the length corresponding to the branch
  from the node to its immediate ancestor (father node). 

  The terminal nodes will hold the name of the taxon's they represent
  (read in from the input file). The non-terminal nodes will have
  \|NULL| pointers here for the name string, this will be used to
  distinguish the terminals from the non-terminal nodes. 

  Each non-terminal node will also hold the renormalized
  $\sigma^2_\infty(i\bar{i})$ and $\sigma^2_\infty(i'\bar{i})$. This
  value will be used in the recursive calculation of the renormalized
  weights.

  */

typedef double LengthT;   /* Type representing the length of a branch */

typedef struct _NodeT {

  char *name;            /* The name of the taxon, \|NULL| for
			    non-terminal */

  int ind;      /* This is used to hold the nodes position in the 
		   distance matrix. */

  LengthT rho;  /* The length of the branch leading back to the nodes
                   father */

  struct _NodeT *child_r, *child_l;  /* Pointers to the childern */

  int cind_r, cind_l; /* index numbers in to the \|nodes| array for the
			 two childern */
  
} NodeT;


/*

  \section{RootNodeT}

  The evolutionary tree is not quite a binary tree. It really is an
  un-rooted tree (graph) of tri-nodes; ie each vertex has 3
  edges. This can be represented as a rooted binary tree but the only
  problem is the root node must be a tri-node. 
  
  \vskip1cm\centerline{\epsfxsize=4in\epsfbox{tree.eps}}

  Note: any non-terminal node in the unrooted tree can serve as the
  root. The root node will be determined by how the joining algorithm
  progresses and will be created by the last three unjoined taxa.

  \|RootNodeT| is the root node data structure and will have 3 child
  pointers. Since it is a non-terminal node it will not have any name
  field. And since it has no ancestor it does not have a length field.

  And once we have the root we do not have to worry about
  renormalizing anymore weights so do not need the \|sigma2inf| values
  either.  

*/
  
typedef struct {

  NodeT *child_r, *child_m, *child_l;

} RootNodeT;


/*

  \section{Function Defs}

   Here are the various function definitions for processing the tree.

   */

NodeT *createNode();
RootNodeT *createRootNode();
void printTree(FILE*, RootNodeT *);
void deleteTree(RootNodeT *);
void deleteNode(NodeT *);

/*

  \section{Tree macros}

  \subsection{Operations on node pointers}

  These macros take as arguments pointers to the \|NodeT| struct.

  */

#define LEAF(n) (!((n)->child_r))   /* test if node \|n| is a leaf */
#define RHO(n)  ( (n)->rho )

#define RIND(n) ( (n)->child_r->ind ) /* return the left and right  */
#define LIND(n) ( (n)->child_l->ind ) /* child node indicies        */ 

#define RIGHT(n) ( (n)->child_r )
#define LEFT(n) ( (n)->child_l )

/*

  \subsection{Operations on node index numbers}

  These macros take as arguments index numbers into the global
  \|nodes| array.

  */


#define SIGMA2INF_R(nidx) ( nodes[(nidx)]->sigma2inf_r )
#define SIGMA2INF_L(nidx) ( nodes[(nidx)]->sigma2inf_l )

#define IND(nidx)  ( nodes[(nidx)]->ind )


/*

  This macro access the $c_{i;j}$ array stored in each node
  $j$=\|nidx| must be a index into the \|nodes| array
  

  */

#endif /* __TREE_H */

/* \endc */





