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
  File: io.h                                    \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{io.h}
  \job{Weighbor}
  
  \synopsis{File I/O Routines}
  

  */

/* $Id: io.h,v 1.2 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: io.h,v $
  Revision 1.2  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.12  1998-06-15 19:56:16-04  nds
  First Release Version!

  Revision 1.11  1998-06-09 18:33:39-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.10  1997-09-18 16:55:42-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.1  1997-08-08 13:56:45-04  nds
  Initial revision


  */

#ifndef __IO_H
#define __IO_H

#include "Weighbor/matrix.h"

FILE * openRead(char *filename);
FILE * openWrite(char *filename);
void printError(char *);
int readPhylipFile(FILE *fp, int *N, MatrixT *D, char **names[]);

#endif /* __IO_H */

/* \endc */
