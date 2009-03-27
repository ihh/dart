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
  File: io.c                                    \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{io.c}
  \job{Weighbor}
  
  \synopsis{File I/O Routines}
  

  */

/* $Id: io.c,v 1.2 2004/09/23 17:00:18 ihh Exp $ */

/*
  
  $Log: io.c,v $
  Revision 1.2  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.13  1998-06-15 19:56:16-04  nds
  First Release Version!

  Revision 1.12  1998-06-09 18:33:39-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.11  1998-05-18 17:49:04-04  nsocci
  Added more information to error messages

  Revision 1.10  1997/09/18 20:55:12  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.1  1997-08-08 13:59:06-04  nds
  Initial revision

  
  */


/*

  \section{Introduction}

  These routines handle the basic file i/o for the weighbor
  program. The the expection of the tree printing routines which are
  in the \|tree| module.


  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Weighbor/matrix.h"

/*

  \subsection{File opening}
  
  Open a file for reading or writing and do a fp pointer
  check. 

  */

FILE *
openRead(char *filename)
{

  FILE *fp;

  fp = fopen(filename, "r");
  if(!fp)
    {
      fprintf(stderr, "%s could not be open for reading\n", filename);
      perror("weighbor::io::openRead");
      exit(1);
    }

  return(fp);

}

FILE *
openWrite(char *filename)
{

  FILE *fp;

  fp = fopen(filename, "w");
  if(!fp)
    {
      fprintf(stderr, "%s could not be open for writing\n", filename);
      perror("weighbor::io::openWrite");
      exit(1);
    }

  return(fp);

}

/*

  \subsection{Read Data File}

  Read in a Phylip data file. 

  The matrix D is reassigned by this routine so it should not be
  pointing to any allocated memory otherwise that memory will be
  lost.
  
  If there is an error in the file the function prints and error and
  exits the program

  \item\ fp - is a filepointer to the file to read from. Had better not
  be null

  \item\ N - Number of taxz (note this is a pointer to an int)

  \item\ D - The distance matrix $N\times N$; it is automatically
  symmetrized by this routine (again a pointer since we change it in
  this function)

  \item\ names - An array of strings holding the taxa
  names. Memory dynamically allocated be the routine

  */


/*
  print an error and die
  */

void
printError(char *s)
{
  perror("weigbor::io");
  fprintf(stderr, "%s\n", s);
  exit(1);
}


int
readPhylipFile(FILE *fp, int *N, MatrixT *D, char **names[])
{
  
  int i, j;
  double tmp;
  char namebuf[81];

  if(fscanf(fp, "%d", N)!=1) /* N already a pointer */
    {
      return(0);
    }

  *D = matrix(*N);

  /* allocate N string pointers, note we index from $0\ldots N-1$ */

  *names = (char **)malloc((*N)*sizeof(char *));
  if(!(*names)) printError("File format error-Out of memory");

  for(i=0;i<(*N);++i) {

    /* read the taxon name */

    if(fscanf(fp, "%s", namebuf)!=1) printError("File format error-1");

    /* allocate enough memory and move string from tmp buf to name
       array */

    (*names)[i] = (char *)malloc((strlen(namebuf)+1)*sizeof(char));
    if(!((*names)[i])) printError("File format error-Out of memory");
    strcpy((*names)[i], namebuf);
    
    /* read in the $i^{\rm th}$ row */

    for(j=0;j<(*N);++j)
      {
	if(fscanf(fp, "%lf", &tmp)!=1) printError("File format error-2");
	(*D)[i][j] = tmp;
      }

  }

  /* Symmetrize Matrix */

  for(i=0;i<(*N);++i) 
    for(j=0;j<(*N);++j)
      if(i!=j) (*D)[i][j] = (*D)[j][i] = .5*((*D)[i][j] + (*D)[j][i]);

  return(1);
  
}



      
/* \endc */ 




