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
  File: weighbor.c                              \hfill\break
  Author: N. D. Socci                           \hfill\break
  
  \title{weighbor.c}
  \job{Weighbor}
  
  \synopsis{\"Weigh"ted neigh\"bor"-joining.
This program will created a evolutionary tree based on distance
between sequences. The tree built attempts to minimize the distances
in it versus the actually distances inputted to the program. The
output is the tree built along with the branch lengths.}


*/

/*

\section{Introduction}

The default input file is: \|infile|

The default output file is \|treefile|

The algorithm used is a weighted neighbor-joining one from W. Bruno
(\'refs'). The weights are computed by the function \|weight(d1,d2)|.

Important note. Array indices are from $0\ldots N-1$ in contrast to
the first version which went from $1\ldots N$

*/

/*
  
  $Log: weighbor.c,v $
  Revision 1.3  2004/09/23 17:00:18  ihh
  fixed paths

  Revision 1.2  2003/12/15 12:28:06  ihh
  Fixed so gcc can compile with -Wall.

  Revision 1.1.1.1  2003/11/10 19:07:59  ihh
  Initial import.

  Revision 1.1  2003/05/16 16:24:22  ihh
  weighbor

  Revision 1.32  1998-06-15 20:02:28-04  nds
  Changes warnings to default off (use -Xw to turn on)

  Revision 1.31  1998-06-15 19:56:18-04  nds
  First Release Version!

  Revision 1.30  1998-06-09 18:33:41-04  nds
  Version 0.2.13
      Added -Xx option
      Changed the normalization on the -Xn option
      Copyright notice added to all files

  Revision 1.29  1998-06-04 18:35:22-04  nds
  Update for version 2.12
      Changed z score
      Change renormalization if z-score negative
      non-deterministic warning to stderr
      -Xn flag to change R(i,j,j') normalization

  Revision 1.28  1998-05-20 17:52:38-04  nsocci
  Bug fix in usage statements

  Revision 1.27  1998/05/20 21:14:25  nsocci
  Initial release 1.0.0 alpha

  Revision 1.26  1998/05/20 19:59:31  nsocci
  Added -b option to change the value of B and changed the default
  setting of the -S option (now default is off, -S turn on barred
  values).

  Revision 1.25  1998/05/20 01:43:45  nsocci
  Added the extended tournament option (-e)

  Revision 1.24  1998/05/19 17:11:37  nsocci
  Added -S and -T options do control the use of the average (barred)
  quanitites in the calculation of R(i,j,j')

  Revision 1.23  1998/05/08 19:26:33  nds
  Added -V flag to print version information

  Revision 1.22  1998-05-01 16:57:20-04  nds
  Made numerous changes to log file and adjusted what was printed at
  various verbose levels. To print the new LLR value of R(q(i),i,j) the
  calcR funtion need to be called from buildTree which require the
  moving of the deltaB and associated matrices to global level. This
  should improve performance since now they are allocated only once per
  tree rather than a every iteration

  Revision 1.21  1998-04-29 19:02:49-04  nds
  Added -i and -o option to set both the input and output. Default is
  still stdin and stdout. Also print a warning if the -L option is not
  used to explicitly set the length.

  Revision 1.20  1998-04-14 16:54:38-04  nds
  Added -z command line switch to use toe old z score routine.
  Changed the meaning of the -r flag. Default is not to do the Delta b
  recalculation, -r now turns this on.

  Revision 1.19  1998-04-08 19:23:43-04  nds
  $\Delta b$ optimization: calculation in $N^2$ time.

  Revision 1.18  1998-03-23 23:51:52-05  nds
  Fixed bug with -L option in case user types -Lnum
  which will set length to zero

  Revision 1.17  1998-03-20 00:02:17-05  nds
  Added -r option to turn of the recalculation of $\Delta b$ and
  associated variables.

  Revision 1.16  1998-03-13 16:52:44-05  nds
  Added -v command line option
  BETA release version

  Revision 1.15  1998-02-06 16:08:52-05  nds
  Fixed syntax error

  Revision 1.14  1998-02-06 16:04:51-05  nds
  Added the constant \|MINB| used in the calculation of the b's (eq
  0.13-0.15, 0.22-0.24)

  Revision 1.13  1997-12-16 15:20:49-05  nds
  Added the checkQQI boolean variable and the -q command line switch to
  turn of the checking of q[q[i]]=i in selection the best pair.

  Revision 1.12  1997-11-15 00:32:06-05  nds
  Changes for weighor version 0.2.7.1

  Revision 1.11  1997-09-19 18:47:03-04  nds
  Added LOGFILE compile flag to switch the logging of extra info to the
  weighbor.out file on and off.

  Revision 1.10  1997-09-18 16:57:59-04  nds
  Final 1.x version of weighbor. Uses approximate residual computation
  to calculate the pairs to join and averages the various branch lengths
  estimates to get the new branch lengths for the joinded pairs.

  Revision 1.1  1997-08-08 13:59:06-04  nds
  Initial revision


  */


#define __WEIGHBOR_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Weighbor/tree.h"
#include "Weighbor/matrix.h"
#include "Weighbor/io.h"
#include "Weighbor/weighbor.h"

BooleanT expertMode=False;

void
usage()
{
  fprintf(stderr, VERSION_INFO);
  fprintf(stderr, "usage: weighbor [-L num] [-b num] [-i name] [-o name] [-v[v[v]]] [-V]\n");
  fprintf(stderr, "         -L num   set the sequence length to num (default %g)\n", DEFAULT_LENGTH);
  fprintf(stderr, "         -b num   set number of bases (default %g)\n", DEFAULT_B);
  fprintf(stderr, "         -i name  read input from file\n");
  fprintf(stderr, "         -o name  write output to file\n");
  fprintf(stderr, "         -v       turn on logfile and dump level 1 info\n");
  fprintf(stderr, "         -vv      turn on logfile and dump level 2 info\n");
  fprintf(stderr, "         -vvv     turn on logfile and dump level 3 info\n");
  fprintf(stderr, "         -V       print version information\n");
  if(expertMode) {
    fprintf(stderr, "\n  Advance (expert) options [-XqrzSTe]\n");
    fprintf(stderr, "         -X       turn on expert mode (allows use of these options)\n");
    fprintf(stderr, "           n      change normalization of R(i,j,j')\n");
    fprintf(stderr, "           q      turn off q[q[i]]==i checking\n");
    fprintf(stderr, "           r      turn on recalculation of Delta B\n");
    fprintf(stderr, "           z      use the old (N^3) method to calculate z(i,j)\n");
    fprintf(stderr, "           S      use barred version of sigma for R(i,j,j')\n");
    fprintf(stderr, "           T      use unbarred versions of d, Delta B and S for R(i,j,j')\n");
    fprintf(stderr, "           e      extended tournament option\n");
    fprintf(stderr, "           w      turn on non-deterministic warnings\n");
    fprintf(stderr, "           x      change normaliztion of R(i,j)\n");
  }
  exit(0);
}

int
main(int argc, char *argv[])
{

  int i;
  FILE *fp=stdout, *inputfp=stdin;

  int N;            /* Number of taxa  */
  MatrixT D;        /* Distance Matrix */
  char **taxon;     /* names of taxa   */
  BooleanT lengthFlag = False; /* set if the -L option was given */

  RootNodeT *tree;

  RootNodeT *buildTree(int, MatrixT, char **);

  for(i=1;i<argc;++i)
    {
      char option, optidx;

      if(argv[i][0]=='-')
	option = argv[i][1];
      else
	usage();

      optidx=1;
      while(optidx && (option=argv[i][(int) optidx++])) 
	{

	  switch(option)
	  {
	      
	    case 'V':
	      fprintf(stderr, VERSION_INFO);
	      break;
	      
	    case 'i':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		inputfp = fopen(argv[i], "r");
		if(!inputfp)
		  {
		    perror(argv[i]);
		    exit(1);
		  }
	      }
	      else
		usage();
	      
	      break;
	      
	    case 'o':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		fp = fopen(argv[i], "w");
		if(!fp)
		  {
		    perror(argv[i]);
		    exit(1);
		  }
	      }
	      else
		usage();
	      
	      break;
	      
	    case 'L':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		L = atof(argv[i]);
		if(L<=0.0)
		  {
		    fprintf(stderr, "Must specify a length eg: -L 1000\n");
		    fprintf(stderr, "Note the space after the `L'\n\n");
		    usage();
		  }
		EPSILON = (EPS_BARE)/L;
		MINB    = (MINB_BARE)/L;
		lengthFlag = True;
	      }
	      else 
		usage();
	      
	      break;
	      
	    case 'b':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		B = atof(argv[i]);
		if(L<=0.0)
		  {
		    fprintf(stderr, 
			    "Must specify a number greater than zero eg: -b 2.3\n");
		    fprintf(stderr, "Note the space after the `b'\n\n");
		    usage();
		  }
	      }
	      else 
		usage();
	      
	      break;
	      
	    case 'v':
	      ++printLevel;
	      break;

	    case 'X':
	      expertMode = True;
	      optidx=2;
	      while((option=argv[i][(int) optidx++]))
		{
		  printf("Expert Option [%c]\n", option);
		  switch(option)
		    {
		      
		    case 'e':
		      extendedTourn=True;
		      break;
		      
		    case 'S':
		      useSigmaBar=True;
		      break;
		      
		    case 'T':
		      useBarValues=False;
		      break;

		    case 'n':
		      n_Flag = True;
		      break;
		      
		    case 'q':
		      checkQQI = False;
		      break;
		      
		    case 'r':
		      recalcB = True;
		      break;
		      
		    case 'w':
		      w_Flag = True;
		      break;

		    case 'x':
		      x_Flag = True;
		      break;

		    case 'z':
		      oldZflag = True;
		      break;
		      
		    default:
		      usage();
		    }
		}
	      optidx=0;
	      break;
	      
	    default:
	      usage();
	}
      }
    }

  if(!lengthFlag)
    fprintf(stderr, "No length was set. Using the default L=%g\n\n", L);

  if(recalcB && n_Flag) {
    fprintf(stderr, "WARNING: Not recommended that you use the -Xr and -Xn flags together\n");
  }

  /*
    
    Read in the Phylip data file. 

    */

  /* open the auxillary output file */
  
  if(printLevel>0)
    {    
      outfile = fopen("weighbor.out", "w"); 
      if(!outfile)
	printError("weighbor::main::open(outfile)");
    }
  else 
    outfile = NULL;


  while(readPhylipFile(inputfp, &N, &D, &taxon))
    {

      int i;
      
      if(printLevel>0) {
	fprintf(outfile, "\n============================================");
	if(printLevel>1) {
	  fprintf(outfile, "\nL = %g, b = %g\n", L, B);
	  fprintf(outfile,"\nInput Distance Matrix = \n");
	  for(i=0;i<N;++i) 
	    {
	      int j;
	      for(j=0;j<N;++j)
		fprintf(outfile, "%7.4lg ", D[i][j]);
	      fprintf(outfile, "\n");
	    }
	  fprintf(outfile, "\n");
	}
      }

      /*
	Reset any global switch and variables
	*/

      warnFlag=!(w_Flag);

      switch(N) {
	
	/* Handle trival cases separately */
	
      case 1: /* 1 taxon; just print the name*/
	
	fprintf(fp, "(%s)\n", taxon[0]);
	if(printLevel>0)
	  fprintf(outfile, "\n(%s)\n\n", taxon[0]);
	  
	break;
	
      case 2: /* 2 taxon; trivial tree with each branch $=D_{0,1}/2$ */
	
	fprintf(fp, "(%s:%f,%s:%f)\n", 
		taxon[0], .5*D[0][1], taxon[1], .5*D[0][1]);
	if(printLevel>0)
	  fprintf(outfile, "\n(%s:%f,%s:%f)\n\n", 
		  taxon[0], .5*D[0][1], taxon[1], .5*D[0][1]);
	break;
	
	
      default: /* 3 or more need to build the tree */

#ifdef HASH
	hashsize = 2*N-3;
	hash = (double *)malloc(hashsize*hashsize*hashsize*sizeof(double));
	if(!hash)
	  {
	    fprintf(stderr, "Can not allocate hash table\n");
	    exit(1);
	  }

	/*
	  fprintf(stderr, "Hashsize = %d %d\n", hashsize,
		hashsize*hashsize*hashsize);
		*/

	for(i=0;i<hashsize*hashsize*hashsize;++i)
	  hash[i] = -1.0;
#endif

	tree = buildTree(N, D, taxon);
	printTree(fp, tree);
	if(printLevel>0) {
	  fprintf(outfile, "\n");
	  printTree(outfile, tree);
	}
	deleteTree(tree);

#ifdef HASH
	free(hash);
	/*
	  fprintf(stderr, "total=%d hash hits=%d percent=%lf%%\n",
		total, hashcnt, 100*((double)hashcnt)/((double)total));
		*/

#endif

	/*
	  Now free the memory used by the name strings
	  */

	for(i=0;i<N;++i)
	  free(taxon[i]);

	break;

      }
    }


  if(printLevel>0)
    fclose(outfile);
  fclose(fp);

  return 1;
}


/* \endc */ 

