#ifndef AMAP_INCLUDED
#define AMAP_INCLUDED

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>
#include <cctype>

#include "amap/SafeVector.h"
#include "amap/MultiSequence.h"
#include "amap/SparseMatrix.h"

#include "amap/dotplot.h"

typedef map<int,map<int,Dotplot> > Dotplot_map;
struct AMAP_adapter {

  /// sequence database
  const FASTA_sequence_database& seq_db;

  /// number of seqs
  int num_seqs;

  /// posterior probability matrices
  const Dotplot_map& dotplots;

  /// constructor
  AMAP_adapter (const FASTA_sequence_database& seq_db, const Dotplot_map& dotplots);

  /// Return the Stockholm alignment produced by AMAP.
  // Must be a Stockade object so that it stores its own Named_profile's.
  Stockade get_alignment (float gapFactor = 1, float edgeWeightThreshold = 0, int numConsistencyReps = 0, bool show_intermediates = false, bool output_for_gui = false);


private:

  /// controls verbosity of AMAP's internal routines
  bool enableVerbose;

  /* Data for AMAP's internal use */
  // sequences (for AMAP)
  MultiSequence sequences;
  // sparse post. prob. matrices (for AMAP)
  SafeVector<SafeVector<SparseMatrix *> > sparseMatrices;

  /* Convert to AMAP's preferred format */
  /// Convert from a dotplot to a posteriorVector.
  VF* posteriorVector_from_dotplot (const Dotplot& dotplot);

  /// consistency checking
  SafeVector<SafeVector<SparseMatrix *> > DoRelaxation ();
  void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior);

};


#endif /* AMAP_INCLUDED */
