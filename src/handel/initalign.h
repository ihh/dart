#ifndef HANDEL_INIT_ALIGN_INCLUDED
#define HANDEL_INIT_ALIGN_INCLUDED

// Code for quick initial estimation of alignment & tree in Handel

#include "handel/alitrans.h"
#include "amap/dotplot.h"
#include "amap/amap_adapter.h"
#include "hmm/pairhmm.h"
#include "hmm/postpairhmm.h"
#include "tree/nj.h"
#include "tree/subdistmat.h"

struct Stockade_initializer {
  // Alignment algorithm
  static Stockade align (FASTA_sequence_database& seq_db, const Alphabet& alphabet, Pair_transducer_factory& trans_fac, double branch_len);

  // Tree algorithm
  static void build_tree (Stockholm& stock, Substitution_matrix_factory& submat_factory);
};

#endif /* HANDEL_INIT_ALIGN_INCLUDED */
