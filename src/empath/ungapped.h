#ifndef UNGAPPED_INCLUDED
#define UNGAPPED_INCLUDED

#include "hmm/singlehmm.h"
#include "empath/trainer.h"
#include "seq/dirichlet.h"
#include "seq/gff.h"

// A simple ungapped motif
struct Ungapped_model : Trainable
{
  int min_size;
  int max_size;
  bool revcomp;

  vector<Alphabet_group> match_emit;  // the motif
  PGroup skip_start;  // start transitions, i.e. no. of motif states to skip at start
  PGroup skip_end;    // end transitions, i.e. no. of motif states to skip at end

  // constructor
  Ungapped_model (int min_size, int max_size, bool revcomp, PScores& pscore, bool mask = 1);
};

// Restricted subclass for fixed-length motifs
struct Fixed_model : Ungapped_model
{
  Fixed_model (int size, bool revcomp, PScores& pscore, bool mask = 1);
  // accessors
  inline int size() const { return max_size; }
  // fast Viterbi/forward method (only difference is in merging of fwd/rev strands)
  void get_local_scores (const Named_profile& np, Metascore& result, bool viterbi = 0);
};

#endif  /* UNGAPPED_INCLUDED */
