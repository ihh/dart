#ifndef FASTGIBBS_INCLUDED
#define FASTGIBBS_INCLUDED

// This module is flotsam, not really used by anything, but a half-finished attempt to make Gibbs sampling work fast.

#include "seq/biosequence.h"
#include "util/score.h"

struct Fast_Gibbs_align
{
  const Alphabet&                      alph;   // the alphabet
  const unsigned int                   mul;    // the multiplier for residue counts (12 for RNA)
  unsigned int                         len;    // motif length
  double                               pseud;  // pseudocount multiplier (pseudocount for symbol S = null[S] * pseud)
  vector<const Digitized_biosequence*> dsq;    // dsq[R] = sequence for row R
  vector<vector<unsigned int> >        pos;    // pos[R][N] = index of N'th motif instance in row R
  vector<vector<unsigned int> >        count;  // count[C][S] = mul * count for symbol S in motif column C
  vector<double>                       null;   // null[S] = null model probability for symbol S
  vector<double>                       pcopy;  // pcopy[N] = probability of emitting N or more copies of motif in any sequence

  // constructor
  // ensure Sequence_database has been converted to dsq's before calling
  //
  Fast_Gibbs_align (const Alphabet& alph, const Sequence_database& seqdb, unsigned int len, double pseud);
};

class Fast_Gibbs_matrix
{
public:
  Fast_Gibbs_align&            align;   // the alignment
  const unsigned int           row;     // index of sampling row
  const Digitized_biosequence& dsq;     // the sequence
  vector<vector<Score> >       motif;   // motif[C][S] = score for symbol S in column C
  vector<Score>                matrix;  // matrix[P] = likelihood of motif pos P

  Fast_Gibbs_matrix (Fast_Gibbs_align& align, unsigned int row);
private:
  void init_motif();
  void fill_matrix();
};


#endif
