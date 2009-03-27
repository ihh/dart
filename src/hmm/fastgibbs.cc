#include "hmm/fastgibbs.h"

Fast_Gibbs_align::Fast_Gibbs_align (const Alphabet& alph, const Sequence_database& seqdb, unsigned int len, double pseud)
  : alph (alph),
    mul (alph.degenerate_lcd),
    len (len),
    pseud (pseud),
    count (len, vector<unsigned int> (alph.size(), (unsigned int) 0)),
    null (alph.size(), (double) 0),
    pcopy (3)
{
  // copy dsq's into align structure, & initialise pos's
  for_const_contents (Sequence_database, seqdb, np) {
    dsq.push_back (&(*np).dsq);
    pos.push_back (vector<unsigned int>());
  }
  // calculate composition
  seqdb.count_symbol_frequencies (null);
  NormalisePr (null);
  // set up pcopy
  pcopy[0] = 1;
  pcopy[1] = .5;
  pcopy[2] = 0;
}

Fast_Gibbs_matrix::Fast_Gibbs_matrix (Fast_Gibbs_align& align, unsigned int row)
  : align(align),
    row(row),
    dsq(*align.dsq[row])
{
  init_motif();
  fill_matrix();
}


void Fast_Gibbs_matrix::init_motif()
{
  motif = vector<vector<Score> > (align.len, vector<Score> (align.alph.size(), -InfinityScore));
  for (int m = 0; m < (int) motif.size(); ++m) {
    vector<Prob> weight (align.alph.size());
    for (int s = 0; s < align.alph.size(); ++s)
      weight[s] = align.pseud * align.null[s] + ((double) align.count[m][s]) / ((double) align.mul);
    NormalisePr (weight);
    for (int s = 0; s < align.alph.size(); ++s) motif[m][s] = Prob2Score (weight[s] / align.null[s]);
  }
}

void Fast_Gibbs_matrix::fill_matrix()
{
  const int endpos = dsq.size() - motif.size();
  matrix = vector<Score> (endpos);
  for (int p = 0; p < endpos; ++p) {
    Score sc = 0;
    for (int i = 0; i < (int) motif.size(); ++i) ScorePMulAcc (sc, motif[i][dsq[p+i]]);
    matrix[p] = sc;
  }
}
