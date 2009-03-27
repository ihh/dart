#ifndef FOLDPIC_INCLUDED
#define FOLDPIC_INCLUDED

#include <gd.h>

#include "scfg/paircfgdp.h"

// fold picture structure
struct Fold_pic
{
  static const int maxcol;  // range of red, green indices
  gdImagePtr image;
  array2d<int> rbcol;  // color indexed by red, blue
  int blackcol;  // black color (background)
  int greencol;  // single green color
  int yellowcol;  // single yellow color
  vector<int> seqcol;  // red, green, yellow, blue = A,C,G,U
  double sigmoid_offset;
  double sigmoid_scale;
  // methods
  Fold_pic (const Named_profile& npx, const Pair_inside_outside_matrix& xmat, const Fold_envelope& xenv,
	    const Named_profile& npy, const Pair_inside_outside_matrix& ymat, const Fold_envelope& yenv,
	    const Pair_CFG_scores& pair_cfg, const Pair_HMM_scores& pair_hmm,
	    double sig_offset = 0., double sig_scale = 1.);
  ~Fold_pic();
  void save (const char* filename) const;
  int color (double rprob, double bprob) const;
  void draw_foldenv (const Named_profile& np, const Pair_inside_outside_matrix& mat, const Fold_envelope& env,
		     int xoffset, int yoffset, int xvec, int yvec);
  double sigmoid (double x) const;
};

// HMM for probabilistic Smith-Waterman
struct Foldpic_HMM : Pair_HMM_scores
{
  enum { LeftPadX = 0, LeftPadY = 1, MatchXY = 2, InsertY = 3, DeleteX = 4, RightPadX = 5, RightPadY = 6, TotalStates = 7 };
  Foldpic_HMM (const array2d<Score>& submat, Prob gap_open_prob, Prob gap_extend_prob, const Alphabet& alph);
};

// simple substitution matrix
struct Simple_submat : array2d<Score>
{
  Simple_submat (int alph_sz, Prob match_prob, const vector<Prob>& null_prob);
};

#endif /* FOLDPIC_INCLUDED */
