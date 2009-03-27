#include <errno.h>
#include "gfx/foldpic.h"
#include "hmm/postpairhmm.h"

const int Fold_pic::maxcol = 15;

Fold_pic::Fold_pic (const Named_profile& npx, const Pair_inside_outside_matrix& xmat, const Fold_envelope& xenv,
		    const Named_profile& npy, const Pair_inside_outside_matrix& ymat, const Fold_envelope& yenv,
		    const Pair_CFG_scores& pair_cfg, const Pair_HMM_scores& pair_hmm,
		    double sig_offset, double sig_scale)
  : rbcol (maxcol, maxcol),
    seqcol (4),
    sigmoid_offset (sig_offset),
    sigmoid_scale (sig_scale)
{
  // print log message
  CTAG(4,FOLDPIC) << "Creating fold pic for '" << npx.name << "' and '" << npy.name << "'\n";

  // create GD image
  CTAG(3,FOLDPIC) << "Allocating GD image\n";
  const int size = npx.size() + npy.size() + 1;
  image = gdImageCreate (size, size);
  
  // allocate colors
  blackcol = gdImageColorAllocate (image, 0, 0, 0);
  greencol = gdImageColorAllocate (image, 0, 255, 0);
  yellowcol = gdImageColorAllocate (image, 255, 255, 0);

  for (int r = 0; r < maxcol; ++r)
    for (int b = 0; b < maxcol; ++b)
      rbcol(r,b) = gdImageColorAllocate (image, r * 255 / (maxcol-1), 0, b * 255 / (maxcol-1));

  seqcol[0] = rbcol (maxcol-1, 0);
  seqcol[1] = greencol;
  seqcol[2] = yellowcol;
  seqcol[3] = rbcol (0, maxcol-1);

  // draw X fold envelope
  CTAG(3,FOLDPIC) << "Drawing X fold envelope\n";
  draw_foldenv (npx, xmat, xenv, 0, npx.size() - 1, 1, 0);

  // draw Y fold envelope
  CTAG(3,FOLDPIC) << "Drawing Y fold envelope\n";
  draw_foldenv (npy, ymat, yenv, npx.size() + 1, npx.size() + 1, 0, 1);

  // do Pair SCFG inside-outside to determine secondary structure alignment coloring
  CTAG(3,FOLDPIC) << "Doing inside-outside algorithm for structural alignment posterior probabilities\n";
  array2d<double> cfgpost (npx.size(), npy.size(), 0.);
  const Pair_inside_outside_matrix io_matrix (npx, npy, xenv, yenv, pair_cfg, TRUE);  // assumes local DP
  for (int i = 0; i < (int) xenv.subseq.size(); ++i)
    {
      const Subseq& xsubseq = xenv.subseq[i];
      for (int j = 0; j < (int) yenv.subseq.size(); ++j)
	{
	  const Subseq& ysubseq = yenv.subseq[j];
	  Score sc = -InfinityScore;
	  for (int state = 0; state < pair_cfg.states(); ++state)
	    ScorePSumAcc (sc, io_matrix.post_state_sc (state, i, j));
	  const Prob p = Score2Prob (sc);
	  cfgpost (xsubseq.start, ysubseq.start) += p;
	  cfgpost (xsubseq.start + xsubseq.len - 1, ysubseq.start + ysubseq.len - 1) += p;
	}
    }
  
  // do Pair HMM forward-backward to determine primary sequence alignment coloring
  CTAG(3,FOLDPIC) << "Doing forward-backward algorithm for primary alignment posterior probabilities\n";
  const Pair_forward_backward_DP_matrix fb (pair_hmm, npx.prof_sc, npy.prof_sc);
  const Post_pair_HMM post_matrix (fb);  // get the posterior match probabilities

  // paint the foldpic
  for (int i = 0; i < npx.size(); ++i)
    for (int j = 0; j < npy.size(); ++j)
      gdImageSetPixel (image, i, j + npx.size() + 1,
		       color (sigmoid(cfgpost(i,j)), sigmoid(post_matrix(i,j))));
}

Fold_pic::~Fold_pic()
{
  gdImageDestroy (image);
}

void Fold_pic::save (const char* filename) const
{
  CTAG(3,FOLDPIC) << "Saving fold pic\n";
  FILE* file = fopen (filename, "w");
  if (file == NULL)
    THROWEXPR ("Couldn't open file '" << filename << "' for writing: " << strerror(errno) << "\n");
  gdImagePng (image, file);
  if (fclose (file) != 0)
    THROWEXPR ("Couldn't close file '" << filename << "': " << strerror(errno) << "\n");
}

int Fold_pic::color (double rprob, double bprob) const
{
  return rbcol (minmax ((int) (.5 + rprob * (double) maxcol), 0, maxcol-1),
		minmax ((int) (.5 + bprob * (double) maxcol), 0, maxcol-1));
}

void Fold_pic::draw_foldenv (const Named_profile& np, const Pair_inside_outside_matrix& mat, const Fold_envelope& env,
			     int xoffset, int yoffset, int xvec, int yvec)
{
  // draw sequence
  for (int i = 0; i < np.size(); ++i)
    gdImageSetPixel (image, xoffset - yvec + xvec * i, yoffset + xvec + yvec * i, seqcol[np.dsq[i]]);

  // draw fold envelope
  array2d<int> env_flag (np.size(), np.size(), 0);
  for_const_contents (vector<Subseq>, env.subseq, s)
    env_flag (s->start, s->len) = 1;

  for (int i = 0; i < (int) mat.inside.xenv.subseq.size(); ++i)
    {
      const Subseq& subseq = mat.inside.xenv.subseq[i];
      Score sc = -InfinityScore;
      for (int state = 0; state < mat.inside.cfg.states(); ++state)
	ScorePSumAcc (sc, mat.post_state_sc (state, i, 0));
      
      const int u = subseq.len / 2 + subseq.start;
      const int v = subseq.len;

      gdImageSetPixel (image, xoffset + u * xvec + v * yvec, yoffset + u * yvec - v * xvec,
		       color ((double) env_flag (subseq.start, subseq.len), sigmoid (Score2Prob (sc))));
    }
}

double Fold_pic::sigmoid (double x) const
{
  double sig = 1. / (1. + exp (-sigmoid_scale * (x + sigmoid_offset)));
  return minmax (sig, 0., 1.);
}

Foldpic_HMM::Foldpic_HMM (const array2d<Score>& submat, double gap_open_prob, double gap_extend_prob, const Alphabet& alph)
  : Pair_HMM_scores (TotalStates, &alph)
{
  // set state types
  state_type[LeftPadX] = state_type[DeleteX] = state_type[RightPadX] = EmitX;
  state_type[LeftPadY] = state_type[InsertY] = state_type[RightPadY] = EmitY;
  state_type[MatchXY] = EmitXY;
  // set transitions
  // left pad
  start[LeftPadX] = start[LeftPadY] = start[MatchXY] = 0;
  transition (LeftPadX, LeftPadY) = 0;
  transition (LeftPadX, LeftPadX) = transition (LeftPadY, LeftPadY) = 0;
  transition (LeftPadX, MatchXY) = transition (LeftPadY, MatchXY) = 0;
  // match/insert/delete
  transition (MatchXY, MatchXY) = Prob2Score (1 - gap_open_prob);
  transition (InsertY, MatchXY) = transition (DeleteX, MatchXY) = Prob2Score (1 - gap_extend_prob);
  transition (MatchXY, InsertY) = transition (MatchXY, DeleteX) = Prob2Score (gap_open_prob / 2.);
  transition (InsertY, InsertY) = transition (DeleteX, DeleteX) = Prob2Score (gap_extend_prob / 2.);
  transition (InsertY, DeleteX) = transition (DeleteX, InsertY) = Prob2Score (gap_extend_prob / 2.);
  // right pad
  end[RightPadX] = end[RightPadY] = end[MatchXY] = 0;
  transition (RightPadX, RightPadY) = 0;
  transition (RightPadX, RightPadX) = transition (RightPadY, RightPadY) = 0;
  transition (MatchXY, RightPadX) = transition (MatchXY, RightPadY) = 0;
  // emissions
  single_emit[LeftPadX] = single_emit[LeftPadY] = single_emit[InsertY] = single_emit[DeleteX] = single_emit[RightPadX] = single_emit[RightPadY] = vector<Score> (alph.size(), (Score) 0);
  pair_emit[MatchXY] = submat;
}

Simple_submat::Simple_submat (int alph_sz, Prob match_prob, const vector<Prob>& null_prob)
  : array2d<Score> (alph_sz, alph_sz)
{
  const Prob mismatch_prob = 1. - match_prob;
  for (int i = 0; i < alph_sz; ++i)
    for (int j = 0; j < alph_sz; ++j)
      (*this) (i, j) = Prob2Score (i==j ? (match_prob / null_prob[j] + mismatch_prob) : mismatch_prob);
}
