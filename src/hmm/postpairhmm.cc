#include "hmm/postpairhmm.h"
#include "hmm/pairenv.h"

Post_pair_HMM::Post_pair_HMM (const Pair_forward_backward_DP_matrix& fb)
  : array2d<Prob> (fb.xsize - 1, fb.ysize - 1)
{
  // make a list of match states in the HMM
  vector<int> match_states;
  for (int s = 0; s < fb.hmm.states(); ++s)
    if (fb.hmm.state_type[s] == Pair_HMM_scores::EmitXY)
      match_states.push_back (s);
  // fill matrix
  for (int x = 0; x < xsize(); ++x)
    for (int y = 0; y < ysize(); ++y)
      {
	Prob p = 0.;
	for_const_contents (vector<int>, match_states, s)
	  p += Score2Prob (fb.post_state_score (*s, x + 1, y + 1));
	if (p > 1)
	  {
	    if (p > 2)
	      CLOGERR << "Warning: posterior probability at (" << x << "," << y << ") is " << p << "\n";
	    p = 1;
	  }
	(*this)(x,y) = p;
      }
}

Submat_ProbArray2d::Submat_ProbArray2d (const Digitized_biosequence& xdsq, const Digitized_biosequence& ydsq, const array2d<Prob>& joint_submat)
  : array2d<Prob> (xdsq.size(), ydsq.size())
{
  const int alph_sz = joint_submat.xsize();
  // find null emission probabilities
  vector<Prob> null_emit (alph_sz, 0.);
  for (int i = 0; i < alph_sz; ++i)
    for (int j = 0; j < alph_sz; ++j)
      null_emit[i] += joint_submat(i,j);
  // fill the matrix
  for (int x = 0; x < xsize(); ++x)
    for (int y = 0; y < ysize(); ++y)
      {
	const Prob aligned_prob = joint_submat (xdsq[x], ydsq[y]);
	const Prob null_prob = null_emit[xdsq[x]] * null_emit[ydsq[y]];
	(*this)(x,y) = Prob2Nats (aligned_prob / null_prob);
      }
}

Smoothed_ProbArray2d::Smoothed_ProbArray2d (const array2d<double>& arr, int neighbourhood_size)
  : array2d<Prob> (arr.xsize(), arr.ysize())
{
  const int xs = xsize();
  const int ys = ysize();
  for (int x = 0; x < xs; ++x)
    for (int y = 0; y < ys; ++y)
      {
	const int xmin = max (x - neighbourhood_size, 0);
	const int xmax = min (x + neighbourhood_size, xs - 1);
	const int ymin = max (y - neighbourhood_size, 0);
	const int ymax = min (y + neighbourhood_size, ys - 1);
	int n = 0;
	double tot = 0;
	for (int i = xmin; i <= xmax; ++i)
	  for (int j = ymin; j <= ymax; ++j)
	    {
	      tot += arr(i,j);
	      ++n;
	    }
	(*this)(x,y) = tot / (double) n;
      }
}

Affine_ProbArray2d::Affine_ProbArray2d (int xsize, int ysize) : array2d<Prob> (xsize, ysize, 0.) { }
void Affine_ProbArray2d::fill_rescaled (const array2d<Prob>& arr, double a, double b)
{
  if (xsize() != arr.xsize() || ysize() != arr.ysize())
    resize (arr.xsize(), arr.ysize(), 0.);
  for (int x = 0; x < xsize(); ++x)
    for (int y = 0; y < ysize(); ++y)
      (*this)(x,y) = a * arr(x,y) + b;
}
Affine_ProbArray2d::Affine_ProbArray2d (const array2d<Prob>& arr, double a, double b)
  : array2d<Prob> (arr.xsize(), arr.ysize())
{
  fill_rescaled (arr, a, b);
}

Normalised_ProbArray2d::Normalised_ProbArray2d (const array2d<Prob>& arr)
  : Affine_ProbArray2d (arr.xsize(), arr.ysize())
{
  const Prob max_prob = *max_element (arr.begin(), arr.end());
  const Prob min_prob = *min_element (arr.begin(), arr.end());
  double a, b;
  if (max_prob > min_prob)
    {
      a = 1. / (max_prob - min_prob);
      b = -min_prob * a;
    }
  else
    {
      a = 1;
      b = 0;
    }
  fill_rescaled (arr, a, b);
}

Optimal_accuracy_DP_matrix::Optimal_accuracy_DP_matrix (const array2d<double>& reward)
  : reward (reward),
    cell (reward.xsize() + 1, reward.ysize() + 1)
{
  fill();
}

void Optimal_accuracy_DP_matrix::fill()
{
  // print reward matrix
  if (CTAGGING(3,OPTACC))
    {
      const double max_reward = *max_element (reward.begin(), reward.end());
      Pair_envelope pair_env;
      pair_env.initialise_from_posterior_matrix (reward, max_reward/8.);  // use min prob of 1/8 for max ANSI color contrast
      CTAG(3,OPTACC) << "Reward matrix for optimal accuracy alignment:\n";  // refresh CTAG() flags
      const Biosequence xseq (reward.xsize(), '*');
      const Biosequence yseq (reward.ysize(), '*');
      pair_env.render_dotplot (CL, xseq, yseq, 8);
    }

  // fill optacc matrix
  for (int i = 0; i < xsize(); ++i)
    for (int j = 0; j < ysize(); ++j)
      best_dir (i, j, cell(i,j));
}

Pairwise_path Optimal_accuracy_DP_matrix::traceback() const
{
  int i = xsize() - 1;
  int j = ysize() - 1;
  double temp;
  vector<EmitDir> emit_path;
  while (i > 0 || j > 0)
    {
      const EmitDir dir = best_dir (i, j, temp);
      emit_path.push_back (dir);
      if (dir & EmitX) --i;
      if (dir & EmitY) --j;
    }
  reverse (emit_path.begin(), emit_path.end());
  Pairwise_path path;
  for_const_contents (vector<EmitDir>, emit_path, e)
    path.append_column (*e & EmitX, *e & EmitY);
  if (path.count_steps_in_row(0) != reward.xsize() || path.count_steps_in_row(1) != reward.ysize())
    THROWEXPR ("Posterior decoding traceback path doesn't match dimensions of reward matrix");
  return path;
}
