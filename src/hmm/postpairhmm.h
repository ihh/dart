#ifndef POST_PAIRHMM_INCLUDED
#define POST_PAIRHMM_INCLUDED

#include "hmm/pairhmm.h"

/** Post_pair_HMM is an array of posterior probabilities for a Pair HMM.
    Entry (i,j) is probability that residues (i,j) of sequences (X,Y) are aligned via a match state of the Pair HMM.
*/
struct Post_pair_HMM : array2d<Prob>
{
  Post_pair_HMM (const Pair_forward_backward_DP_matrix& fb);
};

/** Submat_ProbArray2d is an array of log-odds-ratios of residues being aligned
    for a pair of sequences (X,Y)
    NB really the entries are Loge's not Prob's.... this is unsatisfying as they're supposed to be distinct types
    Should really use a base class of "array2d<double>" throughout this file, not "array2d<Prob>"
 */
struct Submat_ProbArray2d : array2d<Prob>
{
  Submat_ProbArray2d (const Digitized_biosequence& xdsq, const Digitized_biosequence& ydsq, const array2d<Prob>& joint_submat);
};

/** Smoothed_ProbArray2d is a smoothed (i.e. neighbourhood-averaged) array of probabilities.
    A neighbourhood size of 0 leaves the source array unchanged.
 */
struct Smoothed_ProbArray2d : array2d<Prob>
{
  Smoothed_ProbArray2d (const array2d<double>& arr, int neighbourhood_size);
};

/** Affine_ProbArray2d is an affine-rescaled array of probabilities
 */
class Affine_ProbArray2d : public array2d<Prob>
{
  // constructor
protected:
  Affine_ProbArray2d (int xsize, int ysize);  // call fill_rescaled manually
  void fill_rescaled (const array2d<Prob>& arr, double a, double b);
public:
  Affine_ProbArray2d (const array2d<Prob>& arr, double a, double b);  // rescales to via affine transformation x <-- ax+b
};

/** Normalised_ProbArray2d is an affine-renormalised array of probabilities
 */
struct Normalised_ProbArray2d : Affine_ProbArray2d
{
  Normalised_ProbArray2d (const array2d<Prob>& arr);  // rescales to [0,1] via affine transformation
};

/** Optimal_accuracy_DP_matrix applies the "optimal accuracy" posterior decoding algorithm
    to a matrix of posterior probabilities
 */
struct Optimal_accuracy_DP_matrix
{
  const array2d<double>& reward;  // user-supplied matrix of non-negative rewards
  array2d<double> cell;   // cell(i,j) = max path for subseq 1..i of X and 1..j of Y
  // helpers
  inline int xsize() const { return cell.xsize(); }   // == reward.xsize() + 1
  inline int ysize() const { return cell.ysize(); }   // == reward.ysize() + 1
  inline double final_score() const { return cell(xsize()-1,ysize()-1); }
  // method for finding best way into a cell
  enum EmitDir { EmitX = 1, EmitY = 2, EmitXY = 3 };
  inline EmitDir best_dir (int i, int j, double& best) const;
  // constructor
  Optimal_accuracy_DP_matrix (const array2d<double>& reward);
  // fill method
  void fill();
  // traceback method
  Pairwise_path traceback() const;
};

inline Optimal_accuracy_DP_matrix::EmitDir Optimal_accuracy_DP_matrix::best_dir (int i, int j, double& best) const
{
  const double x_reward = i==0 ? -1. : cell(i-1,j);
  const double y_reward = j==0 ? -1. : cell(i,j-1);
  const double xy_reward = i==0 ? (j==0 ? 0. : -1.) : (j==0 ? -1. : cell(i-1,j-1) + reward(i-1,j-1));
  if (xy_reward >= x_reward)
    {
      if (xy_reward >= y_reward)
	{
	  best = xy_reward;
	  return EmitXY;
	}
    }
  else if (x_reward > y_reward)
    {
      best = x_reward;
      return EmitX;
    }
  best = y_reward;
  return EmitY;
}

#endif /* POST_PAIRHMM_INCLUDED */
