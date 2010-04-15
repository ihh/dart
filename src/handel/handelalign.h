#include <fstream>
#include <deque>
#include "handel/talike.h"
#include "handel/handelmodel.h"
#include "handel/handelbase.h"

// main class for Handel model tree, alignment & scores
struct Handel_alignment : Handel_base
{
  // Model
  Handel_model* model;

  // constructors
  Handel_alignment();
  Handel_alignment (const Handel_alignment& hand);

  // virtual assignment operator
  Handel_alignment& operator= (const Handel_alignment& hand);

  // virtual clone method
  Handel_alignment* clone();

  /** Optimising methods
   */
  
  /** optimise_branch (branch)
   // Obtains the Optimal Accuracy pairwise alignment along a branch of the
   // tree, conditioned on the sequences at each end of the branch.
   // (See Holmes & Durbin, JCB, 1998.)
   // Returns true if any optimisation was found.
   */
  bool propose_optimise_branch (const Undirected_pair& branch);

  /** propose_optimise_node (node)
   // unimplemented - should trim unnecessary indels
   */
  bool propose_optimise_node (Node node);

  /** optimise_branch_length
   // optimises length of branch to resolution tres, up to tmax.
   */
  void propose_optimise_branch_length (const Undirected_pair& branch,
				       double tmax = DART_MAX_BRANCH_LENGTH, double tres = .01);

  /** MCMC sampling methods.
   */

  /**
  // If target_loglike!=0, it is used as a target distribution
  // for Metropolis propose/accept/reject sampling.
  // In all these methods, the posterior distribution is effectively
  // raise to the power (1/kT)=beta, to facilitate simulated annealing.
  // Use kT=1 for regular MCMC sampling.
  */

  /** DP methods
   */

  /** sample_node (node, ...)
  // unimplemented
  */
  void propose_sample_node (Node node, double kT = 1.);

  /** sample_branch (branch, ...)
  // Samples the pairwise alignment along an ancestral branch,
  // conditioned on the length of the sequences at each end of the branch.
  */
  void propose_sample_branch (const Undirected_pair& branch, double kT = 1.);

  /** sample_seq (node, ...)
  // Samples an ancestral sequence, conditioned on its length & alignment
  // to adjacent nodes.
  */
  void propose_sample_sequence (Node node, double kT = 1.);

  /** Tree-sampling methods
   */

  /** sample_branch_swap
  // swaps aunt & nephew
  */
  bool propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT = 1);

  /** sample_branch_length
  // samples length of branch to resolution tres, up to tmax
  */
  void propose_sample_branch_length (const Undirected_pair& branch, double kT = 1., double tmax = DART_MAX_BRANCH_LENGTH, int sample_points = 4);

  /** align_and_infer_parent
   // aligns two siblings and infers their common ancestor
   // Three-way (ancestor,x,y) alignment returned in axy_path
   */
  void align_and_infer_parent (const Score_profile& xseq,
			       const Score_profile& yseq,
			       Node ancestor_node,
			       Alignment_path& axy_path);

  // conditioned_branch_path_score
  // Returns score of a pairwise branch alignment,
  // conditioned on the parent sequence's length.
  Score conditioned_branch_path_score (const Node_pair& branch) const;

  /** gamma
   // Returns the geometric probability distribution parameter, gamma.
   //   P(sequence length at equilibrium = L) = (1-gamma) * gamma^L
   */
  Prob gamma() const;

  /** submat_factory
   // Returns the substitution matrix factory
   */
  Substitution_matrix_factory& submat_factory() const;
};

