#include <fstream>
#include <deque>
#include "handel/talike.h"
#include "handel/treeshuffler.h"
#include "tree/substitution_matrix_factory.h"

// Stockholm tag for alignment method
#define Stockholm_alignment_type_tag     "TYPE"
#define Stockholm_alignment_type_sampled "Sampled"
#define Stockholm_alignment_type_refined "Refined"
#define Stockholm_alignment_type_final   "Final"
#define Stockholm_indel_parameter_tag	 "PARAMS"

// Stockholm tag for alignment sample number
#define Stockholm_alignment_step_tag     "STEP"

// memoizations of sequence & branch scores
struct Handel_seq_scores
{
  Score gamma, not_gamma;  // Prob2Score(gamma) & Prob2Score(1-gamma)
  vector<Score> prior;  // prior distribution over residues

  // helpers
  Score seq_len_score (int len) const { return ScorePMul (len * gamma, not_gamma); }
};

struct Handel_branch_scores
{
  array2d<Score> cond_submat;
};

// virtual base class for Handel model tree, alignment & scores
struct Handel_base : Tree_alignment
{
  // typedefs
  typedef Phylogeny::Node Node;
  typedef Phylogeny::Node_pair Node_pair;
  typedef Phylogeny::Undirected_pair Undirected_pair;

  typedef map<Undirected_pair,Handel_branch_scores> Handel_branch_scores_map;

  // score cache
  Handel_seq_scores        handel_seq_scores;
  Handel_branch_scores_map handel_branch_scores;

  // MCMC sampling
  Tree_alignment_likelihood* target_loglike;
  ostream* sample_stream;

  // constructors
  // NB constructors do NOT call update_seq_scores, tree_changed, etc; subclasses must do that.
  Handel_base();
  Handel_base (const Handel_base& hand);

  // assignment operator
  Handel_base& operator= (const Handel_base& hand);

  // virtual clone method
  virtual Handel_base* clone() = 0;

  // general housekeeping methods
  void tree_changed();  // virtual observer method inherited from Tree_alignment; calls update_branch_scores

  /** Optimising methods
   */

  /** optimise_node (node, ...)
   // Obtains the Viterbi inference of an ancestral sequence,
   // conditioned on its alignment to adjacent nodes.
   // Returns true if any optimisation was found.
   // If sample_seq==0, Score_profile can contain Felsenstein
   // wildcards (Holmes & Bruno, Bioinformatics, 2001).
   */
  virtual bool propose_optimise_node (Node node) = 0;
  bool optimise_node (Node node, bool sample_seq = 0);
  
  /** optimise_branch (branch)
   // Obtains the Optimal Accuracy pairwise alignment along a branch of the
   // tree, conditioned on the sequences at each end of the branch.
   // (See Holmes & Durbin, JCB, 1998.)
   // Returns true if any optimisation was found.
   */
  virtual bool propose_optimise_branch (const Undirected_pair& branch) = 0;
  bool optimise_branch (const Undirected_pair& branch);

  /** optimise_sequence (node)
   // Obtains the Viterbi inference of an ancestral sequence,
   // conditioned on its length & its alignment to adjacent nodes.
   */
  void propose_optimise_sequence (Node node);
  void optimise_sequence (Node node);

  /** optimise_branch_length
   // optimises length of branch to resolution tres, up to tmax.
   */
  virtual void propose_optimise_branch_length (const Undirected_pair& branch,
					       double tmax = 10., double tres = .01) = 0;
  void optimise_branch_length (const Undirected_pair& branch,
			       double tmax = 10., double tres = .01);

  /** align_and_infer_parent
   // aligns two siblings and infers their common ancestor
   // Three-way (ancestor,x,y) alignment returned in axy_path
   */
  virtual void align_and_infer_parent (const Score_profile& xseq,
				       const Score_profile& yseq,
				       Node ancestor_node,
				       Alignment_path& axy_path) = 0;

  /** MCMC sampling methods.
   */

  /**
   // General notes:
   //
   // (1) Metropolis-Hastings sampling.
   // The virtual methods (propose_*) are called by the nonvirtual methods
   // If target_loglike!=0, it is then used as a target distribution
   // for Metropolis propose/accept/reject sampling.
   //
   // (2) Simulated annealing.
   // In all these methods, the posterior distribution is effectively
   // raise to the power (1/kT)=beta, to facilitate simulated annealing.
   // Use kT=1 for regular MCMC sampling.
  */

  /** DP methods
   */

  /** sample_node (node, ...)
  // Samples an ancestral sequence, conditioned on its alignment
  // to adjacent nodes.
  // If sample_seq==0, Score_profile can contain Felsenstein wildcards.
  */

  virtual void propose_sample_node (Node node, double kT = 1.) = 0;
  void sample_node (Node node, double kT = 1., bool sample_seq = 0);

  /** sample_branch (branch, ...)
  // Samples the pairwise alignment along an ancestral branch,
  // conditioned on the length of the sequences at each end of the branch.
  */
  virtual void propose_sample_branch (const Undirected_pair& branch, double kT = 1.) = 0;
  void sample_branch (const Undirected_pair& branch, double kT = 1.);

  /** sample_seq (node, ...)
  // Samples an ancestral sequence, conditioned on its length & alignment
  // to adjacent nodes.
  */
  void propose_sample_sequence (Node node, double kT = 1.);
  void sample_sequence (Node node, double kT = 1.);

  /** Tree-sampling methods
   */

  /** sample_branch_slide
   // slides dad between grumpa & son
   */
  virtual void propose_sample_branch_slide (Node grumpa, Node dad, Node son, double kT = 1., int sample_points = 4);
  void sample_branch_slide (Node grumpa, Node dad, Node son, double kT = 1., int sample_points = 4);

  /** sample_branch_swap
  // swaps aunt & nephew
  // returns true if succeeded in swap
  */
  virtual bool propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT = 1) = 0;
  bool sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT = 1);

  /** sample_branch_length
   // samples length of branch to resolution tres, up to tmax
   */
  virtual void propose_sample_branch_length (const Undirected_pair& branch, double kT = 1., double tmax = 10., int sample_points = 4);
  void sample_branch_length (const Undirected_pair& branch, double kT = 1., double tmax = 10., int sample_points = 4);
	
  /** indel_parameters
   // dummy functions - samples indel parameters in derived classes
   */
   virtual void sample_indel_params();
   virtual sstring indel_parameter_string() const; 

  /** proposal_accept
   // Metropolis-Hastings sampling method.
   // Calls target_loglike (if nonzero) to get new tree/alignment score.
   // DOES NOT delete old_align afterwards.
   // Returns true if accepted.
   */
  bool proposal_accept (Handel_base* old_align, double kT = 1.);

  /** Methods that do several rounds of DP
   */

  /** anneal(...)
   // Main MCMC sampling method.
   // Can also (as the name implies) be used for ad-hoc simulated annealing
   // If sample_seq==0, Felsenstein wildcards are used
   // If target_loglike is non-null, it's used to accept/reject moves
   */
  Score anneal (double kT_start, double kT_end, int annealing_steps,
		Tree_shuffler& shuffler, vector<int>& scores,
		bool sample_seq = 0, bool use_best = 0,
		bool refine = 0, int refine_period = 0);

  /** optimise_missing_nodes
  // call optimise_node() on each missing node,
  // within constraints of existing alignment
  */
  void optimise_missing_nodes (bool sample_seq = 0);

  /** viterbi_progressive_alignment
  // discard existing alignment, then do Viterbi
  // child-pair DP for each missing node
  */
  void viterbi_progressive_alignment (bool sample_seq = 0);

  /** refine_nodes_or_branches
  // call optimise_branch() or optimise_node() [depending on node_flag] until every last bit squeezed from alignment
  // return TRUE if anything changed
  */
  bool refine_nodes_or_branches (bool node_flag, bool sample_seq = false);

  // Score calculation methods
  // alignment_path_score
  // calls assert_nodes_equal_rows(), returns 0
  Score alignment_path_score() const;

  // column_emit_score -- implementation of Felsenstein's pruning algorithm
  // assumes that nodes_equal_rows() == true, but doesn't check
  Score column_emit_score (int col, const vector<int>& seq_coords) const;

  // alignment_emit_score
  // calls assert_nodes_equal_rows()
  // returns sum of column_emit_score's
  Score alignment_emit_score() const;

  // alignment_score
  // Returns alignment_path_score(indels) + alignment_emit_score(Felsenstein)
  Score alignment_score() const;

  // conditioned_alignment_path_score
  // Returns score of subtree rooted at dad, conditioned on length of dad,
  // excluding nodes in grumpa's direction.
  Score conditioned_alignment_path_score (Node dad, Node grumpa = -1) const;

  // equilibrium_node_length_score
  // Returns score for observing a given node's length at time t=0
  Score equilibrium_node_length_score (Node node) const;

  // Null scoring methods
  Score null_emit_score() const;
  Score null_length_score() const;
  Score null_score() const;  // null_emit_score + null_length_score

  /** Virtual scoring methods
   */

  /** By overriding the following methods, subclasses can provide their own scoring cache.
   // Remember to call the super method!
   */
  virtual void update_seq_scores();
  virtual void update_branch_scores();
  virtual void clear_branch_scores();

  /* conditioned_branch_path_score
  // Returns score of a pairwise branch alignment,
  // conditioned on the parent sequence's length.
  */
  virtual Score conditioned_branch_path_score (const Node_pair& branch) const = 0;

  /** gamma
   // Returns the geometric probability distribution parameter, gamma.
   //   P(sequence length at equilibrium = L) = (1-gamma) * gamma^L
   */
  virtual Prob gamma() const = 0;

  /** submat_factory
   // Returns the substitution matrix factory
   */
  virtual Substitution_matrix_factory& submat_factory() const = 0;

  // Logging methods
  void show_branch_score_breakdown (ostream& o, const Node_pair& branch);
  void show_branch_emit_breakdown_by_column (ostream& o, const Node_pair& branch);
  void show_node_score_breakdown (ostream& o, Node node);

  void show_column_emit_score_breakdown (ostream& o, int col, const vector<int>& seq_coords) const;

  // Stockholm output method
  void write_Stockholm_with_score (ostream& out, const char* alignment_type, int alignment_step) const;

  // DP helper methods

  // calculate_conditional_score_profile()
  // Calculates a score profile of Felsenstein sums for node #dad,
  // conditional on alignment to child nodes (i.e. all nodes except #grumpa).
  // if with_prior == TRUE, then the likelihoods are multiplied by the prior probabilities as well.
  // if normalise == TRUE, then the score profile will be normalised to become a posterior probability profile.
  void calculate_conditional_score_profile (Score_profile& profile, Node dad, Node grumpa = -1, bool with_prior = 0, bool normalise = 0);

  // more general version of calculate_conditional_score_profile() that takes an explicit list of neighbors to include in the sum
  // (instead of a single node, grumpa, that should be avoided)
  void calculate_conditional_score_profile (Score_profile& profile, Node dad, const vector<Node>& neighbors_to_include, bool with_prior = 0, bool normalise = 0);

  // alphabet
  const Alphabet& alphabet() const { return submat_factory().alphabet(); }
};
