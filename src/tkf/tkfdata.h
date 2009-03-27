#ifndef TKFDATA_INCLUDED
#define TKFDATA_INCLUDED

#include <fstream>
#include <deque>
#include "handel/handelbase.h"
#include "tkf/tkfparams.h"
#include "tkf/tkfhmm.h"

// main class for TKF model tree, alignment & scores
struct TKF_align : Handel_base
{
  // typedefs
  typedef map<Phylogeny::Undirected_pair,TKF_branch_scores*> Branch_scores_map;
  
  // parameters, scores by sequence, scores by branch
  const TKF_params&  params;
  TKF_seq_scores     seq_scores;
  Branch_scores_map  branch_scores;

  // constructors, destructor, assignment operator, clone method
  TKF_align (const TKF_params& params);
  TKF_align (const TKF_align& tkf);
  ~TKF_align();

  TKF_align& operator= (const TKF_align& tkf);

  TKF_align* clone();

  // general housekeeping methods
  void update_seq_scores();
  void clear_branch_scores();
  void update_branch_scores();

  // virtual DP/alignment methods
  bool propose_optimise_node (Node node);
  bool propose_optimise_branch (const Undirected_pair& branch);
  void propose_sample_node (Node node, double kT);
  void propose_sample_branch (const Undirected_pair& branch, double kT);
  void propose_sample_sequence (Node node, double kT);

  void align_and_infer_parent (const Score_profile& xseq,
			       const Score_profile& yseq,
			       Node ancestor_node,
			       Alignment_path& axy_path);

  // sample_progressive_alignment:
  // Discard existing alignment, then sample a child-pair DP path for each missing node
  void sample_progressive_alignment (double kT, bool sample_seq = 0);

  // tree-sampling methods
  void propose_sample_branch_slide (Node grumpa, Node dad, Node son, double kT = 1., int sample_points = 4);
  bool propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT = 1);
  void propose_sample_branch_length (const Undirected_pair& branch, double kT = 1., double tmax = 10., int sample_points = 4);
  void propose_optimise_branch_length (const Undirected_pair& branch, double tmax = 10., double tres = .01);

  // scoring methods
  Score conditioned_branch_path_score (const Node_pair& branch) const;
  Prob gamma() const;
  Substitution_matrix_factory& submat_factory() const;
};

// TKF_functions functors for branch length likelihood
struct TKF_branch_length_funcs
{
  // data
  Score_profile prof1, prof2;
  Pairwise_path path;
  TKF_aligned_counts_function *counts;
  TKF_functions *funcs;

  // constructor
  TKF_branch_length_funcs (TKF_align& tkf, const Phylogeny::Undirected_pair& branch, double tres = .01);

  // destructor
  ~TKF_branch_length_funcs();
};

#endif
