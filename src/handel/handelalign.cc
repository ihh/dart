#include <fstream>
#include <queue>
#include <algorithm>
#include <functional>
#include "handel/handelalign.h"
#include "newmat/newmatio.h"
#include "util/rnd.h"
#include "tree/hasegawa.h"
#include "tree/rate_matrix.h"
#include "util/vector_output.h"
#include "util/logfile.h"
#include "util/vector_output.h"
#include "util/math_fn.h"

Handel_alignment* Handel_alignment::clone()
{
  return new Handel_alignment (*this);
}

Handel_alignment& Handel_alignment::operator= (const Handel_alignment& hand)
{
  model = hand.model;
  Handel_base::operator= (hand); // this will call the virtual method tree_changed(), which then updates handel_branch_scores
  return *this;
}

Handel_alignment::Handel_alignment()
  : Handel_base(),
    model (0)
{ }

Handel_alignment::Handel_alignment (const Handel_alignment& hand)
  : Handel_base (hand),
    model (hand.model)
{
    tree_changed();  // updates branch scores
    align_changed();  // does nothing
}

bool Handel_alignment::propose_optimise_node (Phylogeny::Node node)
{
  CLOGERR << "Warning: general node-optimisation unimplemented (we should trim unnecessary ancestral indels here)\n";
  return false;
}

bool Handel_alignment::propose_optimise_branch (const Phylogeny::Undirected_pair& branch)
{
  Pair_HMM_scores hmm = model->joint_pair_HMM (tree.branch_length (branch));
  Score_profile x_prof;
  Score_profile y_prof;
  calculate_conditional_score_profile (x_prof, branch.first, branch.second, 0, 0);
  calculate_conditional_score_profile (y_prof, branch.second, branch.first, 0, 0);

  Pairwise_path old_path (align.path, node2row[branch.first], node2row[branch.second], (bool) 1);
  
  // Find the Viterbi path through the branch pair HMM
  // TODO: use optimal accuracy instead of Viterbi here
  Pair_Viterbi_DP_matrix matrix (hmm, x_prof, y_prof);
  vector<int> state_path = matrix.optimal_state_path();
  
  Pairwise_path new_path = hmm.convert_state_path_to_alignment (state_path);

  if (new_path == old_path) { CLOG(4) << "New path is same as old path - why change a successful formula?\n"; return 0; }

  realign_pair (branch, new_path);
  return 1;
}

void Handel_alignment::propose_sample_node (Phylogeny::Node node, double kT)
{
  if (tree.is_leaf(node))
    THROWEXPR ("Attempted to optimise sequence at leaf node '" << tree.node_specifier(node) << "'");
  CLOGERR << "Warning: general node-sampling unimplemented (should add/remove residues here; possibly Single HMM traceback)\n";
  return;
}

void Handel_alignment::propose_sample_branch (const Phylogeny::Undirected_pair& branch, double kT)
{
  Pair_HMM_scores hmm = model->joint_pair_HMM (tree.branch_length (branch));
  hmm.scale_all_scores (1/kT);
  Score_profile x_prof;
  Score_profile y_prof;
  calculate_conditional_score_profile (x_prof, branch.first, branch.second, 1, 0);
  calculate_conditional_score_profile (y_prof, branch.second, branch.first, 0, 0);

  Pairwise_path old_path (align.path, node2row[branch.first], node2row[branch.second], (bool) 1);

  // Sample a path through the branch pair HMM using the Forward matrix
  Pair_forward_DP_matrix fwd (hmm, x_prof, y_prof);
  vector<int> state_path = fwd.sample_state_path();
  
  Pairwise_path new_path = hmm.convert_state_path_to_alignment (state_path);

  if (new_path == old_path) { CLOG(4) << "New path is same as old path - why change a successful formula?\n"; return; }

  realign_pair (branch, new_path);
  return;
}


bool Handel_alignment::propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT)
{
  CTAG(5,MCMC BRANCH_SWAP) << "Swapping aunt node '" << tree.node_specifier(aunt) << "' (parent '" << tree.node_specifier(grumpa) << "') with nephew node '" << tree.node_specifier(nephew) << "' (parent '" << tree.node_specifier(dad) << "')\n";

  // Old tree:          New tree:
  //
  //      grumpa           grumpa
  //      /  |             /  |
  //    dad  aunt        dad  nephew
  //    / |              / |
  // ...  nephew      ...  aunt

  // Bail if grumpa-dad branch has gaps
  // (TODO: simultaneously resample grandparent-nephew and parent-aunt alignments)
  const Pairwise_path grumpa_dad_path = subpath (Phylogeny::Node_pair (grumpa, dad));
  if (!grumpa_dad_path.is_ungapped())
    {
      CLOG(4) << "Intermediate branch contains gaps; topology unchanged\n";
      return false;
    }

  Handel_alignment* old_handalign = new Handel_alignment (*this);

  const Score old_score = alignment_score();

  const double aunt_grumpa_len = tree.branch_length (aunt, grumpa);
  const double dad_nephew_len = tree.branch_length (dad, nephew);
  tree.remove_branch (aunt, grumpa);
  tree.remove_branch (dad, nephew);
  tree.add_branch (aunt, dad, aunt_grumpa_len);
  tree.add_branch (nephew, grumpa, dad_nephew_len);
  tree.setup_parents_vector();
  tree_changed();

  const Score new_score = alignment_score();

  const bool accept = Rnd::decide (1. / (1. + Score2Prob (old_score - new_score)));  // new/(new+old) = 1/(1+old/new)
  if (!accept)
    *this = *old_handalign;

  delete old_handalign;
  return accept;
}

void Handel_alignment::propose_sample_branch_length (const Undirected_pair& branch, double kT, double tmax, int sample_points)
{
  CTAG(5,MCMC BRANCH_LENGTH) << "Sampling length of branch from node '" << tree.node_specifier(branch.first) << "' to node '" << tree.node_specifier(branch.second) << "'\n";

  CLOGERR << "Warning: sample_branch_length unimplemented\n";
}

void Handel_alignment::propose_optimise_branch_length (const Undirected_pair& branch, double tmax, double tres)
{
  CTAG(5,MCMC BRANCH_LENGTH) << "Optimising length of branch from node '" << tree.node_specifier(branch.first) << "' to node '" << tree.node_specifier(branch.second) << "'\n";

  THROWEXPR ("optimise_branch_length unimplemented");
}

Score Handel_alignment::conditioned_branch_path_score (const Node_pair& branch) const {
  Wildcard_transducer trans = model->transducer (tree.branch_length (branch));
  Pairwise_path branch_path (align.path, node2row[branch.first], node2row[branch.second], (bool) 1);
  Symbol_score_map dummy_score_map = Dummy_alphabet.flat_score_map (0);
  Score_profile xseq (branch_path.count_steps_in_row(0), dummy_score_map);
  Score_profile yseq (branch_path.count_steps_in_row(1), dummy_score_map);
  return trans.pairwise_path_score (xseq, yseq, branch_path);
}

Prob Handel_alignment::gamma() const { if (model==0) THROWEXPR("no model"); return model->gamma(); }
Substitution_matrix_factory& Handel_alignment::submat_factory() const { if (model==0) THROWEXPR("no model"); return *model; }

void Handel_alignment::align_and_infer_parent (const Score_profile& xseq,
					       const Score_profile& yseq,
					       Node ancestor_node,
					       Alignment_path& axy_path)
{
  const Node xchild = tree.children (ancestor_node, tree.parent[ancestor_node]) [0];
  const Node ychild = tree.children (ancestor_node, tree.parent[ancestor_node]) [1];

  const double xdist = tree.branch_length (ancestor_node, xchild);
  const double ydist = tree.branch_length (ancestor_node, ychild);

  Pair_HMM_scores hmm = model->joint_pair_HMM (xdist + ydist);

  // TODO: use optimal accuracy instead of Viterbi here
  Pair_Viterbi_DP_matrix matrix (hmm, xseq, yseq);
  vector<int> state_path = matrix.optimal_state_path();
  Pairwise_path xy_path = hmm.convert_state_path_to_alignment (state_path);

  // give the ancestor the same length as the closest child
  axy_path = xy_path;
  if (xdist < ydist)
    axy_path.insert_rows (0, xy_path.parent());
  else
    axy_path.insert_rows (0, xy_path.child());
}
