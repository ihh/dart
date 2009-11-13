#include <fstream>
#include <queue>
#include <algorithm>
#include <functional>
#include "tkf/tkfdata.h"
#include "tkf/tkfnodedp.h"
#include "tkf/tkfchildpairdp.h"
#include "newmat/newmatio.h"
#include "util/rnd.h"
#include "tree/hasegawa.h"
#include "tree/rate_matrix.h"
#include "util/vector_output.h"
#include "util/logfile.h"
#include "tkf/tkfhmm.h"
#include "util/vector_output.h"
#include "util/math_fn.h"
#include "util/maximise.h"

TKF_align::TKF_align (const TKF_params& params)
  : Handel_base(),
    params (params),
    seq_scores (params)
{
  update_seq_scores();
  tree_changed();  // updates branch scores
  align_changed();  // does nothing
}

TKF_align::TKF_align (const TKF_align& tkf)
  : Handel_base (tkf),
    params (tkf.params),
    seq_scores (tkf.seq_scores)
{
  update_seq_scores();
  tree_changed();  // updates branch scores
  align_changed();  // does nothing
}

TKF_align::~TKF_align()
{
  clear_branch_scores();
}

TKF_align& TKF_align::operator= (const TKF_align& tkf)
{
  if (&params != &tkf.params)
    THROW Standard_exception ("Can't assign TKF_align when parameter sets don't match, foole");
  Handel_base::operator= (tkf);   // this will call the virtual method tree_changed(), which updates branch_scores
  seq_scores = tkf.seq_scores;
  return *this;
}

TKF_align* TKF_align::clone()
{
  return new TKF_align (*this);
}

void TKF_align::update_seq_scores()
{
  Handel_base::update_seq_scores();
  seq_scores = TKF_seq_scores (params);
}

void TKF_align::clear_branch_scores()
{
  Handel_base::clear_branch_scores();
  for_contents (Branch_scores_map, branch_scores, i)
    delete (*i).second;
  branch_scores.clear();
}

void TKF_align::update_branch_scores()
{
  Handel_base::update_branch_scores();
  if (tree.nodes())
    {
      if (CLOGGING(3)) { CL << "(TKF_align) Updating branch scores. "; TKF_branch_scores::explain_labels (CL); }
      for_iterator (Phylogeny::Branch_iter, b,
		    tree.branches_begin (tree.root, -1),
		    tree.branches_end())
	{
	  TKF_branch_scores* bs = branch_scores[Phylogeny::Undirected_pair(*b)] = new TKF_branch_scores (params, (*b).length);
	  if (CLOGGING(3)) { CL << "Branch " << tree.branch_specifier(*b) << ":\t"; bs->show (CL); }
	}
    }
}

bool TKF_align::propose_optimise_node (Phylogeny::Node node)
{
  TKF_node_DP npath (*this, node);

  if (CLOGGING(3))
    {
      Pairwise_path align_to_row_map (align.path, node2row[node], node2row[node], 0);
      align_to_row_map[0] = align.path.create_full_row();
      int p = npath.path_transition_score (align_to_row_map);
      int e = npath.path_emit_score (align_to_row_map);
      CLOG(3) << "Node HMM log-likelihood (before) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits\n";
    }

  npath.do_max_dp();
  Pairwise_path align_to_new_row_map = npath.max_traceback();
  if (align_to_new_row_map.row(1) == align.path.row(node2row[node]))
    { CLOG(4) << "New row is same as old row; \"nothing new under the sun\" (found written in hieroglyphs on a pyramid wall)\n"; return 0; }

  if (CLOGGING(3))
    {
      int p = npath.path_transition_score (align_to_new_row_map);
      int e = npath.path_emit_score (align_to_new_row_map);
      CLOG(3) << "Node HMM log-likelihood (after) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits\n";
    }

  realign_node (node, align_to_new_row_map);
  
  return 1;
}

bool TKF_align::propose_optimise_branch (const Phylogeny::Undirected_pair& branch)
{
  TKF_joint_pair_HMM_scores hmm (params, tree.branch_length (branch));
  Score_profile x_prof;
  Score_profile y_prof;
  calculate_conditional_score_profile (x_prof, branch.first, branch.second, 0, 0);
  calculate_conditional_score_profile (y_prof, branch.second, branch.first, 0, 0);

  Pairwise_path old_path (align.path, node2row[branch.first], node2row[branch.second], (bool) 1);

  if (CTAGGING(3,BRANCHLIKELIHOOD STATEPATH))
    {
      vector<int> state_path = hmm.find_optimal_state_path (old_path);
      int p = hmm.path_transition_score (state_path);
      int e = hmm.path_emit_score (state_path, x_prof, y_prof);
      if (CTAGGING(0,STATEPATH)) CL << "State path (before): [" << state_path << "]\n";
      CTAG(3,BRANCHLIKELIHOOD) << "Pair HMM log-likelihood (before) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits\n";
    }
  
  Pair_Viterbi_DP_matrix matrix (hmm, x_prof, y_prof);
  vector<int> state_path = matrix.optimal_state_path();

  if (CTAGGING(3,BRANCHLIKELIHOOD STATEPATH))
    {
      if (CTAGGING(0,STATEPATH)) CL << "State path (after): [" << state_path << "]\n";
      int p = hmm.path_transition_score (state_path);
      int e = hmm.path_emit_score (state_path, x_prof, y_prof);
      CTAG(3,BRANCHLIKELIHOOD) << "Pair HMM log-likelihood (after) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits\n";
    }
  
  Pairwise_path new_path = hmm.convert_state_path_to_alignment (state_path);
  if (new_path == old_path) { CLOG(4) << "New path is same as old path - why change a successful formula?\n"; return 0; }
  realign_pair (branch, new_path);

  return 1;
}

void TKF_align::propose_sample_node (Phylogeny::Node node, double kT)
{
  TKF_node_DP npath (*this, node, kT);

  if (CLOGGING(3))
    {
      Pairwise_path align_to_row_map (align.path, node2row[node], node2row[node], 0);
      align_to_row_map[0] = align.path.create_full_row();
      int p, e;
      ScaleScore (p = npath.path_transition_score (align_to_row_map), kT);
      ScaleScore (e = npath.path_emit_score (align_to_row_map), kT);
      CLOG(3) << "Node HMM log-likelihood (before) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits\n";
    }

  npath.do_sum_dp();
  Pairwise_path align_to_new_row_map = npath.sample_traceback();
  if (align_to_new_row_map.row(1) == align.path.row(node2row[node]))
    {
      CLOG(4) << "New row is same as old row; \"nothing new under the sun\" (ancient Egyptian proverb)\n";
      CTAG(5,MCMC SAMPLE_NODE) << "Node-sampling left alignment unchanged\n";
      return;
    }

  if (CLOGGING(3))
    {
      int p, e;
      ScaleScore (p = npath.path_transition_score (align_to_new_row_map), kT);
      ScaleScore (e = npath.path_emit_score (align_to_new_row_map), kT);
      CLOG(3) << "Node HMM log-likelihood (after) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits\n";
    }

  realign_node (node, align_to_new_row_map);
}

void TKF_align::propose_sample_branch (const Phylogeny::Undirected_pair& branch, double kT)
{
  const int over_relaxation_trials = 1;  // HACK, can't be bothered to delete over-relaxation code

  TKF_joint_pair_HMM_scores hmm (params, tree.branch_length (branch));
  hmm.scale_all_scores (1/kT);
  Score_profile x_prof;
  Score_profile y_prof;
  calculate_conditional_score_profile (x_prof, branch.first, branch.second, 1, 0);
  calculate_conditional_score_profile (y_prof, branch.second, branch.first, 0, 0);

  const int total_paths = max (over_relaxation_trials + 1, 2);
  vector<Pairwise_path> trial_paths (total_paths);                  // trial_paths[0] will be the old path
  vector<double>        mean_displacement (total_paths);

  trial_paths[0] = Pairwise_path (align.path, node2row[branch.first], node2row[branch.second], (bool) 1);
  if (over_relaxation_trials > 1)
    mean_displacement[0] = trial_paths[0].mean_match_displacement();

  const int matches = trial_paths[0].match_columns();

  if (CLOGGING(3))
    {
      vector<int> state_path = hmm.find_optimal_state_path (trial_paths[0]);
      int p, e;
      ScaleScore (p = hmm.path_transition_score (state_path), kT);
      ScaleScore (e = hmm.path_emit_score (state_path, x_prof, y_prof), kT);
      CL << "Pair HMM log-likelihood (before) is ";
      CL << Score2Bits (p) << " (indels) + ";
      CL << Score2Bits (e) << " (substitutions)  = ";
      CL << Score2Bits (p+e) << " bits";
      CL << " (mean displacement = " << mean_displacement[0] << ")\n";
    }
  
  Pair_forward_DP_matrix matrix (hmm, x_prof, y_prof);

  int trial = 0;
  do
    {
      vector<int> state_path = matrix.sample_state_path();
      trial_paths [trial+1] = hmm.convert_state_path_to_alignment (state_path);
      if (over_relaxation_trials > 1)
	mean_displacement [trial+1] = trial_paths [trial+1] . mean_match_displacement();
      if (CLOGGING(3))
	{
	  if (over_relaxation_trials > 1)
	    CL << "Trial path #" << trial+1 << ": Pair HMM log-likelihood is ";
	  else
	    CL << "Pair HMM log-likelihood (after) is ";
	  int p, e;
	  ScaleScore (p = hmm.path_transition_score (state_path), kT);
	  ScaleScore (e = hmm.path_emit_score (state_path, x_prof, y_prof), kT);
	  CL << Score2Bits (p) << " (indels) + ";
	  CL << Score2Bits (e) << " (substitutions)  = ";
	  CL << Score2Bits (p+e) << " bits";
	  if (over_relaxation_trials > 1)
	    {
	      CL << " (mean displacement = " << mean_displacement[trial+1] << "; ";
	      CL << "path has " << (100 * trial_paths[trial+1].match_overlap (trial_paths[0]) / matches) << "% of original's match columns)";
	    }
	  CL << "\n";
	}
    }
  while (++trial < over_relaxation_trials);
  
  const Pairwise_path* new_path;
  if (over_relaxation_trials <= 1)
    new_path = &trial_paths[1];
  else
    {
      vector<int> ordered_trial_index (total_paths);
      for (int i = 0; i < total_paths; ++i) ordered_trial_index[i] = i;

      Schwartzian<double> by_mean_displacement (mean_displacement);
      sort (ordered_trial_index.begin(),
	    ordered_trial_index.end(),
	    by_mean_displacement);

      const int old_path_rank = find (ordered_trial_index.begin(), ordered_trial_index.end(), 0) - ordered_trial_index.begin();
      const int new_path_rank = total_paths - 1 - old_path_rank;
      const int new_path_index = ordered_trial_index [new_path_rank];
      new_path = &trial_paths [new_path_index];
      if (CLOGGING(3))
	{
	  CL << "Old path has displacement rank " << old_path_rank+1 << " out of " << total_paths;
	  if (new_path_index == 0) CL << ", i.e. midway; retaining (should probably choose an odd number of trials)\n";
	  else CL << "; picking path at rank " << new_path_rank+1 << " (trial #" << new_path_index+1 << ")\n";
	}
    }
  
  if (*new_path == trial_paths[0])
    {
      CLOG(4) << "New path is same as old path - if it's not broke don't fix it\n";
      CTAG(5,MCMC SAMPLE_NODE) << "Branch-sampling left alignment unchanged\n";
      return;
    }
  realign_pair (branch, *new_path);
}


void TKF_align::align_and_infer_parent (const Score_profile& xseq,
					const Score_profile& yseq,
					Node node,
					Alignment_path& axy_path)
{
  TKF_child_pair_DP cpdp (*this, node, xseq, yseq, 1);
  cpdp.do_max_dp();
  Score_profile tmp;
  axy_path = cpdp.max_traceback (tmp, true);
}

void TKF_align::sample_progressive_alignment (double kT, bool sample_seq)
{
  tree.assert_tree_is_binary();

  // for each internal node, (1) count the number of missing kids; (2) add to queue if all kids are present
  //
  Phylogeny::Node_vector   first_son (tree.nodes(), -1);
  vector<int>              missing_kids (tree.nodes(), 0);
  Phylogeny::Branch_vector branches;
  for_iterator (Phylogeny::Node_const_iter, i, tree.internals_begin(), tree.internals_end())
    {
      int& missing = missing_kids[*i];
      for_iterator (Phylogeny::Child_iter, c, tree.children_begin(*i,tree.parent[*i]), tree.children_end(*i,tree.parent[*i]))
	{
	  if (first_son[*i] == -1) first_son[*i] = *c;
	  if (!tree.is_leaf(*c)) ++missing;
	}
      if (missing == 0) branches.push_back (Phylogeny::Branch (*i, first_son[*i], tree));
    }

  // do the progressive alignment
  //
  Alignment_path::Decomposition decomp;
  vector<Score_profile>         felsenstein (tree.nodes(), Score_profile());

  while (branches.size())
    {
      // choose a random branch
      //
      int b = Rnd::rnd_int (branches.size());
      swap (branches[b], branches.back());

      Phylogeny::Node        dad    = branches.back().first;
      Phylogeny::Node        grumpa = tree.parent[dad];
      Phylogeny::Node_vector kids   = tree.children (dad, grumpa);

      vector<int> row_set (3);
      row_set[0] = node2row[dad];       // parent
      row_set[1] = node2row[kids[0]];   // child x
      row_set[2] = node2row[kids[1]];   // child y
      
      Alignment axy (align, row_set, 1);
      
      CLOG(5) << "Aligning nodes " << axy.row_name[1] << " and " << axy.row_name[2] << "\n";
      
      const Score_profile& x_prof = tree.is_leaf(kids[0]) ? *align.prof[node2row[kids[0]]] : felsenstein[kids[0]];
      const Score_profile& y_prof = tree.is_leaf(kids[1]) ? *align.prof[node2row[kids[1]]] : felsenstein[kids[1]];
      TKF_child_pair_DP cpdp (*this, dad, x_prof, y_prof, kT);
      cpdp.do_sum_dp();
      axy.path = cpdp.sample_traceback (felsenstein[dad], !sample_seq);
      if (sample_seq)
	set_node_profile (dad, new Score_profile (felsenstein[dad]));
      
      if (CTAGGING(4,ALIGN_SIBLINGS))
	{
	  axy.prof[0] = &felsenstein[dad];
	  axy.prof[1] = &x_prof;
	  axy.prof[2] = &y_prof;
	  axy.write_MUL (CL, alphabet());
	}
      
      decomp [Alignment_path::Row_pair (row_set[0], row_set[1])] = Pairwise_path (axy.path, 0, 1, 1);
      decomp [Alignment_path::Row_pair (row_set[0], row_set[2])] = Pairwise_path (axy.path, 0, 2, 1);
      
      // erase child branch from the queue; add branch to grumpa iff we've got all the uncles
      //
      branches.pop_back();
      if (grumpa != -1)
	if (--missing_kids[grumpa] == 0)
	  branches.push_back (Phylogeny::Branch (grumpa, first_son[grumpa], tree));
    }
  align.path.compose_and_log (decomp, 1);

  if (CLOGGING(7))
    {
      int path_score = alignment_path_score();
      int emit_score = alignment_emit_score();
      CL << "Alignment found with log-likelihood " << Score2Bits (path_score + emit_score) << " bits\n";

      if (CLOGGING(5))
	CL << "Log path likelihood = " << Score2Bits (path_score) << " bits, log emission likelihood = " << Score2Bits (emit_score) << " bits\n";
    }
}

void TKF_align::propose_sample_branch_slide (Node grumpa, Node dad, Node son, double kT, int sample_points)
{
  // TODO: accelerate this by constructing a mini TKF_align for the immediate neighborhood of dad
  const double grumpa_dad_len = tree.branch_length (grumpa, dad);
  const double dad_son_len = tree.branch_length (dad, son);
  const double grumpa_son_len = grumpa_dad_len + dad_son_len;
  vector<double> x (sample_points + 1), p (sample_points + 1);
  x[0] = grumpa_dad_len / grumpa_son_len;
  for (int i = 1; i <= sample_points; ++i)
    x[i] = Rnd::prob();
  for (int i = 0; i <= sample_points; ++i)
    {
      tree.branch_length (grumpa, dad) = grumpa_son_len * x[i];
      tree.branch_length (dad, son) = grumpa_son_len * (1. - x[i]);
      tree_changed();
      p[i] = Score2Prob ((Score) (kT * (double) alignment_score()));
    }

  const int i = Rnd::choose (p);

  tree.branch_length (grumpa, dad) = grumpa_son_len * x[i];
  tree.branch_length (dad, son) = grumpa_son_len * (1. - x[i]);
}

bool TKF_align::propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT)
{
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

  TKF_align* old_tkfalign = new TKF_align (*this);

  const Score old_score = alignment_score();

  const double aunt_grumpa_len = tree.branch_length (aunt, grumpa);
  const double dad_nephew_len = tree.branch_length (dad, nephew);
  tree.remove_branch (aunt, grumpa);
  tree.remove_branch (dad, nephew);
  tree.add_branch (aunt, dad, aunt_grumpa_len);
  tree.add_branch (nephew, grumpa, dad_nephew_len);
  tree.rebuild_parents();
  tree_changed();

  const Score new_score = alignment_score();

  const bool accept = Rnd::decide (1. / (1. + Score2Prob (old_score - new_score)));  // new/(new+old) = 1/(1+old/new)
  if (!accept)
    *this = *old_tkfalign;

  delete old_tkfalign;
  return accept;
}

void TKF_align::propose_sample_branch_length (const Undirected_pair& branch, double kT, double tmax, int sample_points)
{
  const double tres = .000000001;  // choose arbitrarily small tres (HACK)
  TKF_branch_length_funcs blf (*this, branch, tres);

  vector<double> t (sample_points + 1), p (sample_points + 1);
  t[0] = tree.branch_length (branch);
  for (int i = 1; i <= sample_points; ++i)
    t[i] = Rnd::prob() * tmax;
  for (int i = 0; i <= sample_points; ++i)
    p[i] = Nats2Prob (kT * blf.funcs->log_like (t[i]));

  const int i = Rnd::choose (p);
  // TODO: check target_loglike
  tree.branch_length (branch) = t[i];
}

void TKF_align::propose_optimise_branch_length (const Undirected_pair& branch, double tmax, double tres)
{
  TKF_branch_length_funcs blf (*this, branch, tres);

  double t1, t2, t3, tbest, fbest;
  bracket_maximum (blf.funcs->log_like, t1, t2, t3, 0., tmax);
  brent_deriv (blf.funcs->log_like, blf.funcs->log_like_dt, t1, t2, t3, tres, tbest, fbest);

  tree.branch_length (branch) = tbest;
}

Score TKF_align::conditioned_branch_path_score (const Phylogeny::Node_pair& branch) const
{
  Pairwise_path path = subpath (branch, 0);
  return ((Branch_scores_map&)branch_scores) [Phylogeny::Undirected_pair(branch)]->conditional_path_score (path);
}

Prob TKF_align::gamma() const
{
  return params.lambda / params.mu;
}

Substitution_matrix_factory& TKF_align::submat_factory() const
{
  return params.submat_factory;
}

TKF_branch_length_funcs::TKF_branch_length_funcs (TKF_align& tkf, const Phylogeny::Undirected_pair& branch, double tres)
{
  tkf.calculate_conditional_score_profile (prof1, branch.first, branch.second, false, false);
  tkf.calculate_conditional_score_profile (prof2, branch.second, branch.first, true, false);
  path = tkf.subpath (branch);
  counts = new TKF_aligned_counts_function (prof1, prof2, path, &tkf.alphabet());
  funcs = new TKF_functions (tkf.params, *counts, true, true, tres);
}

TKF_branch_length_funcs::~TKF_branch_length_funcs()
{
  delete funcs;
  delete counts;
}
