#include <queue>
#include "handelbase.h"
#include "util/vector_output.h"

Handel_base::Handel_base (const Handel_base& hand)
  : Tree_alignment (hand),
    handel_seq_scores (hand.handel_seq_scores),
    handel_branch_scores (hand.handel_branch_scores),
    target_loglike (hand.target_loglike),
    sample_stream (hand.sample_stream)
{ }

Handel_base::Handel_base()
  : target_loglike (0),
    sample_stream (0)
{ }

Handel_base& Handel_base::operator= (const Handel_base& hand)
{
  handel_seq_scores = hand.handel_seq_scores;
  handel_branch_scores = hand.handel_branch_scores;
  target_loglike = hand.target_loglike;
  sample_stream = hand.sample_stream;
  Tree_alignment::operator= (hand);   // this will call the virtual method tree_changed(), which then updates handel_branch_scores
  return *this;
}

void Handel_base::update_seq_scores()
{
  const Prob g = gamma();
  handel_seq_scores.gamma = Prob2Score (g);
  handel_seq_scores.not_gamma = Prob2Score (1. - g);
  handel_seq_scores.prior = Prob2ScoreVec (submat_factory().create_prior());
}

void Handel_base::clear_branch_scores()
{
  handel_branch_scores.clear();
}

void Handel_base::update_branch_scores()
{
  handel_branch_scores.clear();
  if (tree.nodes())
    {
      CTAG(3,HANDEL) << "(Handel_base) Updating branch scores.\n";
      Substitution_matrix_factory& sm = submat_factory();
      for_iterator (Phylogeny::Branch_iter, b,
		    tree.branches_begin (tree.root, -1),
		    tree.branches_end())
	{
	  Handel_branch_scores bs;
	  bs.cond_submat = Prob2ScoreArray2d (sm.create_conditional_substitution_matrix ((*b).length));
	  handel_branch_scores[Phylogeny::Undirected_pair(*b)] = bs;
	  if (CTAGGING(3,HANDEL))
	    CL << "Branch " << tree.branch_specifier(*b) << "\n" << bs.cond_submat;
	}
    }
}

void Handel_base::tree_changed()
{
  update_branch_scores();
  update_alignment_row_names_from_maps();
}

void Handel_base::sample_node (Node node, double kT, bool sample_seq)
{
  if (tree.is_leaf(node)) { sstring e; e << "Attempted to sample sequence at leaf node '" << tree.node_specifier(node) << "'\n"; THROW Standard_exception(e); }

  bool row_is_new = node2row[node] == -1;
  const Score_profile* old_profile = 0;
  if (!row_is_new)
    {
      old_profile = align.prof[node2row[node]];
      if (old_profile != 0 && !sample_seq)
	{ CTAG(5,HANDEL) << "Sequence at node '" << tree.node_specifier(node) << "' is specified; skipping sampling step\n"; return; }
    }

  CTAG(5,MCMC SAMPLE_NODE) << "Sampling sequence at internal node '" << tree.node_specifier(node) << "'\n";

  Handel_base* old_align = 0;
  if (target_loglike)
    old_align = clone();

  if (nodes_equal_rows())
    if (CTAGGING(4,BREAKDOWN))
      {
	CL << "Log-likelihood node breakdown (before node sampling):\n";
	show_node_score_breakdown (CL, node);
      }

  if (row_is_new) add_empty_alignment_row_for_node (node);
  if (!row_is_new && old_profile != 0 && sample_seq)
    {
      align.prof[node2row[node]] = 0;
      if (CTAGGING(4,HANDEL))
	{
	  Score_profile cond;
	  calculate_conditional_score_profile (cond, node, -1, 1, 1);
	  CL << "Temporarily replacing sequence at node '" << tree.node_specifier(node) << "' with wildcards '*'";
	  CL << " (sequence posterior log-probability was " << Score2Bits ((*old_profile).inner_product (cond)) << " bits)\n";
	}
    }

  propose_sample_node (node, kT);

  if (sample_seq)
    propose_sample_sequence (node, kT);

  if (nodes_equal_rows())
    if (CTAGGING(4,BREAKDOWN))
      {
	CL << "Log-likelihood node breakdown (after node sampling):\n";
	show_node_score_breakdown (CL, node);
      }

  proposal_accept (old_align, kT);

  if (target_loglike)
    delete old_align;
}

void Handel_base::sample_branch (const Undirected_pair& branch, double kT)
{
  CTAG(5,MCMC SAMPLE_BRANCH) << "Sampling alignment of branch " << tree.branch_specifier(branch) << ", kT=" << kT << "\n";
  if (align.path.count_steps_in_row(node2row[branch.first]) == 0 || align.path.count_steps_in_row(node2row[branch.second]) == 0)
    {
      CTAG(5,MCMC SAMPLE_BRANCH) << "[aborted - one of the nodes had a sequence of zero length]\n";
      return;
    }
  
  if (CTAGGING(4,BREAKDOWN))
    {
      CL << "Log-likelihood branch breakdown (before branch sampling):\n";
      show_branch_score_breakdown (CL, branch);
    }

  if (CTAGGING(1,EMIT_BREAKDOWN))
    {
      CL << "Emission log-likelihood breakdown by column (before branch sampling):\n";
      show_branch_emit_breakdown_by_column (CL, branch);
    }

  propose_sample_branch (branch, kT);

  if (CTAGGING(4,BREAKDOWN))
    {
      CL << "Log-likelihood branch breakdown (after branch sampling):\n";
      show_branch_score_breakdown (CL, branch);
    }

  if (CTAGGING(1,EMIT_BREAKDOWN))
    {
      CL << "Emission log-likelihood breakdown by column (after branch sampling):\n";
      show_branch_emit_breakdown_by_column (CL, branch);
    }

  Handel_base* old_align = 0;
  if (target_loglike)
    old_align = clone();

  proposal_accept (old_align, kT);

  if (target_loglike)
    delete old_align;
}

void Handel_base::sample_sequence (Node node, double kT)
{
  Handel_base* old_align = 0;
  if (target_loglike)
    old_align = clone();

  propose_sample_sequence (node, kT);

  proposal_accept (old_align, kT);

  if (target_loglike)
    delete old_align;
}

void Handel_base::propose_sample_branch_slide (Node grumpa, Node dad, Node son, double kT, int sample_points)
{
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

void Handel_base::sample_branch_slide (Node grumpa, Node dad, Node son, double kT, int sample_points)
{
  CTAG(5,MCMC BRANCH_SLIDE) << "Sliding node '" << tree.node_specifier(dad) << "' between '" << tree.node_specifier(grumpa) << "' and '" << tree.node_specifier(son) << "'\n";

  Handel_base* old_align = 0;
  if (target_loglike)
    old_align = clone();

  propose_sample_branch_slide (grumpa, dad, son, kT, sample_points);

  proposal_accept (old_align, kT);

  if (target_loglike)
    delete old_align;
}

bool Handel_base::sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT)
{
  CTAG(5,MCMC BRANCH_SWAP) << "Contemplating swapping aunt node '" << tree.node_specifier(aunt) << "' (parent '" << tree.node_specifier(grumpa) << "') with nephew node '" << tree.node_specifier(nephew) << "' (parent '" << tree.node_specifier(dad) << "')\n";

  Handel_base* old_align = 0;
  if (target_loglike)
    old_align = clone();

  propose_sample_branch_swap (aunt, nephew, grumpa, dad, kT);

  const bool accept = proposal_accept (old_align, kT);

  if (target_loglike)
    delete old_align;

  return accept;
}

void Handel_base::propose_sample_branch_length (const Undirected_pair& branch, double kT, double tmax, int sample_points)
{
  vector<double> x (sample_points + 1), p (sample_points + 1);
  x[0] = tree.branch_length (branch);
  for (int i = 1; i <= sample_points; ++i)
    x[i] = Rnd::prob() * tmax;
  for (int i = 0; i <= sample_points; ++i)
    {
      tree.branch_length (branch) = x[i];
      tree_changed();
      p[i] = Score2Prob ((Score) (kT * (double) alignment_score()));
    }

  const int i = Rnd::choose (p);
  tree.branch_length (branch) = x[i];
}

void Handel_base::sample_branch_length (const Undirected_pair& branch, double kT, double tmax, int sample_points)
{
  CTAG(5,MCMC BRANCH_LENGTH) << "Sampling length of branch from node '" << tree.node_specifier(branch.first) << "' to node '" << tree.node_specifier(branch.second) << "'\n";

  Handel_base* old_align = 0;
  if (target_loglike)
    old_align = clone();

  propose_sample_branch_length (branch, kT, tmax, sample_points);

  proposal_accept (old_align, kT);

  if (target_loglike)
    delete old_align;
}

void Handel_base::sample_indel_params()
{
  // base class; do nothing
  return;
}

sstring Handel_base::indel_parameter_string() const
{
   CTAG(1, PARAM_SAMPLE) << "Called Handel_base::indel_parameter_string()\n";
   sstring dummy;
   return dummy;
}

bool Handel_base::optimise_node (Node node, bool sample_seq)
{
  if (target_loglike)
    THROWEXPR ("Can't use Metropolis-Hastings sampling with optimization methods");

  if (tree.is_leaf(node)) { sstring e; e << "Attempted to optimise sequence at leaf node '" << tree.node_specifier(node) << "'\n"; THROW Standard_exception(e); }

  bool row_is_new = node2row[node] == -1;
  const Score_profile* old_profile = 0;
  if (!row_is_new)
    {
      old_profile = align.prof[node2row[node]];
      if (old_profile != 0 && !sample_seq)
	{ CTAG(5,HANDEL) << "Sequence at node '" << tree.node_specifier(node) << "' is specified; skipping optimisation step\n"; return 0; }
    }

  CTAG(5,HANDEL) << "Optimising internal node '" << tree.node_specifier(node) << "'\n";

  if (nodes_equal_rows())
    if (CTAGGING(4,BREAKDOWN))
      {
	CL << "Log-likelihood node breakdown (before node optimisation):\n";
	show_node_score_breakdown (CL, node);
      }
  
  if (row_is_new) add_empty_alignment_row_for_node (node);
  if (!row_is_new && old_profile != 0 && sample_seq)
    {
      align.prof[node2row[node]] = 0;
      if (CTAGGING(4,HANDEL))
	{
	  Score_profile cond;
	  calculate_conditional_score_profile (cond, node, -1, 1, 1);
	  CL << "Temporarily replacing sequence at node '" << tree.node_specifier(node) << "' with wildcards '*'";
	  CL << " (sequence posterior log-probability was " << Score2Bits ((*old_profile).inner_product (cond)) << " bits)\n";
	}
    }

  bool changed = propose_optimise_node (node);

  if (sample_seq)
    {
      propose_optimise_sequence (node);
      changed = true;  // HACK; probably incorrect; should test to see if new sequence is same as old
    }

  if (nodes_equal_rows())
    if (CTAGGING(4,BREAKDOWN))
      {
	CL << "Log-likelihood node breakdown (after node optimisation):\n";
	show_node_score_breakdown (CL, node);
      }

  return changed;
}

bool Handel_base::optimise_branch (const Undirected_pair& branch)
{
  if (target_loglike)
    THROWEXPR ("Can't use Metropolis-Hastings sampling with optimization methods");

  CTAG(5,HANDEL) << "Finding optimal alignment of branch " << tree.branch_specifier(branch) << "\n";
  if (align.path.count_steps_in_row(node2row[branch.first]) == 0 || align.path.count_steps_in_row(node2row[branch.second]) == 0)
    {
      CTAG(5,HANDEL) << "[aborted - one of the nodes had a sequence of zero length]\n";
      return 0;
    }

  if (CTAGGING(4,BREAKDOWN EMIT_BREAKDOWN))
    {
      CL << "Log-likelihood branch breakdown (before branch optimisation):\n";
      show_branch_score_breakdown (CL, branch);
    }

  if (CTAGGING(1,EMIT_BREAKDOWN))
    {
      CL << "Emission log-likelihood breakdown by column (before branch optimisation):\n";
      show_branch_emit_breakdown_by_column (CL, branch);
    }

  const bool changed = propose_optimise_branch (branch);

  if (CTAGGING(4,BREAKDOWN))
    {
      CL << "Log-likelihood branch breakdown (after branch optimisation):\n";
      show_branch_score_breakdown (CL, branch);
    }

  if (CTAGGING(1,EMIT_BREAKDOWN))
    {
      CL << "Emission log-likelihood breakdown by column (after branch optimisation):\n";
      show_branch_emit_breakdown_by_column (CL, branch);
    }

  return changed;
}

void Handel_base::optimise_sequence (Node node)
{
  if (target_loglike)
    THROWEXPR ("Can't use Metropolis-Hastings sampling with optimization methods");

  CTAG(5,HANDEL) << "Optimising sequence at internal node '" << tree.node_specifier(node) << "'\n";
  propose_optimise_sequence (node);
}

void Handel_base::optimise_branch_length (const Undirected_pair& branch,
					  double tmax, double tres)
{
  if (target_loglike)
    THROWEXPR ("Can't use Metropolis-Hastings sampling with optimization methods");

  CTAG(5,HANDEL) << "Optimising length of branch " << tree.branch_specifier(branch) << "\n";
  optimise_branch_length (branch, tmax, tres);
}

bool Handel_base::proposal_accept (Handel_base* old_align, double kT)
{
  if (target_loglike != 0)
    {
      const Loge old_proposal_loglike = Score2Nats (old_align->alignment_score());
      const Loge old_target_loglike = target_loglike->loglike (*old_align);
      
      const Loge new_proposal_loglike = Score2Nats (alignment_score());
      const Loge new_target_loglike = target_loglike->loglike (*this);

      const Prob accept_prob = min (1., Nats2Prob (kT * (new_target_loglike - old_target_loglike - new_proposal_loglike + old_proposal_loglike)));

      if (CTAGGING(5,MCMC BREAKDOWN))
	{
	  CL << "Pre-sampling log-likelihoods: " << Nats2Bits(old_proposal_loglike) << " bits (proposal), " << Nats2Bits(old_target_loglike) << " bits (target)\n";
	  CL << "Post-sampling log-likelihoods: " << Nats2Bits(new_proposal_loglike) << " bits (proposal), " << Nats2Bits(new_target_loglike) << " bits (target)\n";
	  CL << "P(accept) = " << accept_prob << "\n";
	}

      if (Rnd::decide (accept_prob))
	{
	  CL << "Accepting move\n";
	  return true;
	}
      else
	{
	  CL << "Rejecting move\n";
	  *this = *old_align;
	  return false;
	}
    }
  else
    return true;
}

Score Handel_base::anneal (double kT_start, double kT_end, int annealing_steps,
			   Tree_shuffler& shuffler, vector<int>& scores,
			   bool sample_seq, bool use_best,
			   bool refine, int refine_period,
			   bool refine_node_triplets)
{
  if (refine && refine_period > 0 && !use_best && sample_stream == 0)
    CLOGERR << "Warning: these annealing parameters may discard alignments!\n";

  Handel_base* best = clone();
  Score best_score = alignment_score();

  double kT_delta = (kT_end - kT_start) / (min ((double) annealing_steps, 1.0));
  double kT = kT_start;
  Score current_score = alignment_score();
  for (int step = 0;; step++)
    {
      if (sample_stream)
	{
	  write_Stockholm_with_score (*sample_stream, Stockholm_alignment_type_sampled, step);
	  sample_stream->flush();
	}
      scores.push_back (current_score);

      bool new_hiscore = 0;
      if (current_score > best_score)
	{
	  *best = *this;
	  best_score = current_score;
	  new_hiscore = 1;
	}

      if (refine && (step >= annealing_steps || (step > 0 && refine_period > 0 && step % refine_period == 0)))
	{
	  CTAG(7,HANDEL) << "Refining alignment\n";

	  Handel_base* saved_align = clone();

	  Score ref_score = alignment_score();
	  while (1)
	    {
	      const Score old_ref_score = ref_score;
	      const bool changed = refine_nodes_or_branches (refine_node_triplets && shuffler.node_realign_rate > 0., sample_seq);
	      if (!changed)
		break;
	      ref_score = alignment_score();
	      if (ref_score < old_ref_score)
		{
		  CTAG(6,HANDEL) << "Score dropped from " << Score2Bits(old_ref_score) << " to " << Score2Bits(ref_score) << "; ditching\n";
		  break;
		}
	      if (ref_score == old_ref_score) break;
	    }

	  if (sample_stream)
	    {
	      write_Stockholm_with_score (*sample_stream, Stockholm_alignment_type_refined, step);
	      sample_stream->flush();
	    }

	  log_benchmark_results ("Refined alignment");

	  if (ref_score > best_score)
	    {
	      *best = *this;
	      best_score = ref_score;
	      new_hiscore = 1;
	    }

	  *this = *saved_align;
	  delete saved_align;
	}

      if (new_hiscore)
	if (CTAGGING(6,HANDEL))
	  {
	    CL << "New best alignment:\n";
	    best->write_Stockholm (CL, alphabet());
	  }
      
      if (step >= annealing_steps) break;

      Tree_shuffler::Action action = shuffler.next_action();
      if (action == Tree_shuffler::Sample_node)
	sample_node (shuffler.next_node(), kT, sample_seq);
      else if (action == Tree_shuffler::Sample_branch)
	sample_branch (shuffler.next_branch(), kT);
      else if (action == Tree_shuffler::Scale_branch)
	sample_branch_length (shuffler.next_branch(), kT);
      else if (action == Tree_shuffler::Slide_node)
	{
	  Node grumpa, dad, son;
	  shuffler.get_next_slide (grumpa, dad, son);
	  sample_branch_slide (grumpa, dad, son, kT);
	}
      else if (action == Tree_shuffler::Swap_branches)
	{
	  Node aunt, nephew, grumpa, dad;
	  shuffler.get_next_swap (aunt, nephew, grumpa, dad);
	  if (sample_branch_swap (aunt, nephew, grumpa, dad, kT))
	    shuffler.reset();
	}
       else if (action == Tree_shuffler::Sample_indel_params)
	 {
	   // print the log message here, because the sample_indel_params() method is virtual
	   // (HACKY: dual implementations in Affine_transducer_factory & Convex_transducer_factory subclasses duplicate a lot of code.)
	   CTAG(5,MCMC PARAM_SAMPLE) << "Sampling indel parameters\n";
	   sample_indel_params();
	 }
      
      tree_changed();
      align_changed();
      current_score = alignment_score();                   // re-calculate the current score

      if (CTAGGING(7,HANDEL HANDEL_ANNEAL))
	{
	  CL << "Annealing iteration #" << step+1 << ": ";
	  CL << "score = " << Score2Bits(current_score) << " bits, ";
	  CL << "kT = " << kT << "\n";
	}
      
      if (CTAGGING(5,HANDSAMPLE))
	{
	  CL << "Sampled alignment:\n";
	  write_Stockholm (CL, alphabet());
	}
      
      log_benchmark_results ("Sampled alignment");

      kT += kT_delta;
    }

  if (use_best) { *this = *best; current_score = best_score; }
  delete best;
  return current_score;
}

void Handel_base::optimise_missing_nodes (bool sample_seq)
{
  if (target_loglike)
    THROWEXPR ("Can't use Metropolis-Hastings sampling with optimization methods");

  if (nodes_equal_rows())
    {
      bool none_missing = 1;
      if (sample_seq)
	for_iterator (Phylogeny::Node_const_iter, n, tree.internals_begin(), tree.internals_end())
	  if (align.prof[node2row[*n]] == 0)
	    {
	      optimise_sequence (*n);
	      none_missing = 0;
	    }
      if (none_missing)
	CTAG(8,HANDEL) << "Didn't need to estimate optimal paths or sequences at missing nodes\n";
      else
	{
	  CTAG(8,HANDEL) << "Didn't need to estimate optimal paths at missing nodes - only sequences\n";
	  if (CTAGGING(4,HANDEL))
	    {
	      CL << "Full alignment including estimated internal nodes:\n";
	      write_Stockholm (CL, alphabet());
	    }
	}
    }
  else
    {
      assert_leaves_equal_rows();

      vector<int> neighbours_present (tree.nodes(), 0);
      queue<Phylogeny::Node> optimise_queue;
      for (int r = 0; r < align.rows(); ++r)
	{
	  Phylogeny::Node p = tree.parent[row2node[r]];
	  if (++neighbours_present[p] == tree.neighbours(p) - 1)
	    optimise_queue.push (p);
	}

      while (!optimise_queue.empty())
	{
	  Phylogeny::Node n = optimise_queue.front();
	  Phylogeny::Node p = tree.parent[n];
	  optimise_queue.pop();
	  optimise_node (n, sample_seq);
	  if (p != -1)
	    if (++neighbours_present[p] == tree.neighbours(p) - 1)
	      optimise_queue.push (p);
	}
      
      CTAG(8,HANDEL) << "Estimation of optimal paths at missing nodes is complete\n";
      if (CTAGGING(4,HANDEL))
	{
	  CL << "Full alignment including estimated internal nodes:\n";
	  write_Stockholm (CL, alphabet());
	}
      assert_nodes_equal_rows();
    }

  if (CTAGGING(5,BREAKDOWN))
    {
      int path_score = alignment_path_score();
      int emit_score = alignment_emit_score();
      CL << "Alignment log-likelihood breakdown (total " << Score2Bits (path_score + emit_score) << " bits):\n";
      CL << "Log path likelihood = " << Score2Bits (path_score) << " bits, log emission likelihood = " << Score2Bits (emit_score) << " bits\n";
    }
}

void Handel_base::viterbi_progressive_alignment (bool sample_seq)
{
  tree.assert_tree_is_binary();

  // for each internal node, (1) find the closest child; (2) count the number of missing kids; (3) add to queue if all kids are present
  //
  Phylogeny::Node_vector fave_son (tree.nodes(), -1);
  vector<int>            missing_kids (tree.nodes(), 0);
  Phylogeny::Branch_set  branches;
  for_iterator (Phylogeny::Node_const_iter, i, tree.internals_begin(), tree.internals_end())
    {
      double fave_len = 0;
      Phylogeny::Node& fave = fave_son[*i];
      int& missing = missing_kids[*i];
      for_iterator (Phylogeny::Child_iter, c, tree.children_begin(*i,tree.parent[*i]), tree.children_end(*i,tree.parent[*i]))
	{
	  if (fave == -1 ? 1 : tree.branch_length(*i,*c) < fave_len) fave_len = tree.branch_length (*i, fave = *c);
	  if (!tree.is_leaf(*c)) ++missing;
	}
      if (missing == 0) branches.insert (Phylogeny::Branch (*i, fave, tree));
    }

  // do the progressive alignment
  //
  Alignment_path::Decomposition decomp;
  vector<Score_profile>         felsenstein (tree.nodes(), Score_profile());

  while (branches.size())
    {
      Phylogeny::Branch_set::const_iterator b = branches.begin();

      Phylogeny::Node        dad    = (*b).first;
      Phylogeny::Node        grumpa = tree.parent[dad];
      Phylogeny::Node_vector kids   = tree.children (dad, grumpa);

      vector<int> row_set (3);
      row_set[0] = node2row[dad];       // parent
      row_set[1] = node2row[kids[0]];   // child x
      row_set[2] = node2row[kids[1]];   // child y
      
      Alignment axy (align, row_set, 1);
      
      CTAG(5,HANDEL) << "Aligning nodes " << axy.row_name[1] << " and " << axy.row_name[2] << "\n";
      
      const Score_profile& x_prof = tree.is_leaf(kids[0]) ? *align.prof[node2row[kids[0]]] : felsenstein[kids[0]];
      const Score_profile& y_prof = tree.is_leaf(kids[1]) ? *align.prof[node2row[kids[1]]] : felsenstein[kids[1]];

      align_and_infer_parent (x_prof, y_prof, dad, axy.path);

      if (sample_seq)
	{
	  Score_profile cond;
	  calculate_conditional_score_profile (cond, dad, tree.parent[dad], false, false);
	  felsenstein[dad] = cond.sample_dsq (1.);
	  set_node_profile (dad, new Score_profile (felsenstein[dad]));
	}
      else
	{
	  felsenstein[dad] = Score_profile (axy.path.count_steps_in_row (0), alphabet().flat_score_map (0));
	  set_node_profile (dad, 0);
	}

      if (CTAGGING(4,ALIGN_SIBLINGS))
	{
	  axy.prof[0] = &felsenstein[dad];
	  axy.prof[1] = &x_prof;
	  axy.prof[2] = &y_prof;
	  axy.write_MUL (CL, alphabet(), 0, true);
	}
      
      decomp [Alignment_path::Row_pair (row_set[0], row_set[1])] = Pairwise_path (axy.path, 0, 1, 1);
      decomp [Alignment_path::Row_pair (row_set[0], row_set[2])] = Pairwise_path (axy.path, 0, 2, 1);
      
      // erase child branch from the queue; add shortest branch to grumpa iff we've got all the uncles
      //
      branches.erase(b);
      if (grumpa != -1)
	if (--missing_kids[grumpa] == 0)
	  branches.insert (Phylogeny::Branch (grumpa, fave_son[grumpa], tree));
    }
  align.path.compose_and_log (decomp, 1);

  if (CTAGGING(7,HANDEL))
    {
      int path_score = alignment_path_score();
      int emit_score = alignment_emit_score();
      CL << "Alignment found with log-likelihood " << Score2Bits (path_score + emit_score) << " bits\n";

      if (CTAGGING(5,HANDEL))
	CL << "Log path likelihood = " << Score2Bits (path_score) << " bits, log emission likelihood = " << Score2Bits (emit_score) << " bits\n";
    }
}

bool Handel_base::refine_nodes_or_branches (bool node_flag, bool sample_seq)
{
  bool changed = false;
  Score best_score = alignment_score();

  Phylogeny::Node_vector nodes;
  for (int node = 0; node < tree.nodes(); ++node)
    if (node_flag ? tree.is_internal(node) : node != tree.root)
      nodes.push_back (node);

  int iterations_since_last_update = 0;
  int total_iterations = 0;
  for (int i = 0; iterations_since_last_update < (int) nodes.size(); i = (i + 1) % nodes.size())
    {
      CTAG(6,HANDEL_REFINE) << (node_flag ? "Node" : "Branch")
			    << "-refinement iteration #" << total_iterations+1
			    << " (" << iterations_since_last_update << " iterations since last update; stopping at " << nodes.size()
			    << "); score " << Score2Bits(best_score)
			    << " bits, attempting to re-align "
			    << (node_flag ? "neighborhood of" : "branch to")
			    << " node " << tree.node_specifier (nodes[i])
			    << '\n';

      bool changed_this_iteration =
	node_flag
	? optimise_node (nodes[i], sample_seq)
	: optimise_branch (Phylogeny::Undirected_pair (tree.parent[nodes[i]], nodes[i]));

      if (changed_this_iteration)
	{
	  if (CTAGGING(2,HANDEL_REFINE))
	    {
	      CL << "New best alignment:\n";
	      write_Stockholm (CL, alphabet());
	    }

	  const Score new_score = alignment_score();

	  if (new_score < best_score)
	    {
	      CTAG(6,HANDEL_REFINE) << "Score decreased to " << Score2Bits(best_score) << " bits; bailing out\n";
	      break;
	    }

	  if (new_score > best_score)
	    best_score = new_score;
	  else
	    changed_this_iteration = false;
	}

      if (changed_this_iteration)
	{
	  changed = true;
	  iterations_since_last_update = 0;
	}
      else
	++iterations_since_last_update;

      ++total_iterations;
    }


  CTAG(6,HANDEL_REFINE) << "While refining "
			<< (node_flag ? "nodes" : "branches")
			<< ": found "
			<< (changed ? "optimizations, but can't find any more" : "no optimizations")
			<< "; stopping\n";
  return changed;
}

Score Handel_base::alignment_path_score() const
{
  assert_nodes_equal_rows();
  const Score root_sc = equilibrium_node_length_score (tree.root);
  const Score path_sc = conditioned_alignment_path_score (tree.root, -1);
  return ScorePMul (root_sc, path_sc);
}

Score Handel_base::column_emit_score (int col, const vector<int>& seq_coords) const
{
  // do Felsenstein's algorithm with rooted cliques.
  //
  // felsenstein[node2row[N]][X] = log(F(N,X))
  // where F(N,X) is the Felsenstein likelihood
  // - the sum of the likelihoods of all subtrees rooted at node N that have symbol X at the root node.
  //
  // Modified Felsenstein recursion is:
  //      F(N,X) = W(N,X) * product_{children C of N} sum_{symbols Y at C} F(C,Y) * S_{NC} (X,Y)
  //
  // where W(N,X) is the observed weight of symbol X at node N
  // and S_{NC} (X,Y) is the substitution matrix for the N-C branch.
  //
  // To sum over all symbols at hidden nodes N, set W(N,X)=1 for all S.
  // For observed nodes with (e.g.) an equal choice of A ambiguous symbols at node N
  // (e.g. A=2 for purine and pyrimidine symbols 'R' and 'Y' in DNA), set W(N,X) = 1/A.
  //
  // Set felsenstein[node2row[N]][X] array initially equal to log(W(N,X))
  // 
  const Symbol_score_map log1 = alphabet().flat_score_map(0);
  Symbol_score_map dummy_ssm;
  vector<Symbol_score_map> felsenstein;
  felsenstein.reserve (align.rows());

  for (int r = 0; r < align.rows(); r++)
    if (align.prof[r] == 0)
      felsenstein.push_back (log1);
    else if (seq_coords[r] == (int) align.prof[r]->size())
      felsenstein.push_back (dummy_ssm);
    else if (seq_coords[r] >= 0 && seq_coords[r] < (int) align.prof[r]->size())
      {
	const Symbol_score_map& seq_ssm = (*align.prof[r])[seq_coords[r]];
	if ((int) seq_ssm.size() > alphabet().size())
	  THROW Standard_exception ("Sequence corrupted");
	felsenstein.push_back (seq_ssm);
      }
    else
      THROWEXPR ("Alignment corrupted: row " << r << ", column " << col << ", seqpos " << seq_coords[r] << ", sequence length " << align.prof[r]->size());
  
  // the following vectors are used to locate cliques:
  //
  vector<bool> row_eligible_to_be_a_clique_root  (align.rows(), (bool) 1);
  vector<int>  parent_of_row (align.rows(), -1);

  // do the recursion, from the leaves up
  //
  for_iterator (Phylogeny::Branch_iter, b,
		tree.branches_begin (tree.root, -1),
		tree.branches_end())
    {
      const Node parent_node = (*b).first;
      const Node child_node  = (*b).second;
      const int  parent_row  = node2row[parent_node];
      const int  child_row   = node2row[child_node];

      if (parent_row != -1 && child_row != -1)
	{
	  parent_of_row[child_row] = parent_row;
	  if (align.path(parent_row,col) && align.path(child_row,col))
	    {
	      Symbol_score_map&       parent_ssm = felsenstein[parent_row];
	      const Symbol_score_map& child_ssm  = felsenstein[child_row];
	      array2d<Score>&         sub        = ((Handel_branch_scores_map&)handel_branch_scores)[Undirected_pair(*b)].cond_submat;
	    
	      for_contents (Symbol_score_map, parent_ssm, ssp)
		{
		  Score sum_score = -InfinityScore;
		  for_const_contents (Symbol_score_map, child_ssm, ssc)
		    ScorePSumAcc (sum_score, ScorePMul ((*ssc).second, sub ((*ssp).first, (*ssc).first)));
		  ScorePMulAcc ((*ssp).second, sum_score);
		}
	    
	      row_eligible_to_be_a_clique_root[child_row] = false;     // if a row has been visited as a child, then it can't be a clique root
	    }
	}
    }

  // Column likelihood = product_{clique root-nodes R} sum_{symbols X at R} P(X) * F(R,X)
  // where P(X) is the prior for symbol X
  //

  Score column_score = 0;
  map<int,Score> row_clique_score;
  for (int r = 0; r < align.rows(); r++)
    if (align.path(r,col) && row_eligible_to_be_a_clique_root[r])
      {
	Score clique_score = -InfinityScore;
	const Symbol_score_map& clique_root_ssm = felsenstein[r];
	for_const_contents (Symbol_score_map, clique_root_ssm, ssr)
	  ScorePSumAcc (clique_score, ScorePMul (handel_seq_scores.prior[(*ssr).first], (*ssr).second));
	ScorePMulAcc (column_score, clique_score);
	row_clique_score[r] = clique_score;
      }

  // dump the entire fucking mess to the logfile
  if (CTAGGING(-3,HANDEL_COLUMN_EMIT_SCORE))
    {
      // header
      CL << "DP matrix for Felsenstein pruning at column " << col
	 << ", seqCoords {" << seq_coords
	 << "}\n";

      // row-by-row
      for (int r = 0; r < align.rows(); r++)
	{
	  CL << "Row " << r
	     << ". path:" << YES_OR_NO(align.path(r,col))
	     << ", canBeRoot:" << YES_OR_NO(row_eligible_to_be_a_clique_root[r])
	     << ", parent=" << parent_of_row[r]
	     << ", matrixRow{";
	  for_contents (Symbol_score_map, felsenstein[r], ssm)
	    {
	      CL << ' ' << ssm->first << ':';
	      ShowScore (ssm->second, CL);
	    }
	  CL << " }";
	  if (row_clique_score.find (r) != row_clique_score.end())
	    CL << ", cliqueScore=" << row_clique_score[r];
	  CL << '\n';
	}
    }

  return column_score;
}

Score Handel_base::alignment_emit_score() const
{
  assert_nodes_equal_rows();
  Score score = 0;
  vector<int> seq_coords = align.path.create_seq_coords();
  for (int col = 0; col < align.columns(); col++)
    {
      ScorePMulAcc (score, column_emit_score (col, seq_coords));
      align.path.inc_seq_coords (seq_coords, col);
    }
  return score;
}

Score Handel_base::alignment_score() const
{
  return alignment_path_score() + alignment_emit_score();
}

Score Handel_base::conditioned_alignment_path_score (Node dad, Node grumpa) const
{
  Score score = 0;
  for_iterator (Phylogeny::Branch_iter, b,
		tree.branches_begin (dad, grumpa),
		tree.branches_end())
    ScorePMulAcc (score, conditioned_branch_path_score (*b));
  return score;
}

Score Handel_base::equilibrium_node_length_score (Node node) const
{
  return handel_seq_scores.seq_len_score (align.path.count_steps_in_row (node2row[node]));
}

Score Handel_base::null_emit_score() const
{
  return align.null_emit_score (handel_seq_scores.prior);
}

Score Handel_base::null_length_score() const
{
  Score sc = 0;
  for (int row = 0; row < align.rows(); ++row)
    ScorePMulAcc (sc, handel_seq_scores.seq_len_score (align.path.count_steps_in_row (row)));
  return sc;
}

Score Handel_base::null_score() const
{
  return ScorePMul (null_emit_score(), null_length_score());
}

#define PlusScore2Bits(S) (S<0?"":" +") << Score2Bits(S)
  
void Handel_base::show_branch_score_breakdown (ostream& out, const Node_pair& branch)
{
  // preserve logfile state and prepare output string
  CLOGSTREAM.save_logfile_state();
  sstring o;

  // calculate path score
  const Score p = alignment_path_score();
  o << "(";
  o << Score2Bits (conditioned_alignment_path_score (branch.first, branch.second)) << " (...<--'" << tree.node_specifier(branch.first) << "') ";
  o << PlusScore2Bits (equilibrium_node_length_score (branch.first)) << " ('" << tree.node_specifier(branch.first) << "') ";
  o << PlusScore2Bits (conditioned_branch_path_score (branch)) << " (" << tree.directed_branch_specifier(branch) << ") ";
  o << PlusScore2Bits (conditioned_alignment_path_score (branch.second, branch.first)) << " ('" << tree.node_specifier(branch.second) << "'-->...) = ";
  o << Score2Bits (p) << ") (indels)\n";

  // calculate emit score
  const Score e = alignment_emit_score();
  o << "   " << PlusScore2Bits (e) << " (substitutions)   = ";

  // sum'em and show it
  o << Score2Bits (p + e) << " bits\n";

  // restore logfile; print breakdown to output stream (which may be the logfile)
  CLOGSTREAM.restore_logfile_state();
  out << o;
}

void Handel_base::show_branch_emit_breakdown_by_column (ostream& out, const Node_pair& branch)
{
  CLOGSTREAM.save_logfile_state();
  sstring o;  // eventual output

  vector<int> seq_coords = align.path.create_seq_coords();
  for (int col = 0; col < align.columns(); col++)
    {
      show_column_emit_score_breakdown (o, col, seq_coords);
      align.path.inc_seq_coords (seq_coords, col);
    }

  CLOGSTREAM.restore_logfile_state();
  out << o;
}

void Handel_base::show_node_score_breakdown (ostream& out, Node node)
{
  CLOGSTREAM.save_logfile_state();
  sstring o;  // eventual output

  Score p_nbr = 0;
  Score p_rest = 0;
  const Score p = alignment_path_score();
  const Score e = alignment_emit_score();
  vector<int> nbr_rows;

  for_iterator (Phylogeny::Relative_iter, rel, tree.relatives_begin(node), tree.relatives_end(node))
    {
      ScorePMulAcc (p_nbr, conditioned_branch_path_score (Phylogeny::Node_pair (node, *rel)));
      ScorePMulAcc (p_rest, conditioned_alignment_path_score (*rel, node));
      if (node2row[*rel] != -1) nbr_rows.push_back (node2row[*rel]);
    }

  o << "(";
  o << Score2Bits (equilibrium_node_length_score (node)) << " ('" << tree.node_specifier(node) << "') ";
  o << PlusScore2Bits (p_nbr) << " (" << tree.node_specifier(node) << "-->neighbours) ";
  o << PlusScore2Bits (p_rest) << " (neighbours->...) = ";
  o << Score2Bits (p) << ") (indels)\n";

  o << "   " << PlusScore2Bits (e) << " (substitutions) = ";
  o << Score2Bits (ScorePMul(p,e)) << " bits\n";

  CLOGSTREAM.restore_logfile_state();
  out << o;
}

void Handel_base::show_column_emit_score_breakdown (ostream& o, int col, const vector<int>& seq_coords) const
{
  o << "Column " << col << ": " << Score2Bits (column_emit_score (col, seq_coords)) << " bits\n";
}

void Handel_base::write_Stockholm_with_score (ostream& out, const char* alignment_type, int alignment_step) const
{
  const Loge align_ll = Score2Nats (alignment_score());
  const Loge null_ll = Score2Nats (null_score());
  const Loge odds_ratio_ll = NatsPMul (align_ll, -null_ll);
  const sstring indel_param = indel_parameter_string();

  out << Stockholm_header;

  out << Stockholm_file_annotation << ' ' << Stockholm_alignment_type_tag << ' ' << alignment_type << '\n';
  out << Stockholm_file_annotation << ' ' << Stockholm_alignment_step_tag << ' ' << alignment_step + 1 << '\n';
  if (!indel_param.empty())
    out << Stockholm_file_annotation << ' ' << Stockholm_indel_parameter_tag << ' '<< indel_param << '\n';
  out << Stockholm_file_annotation << ' ' << Stockholm_bit_score_tag << ' ' << Nats2Bits (odds_ratio_ll)  << '\n';

  tree.write_Stockholm (out);
  align.write_MUL (out, alphabet(), 0, true);

  out << Stockholm_file_annotation << ' ' << Stockholm_comment_tag << ' ';
  out << "LgP(alignment|tree)= " << Nats2Bits (align_ll) << " bits; ";
  out << "LgP(sequences|unrelated)= " << Nats2Bits (null_ll) << " bits; ";
  out << "LgOddsRatio= " << Nats2Bits (odds_ratio_ll) << " bits\n";

  out << Stockholm_alignment_separator << '\n';
}

void Handel_base::calculate_conditional_score_profile (Score_profile& profile, Node dad, Node grumpa, bool with_prior, bool normalise)
{
  vector<Node> neighbors_to_include;
  for_children (tree, dad, grumpa, c)
    neighbors_to_include.push_back (*c);
  calculate_conditional_score_profile (profile, dad, neighbors_to_include, with_prior, normalise);
}

void Handel_base::calculate_conditional_score_profile (Score_profile& profile, Node dad, const vector<Node>& neighbors_to_include, bool with_prior, bool normalise)
{
  assert_nodes_equal_rows();
  int dads_row = node2row[dad];

  // save a bit of time for leaf nodes: just copy the leaf profile
  if (tree.is_leaf(dad) && neighbors_to_include.size() == 0)
    profile = *align.prof[dads_row];
  else
    {

      const Symbol_score_map wildcard_score_map = alphabet().flat_score_map (0);
      if (align.prof[dads_row] != 0)
	profile = *align.prof[dads_row];
      else
	profile = Score_profile (align.path.count_steps_in_row(dads_row), wildcard_score_map);
  
      // do Felsenstein's algorithm for each column in which dad has a wildcard
      //
      vector<int> seq_coords = align.path.create_seq_coords();
      for (int col = 0; col < align.columns(); col++)
	{
	  if (align.path(dads_row,col) && profile[seq_coords[dads_row]] == wildcard_score_map)
	    {
	      // logF[N][X] = log(F(N,X))
	      // where F(N,X) is the Felsenstein likelihood
	      // - the sum of the likelihoods of all subtrees rooted at node N that have symbol X at the root node.
	      //
	      // Modified Felsenstein recursion is:
	      //      F(N,X) = W(N,X) * product_{children C of N} sum_{symbols Y at C} F(C,Y) * S_{NC} (X,Y)
	      //
	      // where W(N,X) is the observed weight of symbol X at node N
	      // and S_{NC} (X,Y) is the substitution matrix for the N-C branch.

	      // identify branches in dad's clique
	      //
	      Phylogeny::Branch_vector dads_clique;
	      for_const_contents (vector<Node>, neighbors_to_include, child)
		for_nodes_pre (tree, *child, dad, b)
		{
		  const int r = node2row[(*b).second];
		  if (align.path(r,col))
		    dads_clique.push_back (*b);
		  else
		    b.skip_children();
		}

	      // get nodes in dad's clique
	      vector<Node> nodes_in_dads_clique (1, dad);
	      for_const_contents (Phylogeny::Branch_vector, dads_clique, b)
		nodes_in_dads_clique.push_back ((*b).second);

	      // set up initial conditions for recursion
	      map<Node,Symbol_score_map> logF;
	      const Symbol_score_map log1 = alphabet().flat_score_map(0);
	      for_const_contents (vector<Node>, nodes_in_dads_clique, n)
		{
		  const int r = node2row[*n];
		  // initially,  logF[N][X] = log(W(N,X))
		  logF[*n] =
		    align.prof[r]
		    ? (*align.prof[r])[seq_coords[r]]
		    : log1;
		}

	      // do the recursion, from the leaves up
	      //
	      for_iterator (Phylogeny::Branch_vector::reverse_iterator, b,
			    dads_clique.rbegin(),
			    dads_clique.rend())
		{
		  const Phylogeny::Node parent_node = (*b).first;
		  const Phylogeny::Node child_node  = (*b).second;

		  Symbol_score_map&       parent_ssm = logF[parent_node];
		  const Symbol_score_map& child_ssm  = logF[child_node];
		  array2d<int>&           sub        = ((Handel_branch_scores_map&)handel_branch_scores)[Undirected_pair(*b)].cond_submat;

		  for_contents (Symbol_score_map, parent_ssm, ssp)
		    {
		      int sum_score = -InfinityScore;
		      for_const_contents (Symbol_score_map, child_ssm, ssc)
			ScorePSumAcc (sum_score, ScorePMul ((*ssc).second, sub ((*ssp).first, (*ssc).first)));
		      ScorePMulAcc ((*ssp).second, sum_score);
		    }
		}

	      // add prior residue distribution
	      //
	      if (with_prior)
		for_contents (Symbol_score_map, logF[dad], ssd)
		  ScorePMulAcc ((*ssd).second, handel_seq_scores.prior [(*ssd).first]);
	  
	      // normalise to get posterior probability distribution for dad
	      //
	      if (normalise)
		{
		  int norm = -InfinityScore;
		  for_contents (Symbol_score_map, logF[dad], ssd)
		    ScorePSumAcc (norm, (*ssd).second);
		  for_contents (Symbol_score_map, logF[dad], ssd)
		    (*ssd).second -= norm;
		}

	      profile[seq_coords[dads_row]] = logF[dad];
	    }
      
	  align.path.inc_seq_coords (seq_coords, col);
	}
    }

  if (CTAGGING(-1,CONDITIONAL_SCORE_PROFILE))
    {
      CL << "Conditional score profile for node " << tree.node_specifier(dad) << " including neighbors (";
      for (int i = 0; i < (int) neighbors_to_include.size(); ++i)
	CL << (i==0 ? "" : " ") << tree.node_specifier (neighbors_to_include[i]);
      CL << "), with" << (with_prior ? "" : "out") << " prior, with" << (normalise ? "" : "out") << " normalisation\n";
      for (int col = 0; col < (int) profile.size(); ++col)
	{
	  CL << "Position " << col << ":";
	  for_const_contents (Symbol_score_map, profile[col], ss)
	    {
	      CL << " ";
	      CL.width (10);
	      ShowScore (ss->second, CL);
	      CL << "(" << alphabet().int2char (ss->first) << ")";
	    }
	  CL << "\n";
	}
    }
}

void Handel_base::propose_optimise_sequence (Phylogeny::Node node)
{
  if (target_loglike)
    THROWEXPR ("Can't use Metropolis-Hastings sampling with optimization methods");

  Score_profile cond;
  align.prof[node2row[node]] = 0;
  calculate_conditional_score_profile (cond, node, -1, 1, 1);
  Score_profile* new_profile = new Score_profile (cond.consensus_dsq());
  set_node_profile (node, new_profile);
  if (CTAGGING(4,HANDEL))
    CL << "Optimal sequence at node '" << tree.node_specifier(node) << "' has posterior log-probability " << Score2Bits (new_profile->inner_product (cond)) << " bits\n";
}

void Handel_base::propose_sample_sequence (Phylogeny::Node node, double kT)
{
  Score_profile cond;
  align.prof[node2row[node]] = 0;
  calculate_conditional_score_profile (cond, node, -1, 1, 1);
  Score_profile* new_profile = new Score_profile (cond.sample_dsq (kT));
  set_node_profile (node, new_profile);
  if (CTAGGING(4,HANDEL))
    CL << "Sampled sequence at node '" << tree.node_specifier(node) << "' has posterior log-probability " << Score2Bits (new_profile->inner_product (cond)) << " bits\n";
}
