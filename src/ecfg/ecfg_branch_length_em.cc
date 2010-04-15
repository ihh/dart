#include "ecfg/ecfg_branch_length_em.h"

ECFG_branch_state_counts_map::ECFG_branch_state_counts_map (ECFG_scores& ecfg, Stockholm& stock, Tree_alignment& tree_align, double prior_param) :
  ecfg (ecfg),
  stock (stock),
  tree_align (tree_align),
  tree (tree_align.tree),
  prior_param (prior_param),
  em_max_iter (-1),
  forgive (0),
  em_min_inc (.001)
{
  clear();
}

void ECFG_branch_state_counts_map::clear()
{
  for_rooted_branches_post (tree, b) {
    branch_state_counts[*b] = vector<Branch_state_counts> (ecfg.matrix_set.chain.size());
    for (int c = 0; c < (int) ecfg.matrix_set.chain.size(); ++c) {
      const int states = ecfg.matrix_set.total_states(c);
      branch_state_counts[*b][c].resize (states, states, 0.);
    }
  }
}

Loge ECFG_branch_state_counts_map::collect_branch_counts (ECFG_EM_matrix& em_matrix, const ECFG_cell_score_map& annot, double weight)
{
  // final log-likelihood
  Loge final_ll = 0.;

  // create dummy counts
  ECFG_counts dummy_counts (ecfg);

  // loop over nonterminals in parse tree, using reverse order for easy comparison with hsm/branch_length_em.cc
  for_const_reverse_contents (ECFG_cell_score_map, annot, ecsm_ptr)
    {
      const Subseq_coords& coords = ecsm_ptr->first.first;
      const int ecfg_state = ecsm_ptr->first.second;

      // get state info
      const ECFG_state_info& info = ecfg.state_info[ecfg_state];
      if (info.emit_size())
	{
	  // get chain
	  const int chain_idx = info.matrix;
	  const ECFG_chain& chain = ecfg.matrix_set.chain[chain_idx];

	  // call fill_down
	  em_matrix.fill_down (dummy_counts, em_matrix.env.find_subseq_idx (coords.start, coords.len), ecfg_state, 1.);

	  // accumulate log-likelihood
	  EM_matrix_base::Column_matrix colmat = em_matrix.colmat[ecfg_state];
	  const Loge colmat_ll = colmat.total_log_likelihood();
	  NatsPMulAcc (final_ll, colmat_ll);

	  // log
	  CTAG(4,TREE_EM) << "Getting counts for nonterminal " << info.name << " at subsequence (" << coords.start << ".." << coords.end()
			  << "); log-likelihood " << colmat_ll << "\n";

	  // loop over branches
	  for_rooted_branches_post (tree, b)
	    {
	      const int p = (*b).first;
	      const int n = (*b).second;

	      if (!colmat.gapped[p] && !colmat.gapped[n])
		{
		  Branch_state_counts& bcounts = branch_state_counts[*b][chain_idx];

		  for_const_contents (vector<int>, colmat.allowed[p], i)
		    for_const_contents (vector<int>, colmat.allowed[n], j)
		    {
		      const double branch_pp = colmat.branch_post_prob (n, *i, *j, tree, *chain.matrix);
		      bcounts(*i,*j) += weight * branch_pp;
		    }
		}
	    }
	}
    }

  // return
  return final_ll;
}

void ECFG_branch_state_counts_map::update_branch_lengths (double resolution, double tmax, double tmin)
{
  for_rooted_branches_post (tree, b)
    {
      ECFG_bell_funcs bell (*this, *b, resolution, tmax, tmin);
      tree.branch_length (*b) = bell.bell_max();
    }
}

Loge ECFG_branch_state_counts_map::do_EM (double resolution, double tmax, double tmin)
{
  // log
  if (CTAGGING(5,TREE_EM))
    {
      CL << "Optimizing tree branch lengths by EM.\nTree before optimization:\n";
      tree.write (CL);
    }

  Loge best_log_likelihood = 0;
  PHYLIP_tree best_tree;

  // get CYK matrix & trace
  const int seqlen = tree_align.align.columns();
  const int max_subseq_len = -1;   // unlimit subsequence lengths
  ECFG_auto_envelope env (seqlen, ecfg, max_subseq_len);  // create fold envelope

  Aligned_score_profile asp;
  asp.init (tree_align.align, stock.np, ecfg.alphabet);

  const bool try_fast_prune = false;
  ECFG_CYK_matrix cyk_mx (ecfg, stock, asp, env, try_fast_prune);   // create CYK matrix
  cyk_mx.fill();
  ECFG_cell_score_map cyk_trace = cyk_mx.traceback();  // get traceback

  int dec = 0;
  for (int iter = 0; ; ++iter)
    {
      // check for max iterations
      if (iter >= em_max_iter && em_max_iter >= 0)
	{
	  CTAG(7,TREE_EM) << "EM hit " << em_max_iter << " iterations; stopping\n";
	  break;
	}

      // collect update statistics
      clear();
      Loge current = collect_branch_counts (cyk_mx, cyk_trace, 1.);

      // check for best likelihood
      const Loge prev_best = best_log_likelihood;
      if (iter == 0 || current > best_log_likelihood)
	{
	  best_log_likelihood = current;
	  best_tree = tree;
	}
      CTAG(6,TREE_EM) << "EM iteration #" << iter+1 << ": log-likelihood = " << Nats2Bits(current) << " bits\n";

      // check for likelihood decrease
      if (iter > 0)
	{
	  const double inc = (current - prev_best) / (abs(prev_best) < TINY ? 1. : -prev_best);  // IH, 4/20/2005
	  CTAG(3,TREE_EM) << "(previous best = " << Nats2Bits(prev_best) << " bits; fractional improvement = " << inc << ")\n";
	  if (inc < em_min_inc)
	    {
	      if (current < prev_best)
		CTAG(7,TREE_EM) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_best) << " to " << Nats2Bits(current) << " bits during EM\n";
	      if (++dec > forgive)
		{
		  CTAG(7,TREE_EM) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
		  break;
		}
	    }
	  else
	    dec = 0;
	}

      // update the tree
      update_branch_lengths (resolution, tmax, tmin);
      tree_align.tree_changed();
      copy_branch_lengths (tree, cyk_mx.tree);

      // log
      if (CTAGGING(5,TREE_EM))
	{
	  CL << "Optimized tree after step #" << iter+1 << ":\n";
	  tree.write (CL);
	}
    }
  
  // since we only update the best log-likelihood the round *after* the new params are calculated, we need to check for an increase one last time
  CTAG(6,TREE_EM) << "Checking post-iteration log-likelihood\n";
  // collect update statistics
  const Loge final_log_likelihood = collect_branch_counts (cyk_mx, cyk_trace, 1.);
  CTAG(6,TREE_EM) << "Post-iteration log-likelihood = " << Nats2Bits(final_log_likelihood) << " bits\n";
  if (em_max_iter != 0 && final_log_likelihood < best_log_likelihood)
  {
    CTAG(6,TREE_EM) << "Restoring previous best branch lengths\n";
    tree = best_tree;
    tree_align.tree_changed();
  }
  else
    best_log_likelihood = final_log_likelihood;

  // log
  if (CTAGGING(5,TREE_EM))
    {
      CL << "Tree after optimization:\n";
      tree.write (CL);
    }

  // and return
  return best_log_likelihood;
}

void ECFG_branch_state_counts_map::copy_branch_lengths (const PHYLIP_tree& t1, PHYLIP_tree& t2)
{
  for_rooted_branches_pre (t1, b)
    t2.branch_length(*b) = t1.branch_length(*b);
}

ECFG_branch_expected_loglike_base::ECFG_branch_expected_loglike_base (const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param)
  : counts (&counts), ecfg (&ecfg), prior_param (prior_param)
{
  if (counts.size() != ecfg.matrix_set.chain.size())
    THROWEXPR ("Size of ECFG_branch_state_counts doesn't match number of chains in ECFG");
}

ECFG_branch_expected_loglike::ECFG_branch_expected_loglike (const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param)
  : ECFG_branch_expected_loglike_base (counts, ecfg, prior_param)
{
  for (int c = 0; c < (int) counts.size(); ++c)
    chain_bell.push_back (Branch_expected_loglike (counts[c], *ecfg.matrix_set.chain[c].matrix, 0.));
}

double ECFG_branch_expected_loglike::operator() (double t)
{
  Loge val = -prior_param * t;
  for_contents (vector<Branch_expected_loglike>, chain_bell, cb)
    NatsPMulAcc (val, (*cb)(t));
  return val;
}

ECFG_branch_expected_loglike_deriv::ECFG_branch_expected_loglike_deriv (const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param)
  : ECFG_branch_expected_loglike_base (counts, ecfg, prior_param)
{
  for (int c = 0; c < (int) counts.size(); ++c)
    chain_bell_deriv.push_back (Branch_expected_loglike_deriv (counts[c], *ecfg.matrix_set.chain[c].matrix, 0.));
}

double ECFG_branch_expected_loglike_deriv::operator() (double t)
{
  double val = -prior_param;
  for_contents (vector<Branch_expected_loglike_deriv>, chain_bell_deriv, cb)
    NatsPMulAcc (val, (*cb)(t));
  return val;
}

ECFG_bell_funcs::ECFG_bell_funcs (const ECFG_branch_state_counts_map& tree_counts, const Phylogeny::Undirected_pair& branch, double tres, double tmax, double tmin)
  : Cached_function <ECFG_branch_expected_loglike, ECFG_branch_expected_loglike_deriv> (func, deriv, tmin, tmax, tres)
{
  const ECFG_branch_state_counts& bcounts = tree_counts.branch_state_counts.find(branch)->second;

  func = ECFG_branch_expected_loglike (bcounts, tree_counts.ecfg, tree_counts.prior_param);
  deriv = ECFG_branch_expected_loglike_deriv (bcounts, tree_counts.ecfg, tree_counts.prior_param);

  CTAG(2,TREE_EM) << "Branch counts for branch " << tree_counts.tree.branch_specifier (branch)
		  << ":\n" << bcounts;
}

double ECFG_bell_funcs::bell_max()
{
  // protect against small counts
  double total_count = 0.;
  for (int chain = 0; chain < (int) func.counts->size(); ++chain)
    for_const_contents (Branch_state_counts, (*func.counts)[chain], c)
      total_count += *c;
  if (total_count < TINY)
    {
      CTAG(2,TREE_EM) << "Branch counts negligible (" << total_count << " residues observed); setting branch length to " << min_arg << "\n";
      return min_arg;
    }
  // delegate to base class
  return find_max();
}

Loge ECFG_EM_tree_alignment_database::optimise_branch_lengths_by_ECFG_EM (ECFG_scores& ecfg, double prior_param, int em_max_iter, int forgive, double em_min_inc, double resolution, double tmax, double tmin)
{
  CTAG(7,TREE_EM) << "Optimizing branch lengths of all trees in alignment database, using EM\n";

  if (stock_db == 0)
    THROWEXPR("ECFG_EM_tree_alignment_database must be initialized from a Stockholm database in order to do branch-length EM");

  Loge loglike = 0.;
  list<Stockholm>::iterator stock_iter = stock_db->align.begin();
  for_contents (list<Tree_alignment>, tree_align_list, tree_align)
    {
      ECFG_branch_state_counts_map tree_counts (ecfg, *stock_iter, *tree_align, prior_param);
      tree_counts.em_max_iter = em_max_iter;
      tree_counts.forgive = forgive;
      tree_counts.em_min_inc = em_min_inc;

      const Loge ta_ll = tree_counts.do_EM (resolution, tmax, tmin);
      NatsPMulAcc (loglike, ta_ll);

      ++stock_iter;
    }

  CTAG(6,TREE_EM) << "Alignment database log-likelihood: " << Nats2Bits(loglike) << " bits\n";
  return loglike;
}
