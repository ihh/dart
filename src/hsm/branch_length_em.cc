#include "hsm/branch_length_em.h"

Branch_state_counts_map::Branch_state_counts_map (EM_matrix_base& hsm, Tree_alignment& tree_align, double prior_param) :
  hsm (hsm),
  tree_align (tree_align),
  tree (tree_align.tree),
  prior_param (prior_param),
  em_max_iter (-1),
  forgive (0),
  em_min_inc (.001)
{
  clear();
}

void Branch_state_counts_map::clear()
{
  for_rooted_branches_post (tree, b)
    branch_state_counts[*b].resize (hsm.m(), hsm.m(), 0.);
}

void Branch_state_counts_map::collect_branch_counts (const Column_matrix& colmat, double weight)
{
  for_rooted_branches_post (tree, b)
    {
      const int p = (*b).first;
      const int n = (*b).second;

      if (!colmat.gapped[p] && !colmat.gapped[n])
	{
	  Branch_state_counts& bcounts = branch_state_counts[*b];

	  for_const_contents (vector<int>, colmat.allowed[p], i)
	    for_const_contents (vector<int>, colmat.allowed[n], j)
	    {
	      const double branch_pp = colmat.branch_post_prob (n, *i, *j, tree, hsm);
	      bcounts(*i,*j) += weight * branch_pp;
	      // CTAG(1,TREE_EM_DEBUG) << "p="<<tree.node_name[p]<<" n="<<tree.node_name[n]<<" i="<<*i<<" j="<<*j<<" prob="<<branch_pp<<"\n";
	    }
	}
    }
}

void Branch_state_counts_map::update_branch_lengths (double resolution, double tmax, double tmin)
{
  for_rooted_branches_post (tree, b)
    {
      Bell_funcs bell (*this, *b, resolution, tmax, tmin);
      tree.branch_length (*b) = bell.bell_max();
    }
}

Loge Branch_state_counts_map::do_EM (double resolution, double tmax, double tmin)
{
  // log
  if (CTAGGING(5,TREE_EM))
    {
      CL << "Optimizing tree branch lengths by EM.\nTree before optimization:\n";
      tree.write (CL);
    }

  Loge best_log_likelihood = 0;
  PHYLIP_tree best_tree;

  EM_matrix_base::Update_statistics dummy_stats (0);
  EM_matrix_base::Column_matrix colmat;
  colmat.alloc (tree.nodes(), hsm.m(), false);
  const Symbol_score_map* wildcard = &hsm.alphabet().wild_ssm;

  const Alignment_path& path = tree_align.align.path;

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
      dummy_stats.clear();
      Alignment_path::Sequence_coords coords = path.create_seq_coords();
      Loge current = 0.;
      for (int col = 0; col < tree_align.align.columns(); ++col)
	{
	  colmat.initialise (tree_align, col, coords, wildcard);
	  colmat.fill_up (hsm, tree, col);
	  colmat.fill_down (hsm, tree, dummy_stats, col, 1.);
	  const Loge colmat_ll = colmat.total_log_likelihood();
	  NatsPMulAcc (current, colmat_ll);
	  CTAG(4,TREE_EM) << "Getting counts for column " << col << ", log-likelihood " << colmat_ll << "\n";
	  collect_branch_counts (colmat, 1.);
	  path.inc_seq_coords (coords, col);
	}

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
  EM_matrix_base::Alignment_matrix alnmat (hsm, tree_align, false);
  alnmat.fill_up();
  const Loge final_log_likelihood = alnmat.total_loglike;
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

Branch_expected_loglike_base::Branch_expected_loglike_base (const Branch_state_counts& counts, Substitution_matrix_factory& submat, double prior_param)
  : counts (&counts), submat (&submat), prior_param (prior_param)
{ }

Branch_expected_loglike::Branch_expected_loglike (const Branch_state_counts& counts, Substitution_matrix_factory& submat, double prior_param)
  : Branch_expected_loglike_base (counts, submat, prior_param)
{ }

double Branch_expected_loglike::operator() (double t)
{
  const array2d<double> cond = submat->create_conditional_substitution_matrix (t);
  const int states = counts->xsize();
  Loge val = -prior_param * t;
  for (int i = 0; i < states; ++i)
    for (int j = 0; j < states; ++j) {
      const Loge ij_loglike = (*counts) (i,j) * Prob2Nats (cond (i, j));
      NatsPMulAcc (val, ij_loglike);
    }
  return val;
}

Branch_expected_loglike_deriv::Branch_expected_loglike_deriv (const Branch_state_counts& counts, Substitution_matrix_factory& submat, double prior_param)
  : Branch_expected_loglike_base (counts, submat, prior_param)
{ }

double Branch_expected_loglike_deriv::operator() (double t)
{
  const array2d<double> cond = submat->create_conditional_substitution_matrix (t);
  const array2d<double> deriv = submat->differentiate_conditional_substitution_matrix (t);
  const int states = counts->xsize();
  double val = -prior_param;   // type of val is d/dt(Loge), we currently don't have a typedef for such a thing, but double will do.... IH 4/13/2010
  for (int i = 0; i < states; ++i)
    for (int j = 0; j < states; ++j) {
      const Loge ij_loglike_deriv = (*counts) (i,j) * deriv (i, j) / cond (i, j);
      NatsPMulAcc (val, ij_loglike_deriv);
    }
  return val;
}

Bell_funcs::Bell_funcs (const Branch_state_counts_map& tree_counts, const Phylogeny::Undirected_pair& branch, double tres, double tmax, double tmin)
  : Cached_function <Branch_expected_loglike, Branch_expected_loglike_deriv> (func, deriv, tmin, tmax, tres)
{
  const Branch_state_counts& bcounts = tree_counts.branch_state_counts.find(branch)->second;
  init (bcounts, tree_counts.hsm, tree_counts.prior_param);

  CTAG(2,TREE_EM) << "Branch counts for branch " << tree_counts.tree.branch_specifier (branch)
		  << ":\n" << bcounts;
}

Bell_funcs::Bell_funcs (const Branch_state_counts& branch_counts, Substitution_matrix_factory& submat, double prior_param, double tres, double tmax, double tmin)
  : Cached_function <Branch_expected_loglike, Branch_expected_loglike_deriv> (func, deriv, tmin, tmax, tres)
{
  init (branch_counts, submat, prior_param);
  CTAG(2,TREE_EM) << "Branch counts:\n" << branch_counts;
}

void Bell_funcs::init (const Branch_state_counts& bcounts, Substitution_matrix_factory& submat, double prior_param)
{
  func = Branch_expected_loglike (bcounts, submat, prior_param);
  deriv = Branch_expected_loglike_deriv (bcounts, submat, prior_param);
}

double Bell_funcs::bell_max()
{
  // protect against small counts
  double total_count = 0.;
  for_const_contents (Branch_state_counts, *func.counts, c)
    total_count += *c;
  if (total_count < TINY)
    {
      CTAG(2,TREE_EM) << "Branch counts negligible (" << total_count << " residues observed); setting branch length to " << min_arg << "\n";
      return min_arg;
    }
  // delegate to base class
  return find_max();
}

Loge EM_tree_alignment_database::optimise_branch_lengths_by_EM (EM_matrix_base& hsm, double prior_param, int em_max_iter, int forgive, double em_min_inc, double resolution, double tmax, double tmin)
{
  CTAG(7,TREE_EM) << "Optimizing branch lengths of all trees in alignment database, using EM\n";

  Loge loglike = 0.;
  for_contents (list<Tree_alignment>, tree_align_list, tree_align)
    {
      Branch_state_counts_map tree_counts (hsm, *tree_align, prior_param);
      tree_counts.em_max_iter = em_max_iter;
      tree_counts.forgive = forgive;
      tree_counts.em_min_inc = em_min_inc;

      const Loge ta_ll = tree_counts.do_EM (resolution, tmax, tmin);
      NatsPMulAcc (loglike, ta_ll);
    }

  CTAG(6,TREE_EM) << "Alignment database log-likelihood: " << Nats2Bits(loglike) << " bits\n";
  return loglike;
}
