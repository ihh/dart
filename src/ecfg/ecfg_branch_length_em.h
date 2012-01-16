#ifndef ECFG_BRANCH_LENGTH_EM_INCLUDED
#define ECFG_BRANCH_LENGTH_EM_INCLUDED

#include "hsm/branch_length_em.h"
#include "ecfg/ecfgdp.h"

// branch update statistics, indexed by ECFG chain index
typedef vector<Branch_state_counts> ECFG_branch_state_counts;

// tree update statistics
struct ECFG_branch_state_counts_map
{
  // typedefs
  typedef Phylogeny::Undirected_pair Undirected_pair;
  typedef EM_matrix_base::Column_matrix Column_matrix;

  // data
  map<Undirected_pair,ECFG_branch_state_counts> branch_state_counts;
  ECFG_scores& ecfg;
  Stockholm& stock;
  Tree_alignment& tree_align;  // this must be constructed from stock
  double prior_param;  // P(branch_length=t) = exp(-prior_param*t)

  // EM params
  int em_max_iter, forgive;
  double em_min_inc;

  // constructor
  ECFG_branch_state_counts_map (ECFG_scores& ecfg, Stockholm& stock, Tree_alignment& tree_align, double prior_param = 0.);

  // reset method
  void clear();

  // method to populate branch_state_counts
  // returns log-likelihood (excluding grammar rule probabilities)
  // At the moment, the counts are conditioned on a particular parse tree; probably would be better to sum over all parse trees
  Loge collect_branch_counts (ECFG_EM_matrix& em_matrix, const ECFG_cell_score_map& annot, double weight = 1.);

  // tree update method
  void update_branch_lengths (ECFG_EM_matrix& em_matrix, double resolution = TINY, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);

  // EM method: optimizes all the branch lengths of the tree
  // At the moment, the counts are conditioned on the CYK parse tree; probably would be better to sum over all parse trees
  Loge do_EM (double resolution = TINY, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);
};

// base class for ECFG_branch_expected_loglike and ECFG_branch_expected_loglike_deriv
struct ECFG_branch_expected_loglike_base
{
  // typedefs
  typedef double argument_type;
  typedef double result_type;

  // data
  const ECFG_branch_state_counts* counts;
  ECFG_scores* ecfg;
  double prior_param;

  // constructors
  ECFG_branch_expected_loglike_base() { }
  ECFG_branch_expected_loglike_base (const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param);
};

// function returning expected log-likelihood, given expected branch counts
struct ECFG_branch_expected_loglike : ECFG_branch_expected_loglike_base
{
  // data
  vector<Branch_expected_loglike> chain_bell;

  // constructors
  ECFG_branch_expected_loglike() { }
  ECFG_branch_expected_loglike (const ECFG_EM_matrix& em_matrix, Phylogeny::Node branch_idx, const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param);
  ECFG_branch_expected_loglike (const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param);

  // evaluation
  Loge operator() (double t);
};

// function returning time derivative of expected log-likelihood, given expected branch counts
struct ECFG_branch_expected_loglike_deriv : ECFG_branch_expected_loglike_base
{
  // data
  vector<Branch_expected_loglike_deriv> chain_bell_deriv;

  // constructors
  ECFG_branch_expected_loglike_deriv() { }
  ECFG_branch_expected_loglike_deriv (const ECFG_EM_matrix& em_matrix, Phylogeny::Node branch_idx, const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param);
  ECFG_branch_expected_loglike_deriv (const ECFG_branch_state_counts& counts, ECFG_scores& ecfg, double prior_param);

  // evaluation
  double operator() (double t);
};

// struct encapsulating time-likelihood function, derivative & caches
// "bell_funcs" stands for "Branch Expected Log-Likelihood FUNCtions"
struct ECFG_bell_funcs : Cached_function <ECFG_branch_expected_loglike, ECFG_branch_expected_loglike_deriv>
{
  // data
  ECFG_branch_expected_loglike func;
  ECFG_branch_expected_loglike_deriv deriv;

  // constructors
  ECFG_bell_funcs (const ECFG_EM_matrix& em_matrix, const ECFG_branch_state_counts_map& tree_counts, const Phylogeny::Undirected_pair& branch, double tres = .01, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);
  ECFG_bell_funcs (const ECFG_branch_state_counts& branch_counts, ECFG_scores& ecfg, double prior_param = 0., double tres = .01, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);

  // find_max wrapper
  double bell_max();

  // helper init methods
  void init (const ECFG_EM_matrix& em_matrix, Phylogeny::Node branch_idx, const ECFG_branch_state_counts& branch_counts, ECFG_scores& ecfg, double prior_param);
  void init (const ECFG_branch_state_counts& branch_counts, ECFG_scores& ecfg, double prior_param);
};

// extension to EM_tree_alignment_database incorporating branch length EM for entire ECFG
struct ECFG_EM_tree_alignment_database : EM_tree_alignment_database
{
  // constructors
  ECFG_EM_tree_alignment_database (Sequence_database& seq_db)
    : EM_tree_alignment_database (seq_db)
  { }

  // ECFG branch-length EM method; returns final log-likelihood
  Loge optimise_branch_lengths_by_ECFG_EM (ECFG_scores& ecfg, double prior_param = 0.,
					   int em_max_iter = 0, int forgive = 0, double em_min_inc = .001,
					   double time_resolution = TINY, double time_max = DART_MAX_BRANCH_LENGTH, double time_min = 0.);
};

#endif /* ECFG_BRANCH_LENGTH_EM_INCLUDED */
