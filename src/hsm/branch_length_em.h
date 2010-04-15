#ifndef BRANCH_LENGTH_EM_INCLUDED
#define BRANCH_LENGTH_EM_INCLUDED

#include "hsm/em_matrix_base.h"
#include "util/maximise.h"

// Sadly, a subset of the functionality in this file is duplicated by tree/subdistmat.*
// e.g.
//  Substitution_counts <--> Branch_state_counts_map
//  Subst_log_like <--> Branch_expected_loglike
//  Subst_log_like_dt <--> Branch_expected_loglike_deriv
// Time to refactor...
// IH, 4/15/2010

// branch update statistics
// these can be used to re-estimate branch lengths by EM, among other things
typedef array2d<double> Branch_state_counts;

// tree update statistics
struct Branch_state_counts_map
{
  // typedefs
  typedef Phylogeny::Undirected_pair Undirected_pair;
  typedef EM_matrix_base::Column_matrix Column_matrix;

  // data
  map<Undirected_pair,Branch_state_counts> branch_state_counts;
  EM_matrix_base& hsm;
  Tree_alignment& tree_align;
  PHYLIP_tree& tree;
  double prior_param;  // P(branch_length=t) = exp(-prior_param*t)

  // EM params
  int em_max_iter, forgive;
  double em_min_inc;

  // constructor
  Branch_state_counts_map (EM_matrix_base& hsm, Tree_alignment& tree_align, double prior_param = 0.);

  // reset method
  void clear();

  // column fill method
  void collect_branch_counts (const Column_matrix& colmat, double weight = 1.);

  // tree update method
  void update_branch_lengths (double resolution = TINY, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);

  // EM method for alignment
  Loge do_EM (double resolution = TINY, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);
};

// base class for Branch_expected_loglike and Branch_expected_loglike_deriv
struct Branch_expected_loglike_base
{
  // typedefs
  typedef double argument_type;
  typedef double result_type;

  // data
  const Branch_state_counts* counts;
  Substitution_matrix_factory* submat;
  double prior_param;

  // constructors
  Branch_expected_loglike_base() { }
  Branch_expected_loglike_base (const Branch_state_counts& counts, Substitution_matrix_factory& submat, double prior_param);
};

// function returning expected log-likelihood, given expected branch counts
struct Branch_expected_loglike : Branch_expected_loglike_base
{
  // constructors
  Branch_expected_loglike() { }
  Branch_expected_loglike (const Branch_state_counts& counts, Substitution_matrix_factory& submat, double prior_param);

  // evaluation
  double operator() (double t);
};

// function returning time derivative of expected log-likelihood, given expected branch counts
struct Branch_expected_loglike_deriv : Branch_expected_loglike_base
{
  // constructors
  Branch_expected_loglike_deriv() { }
  Branch_expected_loglike_deriv (const Branch_state_counts& counts, Substitution_matrix_factory& submat, double prior_param);

  // evaluation
  double operator() (double t);
};

// struct encapsulating time-likelihood function, derivative & caches
// "Bell_funcs" stands for "Branch Expected Log-Likelihood FUNCtions"
struct Bell_funcs : Cached_function <Branch_expected_loglike, Branch_expected_loglike_deriv>
{
  // data
  Branch_expected_loglike func;
  Branch_expected_loglike_deriv deriv;

  // constructors
  Bell_funcs (const Branch_state_counts_map& tree_counts, const Phylogeny::Undirected_pair& branch, double tres = .01, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);
  Bell_funcs (const Branch_state_counts& branch_counts, Substitution_matrix_factory& submat, double prior_param = 0., double tres = .01, double tmax = DART_MAX_BRANCH_LENGTH, double tmin = 0.);

  // find_max wrapper
  double bell_max();

  // helper init method
  void init (const Branch_state_counts& branch_counts, Substitution_matrix_factory& submat, double prior_param);
};

// extension to Tree_alignment_database incorporating branch length EM
struct EM_tree_alignment_database : Tree_alignment_database
{
  // constructors
  EM_tree_alignment_database (Sequence_database& seq_db)
    : Tree_alignment_database (seq_db)
  { }

  EM_tree_alignment_database (Sequence_database& seq_db, const char* index_filename)
    : Tree_alignment_database (seq_db, index_filename)
  { }

  // Branch-length EM method; returns final log-likelihood
  Loge optimise_branch_lengths_by_EM (EM_matrix_base& hsm, double prior_param = 0., int em_max_iter = 0, int forgive = 0, double em_min_inc = .001, double time_resolution = TINY, double time_max = DART_MAX_BRANCH_LENGTH, double time_min = 0.);
};

#endif /* BRANCH_LENGTH_EM_INCLUDED */
