#ifndef EM_MATRIX
#define EM_MATRIX

#include "hsm/em_matrix_base.h"
#include "hsm/newtonraphson.h"
#include "util/vector_output.h"
#include "ecfg/ecfgkeywords.h"

typedef EM_matrix_base::Update_statistics Update_statistics;
typedef EM_matrix_base::Column_matrix Column_matrix;

// Reversible EM algorithm
// see em_matrix_base.h
struct EM_matrix : EM_matrix_base
{
  // Newton-Raphson params. NB Newton-Raphson methods are obsolete (see below).
  double nr_tol;  // tolerance
  int nr_max_iter;  // max iterations
  vector<sstring> param_label;  // param labels

  // do_equilibrate flag added by IH, 4/20/2005
  // If true, call equilibrate() to incorporate start counts into
  // wait times & transition usage counts.
  // If false (default), discard start counts for reversible process.
  bool do_equilibrate;

  // methods
  // constructor
  EM_matrix (int C,
	     int A,
	     int max_fork = 1,
	     const Tree_alignment_database* align_db = 0,
	     double timepoint_res = DEFAULT_TIMEPOINT_RES);

  // descriptor
  const char* update_policy() const { return EG_POLICY_REV; }

  // initialisers
  void init_nr();	// Newton-Raphson specific initialisation
  void init_param_labels (const Alphabet& base_alphabet);
  void init_alphabet (const Alphabet& base_alphabet);

  // reversible diagonalize() method
  void diagonalize();

  // reversible transform_symmetrised... method
  void transform_symmetrised_waits_transitions (Update_statistics& stats, bool symmetrise) const;

  // up/down algorithm
  Update_statistics get_stats (bool infer_class_labels = FALSE) const;

  // reversible M-step hacks
  void equilibrate (Update_statistics& stats) const;  // adds s[] counts into w[] and u()
  void update_pi (const Update_statistics& stats);	// pi M-step, specific to reversible matrices

  // Newton-Raphson EM
  Update_statistics single_EM (bool infer_class_labels = FALSE);
  Loge iterate_EM (bool infer_class_labels = FALSE);  // returns best log-likelihood

  // "quick" EM
  Update_statistics single_quick_EM (bool intra = TRUE, bool inter = TRUE, bool infer_class_labels = FALSE); // does partial M-step only
  void quick_M (const Update_statistics& stats, bool intra = TRUE, bool inter = TRUE);  // the M-step

  // interface to phylo_em.h
  array2d<double> clean_phylo_EM (int src, int dest, int T) const;

  //////////////////////////////////
  // BEGIN NEWTON-RAPHSON METHODS //
  //////////////////////////////////
  //
  // These methods are very outdated, since we no longer attempt to maximise a tricky expected log-likelihood
  // (e.g. including reversibility constraints etc) but instead use a simple, constraint-free setup
  // and impose reversibility by other means (forcibly symmetrising the dataset)
  //
  // methods to get params from a single vector (including Lagrange multipliers), for Newton-Raphson
  void set_params_from_vector (const vector<double>& params);

  // indices
  inline int pi_idx (int c, int a) const         { return ca(c,a); }
  inline int X_idx (int c, int ai, int aj) const
  	{ return pi_idx(C,0) + c*A*(A-1)/2 + ai*A - (ai*(ai+3))/2 + aj - 1; }
  inline int Y_idx (int ci, int cj, int a) const
  	{ return X_idx(C,0,1) + a*C*(C-1)/2 + ci*C - (ci*(ci+3))/2 + cj - 1; }
  inline int beta_idx() const                    { return Y_idx(0,1,A); }
  inline int kappa_idx (int c, int a) const      { return beta_idx() + 1 + ca(c,a); }
  inline int n_params() const                    { return kappa_idx(C,0); }

  // method returning constraints for Newton-Raphson
  map<int,double> param_minima() const;
  
  // Diffs
  Diff pi_diff (int c, int a) const;
  Diff X_diff (int c, int ai, int aj) const;
  Diff Y_diff (int ci, int cj, int a) const;

  // Lagrangian ("L" function in Holmes/Rubin paper)
  Diff lagrangian (const Update_statistics& stats);

  // unconstrained optimal params. used as seed for Newton-Raphson
  vector<double> unconstrained_optimum (const Update_statistics& stats) const;

  //
  ////////////////////////////////
  // END NEWTON-RAPHSON METHODS //
  ////////////////////////////////
};

typedef EM_matrix Reversible_EM_matrix;

#endif /* EM_MATRIX */
