#ifndef RIVAS_TRANSDUCER_FACTORY_INCLUDED
#define RIVAS_TRANSDUCER_FACTORY_INCLUDED

#include "handel/alitrans.h"
#include "irrev/irrev_em_matrix.h"
#include "util/rnd.h"
#include "util/score.h"
#include "randlib/randlib.h"

// Rivas_transducer_factory
// Implements a class of time-parametric transducers described in the following article:
// Rivas E. Evolutionary models for insertions and deletions in a probabilistic modeling framework.
// BMC Bioinformatics. 2005 Mar 21;6(1):63. 
struct Rivas_transducer_factory : Transducer_alignment
{
  // data
  Pair_transducer_scores prior;
  Pair_transducer_scores branch;

  // banding coefficient variables
  double del_rate, mean_del_size;

  // substitution model
  Irrev_EM_matrix subst_model;

  // geometric equilibrium length parameter
  const double gamma_val;

  // models for evolution of transition probabilities
  typedef map<int,Irrev_EM_matrix*> Trans_model_map;
  typedef map<int,vector<int> > Dest_states_map;
  Trans_model_map trans_model;  // transition evolution matrices by source state
  Dest_states_map dest_states;  // allowable dest states from each evolving source

  // constructor
  Rivas_transducer_factory (int branch_states, double gamma);

  // destructor
  ~Rivas_transducer_factory();

  // gamma method
  double gamma() const { return gamma_val; }

  // submat_factory method
  Substitution_matrix_factory& submat_factory() const { return (Substitution_matrix_factory&) subst_model; }

  // clone method
  Rivas_transducer_factory* clone();

  // Pair_transducer_factory methods
  Pair_transducer_scores prior_pair_trans_sc();
  Pair_transducer_scores branch_pair_trans_sc (double time);

  double gap_rate() { return del_rate; }
  double mean_gap_size() { return mean_del_size; }

  // helper methods to set state types for prior & branch transducers
  void init_prior_wait   (int state, const char* name) { prior.state_type[state] = 0; prior.transition(state,End) = 0; prior.state_name[state] = name; }
  void init_prior_ins    (int state, const char* name) { prior.state_type[state] = 2; prior.state_name[state] = name; }
  void init_branch_wait  (int state, const char* name) { branch.state_type[state] = 0; branch.transition(state,End) = 0; branch.state_name[state] = name; }
  void init_branch_del   (int state, const char* name) { branch.state_type[state] = 1; branch.state_name[state] = name; }
  void init_branch_ins   (int state, const char* name) { branch.state_type[state] = 2; branch.state_name[state] = name; }
  void init_branch_match (int state, const char* name) { branch.state_type[state] = 3; branch.state_name[state] = name; }

  // helper methods to create & populate Irrev_EM_matrix objects for a state
  void init_trans_model (int src_state, int n_dest_states);
  int  get_dest_state_index (int src_state, int dest_state);
  void set_trans_prob_init (int src_state, int dest_state, Prob prob);
  void set_trans_prob_rate (int src_state, int old_dest_state, int new_dest_state, double rate);
  void update_trans_probs();  // calls update() on the Irrev_EM_matrix objects in trans_model
};


struct Affine_transducer_factory_param_container
{
  double gamma_param;
  double delete_rate_param;
  double delete_extend_prob_param;
};

struct Convex_transducer_factory_param_container
{
  double gamma_param;
  double delete_rate_param;
  vector<double> cpt_weight_param;
  vector<double> cpt_delete_extend;
};

// Affine_transducer_factory
// A Rivas_transducer_factory for simple affine gap models (geometric distribution)
struct Affine_transducer_factory : Affine_transducer_factory_param_container, Rivas_transducer_factory
{
  static double alpha, beta, gamma_alpha, gamma_beta, delete_extend_alpha, delete_extend_beta; 

  Affine_transducer_factory (double gamma, double delete_rate, double delete_extend_prob);
  
  void set_transitions (double gamma, double delete_rate, double delete_extend_prob);

  void sample_indel_params ();

  sstring indel_parameter_string() const;

  Affine_transducer_factory_param_container* propose_indel_params ();
};

// Convex_transducer_factory
// A Rivas_transducer_factory for convex gap models (mixture-of-geometric distribution)
struct Convex_transducer_factory : Convex_transducer_factory_param_container, Rivas_transducer_factory
{
  static double alpha, beta, gamma_alpha, gamma_beta, cpt1_alpha, cpt1_beta, cpt2_alpha, cpt2_beta; 
  // constructor
  Convex_transducer_factory (double gamma, double delete_rate, const vector<double>& cpt_weight, const vector<double>& cpt_delete_extend);

  void set_transitions (double gamma, double delete_rate, const vector<double>& cpt_weight, const vector<double>& cpt_delete_extend);
  
  void sample_indel_params ();

  sstring indel_parameter_string() const;

  Convex_transducer_factory_param_container* propose_indel_params ();

};

#endif /* RIVAS_TRANSDUCER_FACTORY_INCLUDED */
