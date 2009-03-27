#ifndef MULTIHIT_INCLUDED
#define MULTIHIT_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include "util/piper.h"
#include "empath/trainer.h"
#include "hmm/singlefastdp.h"
#include "seq/distmat.h"

// abstract model for sequences plus whatever other data (expression, GO ids...)
struct Hitter
{
  virtual void seed (const Named_profile& np, int pos) = 0;
  virtual void initialise() = 0;  // resets counts, updates params
  virtual Loge calc_loglike (const Named_profile& np) = 0;
  virtual void optimise() = 0;
  virtual vector<double> get_params() = 0;
  virtual void set_params (const vector<double>& params) = 0;
  // parallelisable update read/write methods
  virtual void send_update (ostream& out, const Named_profile& np, Prob weight) = 0;
  virtual void receive_update (istream& in) = 0;
  // display method
  virtual void display (ostream& out) = 0;
  // virtual destructor, to keep gcc happy
  virtual ~Hitter() { }
};

// concrete Hitter implementation for HMMs
struct Sequence_hitter : Hitter
{
  // data
  Local_trainer& trainer;
  const Single_PHMM& hmm;
  Dirichlet_prior& prior;
  PScores& pscore;
  Single_matrix_factory& dp_factory;
  Single_HMM_scores hmm_scores;
  PCounts_like pcount;
  // constructor
  Sequence_hitter (Local_trainer& trainer);
  // virtual methods
  void seed (const Named_profile& np, int pos);
  void initialise();
  Loge calc_loglike (const Named_profile& np);
  void optimise();
  vector<double> get_params();
  void set_params (const vector<double>& params);
  // funky parallel methods
  void send_update (ostream& out, const Named_profile& np, Prob weight);
  void receive_update (istream& in);
  // display method
  void display (ostream& out);
};

// adapter for Sequence_hitter using Local_trainer constructor
struct Local_hitter : Sequence_hitter
{
  Local_hitter (Local_trainer& trainer);
};

// class for multiple hits of a Hitter to a database
class Multihitter : public HMM_state_enum, private Piper
{
private:
  // db_hmm is a pseudo-HMM, containing a "hit" state and a "miss" state for each state in the supplied transition matrix.
  // The hit states in the pseudo-HMM are connected by the supplied transition matrix;
  // thus, this matrix specifies the distribution of hits.
  // The miss states can be entered at zero cost following each hit state; the self-loop cost for these states is also zero.
  // This is because each "residue" in the pseudo-model is actually a whole sequence in the database,
  // and sequences without a hit can occur anywhere.
  // For the same reason, all the hit states are effectively identical.
  Single_PHMM db_hmm;  // pseudo HMM

  inline static int hit_state (int s) { return s * 2 + 1; }
  inline static int miss_state (int s) { return s * 2 + 2; }
  inline static int start_miss_state() { return 0; }

  const int hit_meta_idx;
  const int dummy_sym;

  // posterior hit threshold (divided by size of database)
  double threshold;

  // method to return dummy Named_profile with hit log-likelihoods as metascores
  Named_profile make_hit_profile (const vector<Named_profile*>& db);

public:
  // data
  Hitter& hitter;
  PScores& pscore;  // for Transition_funcs
  Dirichlet_prior prior;  // contains prior for Transition_funcs *only* (not model prior)

  // constructor
  Multihitter (const Transition_funcs& trans, Hitter& hitter, PScores& pscore, int max_fork, double threshold = 1.0);

  // methods
  Loge multihit_forward (const vector<Named_profile*>& db);  // does forward algorithm on single sequence
  Loge multihit_EM (const vector<Named_profile*>& db);  // does one round of EM using multihit model
  void iterate_multihit_EM (const vector<Named_profile*>& db, double min_inc, int forgive);  // does "forgiving" EM until convergence
};

// Range_transitions is a set of constant Transition_funcs for a simple range
class Range_transitions : public Transition_funcs_concrete
{
public:
  Range_transitions (int range_min, int range_max);
};

#endif
