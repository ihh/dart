#ifndef TRAINER_INCLUDED
#define TRAINER_INCLUDED

#include "hmm/singlephmm.h"
#include "hmm/singlefastdp.h"
#include "seq/dirichlet.h"
#include "seq/gff.h"

// A Trainable has a Single_PHMM plus scores, prior, null model & matrix factory
class Trainable : public HMM_state_enum, protected Stream_saver
{
private:
  vector<Alphabet_group> emit_group;  // PGroup's that take scaled null model pseudocounts for a prior
protected:
  Alphabet_group new_emit_group (const char* name = "");  // updates emit_group
public:
  // model, scores, prior
  Single_PHMM     hmm;
  PScores&        pscore;
  Dirichlet_prior prior;
  // index of mask metascore index (-1 if no mask)
  int mask_metascore_idx;
  // forward and reverse start states
  // by default, these are connected to Start & End with transition probability 1
  enum { FwdStart = 0, RevStart = 1, FwdEnd = 2, RevEnd = 3, StateOffset = 4 };
  // list of reverse states in HMM, for setting "strand" field in GFF results
  vector<int> reverse_states;
  // null model
  Alphabet_group  null_emit;
  Boolean_group   null_extend;
  // "seed path": path through model, used to seed model from data in MEME-like algorithms
  vector<int>     seed_path;
  // constructor
  Trainable (int states, const Alphabet& alphabet, PScores& pscore);
  // virtual destructor
  virtual ~Trainable();
  // virtual turtle-friendly output method
  virtual void instruct_turtle (ostream& turtle_stream) const;
  // helpers
  int seed_path_residues() const;  // counts the residues on the seed path
  Score_profile prof_sc() const;  // creates a Score_profile from emit_group
  // null model setup. user MUST call set_null_emit_prob() if emit priors are to be automagically assigned
  void set_null_emit_prob (const vector<Prob>& null_emit_prob, double pseudocount_multiplier);  // assigns emit priors
  void set_null_length (const double length);
  void optimise_null_model (const Sequence_database& db, double pseudocount_multiplier);  // trains the null model on a dataset
  // training methods
  void reset_to_prior();  // resets probabilities (scores) to prior mean
  void optimise_pscore (const Single_HMM_counts& hmm_counts);  // updates pscore from HMM counts
  void seed (const Named_profile& np);   // seed on a particular sequence
  void train (const Sequence_database& db,
	      Single_matrix_factory& dp_factory = Single_fast_matrix_factory::instance);  // repeatedly calls do_EM until convergence
  Loge do_EM (const Sequence_database& db,
	      Single_matrix_factory& dp_factory = Single_fast_matrix_factory::instance);  // does one round of EM using this model
};

// Local_trainer has a Trainable plus a Single_local_PHMM & local training/search algorithms
class Local_trainer
{
private:
  void optimise_pscore (const Single_HMM_counts& local_hmm_counts);
public:
  // data
  Trainable&             model;  // the global model
  Single_matrix_factory& dp_factory;
  Single_local_PHMM      local_hmm;
  Dirichlet_prior        prior;  // = model.prior
  // constructor
  // mask argument is passed to Single_local_PHMM
  Local_trainer (Trainable& model,
		 Single_matrix_factory& dp_factory = Single_fast_matrix_factory::instance);
  // data methods
  void global_train (const Sequence_database& db);  // wrapper for Trainable::train()
  void local_search (const Sequence_database& db, GFF_list& results, const char* source = "local", const char* feature = "match");  // wrapper for Single_local_model::search()
  void local_seed (const Named_profile& np, int pos);  // extracts subsequence from np, calls Trainable::seed()
  Loge local_EM (const Sequence_database& db);  // does one round of EM on full database using local model
  Loge local_EM (const Named_profile& np);  // does one round of EM on single sequence using local model
  Loge local_forward (const Named_profile& np);  // does Forward algorithm for local model on one sequence only
  void local_mask (Sequence_database& db) const;  // masks out hits to local model
  void get_local_mask (const Sequence_database& db, GFF_list& gff_list, Prob min_prob, const char* source = "", const char* feature = "") const;  // gets GFF hits
  // display methods
  void display_vars (ostream& o) const;
  sstring seqlogo() const;
};

#endif
