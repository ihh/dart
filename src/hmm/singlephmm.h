#ifndef SINGLE_PHMM_INCLUDED
#define SINGLE_PHMM_INCLUDED

#include "hmm/singletmpl.h"
#include "hmm/singlefastdp.h"
#include "seq/pfunc.h"
#include "seq/gff.h"

class Single_PHMM : public Single_meta_HMM<PFunc>
{
private:
  const PScope* pscope;  // for Telegraph calibration
  const vector<vector<sstring> >* group_suffix;   // for display only

public:
  Single_PHMM (int states, const Alphabet& alphabet);
  Single_PHMM (int states, const Alphabet& alphabet, const vector<vector<sstring> >& group_suffix);

  Single_HMM_scores eval_hmm_scores (const PScores& var_scores) const;
  Single_HMM_scores eval_hmm_scores_gs (const PScores& var_scores) const;   // same as eval_hmm_scores(), but sets the group_suffix

  void inc_var_counts (PCounts&                 var_counts,
		       const PScores&           var_scores,
		       const Single_HMM_counts& hmm_counts,
		       const Prob               model_count = 1.0) const;

  void set_group_suffix (const vector<vector<sstring> >& g_suffix);
  void set_pscope (const PScope& pscope);

  const char* element_descriptor() const { return "probability functions"; }
  int  element_width() const { return 16; }
  void show_element (const PFunc& element, ostream& o) const { element.show(o,group_suffix); }
};

// Single_local_PHMM is an adapter class that makes local likelihood-ratio models out of global absolute-likelihood models
class Single_local_PHMM : public Single_PHMM
{
public:
  // the global model
  const Single_PHMM& global;

  // special state indices:
  enum { LeftPad = 0, LocalStart = 1, LocalEnd = 2, RightPad = 3, GlobalOffset = 4 };

  // if S is a state in the old model and T is a state in this model, then local_state(S) = T and global_state(T) = S
  static inline int local_state (int global_state) { return GlobalOffset + global_state; }
  static inline int global_state (int local_state) { return local_state - GlobalOffset; }
  
  // null model PGroups
  Alphabet_group null_emit;
  Boolean_group  null_extend;

  // Metascore index of mask (-1 if no mask)
  const int mask_metascore_idx;

  // flag to factor in null model
  const bool factor_in_null_model;

  // update() method -- call this if any of the above PGroups are changed
  void update();

  // constructor
  Single_local_PHMM (const Single_PHMM& global,
		     const Alphabet_group& null_emit,
		     const Boolean_group& null_extend,
		     PScope& pscope,
		     int mask_metascore_idx = -1,
		     bool factor_in_null_model = 0);
  
  // method to add multiple-hits transition
  void allow_multiple_hits();

  // database search method
  void search (const PScores& pscores,
	       const Sequence_database& db,
	       GFF_list& results,
	       const vector<int>& reverse_states,
	       Single_matrix_factory& matrix_factory = Single_fast_matrix_factory::instance,
	       const char* source = "local",
	       const char* feature = "match",
	       const char* pattern = "") const;

  // masking methods
  void reset_mask (Named_profile& np) const;
  void reset_masks (Sequence_database& db) const;
  void mask (Named_profile& np, const vector<Metaprob>& expected_metacounts) const;
  void get_mask_gff (GFF_list& gff_list, const vector<Metaprob>& expected_metacounts, Prob min_prob, const char* seqname, const char* seq = 0, const char* source = "", const char* feature = "") const;
};

typedef Transition_matrix<PFunc> Transition_funcs;

class Transition_funcs_concrete : public Transition_funcs
{
public:
  Transition_funcs_concrete (int size) : Transition_funcs (size, PFunc(0.0)) { }
  // inherited virtuals
  const char* element_descriptor() const { return "probability functions"; }
  int  element_width() const { return 16; }
  void show_element (const PFunc& element, ostream& o) const { element.show(o); }
};

#endif
