#ifndef PAIR_PCFG_INCLUDED
#define PAIR_PCFG_INCLUDED

#include "scfg/paircfgdp.h"
#include "hmm/singlephmm.h"
#include "hmm/pairphmm.h"
#include "seq/dirichlet.h"

class Pair_PCFG : public Pair_CFG<PFunc>
{
private:
  const vector<vector<sstring> >* group_suffix;   // for display only

public:
  Pair_PCFG (int states);
  Pair_PCFG (int states, const vector<vector<sstring> >& group_suffix);
  Pair_PCFG (const Pair_PHMM& hmm);
  Pair_PCFG (const Pair_CFG_scores& cfg);

  Pair_CFG_scores eval_cfg_scores (const PScores& var_scores) const;
  Pair_CFG_scores eval_cfg_scores_gs (const PScores& var_scores) const;   // same as eval_cfg_scores(), but sets the group_suffix

  void inc_var_counts (PCounts&               var_counts,
		       const PScores&         var_scores,
		       const Pair_CFG_counts& cfg_counts,
		       const Prob             model_count = 1.0) const;

  void set_group_suffix (const vector<vector<sstring> >& g_suffix);

  const char* element_descriptor() const { return "probability functions"; }
  int  element_width() const { return 16; }
  void show_element (const PFunc& element, ostream& o) const { element.show(o,group_suffix); }

  // method to factor in null model
  void factor_in_null_model (const Alphabet_group& null_emit);
  void factor_in_null_extend (const Boolean_group& null_extend);
};

// Pair_local_PCFG is an adapter class that makes local likelihood-ratio models out of global absolute-likelihood models
class Pair_local_PCFG : public Pair_PCFG
{
public:
  // the global model
  const Pair_PCFG& global;
  
  // special state indices:
  enum { LFlankBif = 0, RFlankBif = 1, PadBif = 2,
	 PreFlank = 3, PrePad = 4, PreBif = 5,
	 XLFlank = 6, XRFlank = 7, XLRFlank = 8,
	 YLFlank = 9, YRFlank = 10, YLRFlank = 11,
	 XLPad = 12, XRPad = 13, XLRPad = 14,
	 YLPad = 15, YRPad = 16, YLRPad = 17,
	 LocalStart = 18, GlobalOffset = 19 };

  // if S is a state in the old model and T is a state in this model, then local_state(S) = T and global_state(T) = S
  static inline int local_state (int global_state) { return GlobalOffset + global_state; }
  static inline int global_state (int local_state) { return local_state - GlobalOffset; }
  
  // Metascore index of mask (-1 if no mask)
  const int mask_metascore_idx;

  // update() method -- call this if any of the above PGroups are changed
  void update();

  // constructor
  Pair_local_PCFG (const Pair_PCFG& global,
		   PScope& pscope,
		   int mask_metascore_idx = -1);
  
  // method to add multiple-hits transition
  void allow_multiple_hits();

  // method to add LR pad states
  void allow_lr_pad();

  // masking methods
  void reset_mask (Named_profile& np) const;
  void reset_masks (Sequence_database& db) const;
  void mask (Named_profile& np, const vector<Metaprob>& expected_metacounts) const;
  void get_mask_gff (GFF_list& gff_list, const vector<Metaprob>& expected_metacounts, Prob min_prob, const char* seqname, const char* seq = 0, const char* source = "", const char* feature = "") const;
};

struct Odds_PCFG : Pair_PCFG
{
  // null model
  Alphabet_group null_emit;
  Boolean_group null_extend;
  // constructors
  Odds_PCFG (int states, Alphabet_group null_emit, Boolean_group null_extend);
  Odds_PCFG (int states, PScope& pscope);
  Odds_PCFG();
};

class Trainable_PCFG : public Odds_PCFG
{
private:
  // private stuff
  int first_pgroup_idx;  // set by constructor
  vector<PGroup> my_pgroups;
  void init_my_pgroups();  // call after creating PGroups
public:
  // public data
  PScores& pscore;
  // constructor
  Trainable_PCFG (PScores& pscore, int states);
  // virtual prior method
  virtual Dirichlet_prior default_prior() const;
};

#endif /* PAIR_PCFG_INCLUDED */
