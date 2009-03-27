#ifndef PAIR_PHMM_INCLUDED
#define PAIR_PHMM_INCLUDED

#include "hmm/pairhmm.h"
#include "seq/pfunc.h"

class Pair_PHMM : public Pair_HMM<PFunc>
{
private:
  const PScope* pscope;  // for Telegraph calibration
  const vector<vector<sstring> >* group_suffix;   // for display only

public:
  Pair_PHMM (Pair_HMM<PFunc>& hmm);
  Pair_PHMM (int states, const Alphabet& alphabet);
  Pair_PHMM (int states, const Alphabet& alphabet, const vector<vector<sstring> >& group_suffix);

  Pair_HMM_scores eval_hmm_scores (const PScores& var_scores) const;
  Pair_HMM_scores eval_hmm_scores_gs (const PScores& var_scores) const;   // same as eval_hmm_scores(), but sets the group_suffix

  void inc_var_counts (PCounts&               var_counts,
		       const PScores&         var_scores,
		       const Pair_HMM_counts& hmm_counts,
		       const Prob             model_count = 1.0) const;

  void set_group_suffix (const vector<vector<sstring> >& g_suffix);
  void set_pscope (const PScope& pscope);

  const char* element_descriptor() const { return "probability functions"; }
  int  element_width() const { return 16; }
  void show_element (const PFunc& element, ostream& o) const { element.show(o,group_suffix); }

  // method to factor in null model
  void factor_in_null_model (const Alphabet_group& null_emit);
};

#endif /* PAIR_PHMM_INCLUDED */

