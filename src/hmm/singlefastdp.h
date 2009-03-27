#ifndef SINGLE_FAST_DP_INCLUDED
#define SINGLE_FAST_DP_INCLUDED

#include "hmm/singledp.h"

// Single HMM DP classes
// the following implementation works with Digitized_biosequence's and supports Metascores.
//
struct Single_fast_emit_calculator : Metascore_calculator
{
  const Named_profile&         np;
  const Digitized_biosequence& dsq;
  const Alphabet*              alphabet;   // optional, for display purposes only

  Single_fast_emit_calculator (const Single_HMM_scores& hmm, const Named_profile& np) :
    Metascore_calculator (hmm, np), np(np), dsq(np.dsq)
  {
    if (hmm.max_metascore_idx() >= (int) meta_sc.size()) THROW Standard_exception ("Can't resolve metascore");
    for_const_contents (vector<Metascore>, meta_sc, msc)
      if ((*msc).size() != dsq.size()) THROW Standard_exception ("Metascore size mismatch");
  }

  static const char* calc_type() { return "digitised"; }

  inline Score calc_emit_score (int state, int seqpos) const
    {
      Score sc = hmm.emit[state][dsq[seqpos]];
      add_meta (sc, state, seqpos);
      return sc;
    }

  inline Score accum_emit_score (Score posterior_sc, Single_HMM_counts& counts, vector<Metaprob>& metacounts, int state, int seqpos)
    {
      const Prob p = Score2Prob (posterior_sc);
      counts.emit[state][dsq[seqpos]] += p;
      count_meta (p, metacounts, state, seqpos);
      return calc_emit_score (state, seqpos);
    }

  char char_at_pos (int x) const
  {
    int i = dsq[x];
    if (alphabet) return alphabet->int2char (i);
    if (i < 10) return i + '0';
    return i-10 + 'a';
  }

  Score path_score (const vector<int>& state_path) const { Score_profile prof_sc (dsq); return hmm.path_score (state_path, np); }
};



typedef Single_Viterbi_matrix_template <Single_fast_emit_calculator> Single_fast_Viterbi_matrix;
typedef Single_forward_matrix_template <Single_fast_emit_calculator> Single_fast_forward_matrix;
typedef Single_forward_backward_matrix_template <Single_fast_emit_calculator> Single_fast_forward_backward_matrix;


// factory for the above DP classes
// converts Digitized_biosequence's into Score_profile's before passing them to DP routines
// hence this class is somewhat slow, and this is reflected in the name

class Single_fast_matrix_factory : public Single_matrix_factory
{
public:
  Single_fast_matrix_factory() { }
  Single_Viterbi_interface* new_Viterbi (const Single_HMM_scores& hmm, const Named_profile& np);
  Single_forward_interface* new_forward (const Single_HMM_scores& hmm, const Named_profile& np);
  Single_forward_backward_interface* new_forward_backward (const Single_HMM_scores& hmm, const Named_profile& np);

  static Single_fast_matrix_factory instance;  // singleton
};

#endif
