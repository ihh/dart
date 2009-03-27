#ifndef SINGLE_DP_INCLUDED
#define SINGLE_DP_INCLUDED

#include "hmm/singledptmpl.h"

// Metascore calculator
struct Metascore_calculator
{
  const Single_HMM_scores& hmm;
  const vector<Metascore>& meta_sc;

  Metascore_calculator (const Single_HMM_scores& hmm, const Meta_profile& mp) :
    hmm (hmm), meta_sc (mp.meta_sc)
  { }

  inline void add_meta (Score& sc, int state, int seqpos) const
  {
    for_const_contents (vector<int>, hmm.metascore_idx[state], metascore_idx)
      ScorePMulAcc (sc, meta_sc[*metascore_idx][seqpos]);
  }
  
  inline void count_meta (const Prob posterior_prob, vector<Metaprob>& metacounts, int state, int seqpos) const
  {
    for_const_contents (vector<int>, hmm.metascore_idx[state], m_idx)
      metacounts[*m_idx][seqpos] += posterior_prob;
  }
};

// Single HMM DP classes
// the following implementation works with ungapped profiles (Score_profile's) and so may be slow.
//
struct Single_emit_calculator : Metascore_calculator
{
  const Meta_profile&      mp;
  const Score_profile&     seq;
  const Alphabet*          alphabet;   // optional, for display purposes only

  Single_emit_calculator (const Single_HMM_scores& hmm, const Meta_profile& mp) :
    Metascore_calculator (hmm, mp), mp (mp), seq (mp.prof_sc), alphabet (&hmm.alphabet())
  { }
  
  static const char* calc_type() { return "profile"; }

  inline Score calc_emit_without_meta (int state, const Symbol_score_map& ssmap) const
  {
    const vector<Score>& emit = hmm.emit[state];
    
    Score sc = -InfinityScore;
    for_const_contents (Symbol_score_map, ssmap, ss)
      ScorePSumAcc (sc, ScorePMul (emit[(*ss).first], (*ss).second));

    return sc;
  }

  inline Score calc_emit_score (int state, int seqpos) const
    {
      Score sc = calc_emit_without_meta (state, seq[seqpos]);
      add_meta (sc, state, seqpos);
      return sc;
    }

  inline Score accum_emit_score (Score posterior_sc, Single_HMM_counts& counts, vector<Metaprob>& metacounts, int state, int seqpos)
    {
      const Symbol_score_map& ssmap          =  seq[seqpos];
      const vector<Score>&    emit           =  hmm.emit[state];
      vector<Prob>&           emit_counts    =  counts.emit[state];

      const Score emit_without_meta = calc_emit_without_meta (state, ssmap);
      const Score posterior_minus_emit_sc = posterior_sc - emit_without_meta;

      for_const_contents (Symbol_score_map, ssmap, ss)
	{
	  const Prob p = Score2Prob (ScorePMul3 ((*ss).second, emit[(*ss).first], posterior_minus_emit_sc));
	  emit_counts[(*ss).first] += p;
	  count_meta (p, metacounts, state, seqpos);
	}
      
      Score emit_sc = emit_without_meta;
      add_meta (emit_sc, state, seqpos);
      return emit_sc;
    }

  char char_at_pos (int x) const
    {
      const Symbol_score_map& m = seq[x];
      if (alphabet) return alphabet->score2char (m);
      int i = Score_profile::score_map_consensus(m);
      if (i < 10) return i + '0';
      return i-10 + 'a';
    }
  
  Score path_score (const vector<int>& state_path) const { return hmm.path_score (state_path, mp); }
};



typedef Single_Viterbi_matrix_template <Single_emit_calculator> Single_Viterbi_matrix;
typedef Single_forward_matrix_template <Single_emit_calculator> Single_forward_matrix;
typedef Single_forward_backward_matrix_template <Single_emit_calculator> Single_forward_backward_matrix;


// factory for the above DP classes
// converts Digitized_biosequence's into Score_profile's before passing them to DP routines
// hence this class is somewhat slow, and this is reflected in the name

class Single_profile_matrix_factory : public Single_matrix_factory
{
public:
  Single_profile_matrix_factory() { }
  Single_Viterbi_interface* new_Viterbi (const Single_HMM_scores& hmm, const Named_profile& np);
  Single_forward_interface* new_forward (const Single_HMM_scores& hmm, const Named_profile& np);
  Single_forward_backward_interface* new_forward_backward (const Single_HMM_scores& hmm, const Named_profile& np);

  static Single_profile_matrix_factory instance;  // singleton
};

#endif
