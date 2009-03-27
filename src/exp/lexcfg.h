#ifndef LEXCFG_INCLUDED
#define LEXCFG_INCLUDED

#include "scfg/pairpcfg.h"

class Lex_PCFG : public Pair_PCFG
{
private:
  inline PFunc& transition (int src, int dest);
  inline const PFunc& transition (int src, int dest) const;
  inline void verify_ltrans (int src, int dest) const; // throws an exception if transition is illegal
public:
  const int unlexed_states;
  Lex_PCFG (int unlexed_states);
  Lex_PCFG (int unlexed_states, const vector<vector<sstring> >& group_suffix);
  void set_unlexed_state_type (int unlexed_state, State_type type);  // sets delta emit profiles for lexed states
  // state indices
  inline int lstate (int unlexed_state, int xl, int xr, int yl, int yr);  // context indices (xl, xr etc) are nondegenerate
  inline int lstate (int unlexed_state, const char* context);
  // transition accessors
  inline PFunc& ltrans (int lexed_src, int lexed_dest);  // calls verify_ltrans()
  inline PFunc& ltrans (int lexed_src, int lexed_dest, const char* src_context, const char* dest_emit);
  inline const PFunc& ltrans (int lexed_src, int lexed_dest) const;  // calls verify_ltrans()
  inline const PFunc& ltrans (int lexed_src, int lexed_dest, const char* src_context, const char* dest_emit) const;
};

// private transition() method wrappers
PFunc& Lex_PCFG::transition (int src, int dest)
{ return Pair_CFG<PFunc>::transition (src, dest); }

const PFunc& Lex_PCFG::transition (int src, int dest) const
{ return Pair_CFG<PFunc>::transition (src, dest); }

// lexed state indices
int Lex_PCFG::lstate (int unlexed_state, int xl, int xr, int yl, int yr)
{ return unlexed_states * emit_xylr_offset (xl, xr, yl, yr) + (unlexed_state % unlexed_states); }

#endif /* LEXCFG_INCLUDED */
