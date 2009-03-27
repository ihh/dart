#include "scfg/lexcfg.h"

Lex_PCFG::Lex_PCFG (int us)
  : Pair_PCFG (us * CFG_lex_size),
    unlexed_states (us)
{ }

Lex_PCFG::Lex_PCFG (int us, const vector<vector<sstring> >& group_suffix)
  : Pair_PCFG (us * CFG_lex_size, group_suffix),
    unlexed_states (us)
{ }
