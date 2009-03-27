#ifndef WEHITS_INCLUDED
#define WEHITS_INCLUDED

#include "scfg/paircfgdp.h"

// Waterman-Eggert-type single-SCFG hits-getter

struct WE_hits : GFF_list, Fold_char_enum
{
  void get_hits (const Pair_CFG_scores& cfg, const Named_profile& np, bool local, int max_subseq_len, Score min_score, int max_hits = -1);
  void get_hits (const Pair_CFG_scores& cfg, const set<int>& paired_states, const Named_profile& np, bool local, int max_subseq_len, Score min_score, int max_hits = -1);
};

#endif /* WEHITS_INCLUDED */
