#ifndef CFGDOTPLOT_INCLUDED
#define CFGDOTPLOT_INCLUDED

#include "util/array2d.h"
#include "amap/dotplot.h"
#include "scfg/paircfgdp.h"

/// Fold dotplot
struct PairCFG_fold_dotplot : Dotplot, Pair_CFG_state_type_enum {

  // seq_index = 0 for x-sequence, 1 for y-sequence
  PairCFG_fold_dotplot (const Pair_inside_outside_matrix& in_out, int seq_index);

};

/// Alignment dotplot
struct PairCFG_alignment_dotplot : Dotplot, Pair_CFG_state_type_enum {

  /// constructor
  /*
   * If sanitize_probs is true, then probs are rounded down (to 0 or 1) as appropriate.
   */
  PairCFG_alignment_dotplot (const Pair_inside_outside_matrix& in_out);

};

#endif /* CFGDOTPLOT_INCLUDED */
