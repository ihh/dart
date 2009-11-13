#ifndef XRATE_PRIMITIVES_INCLUDED
#define XRATE_PRIMITIVES_INCLUDED

#include "guile/guile-defs.h"
#include "ecfg/ecfg.h"

// ECFG
SCM ecfg_to_scm (const ECFG_scores& ecfg, const ECFG_counts* counts = 0);

// xrate functions
void init_xrate_primitives (void);

#endif /* XRATE_PRIMITIVES_INCLUDED */
