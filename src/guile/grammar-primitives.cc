#include "guile/grammar-primitives.h"

#include "guile/guile-keywords.h"
#include "guile/guile-defs.h"
#include "guile/stockholm-type.h"
#include "guile/newick-type.h"
#include "guile/xrate-primitives.h"

#include "ecfg/ecfg.h"

static void*
register_grammar_functions (void* data)
{
  ECFG_Scheme_evaluator& evaluator = *((ECFG_Scheme_evaluator*) data);

  init_stockholm_type();
  init_newick_type();
  init_xrate_primitives();

  if (evaluator.stock)
    scm_c_define (GUILE_ALIGNMENT, make_stockholm_smob (*evaluator.stock));

  return NULL;
}

ECFG_Scheme_evaluator::ECFG_Scheme_evaluator (const Stockholm* stock)
  : stock(stock)
{
  register_functions = &register_grammar_functions;
  data = (void*) this;
}
