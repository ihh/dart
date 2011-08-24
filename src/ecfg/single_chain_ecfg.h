#ifndef SINGLE_CHAIN_ECFG_INCLUDED
#define SINGLE_CHAIN_ECFG_INCLUDED

#include "ecfg/ecfg.h"

struct Single_chain_ECFG : ECFG_scores
{
  Single_chain_ECFG (const Irrev_EM_matrix& matrix);
};


#endif /* SINGLE_CHAIN_ECFG_INCLUDED */

