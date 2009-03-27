#ifndef TELEGRAPH_COMPILER_INCLUDED
#define TELEGRAPH_COMPILER_INCLUDED

#include <math.h>

/* 80-column formatting
.........1.........2.........3.........4.........5.........6.........7.........
1234567890123456789012345678901234567890123456789012345678901234567890123456789
*/

#ifdef __cplusplus 
extern "C" {
#endif /* defined __cplusplus */

  /* Implemented */
#include "tgvector.h"
#include "tgtypetable.h"
#include "tgscore.h"
#include "tgindices.h"

#include "tgsymbol.h"
#include "tgprob.h"
#include "tgexpression.h"
#include "tgrule.h"
#include "tgparsetree.h"
#include "tggrammar.h"

#include "tgsequence.h"
#include "tgbuildgrammar.h"
#include "tgoutputgrammar.h"
#include "tgprobscore.h"
#include "list_assignment.h"

#include "stkutil.h"

  /* Compile method */
  void tg_compile (tg_grammar* tg);
  
#ifdef __cplusplus 
}
#endif /* defined __cplusplus */
#endif  /* defined TELEGRAPH_COMPILER_INCLUDED */
