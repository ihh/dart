#ifndef QUICKCDSPAIR_INCLUDED
#define QUICKCDSPAIR_INCLUDED

#include "cdspair.h"

struct Quick_CDS_pair_state_enum
{
  enum {
    // 4-lexed states
    XYCOD2 = 0,   // type EmitXY * 256
    XYCOD3 = 256, // type EmitXY * 256
    // 2-lexed states
    XYCOD1 = 512, // type EmitXY * 16
    XCOD2 = 528,  // type EmitX * 16
    YCOD2 = 544,  // type EmitY * 16
    XCOD3 = 560,  // type EmitX * 16
    YCOD3 = 576,  // type EmitY * 16
    // 1-lexed states
    XCOD1 = 592,  // type EmitX * 4
    YCOD1 = 596,  // type EmitY * 4
    // 0-lexed states
    // Pad emit states
    XYMAT = 600,  // type EmitXY
    XDEL,         // type EmitX
    YINS,         // type EmitY
    // Null states
    PAD,          // type Null
    ENDPAD,       // type Null
    GENE,         // type Null
    ENDGENE,      // type Null
    XYCOD,        // type Null
    ENDXYCOD,     // type Null
    XCOD,         // type Null
    ENDXCOD,      // type Null
    YCOD,         // type Null
    ENDYCOD,      // type Null
    TotalStates  // = 613
  };

  // lexed state index helpers
  static inline int lex4xxyy (int state, const vector<int>& x, const vector<int>& y)
  { return state + x[0] + 4 * (x[1] + 4 * (y[0] + 4 * y[1])); }
  static inline int lex2xy (int state, const vector<int>& x, const vector<int>& y)
  { return state + x[0] + 4 * y[0]; }
  static inline int lex2xx (int state, const vector<int>& x)
  { return state + x[0] + 4 * x[1]; }
  static inline int lex1x (int state, const vector<int>& x)
  { return state + x[0]; }
};

struct Quick_CDS_pair : Quick_CDS_pair_state_enum, CDS_pair_params, Pair_PHMM
{
  // constructor
  Quick_CDS_pair (PScores& pscore);
};

#endif /* QUICKCDSPAIR_INCLUDED */
