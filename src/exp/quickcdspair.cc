#include "stemloc/quickcdspair.h"

Quick_CDS_pair::Quick_CDS_pair (PScores& pscore) :
  CDS_pair_params (pscore),
  Pair_PHMM (TotalStates, CFG_alphabet)
{
  // some useful variables
  const PFunc zero = 0.;
  const PFunc one = 1.;

  // initialize 0-lexed states
  init_emit (XYMAT, EmitXY, zero);
  init_emit (XDEL, EmitX, one);
  init_emit (YINS, EmitY, one);

  state_type[PAD] = Null;
  state_type[ENDPAD] = Null;
  state_type[XYCOD] = Null;
  state_type[ENDXYCOD] = Null;
  state_type[XCOD] = Null;
  state_type[ENDXCOD] = Null;
  state_type[YCOD] = Null;
  state_type[ENDYCOD] = Null;

  // transitions for 0-lexed states
  transition (Start, PAD) = padOpen.y;
  transition (Start, GENE) = padOpen.n * geneExtend.y;
  transition (Start, End) = padOpen.n * geneExtend.n;

  transition (PAD, YINS) = padGapOpen.y;
  transition (PAD, XDEL) = padGapOpen.n * padGapOpen.y;
  transition (PAD, XYMAT) = padGapOpen.n * padGapOpen.n;

  transition (XYMAT, YINS) = padGapOpen.y;
  transition (XYMAT, XDEL) = padGapOpen.n * padGapOpen.y;
  transition (XYMAT, XYMAT) = padGapOpen.n * padGapOpen.n * padExtend.y;
  transition (XYMAT, ENDPAD) = padGapOpen.n * padGapOpen.n * padExtend.n;

  transition (YINS, YINS) = padGapExtend.y;
  transition (YINS, XDEL) = padGapExtend.n * padGapOpen.y;
  transition (YINS, XYMAT) = padGapExtend.n * padGapOpen.n * padExtend.y;
  transition (YINS, ENDPAD) = padGapExtend.n * padGapOpen.n * padExtend.n;

  transition (XDEL, XDEL) = padGapExtend.y;
  transition (XDEL, XYMAT) = padGapExtend.n * padExtend.y;
  transition (XDEL, ENDPAD) = padGapExtend.n * padExtend.n;

  transition (ENDPAD, GENE) = geneExtend.y;
  transition (ENDPAD, End) = geneExtend.n;

  transition (GENE, XCOD) = codGapOpen.y;
  transition (GENE, YCOD) = codGapOpen.n * codGapOpen.y;
  transition (GENE, XYCOD) = codGapOpen.n * codGapOpen.n;

  transition (ENDXCOD, XCOD) = codGapExtend.y;
  transition (ENDXCOD, YCOD) = codGapExtend.n * codGapOpen.y;
  transition (ENDXCOD, XYCOD) = codGapExtend.n * codGapOpen.n * geneExtend.y;
  transition (ENDXCOD, ENDGENE) = codGapExtend.n * codGapOpen.n * geneExtend.n;

  transition (ENDYCOD, YCOD) = codGapExtend.y;
  transition (ENDYCOD, XYCOD) = codGapExtend.n * geneExtend.y;
  transition (ENDYCOD, ENDGENE) = codGapExtend.n * geneExtend.n;

  transition (ENDXYCOD, XCOD) = codGapOpen.y;
  transition (ENDXYCOD, YCOD) = codGapOpen.n * codGapOpen.y;
  transition (ENDXYCOD, XYCOD) = codGapOpen.n * codGapOpen.n * geneExtend.y;
  transition (ENDXYCOD, ENDGENE) = codGapOpen.n * codGapOpen.n * geneExtend.n;

  transition (ENDGENE, PAD) = padOpen.y;
  transition (ENDGENE, End) = padOpen.n;

  // emissions for 0-lexed XYMAT state
  for (int x = 0; x < CFG_alphabet_size; ++x)
    for (int y = 0; y < CFG_alphabet_size; ++y)
      pair_emit[XYMAT](x,y) = nucSub[x][y] / nullEmit[y];

  // initialize lexed states
  vector<int> xcod (3);
  vector<int> ycod (3);
  for (xcod[0] = 0; xcod[0] < CFG_alphabet_size; ++xcod[0])
    {
      init_emit (lex1x(XCOD1,xcod), EmitX, zero);
      init_emit (lex1x(YCOD1,xcod), EmitY, zero);

      single_emit[lex1x(XCOD1,xcod)][xcod[0]] = one;
      single_emit[lex1x(YCOD1,xcod)][xcod[0]] = one;

      transition (XCOD, lex1x(XCOD1,xcod)) = one;
      transition (YCOD, lex1x(YCOD1,xcod)) = one;

      for (xcod[1] = 0; xcod[1] < CFG_alphabet_size; ++xcod[1])
	{
	  init_emit (lex2xx(XCOD2,xcod), EmitX, zero);
	  init_emit (lex2xx(YCOD2,xcod), EmitY, zero);
	  init_emit (lex2xx(XCOD3,xcod), EmitX, zero);
	  init_emit (lex2xx(YCOD3,xcod), EmitY, zero);

	  single_emit[lex2xx(XCOD2,xcod)][xcod[1]] = one;
	  single_emit[lex2xx(YCOD2,xcod)][xcod[1]] = one;

	  transition (lex1x(XCOD1,xcod), lex2xx(XCOD2,xcod)) = one;
	  transition (lex2xx(XCOD2,xcod), lex2xx(XCOD3,xcod)) = one;
	  transition (lex2xx(XCOD3,xcod), ENDXCOD) = one;

	  transition (lex1x(YCOD1,xcod), lex2xx(YCOD2,xcod)) = one;
	  transition (lex2xx(YCOD2,xcod), lex2xx(YCOD3,xcod)) = one;
	  transition (lex2xx(YCOD3,xcod), ENDYCOD) = one;

	  for (xcod[2] = 0; xcod[2] < CFG_alphabet_size; ++xcod[2])
	    {
	      const int xc = codEmit.intvec2index (xcod);

	      PFunc xnull = one;
	      for (int i = 0; i < 3; ++i)
		xnull *= nullEmit[xcod[i]];

	      single_emit[lex2xx(XCOD3,xcod)][xcod[2]] = codEmit[xc] / xnull;
	      single_emit[lex2xx(YCOD3,xcod)][xcod[2]] = codEmit[xc] / xnull;
	    }
	}

      for (ycod[0] = 0; ycod[0] < CFG_alphabet_size; ++ycod[0])
	{
	  init_emit (lex2xy(XYCOD1,xcod,ycod), EmitXY, zero);

	  pair_emit[lex2xy(XYCOD1,xcod,ycod)](xcod[0],ycod[0]) = one;

	  transition (XYCOD, lex2xy(XYCOD1,xcod,ycod)) = one;

	  for (xcod[1] = 0; xcod[1] < CFG_alphabet_size; ++xcod[1])
	    for (ycod[1] = 0; ycod[1] < CFG_alphabet_size; ++ycod[1])
	      {
		init_emit (lex4xxyy(XYCOD2,xcod,ycod), EmitXY, zero);
		init_emit (lex4xxyy(XYCOD3,xcod,ycod), EmitXY, zero);

		pair_emit[lex4xxyy(XYCOD2,xcod,ycod)](xcod[1],ycod[1]) = one;

		transition (lex2xy(XYCOD1,xcod,ycod), lex4xxyy(XYCOD2,xcod,ycod)) = one;
		transition (lex4xxyy(XYCOD2,xcod,ycod), lex4xxyy(XYCOD3,xcod,ycod)) = one;
		transition (lex4xxyy(XYCOD3,xcod,ycod), ENDXYCOD) = one;

		for (xcod[2] = 0; xcod[2] < CFG_alphabet_size; ++xcod[2])
		  {
		    const int xc = codEmit.intvec2index (xcod);

		    PFunc xnull = one;
		    for (int i = 0; i < 3; ++i)
		      xnull *= nullEmit[xcod[i]];

		    for (ycod[2] = 0; ycod[2] < CFG_alphabet_size; ++ycod[2])
		      {
			const int yc = codEmit.intvec2index (ycod);

			PFunc xynull = xnull;
			for (int i = 0; i < 3; ++i)
			  xynull *= nullEmit[ycod[i]];

			pair_emit[lex4xxyy(XYCOD3,xcod,ycod)](xcod[2],ycod[2]) = codEmit[xc] * codSub[xc][yc] / xynull;
		      }
		  }
	      }
	}
    }
}
