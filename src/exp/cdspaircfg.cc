#include "stemloc/cdspaircfg.h"

// UNFINISHED

CDS_pair_CFG::CDS_pair_CFG (PScores& pscore) :
  CDS_pair_params (pscore),
  Odds_PCFG (TotalStates, nullEmit, nullExtend)
{
  // initialize emit states
  PFunc zeroFunc = 0.;

  init_emit (COD3, EmitXRYR, zeroFunc);
  init_emit (COD12, EmitXLRYLR, zeroFunc);
  init_emit (INT1COD1, EmitXLYL, zeroFunc);
  init_emit (COD23, EmitXLRYLR, zeroFunc);
  init_emit (INT2COD3, EmitXRYR, zeroFunc);

  init_emit (XCOD3, EmitXR, zeroFunc);
  init_emit (XCOD12, EmitXLR, zeroFunc);
  init_emit (XINT1COD1, EmitXL, zeroFunc);
  init_emit (XCOD23, EmitXLR, zeroFunc);
  init_emit (XINT2COD3, EmitXR, zeroFunc);

  init_emit (YCOD3, EmitYR, zeroFunc);
  init_emit (YCOD12, EmitYLR, zeroFunc);
  init_emit (YINT1COD1, EmitYL, zeroFunc);
  init_emit (YCOD23, EmitYLR, zeroFunc);
  init_emit (YINT2COD3, EmitYR, zeroFunc);

  init_emit (XYINT_GG, EmitXLRYLR, zeroFunc);
  init_emit (XYINT_TA, EmitXLRYLR, zeroFunc);
  init_emit (XYINTMAT, EmitXLYL, zeroFunc);
  init_emit (XYINTDEL, EmitXL, zeroFunc);
  init_emit (XYINTINS, EmitYL, zeroFunc);

  init_emit (XINT_GG, EmitXLR, zeroFunc);
  init_emit (XINT_TA, EmitXLR, zeroFunc);
  init_emit (XINT, EmitXL, zeroFunc);

  init_emit (YINT_GG, EmitYLR, zeroFunc);
  init_emit (YINT_TA, EmitYLR, zeroFunc);
  init_emit (YINT, EmitYL, zeroFunc);

  // initialize null states
  state_type[COD] = Null;
  state_type[INTCOD] = Null;
  state_type[INT1COD] = Null;
  state_type[INT2COD] = Null;

  state_type[XCOD] = Null;
  state_type[XINTCOD] = Null;
  state_type[XINT1COD] = Null;
  state_type[XINT2COD] = Null;

  state_type[YCOD] = Null;
  state_type[YINTCOD] = Null;
  state_type[YINT1COD] = Null;
  state_type[YINT2COD] = Null;

  state_type[INT] = Null;

  state_type[EXMAT] = Null;
  state_type[EXDEL] = Null;
  state_type[EXINS] = Null;

  state_type[CODON] = Null;
  state_type[XCODON] = Null;
  state_type[YCODON] = Null;

  // initialize bifurcation states
  init_bifurc (INT0COD, INT, COD);
  init_bifurc (XINT0COD, XINT, XCOD);
  init_bifurc (YINT0COD, YINT, YCOD);
  init_bifurc (CODEXMAT, CODON, EXMAT);
  init_bifurc (CODEXDEL, XCODON, EXDEL);
  init_bifurc (CODEXINS, YCODON, EXINS);

  // set up pair-lexed transitions, emissions & bifurcations
  for (int i = 0; i < CFG_alphabet_size; ++i)
    for (int j = 0; j < CFG_alphabet_size; ++j)
      {
	// lex-emit indices
	const int xlyl_idx = i*emit_xl_mul(EmitXLYL) + j*emit_yl_mul(EmitXLYL);
	const int xryr_idx = i*emit_xr_mul(EmitXRYR) + j*emit_yr_mul(EmitXRYR);

	// lex-emissions
	emit[COD3][xryr_idx] = 1;
	emit[INT1COD1][xlyl_idx] = 1;
	emit[INT2COD3][xryr_idx] = 1;

	// transitions
	const int idx = xlyl_idx;  // transition index

	transition (COD3 + idx, COD12 + idx) = 1;
	transition (COD12 + idx, End) = 1;
	transition (INT1COD1 + idx, INT1COD23 + idx) = 1;
	transition (COD23 + idx, End) = 1;
	transition (INT2COD3 + idx, INT2COD12 + idx) = 1;

	transition (COD, COD3 + idx) = 1;
	transition (INT1COD, INT1COD1 + idx) = 1;
	transition (INT2COD, INT2COD3 + idx) = 1;

	transition (XCOD, XCOD3 + idx) = 1;
	transition (XINT1COD, XINT1COD1 + idx) = 1;
	transition (XINT2COD, XINT2COD3 + idx) = 1;

	transition (YCOD, YCOD3 + idx) = 1;
	transition (YINT1COD, YINT1COD1 + idx) = 1;
	transition (YINT2COD, YINT2COD3 + idx) = 1;

	// bifurcations
	init_bifurc (INT1COD23 + idx, INT, COD23 + idx);
	init_bifurc (INT2COD12 + idx, COD12 + idx, INT);
      }

  // set up single-lexed transitions & emissions
  for (int i = 0; i < CFG_alphabet_size; ++i)
    {
      // lex-emissions






      // bifurcations
      init_bifurc (XINT1COD23 + i, XINT, XCOD23 + i);
      init_bifurc (XINT2COD12 + i, XCOD12 + i, XINT);
      init_bifurc (YINT1COD23 + i, YINT, YCOD23 + i);
      init_bifurc (YINT2COD12 + i, YCOD12 + i, YINT);
    }
  

  // further initialization of emit states
}

