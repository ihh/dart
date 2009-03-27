#include "empath/bindpaircfg.h"

Bind_pair_CFG::Bind_pair_CFG (const RNA_energy& energy)
  : Pair_CFG_scores (TotalStates),
    Duplex_lookup()
{
  // initialise padding states
  init_emit (xPad5, EmitXL, 0);
  init_emit (yPad3, EmitYR, 0);
  init_emit (xPad3, EmitXL, 0);
  init_emit (yPad5, EmitYR, 0);
  // initialise duplex emit states
  for (int d = 0; d < Duplexes; ++d)
    {
      const int xsym = duplex_sym5[d];
      const int ysym = duplex_sym3[d];
      const int emit_idx = xsym + Symbols * ysym;
      init_emit (xyDuplex + d, EmitXLYR, -InfinityScore);
      emit[xyDuplex + d][emit_idx] = 0;
    }
  // initialise interior loop terminal states
  for (int d = 0; d < Duplexes; ++d)
    {
      const int xsym = duplex_sym5[d];
      const int ysym = duplex_sym3[d];
      init_emit (xyIntBeg + d, EmitXLYR, -InfinityScore);
      init_emit (xyIntEnd + d, EmitXLYR, -InfinityScore);
      for (int xsym_term = 0; xsym_term < Symbols; ++xsym_term)
	for (int ysym_term = 0; ysym_term < Symbols; ++ysym_term)
	  {
	    const int emit_idx_term = xsym_term + Symbols * ysym_term;
	    emit[xyIntBeg + d][emit_idx_term] = energy.interior (xsym, ysym) (xsym_term, ysym_term);
	    emit[xyIntEnd + d][emit_idx_term] = energy.interior (xsym_term, ysym_term) (xsym, ysym);
	  }
    }
  // initialise loop emit state
  init_emit (xyIntLoop, EmitXLYR, 0);
  // duplex state transitions
  for (int d = 0; d < Duplexes; ++d)
    {
      const int xsym = duplex_sym5[d];
      const int ysym = duplex_sym3[d];
      // new duplex
      transition (Start, xyDuplex + d) = 0;
      transition (xPad5, xyDuplex + d) = 0;
      transition (yPad3, xyDuplex + d) = 0;
      // stacking
      for (int d_prev = 0; d_prev < Duplexes; ++d_prev)
	{
	  const int xsym_prev = duplex_sym5[d_prev];
	  const int ysym_prev = duplex_sym3[d_prev];
	  transition (xyDuplex + d_prev, xyDuplex + d) = energy.stacked (xsym_prev, ysym_prev) (xsym, ysym);
	}
      // end of duplex
      transition (xyDuplex + d, xPad3) = 0;
      transition (xyDuplex + d, yPad5) = 0;
      transition (xyDuplex + d, End) = 0;
      // interior loop transitions
      transition (xyDuplex + d, xyIntBeg + d) = 0;
      for (int d_prev = 0; d_prev < Duplexes; ++d_prev)
	transition (xyIntBeg + d_prev, xyIntEnd + d) = energy.interior_loop[2];
      transition (xyIntEnd + d, xyDuplex + d) = 0;
      //      transition (xyIntBeg + d, xyIntLoop) = affine_approximation_to_interior_loop ...
    }
  // padding state transitions
  transition (Start, xPad5) = 0;
  transition (Start, yPad3) = 0;

  transition (xPad5, xPad5) = 0;
  transition (xPad5, yPad3) = 0;
  transition (xPad5, End) = 0;

  transition (yPad3, yPad3) = 0;
  transition (yPad3, End) = 0;

  transition (xPad3, xPad3) = 0;
  transition (xPad3, yPad5) = 0;
  transition (xPad3, End) = 0;

  transition (yPad5, yPad5) = 0;
  transition (yPad5, End) = 0;
  // flag mask states
  for (int s = xyDuplex; s < TotalStates; ++s) mask_states.insert (s);
}
