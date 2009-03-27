#include "stemloc/quickstem.h"

Quick_stem::Quick_stem (PScores& pscore, const Alphabet_group& null_emit)
  : Pair_PCFG (Loop + loop_states),
    pscore (pscore),
    pad_extend_pg (pscore.new_boolean_group ("PadExtend")),
    stem_prior_pg (pscore.new_boolean_group ("StemPrior")),
    stem_extend_pg (pscore.new_boolean_group ("StemExtend")),
    stem_bifurc_pg (pscore.new_boolean_group ("StemBifurc")),
    sym_bulge_open_pg (pscore.new_boolean_group ("SymBulgeOpen")),
    sym_bulge_extend_pg (pscore.new_boolean_group ("SymBulgeExtend")),
    asym_bulge_open_pg (pscore.new_boolean_group ("AsymBulgeOpen")),
    asym_bulge_extend_pg (pscore.new_boolean_group ("AsymBulgeExtend")),
    loop_extend_pg (pscore.new_boolean_group ("LoopExtend")),
    loop_start_pos (pscore.new_group (loop_states, "LoopStartPos")),
    null_emit (null_emit),
    dinuc (pscore.new_alphabet_group (CFG_alphabet, 2, "emit"))
{
  // set state types
  const PFunc one (1.0);
  init_emit (PrePad,  Null,    one);
  init_emit (PadL,    EmitXL,  one);
  init_emit (StemBif, Bifurc,  one);
  init_emit (PreStem, Null,    one);
  init_emit (StemLR,  EmitXLR, one);
  init_emit (SymL,    EmitXL,  one);
  init_emit (SymR,    EmitXR,  one);
  init_emit (AsymL,   EmitXL,  one);
  init_emit (AsymR,   EmitXR,  one);
  for (int i = 0; i < loop_states; ++i)
    init_emit (loop(i), EmitXL, one);

  // set transitions
  transition (Start,   PrePad)  = 1;

  transition (PrePad,  StemBif) = stem_prior_pg.YES;
  transition (PrePad,  PadL)    = stem_prior_pg.NO * pad_extend_pg.YES;
  transition (PrePad,  End)     = stem_prior_pg.NO * pad_extend_pg.NO;

  // the following transition allows stems flush to 3' end, which otherwise would be lost since they would generate a zero-length subseq... 
  // this behaviour is slightly unsatisfactory...
  transition (PrePad,  PreStem) = stem_prior_pg.YES * pad_extend_pg.NO;

  transition (PadL,    PrePad)  = 1;

  transition (PreStem, StemLR)  = 1;

  transition (StemLR,  StemLR)  = stem_extend_pg.YES * sym_bulge_open_pg.NO * asym_bulge_open_pg.NO * stem_bifurc_pg.NO;
  transition (StemLR,  StemBif) = stem_extend_pg.YES * sym_bulge_open_pg.NO * asym_bulge_open_pg.NO * stem_bifurc_pg.YES;
  transition (StemLR,  SymL)    = stem_extend_pg.YES * sym_bulge_open_pg.YES;
  transition (StemLR,  AsymL)   = stem_extend_pg.YES * sym_bulge_open_pg.NO * asym_bulge_open_pg.YES / 2;
  transition (StemLR,  AsymR)   = stem_extend_pg.YES * sym_bulge_open_pg.NO * asym_bulge_open_pg.YES / 2;

  for (int i = 0; i < loop_states; ++i)
    transition (StemLR, Loop + i) = loop_start_pos[i];

  transition (SymL,    SymR)    = 1;

  transition (SymR,    StemLR)  = sym_bulge_extend_pg.NO * asym_bulge_open_pg.NO * stem_bifurc_pg.NO;
  transition (SymR,    StemBif) = sym_bulge_extend_pg.NO * asym_bulge_open_pg.NO * stem_bifurc_pg.YES;
  transition (SymR,    SymL)    = sym_bulge_extend_pg.YES;
  transition (SymR,    AsymL)   = sym_bulge_extend_pg.NO * asym_bulge_open_pg.YES / 2;
  transition (SymR,    AsymR)   = sym_bulge_extend_pg.NO * asym_bulge_open_pg.YES / 2;
  
  transition (AsymL,   StemLR)  = asym_bulge_extend_pg.NO * stem_bifurc_pg.NO;
  transition (AsymL,   StemBif) = asym_bulge_extend_pg.NO * stem_bifurc_pg.YES;
  transition (AsymL,   AsymL)   = asym_bulge_extend_pg.YES;

  transition (AsymR,   StemLR)  = asym_bulge_extend_pg.NO * stem_bifurc_pg.NO;
  transition (AsymR,   StemBif) = asym_bulge_extend_pg.NO * stem_bifurc_pg.YES;
  transition (AsymR,   AsymR)   = asym_bulge_extend_pg.YES;

  for (int i = 0; i < loop_states; ++i)
    {
      const int next_loop = i+1 < loop_states ? loop(i+1) : End;
      transition (loop(i), loop(i))   = loop_extend_pg.YES;
      transition (loop(i), next_loop) = loop_extend_pg.NO;
    }

  // set bifurcations
  bifurc[StemBif] = Bifurcation (PreStem, PrePad);

  // set emissions
  const int A = alphabet().size();
  vector<int> dinuc_i (2);
  int& xl (dinuc_i[0]);
  int& xr (dinuc_i[1]);
  for (xl = 0; xl < A; ++xl)
    for (xr = 0; xr < A; ++xr)
      {
	const int idx = emit_xl_mul(EmitXLR) * xl + emit_xr_mul(EmitXLR) * xr;
	const int i = dinuc.intvec2index (dinuc_i);
	emit[StemLR][idx] = dinuc[i] / (null_emit[xl] * null_emit[xr]);
      }
}

Dirichlet_prior Quick_stem::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  // create the prior, assign idiotically simple pseudocounts to all PGroup's
  // (not quite Laplace, because k<1; but close)
  Dirichlet_prior prior (pscore);
  prior.assign (pad_extend_pg, Laplace_prior (pad_extend_pg, k));
  prior.assign (stem_prior_pg, Laplace_prior (stem_prior_pg, k));
  prior.assign (stem_extend_pg, Laplace_prior (stem_extend_pg, k));
  prior.assign (stem_bifurc_pg, Laplace_prior (stem_bifurc_pg, k));
  prior.assign (sym_bulge_open_pg, Laplace_prior (sym_bulge_open_pg, k));
  prior.assign (sym_bulge_extend_pg, Laplace_prior (sym_bulge_extend_pg, k));
  prior.assign (asym_bulge_open_pg, Laplace_prior (asym_bulge_open_pg, k));
  prior.assign (asym_bulge_extend_pg, Laplace_prior (asym_bulge_extend_pg, k));
  prior.assign (loop_extend_pg, Laplace_prior (loop_extend_pg, k));
  prior.assign (loop_start_pos, Laplace_prior (loop_start_pos, k));
  prior.assign (dinuc, Laplace_prior (dinuc, k));
  // return the prior
  return prior;
}
