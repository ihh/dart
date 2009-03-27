#include "stemloc/stempair.h"

Stem_pair::Stem_pair (PScores& pscore, bool use_bif)
  : Pair_PCFG (use_bif ? StatesWithBif : StatesWithoutBif),
    pscore (pscore),
    stem_LR_pg (pscore.new_group ("stemLR_", "L LR B")),
    stem_L_pg (pscore.new_group ("stemL_", "L R LR")),
    stem_pg (pscore.new_group ("stem_", "stem loop")),
    loop_pg (pscore.new_group ("loop_", "loop end")),
    stem_mat_pg (pscore.new_group ("stemMatch", "Match Ins")),
    stem_ins_pg (pscore.new_group ("stemIns", "Match Ins")),
    loop_mat_pg (pscore.new_group ("loopMatch", "Match Ins")),
    loop_ins_pg (pscore.new_group ("loopIns", "Match Ins")),
    single_nuc (pscore.new_alphabet_group (CFG_alphabet, 1, "emit")),
    pair_nuc (CFG_alphabet.size()),
    single_dinuc (pscore.new_alphabet_group (CFG_alphabet, 2, "emit")),
    pair_dinuc (CFG_alphabet.size() * CFG_alphabet.size()),
    null_extend (pscore.new_boolean_group ("nullExtend")),
    null_emit (pscore.new_alphabet_group (CFG_alphabet, 1, "nullEmit"))
{
  // create substitution PGroups
  const int A = alphabet().size();
  for (int i = 0; i < A; ++i)
    {
      sstring name;
      name << "sub" << single_nuc.index2word(i) << '_';
      pair_nuc[i] = pscore.new_alphabet_group (alphabet(), 1, name.c_str());
    }
  for (int i = 0; i < A*A; ++i)
    {
      sstring name;
      name << "sub" << single_dinuc.index2word(i) << '_';
      pair_dinuc[i] = pscore.new_alphabet_group (alphabet(), 2, name.c_str());
    }

  // set state types
  PFunc dummy;

  init_emit (StemMatLR, EmitXLRYLR, dummy);
  init_emit (StemMatL,  EmitXLYL,   dummy);
  init_emit (StemMatR,  EmitXRYR,   dummy);

  init_emit (StemDelLR, EmitXLR,    dummy);
  init_emit (StemDelL,  EmitXL,     dummy);
  init_emit (StemDelR,  EmitXR,     dummy);

  init_emit (StemInsLR, EmitYLR,    dummy);
  init_emit (StemInsL,  EmitYL,     dummy);
  init_emit (StemInsR,  EmitYR,     dummy);

  init_emit (LoopMat,   EmitXLYL,   dummy);
  init_emit (LoopDel,   EmitXL,     dummy);
  init_emit (LoopIns,   EmitYL,     dummy);

  if (use_bif)
    {
      state_type[Bif]      = Bifurc;
      state_type[BifStart] = Null;
    }

  // set transitions

  // transitions from Start into stem match states

  transition (Start, StemMatLR) = stem_mat_to_mat() * stem_LR_to_LR();
  transition (Start, StemMatL)  = stem_mat_to_mat() * stem_LR_to_L() / 2;
  transition (Start, StemMatR)  = stem_mat_to_mat() * stem_LR_to_L() / 2;

  // allow transitions from Start into non-match states, so training alignments don't have zero likelihood

  transition (Start, StemDelLR) = stem_mat_to_ins() * stem_LR_to_LR() / 2;
  transition (Start, StemDelL)  = stem_mat_to_ins() * stem_LR_to_L() / 4;
  transition (Start, StemDelR)  = stem_mat_to_ins() * stem_LR_to_L() / 4;

  transition (Start, StemInsLR) = stem_mat_to_ins() * stem_LR_to_LR() / 2;
  transition (Start, StemInsL)  = stem_mat_to_ins() * stem_LR_to_L() / 4;
  transition (Start, StemInsR)  = stem_mat_to_ins() * stem_LR_to_L() / 4;

  transition (Start, End) = 1;

  for (int src = StemMatLR; src <= StemInsR; ++src)
    {
      const int src_mdi  = src / 3;
      const int src_lr   = src % 3;
      if (use_bif && src_lr == 0)
	transition (src, Bif) = stem_to_stem() * stem_LR_to_B();

      for (int dest = StemMatLR; dest <= StemInsR; ++dest)
	{
	  const int dest_mdi = dest / 3;
	  const int dest_lr  = dest % 3;
	  PFunc t = 1;
	  if (src_lr  == 0) t *= stem_to_stem();
	  if (src_mdi == 0 && dest_mdi == 0) t *= stem_mat_to_mat();
	  if (src_mdi == 0 && dest_mdi != 0) t *= stem_mat_to_ins() / 2;
	  if (src_mdi != 0 && dest_mdi == 0) t *= stem_ins_to_mat();
	  if (src_mdi != 0 && dest_mdi != 0) t *= stem_ins_to_ins() / 2;
	  if (src_lr == 0 && dest_lr == 0) t *= stem_LR_to_LR();
	  if (src_lr == 0 && dest_lr != 0) t *= stem_LR_to_L() / 2;
	  if (src_lr != 0 && dest_lr == 0) t *= stem_L_to_LR();
	  if (src_lr != 0 && dest_lr != 0 && dest_lr == src_lr) t *= stem_L_to_L();
	  if (src_lr != 0 && dest_lr != 0 && dest_lr != src_lr) t *= stem_L_to_R();
	  transition (src, dest) = t;
	  if (src_lr == 0)
	    {
	      transition (src, LoopMat) = stem_to_loop() * loop_mat_to_mat();
	      transition (src, LoopDel) = stem_to_loop() * loop_mat_to_ins() / 2;
	      transition (src, LoopIns) = stem_to_loop() * loop_mat_to_ins() / 2;
	    }
	}
    }

  transition (LoopMat, LoopMat) = loop_to_loop() * loop_mat_to_mat();
  transition (LoopMat, LoopDel) = loop_to_loop() * loop_mat_to_ins() / 2;
  transition (LoopMat, LoopIns) = loop_to_loop() * loop_mat_to_ins() / 2;
  transition (LoopMat, End)     = loop_to_end();

  transition (LoopDel, LoopMat) = loop_to_loop() * loop_ins_to_mat();
  transition (LoopDel, LoopDel) = loop_to_loop() * loop_ins_to_ins() / 2;
  transition (LoopDel, LoopIns) = loop_to_loop() * loop_ins_to_ins() / 2;
  transition (LoopDel, End)     = loop_to_end();

  transition (LoopIns, LoopMat) = loop_to_loop() * loop_ins_to_mat();
  transition (LoopIns, LoopDel) = loop_to_loop() * loop_ins_to_ins() / 2;
  transition (LoopIns, LoopIns) = loop_to_loop() * loop_ins_to_ins() / 2;
  transition (LoopIns, End)     = loop_to_end();

  if (use_bif)
    {
      bifurc[Bif] = Bifurcation (BifStart, BifStart);
      for (int dest = StemMatLR; dest <= StemInsR; ++dest)
	transition (BifStart, dest) = transition (StemMatLR, dest);
      transition (BifStart, Bif) = transition (StemMatLR, Bif);
      transition (Start, Bif) = stem_mat_to_mat() * stem_LR_to_B();
    }

  // set emissions
  // single nucleotide indels
  for (int i = 0; i < A; ++i)
    emit[StemDelL][i] = emit[StemDelR][i] = emit[StemInsL][i] = emit[StemInsR][i] = emit[LoopIns][i] = emit[LoopDel][i] = PFunc(single_nuc[i]);
  // single nucleotide matches
  for (int i = 0; i < A; ++i)
    for (int j = 0; j < A; ++j)
      {
	const int l_idx = i * emit_xl_mul(EmitXLYL) + j * emit_yl_mul(EmitXLYL);
	const int r_idx = i * emit_xr_mul(EmitXRYR) + j * emit_yr_mul(EmitXRYR);
	emit[StemMatL][l_idx] = emit[StemMatR][r_idx] = emit[LoopMat][l_idx] = (PFunc(single_nuc[i]) * PFunc(pair_nuc[i][j]) + PFunc(single_nuc[j]) * PFunc(pair_nuc[j][i])) / 2;
      }
  // dinucleotide indels
  vector<int> dinuc_i (2);
  vector<int> dinuc_j (2);
  for (dinuc_i[0] = 0; dinuc_i[0] < A; ++dinuc_i[0])
    for (dinuc_i[1] = 0; dinuc_i[1] < A; ++dinuc_i[1])
      {
	const int x_idx = dinuc_i[0] * emit_xl_mul(EmitXLR) + dinuc_i[1] * emit_xr_mul(EmitXLR);
	const int y_idx = dinuc_i[0] * emit_yl_mul(EmitYLR) + dinuc_i[1] * emit_yr_mul(EmitYLR);
	const int dinuc_idx = single_dinuc.intvec2index (dinuc_i);
	emit[StemDelLR][x_idx] = emit[StemInsLR][y_idx] = PFunc(single_dinuc[dinuc_idx]);
      }
  // dinucleotide matches
  for (dinuc_i[0] = 0; dinuc_i[0] < A; ++dinuc_i[0])
    for (dinuc_i[1] = 0; dinuc_i[1] < A; ++dinuc_i[1])
      for (dinuc_j[0] = 0; dinuc_j[0] < A; ++dinuc_j[0])
	for (dinuc_j[1] = 0; dinuc_j[1] < A; ++dinuc_j[1])
	  {
	    const int idx = dinuc_i[0] * emit_xl_mul(EmitXLRYLR) + dinuc_i[1] * emit_xr_mul(EmitXLRYLR) + dinuc_j[0] * emit_yl_mul(EmitXLRYLR) + dinuc_j[1] * emit_yr_mul(EmitXLRYLR);
	    const int i = single_dinuc.intvec2index (dinuc_i);
	    const int j = single_dinuc.intvec2index (dinuc_j);
	    emit[StemMatLR][idx] = (PFunc(single_dinuc[i]) * PFunc(pair_dinuc[i][j]) + PFunc(single_dinuc[j]) * PFunc(pair_dinuc[j][i])) / 2;
	  }
  // subtract off the null emission scores
  factor_in_null_model (null_emit);
  factor_in_null_extend (null_extend);
}

Dirichlet_prior Stem_pair::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  // create the prior, assign idiotically simple pseudocounts to all PGroup's
  // (not quite Laplace, because k<1; but close)
  Dirichlet_prior prior (pscore);
  prior.assign (stem_LR_pg, Laplace_prior (stem_LR_pg, k));
  prior.assign (stem_L_pg, Laplace_prior (stem_L_pg, k));
  prior.assign (stem_mat_pg, Laplace_prior (stem_mat_pg, k));
  prior.assign (stem_ins_pg, Laplace_prior (stem_ins_pg, k));
  prior.assign (loop_mat_pg, Laplace_prior (loop_mat_pg, k));
  prior.assign (loop_ins_pg, Laplace_prior (loop_ins_pg, k));
  prior.assign (stem_pg, Laplace_prior (stem_pg, k));
  prior.assign (loop_pg, Laplace_prior (loop_pg, k));
  prior.assign (single_nuc, Laplace_prior (single_nuc, k));
  for_const_contents (vector<Alphabet_group>, pair_nuc, pn)
    prior.assign (*pn, Laplace_prior (*pn, k));
  prior.assign (single_dinuc, Laplace_prior (single_dinuc, k));
  for_const_contents (vector<Alphabet_group>, pair_dinuc, pd)
    prior.assign (*pd, Laplace_prior (*pd, k));
  // null model
  prior.assign (null_emit, Laplace_prior (null_emit, k));
  prior.assign (null_extend, Laplace_prior (null_extend, k));
  // return the prior
  return prior;
}
