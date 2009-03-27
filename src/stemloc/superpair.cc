#include "stemloc/superpair.h"

Super_pair::Super_pair() : Odds_PCFG(), pscore (0)
{ }

Super_pair::Super_pair (PScores* pscore)
  : Odds_PCFG (TotalStates, *pscore),
    pscore (pscore),
    start_in_stem (pscore->new_boolean_group ("startInStem")),
    loop_extend (pscore->new_boolean_group ("loopExtend")),
    loop_gap_open (pscore->new_boolean_group ("loopGapOpen")),
    loop_gap_extend (pscore->new_boolean_group ("loopGapExtend")),
    loop_gap_swap (pscore->new_boolean_group ("loopGapSwap")),
    stem_extend (pscore->new_boolean_group ("stemExtend")),
    stem_gap_open (pscore->new_boolean_group ("stemGapOpen")),
    stem_gap_extend (pscore->new_boolean_group ("stemGapExtend")),
    stem_gap_swap (pscore->new_boolean_group ("stemGapSwap")),
    multi_extend (pscore->new_boolean_group ("multiExtend")),
    multi_bulge_open (pscore->new_boolean_group ("multiBulgeOpen")),
    post_stem (pscore->new_group ("stem_to_", "loop bulge lrbulge multiloop")),
    single_nuc (pscore->new_alphabet_group (CFG_alphabet, 1, "emit")),
    pair_nuc (CFG_alphabet.size()),
    single_dinuc (pscore->new_alphabet_group (CFG_alphabet, 2, "pair")),
    pair_dinuc (CFG_alphabet.size() * CFG_alphabet.size())
{
  // set grammar name
  name = "Super_pair";

  // set state names
  state_name[Loop] = "Loop";
  state_name[LoopMatch] = "LoopMatch";
  state_name[LoopIns] = "LoopIns";
  state_name[LoopDel] = "LoopDel";
  state_name[LBulge] = "LBulge";
  state_name[LBulgeMatch] = "LBulgeMatch";
  state_name[LBulgeIns] = "LBulgeIns";
  state_name[LBulgeDel] = "LBulgeDel";
  state_name[RBulge] = "RBulge";
  state_name[RBulgeMatch] = "RBulgeMatch";
  state_name[RBulgeIns] = "RBulgeIns";
  state_name[RBulgeDel] = "RBulgeDel";
  state_name[LRBulge] = "LRBulge";
  state_name[LRBulgeMatch] = "LRBulgeMatch";
  state_name[LRBulgeIns] = "LRBulgeIns";
  state_name[LRBulgeDel] = "LRBulgeDel";
  state_name[Stem] = "Stem";
  state_name[StemMatch] = "StemMatch";
  state_name[StemIns] = "StemIns";
  state_name[StemDel] = "StemDel";
  state_name[StemEnd] = "StemEnd";
  state_name[Multi] = "Multi";
  state_name[LMulti] = "LMulti";
  state_name[RMulti] = "RMulti";

  // create substitution PGroups
  const int A = alphabet().size();
  for (int i = 0; i < A; ++i)
    {
      sstring name;
      name << "sub" << single_nuc.index2word(i);
      pair_nuc[i] = pscore->new_alphabet_group (alphabet(), 1, name.c_str());
    }
  for (int i = 0; i < A*A; ++i)
    {
      sstring name;
      name << "sub" << single_dinuc.index2word(i);
      pair_dinuc[i] = pscore->new_alphabet_group (alphabet(), 2, name.c_str());
    }

  // set state types
  PFunc dummy;

  state_type[Stem]    = Null;
  state_type[StemEnd] = Null;

  init_emit (StemMatch, EmitXLRYLR, dummy);
  init_emit (StemIns,   EmitYLR,    dummy);
  init_emit (StemDel,   EmitXLR,    dummy);

  state_type[Multi]   = Bifurc;
  state_type[LMulti]  = Null;
  state_type[RMulti]  = Null;

  // initialise single-stranded regions
  init_unpaired (Loop,    LoopMatch,    LoopIns,    LoopDel,    End,    EmitXL, EmitYL);
  init_unpaired (LBulge,  LBulgeMatch,  LBulgeIns,  LBulgeDel,  Stem,   EmitXL, EmitYL);
  init_unpaired (RBulge,  RBulgeMatch,  RBulgeIns,  RBulgeDel,  Stem,   EmitXR, EmitYR);
  init_unpaired (LRBulge, LRBulgeMatch, LRBulgeIns, LRBulgeDel, RBulge, EmitXL, EmitYL);

  // set transitions
  // stems
  transition (Stem, StemMatch) = stem_gap_open.n;
  transition (Stem, StemIns)   = stem_gap_open.y / 2;
  transition (Stem, StemDel)   = stem_gap_open.y / 2;

  transition (StemMatch, StemMatch) = stem_gap_open.n * stem_extend.y;
  transition (StemMatch, StemIns)   = stem_gap_open.y / 2;
  transition (StemMatch, StemDel)   = stem_gap_open.y / 2;
  transition (StemMatch, StemEnd)   = stem_gap_open.n * stem_extend.n;

  transition (StemIns,   StemMatch) = stem_gap_extend.n * stem_gap_swap.n * stem_extend.y;
  transition (StemIns,   StemIns)   = stem_gap_extend.y;
  transition (StemIns,   StemDel)   = stem_gap_extend.n * stem_gap_swap.y;
  transition (StemIns,   StemEnd)   = stem_gap_extend.n * stem_gap_swap.n * stem_extend.n;

  transition (StemDel,   StemMatch) = stem_gap_extend.n * stem_gap_swap.n * stem_extend.y;
  transition (StemDel,   StemDel)   = stem_gap_extend.y;
  transition (StemDel,   StemIns)   = stem_gap_extend.n * stem_gap_swap.y;
  transition (StemDel,   StemEnd)   = stem_gap_extend.n * stem_gap_swap.n * stem_extend.n;

  // ends of stems
  transition (StemEnd,   Loop)    = post_stem[0];
  transition (StemEnd,   LBulge)  = post_stem[1] / 2;
  transition (StemEnd,   RBulge)  = post_stem[1] / 2;
  transition (StemEnd,   LRBulge) = post_stem[2];
  transition (StemEnd,   Multi)   = post_stem[3];

  // start state
  // disallow transitions from Start -> Loop, so all hits must have at least one basepair
  const PFunc not_loop        = post_stem[1] + post_stem[2] + post_stem[3];
  transition (Start, Stem)    = start_in_stem.y;
  transition (Start, LBulge)  = start_in_stem.n * post_stem[1] / (2 * not_loop);
  transition (Start, RBulge)  = start_in_stem.n * post_stem[1] / (2 * not_loop);
  transition (Start, LRBulge) = start_in_stem.n * post_stem[2] / not_loop;
  transition (Start, Multi)   = start_in_stem.n * post_stem[3] / not_loop;

  // multiloops
  bifurc[Multi] = Bifurcation (LMulti, RMulti);

  transition (LMulti, LBulge) = multi_bulge_open.y;
  transition (LMulti, Stem)   = multi_bulge_open.n;

  transition (RMulti, Multi)   = multi_extend.y;
  transition (RMulti, Stem)    = multi_extend.n * multi_bulge_open.n * multi_bulge_open.n;
  transition (RMulti, LBulge)  = multi_extend.n * multi_bulge_open.y * multi_bulge_open.n;
  transition (RMulti, RBulge)  = multi_extend.n * multi_bulge_open.n * multi_bulge_open.y;
  transition (RMulti, LRBulge) = multi_extend.n * multi_bulge_open.y * multi_bulge_open.y;

  // set emissions
  // dinucleotide indels
  vector<int> dinuc_i (2);
  vector<int> dinuc_j (2);
  for (dinuc_i[0] = 0; dinuc_i[0] < A; ++dinuc_i[0])
    for (dinuc_i[1] = 0; dinuc_i[1] < A; ++dinuc_i[1])
      {
	const int x_idx = dinuc_i[0] * emit_xl_mul(EmitXLR) + dinuc_i[1] * emit_xr_mul(EmitXLR);
	const int y_idx = dinuc_i[0] * emit_yl_mul(EmitYLR) + dinuc_i[1] * emit_yr_mul(EmitYLR);
	const int dinuc_idx = single_dinuc.intvec2index (dinuc_i);
	emit[StemDel][x_idx] = emit[StemIns][y_idx] = PFunc(single_dinuc[dinuc_idx]);
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
	    emit[StemMatch][idx] = (PFunc(single_dinuc[i]) * PFunc(pair_dinuc[i][j]) + PFunc(single_dinuc[j]) * PFunc(pair_dinuc[j][i])) / 2;
	  }
  // subtract off the null emission scores
  factor_in_null_model (null_emit);
  factor_in_null_extend (null_extend);
}

void Super_pair::init_unpaired (int start, int match, int ins, int del, int end, State_type xemit, State_type yemit)
{
  // set state types
  const PFunc dummy;
  const State_type xyemit = (State_type) ((int) xemit | (int) yemit);

  state_type[start] = Null;

  init_emit (match, xyemit, dummy);
  init_emit (ins,   yemit,  dummy);
  init_emit (del,   xemit,  dummy);

  // set transitions
  transition (start, match) = loop_gap_open.n;
  transition (start, ins)   = loop_gap_open.y / 2;
  transition (start, del)   = loop_gap_open.y / 2;

  transition (match, match) = loop_gap_open.n * loop_extend.y;
  transition (match, ins)   = loop_gap_open.y / 2;
  transition (match, del)   = loop_gap_open.y / 2;
  transition (match, end)   = loop_gap_open.n * loop_extend.n;

  transition (ins,   match) = loop_gap_extend.n * loop_gap_swap.n * loop_extend.y;
  transition (ins,   ins)   = loop_gap_extend.y;
  transition (ins,   del)   = loop_gap_extend.n * loop_gap_swap.y;
  transition (ins,   end)   = loop_gap_extend.n * loop_gap_swap.n * loop_extend.n;

  transition (del,   match) = loop_gap_extend.n * loop_gap_swap.n * loop_extend.y;
  transition (del,   del)   = loop_gap_extend.y;
  transition (del,   ins)   = loop_gap_extend.n * loop_gap_swap.y;
  transition (del,   end)   = loop_gap_extend.n * loop_gap_swap.n * loop_extend.n;

  // set emissions
  const int A = alphabet().size();

  // single nucleotide indels
  for (int i = 0; i < A; ++i)
    emit[ins][i] = emit[del][i] = PFunc(single_nuc[i]);

  // single nucleotide matches
  for (int i = 0; i < A; ++i)
    for (int j = 0; j < A; ++j)
      {
	const int idx = i * (xemit == EmitXL ? emit_xl_mul(xyemit) : emit_xr_mul(xyemit)) + j * (yemit == EmitYL ? emit_yl_mul(xyemit) : emit_yr_mul(xyemit));
	emit[match][idx] = (PFunc(single_nuc[i]) * PFunc(pair_nuc[i][j]) + PFunc(single_nuc[j]) * PFunc(pair_nuc[j][i])) / 2;
      }
}

Dirichlet_prior Super_pair::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  // create the prior, assign idiotically simple pseudocounts to all PGroup's
  // (not quite Laplace, because k<1; but close)
  Dirichlet_prior prior (*pscore);
  prior.assign_Laplace (start_in_stem, k);
  prior.assign_Laplace (loop_extend, k);
  prior.assign_Laplace (loop_gap_open, k);
  prior.assign_Laplace (loop_gap_extend, k);
  prior.assign_Laplace (loop_gap_swap, k);
  prior.assign_Laplace (stem_extend, k);
  prior.assign_Laplace (stem_gap_open, k);
  prior.assign_Laplace (stem_gap_extend, k);
  prior.assign_Laplace (stem_gap_swap, k);
  prior.assign_Laplace (multi_extend, k);
  prior.assign_Laplace (multi_bulge_open, k);
  prior.assign_Laplace (post_stem, k);
  prior.assign_Laplace (single_nuc, k);
  for_const_contents (vector<Alphabet_group>, pair_nuc, pn)
    prior.assign_Laplace (*pn, k);
  prior.assign_Laplace (single_dinuc, k);
  for_const_contents (vector<Alphabet_group>, pair_dinuc, pd)
    prior.assign_Laplace (*pd, k);
  // null model
  prior.assign_Laplace (null_emit, k);
  prior.assign_Laplace (null_extend, k);
  // return the prior
  return prior;
}
