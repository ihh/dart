#include "stemloc/superstem.h"

Super_stem::Super_stem (PScores& pscore)
  : Odds_PCFG (TotalStates, pscore),
    pscore (pscore),
    start_in_stem (pscore.new_boolean_group ("startInStem")),
    loop_extend (pscore.new_boolean_group ("loopExtend")),
    bulge_extend (pscore.new_boolean_group ("bulgeExtend")),
    stem_extend (pscore.new_boolean_group ("stemExtend")),
    multi_extend (pscore.new_boolean_group ("multiExtend")),
    multi_bulge_open (pscore.new_boolean_group ("multiBulgeOpen")),
    post_stem (pscore.new_group ("stem_to_", "loop bulge lrbulge multiloop")),
    nuc (pscore.new_alphabet_group (CFG_alphabet, 1, "emit")),
    opening_dinuc (pscore.new_alphabet_group (CFG_alphabet, 2, "open")),
    stacked_dinuc (CFG_alphabet.size() * CFG_alphabet.size())
{
  // set grammar name
  name = "Super_stem";

  // set unlexicalized state names
  state_name[Loop] = "Loop";
  state_name[LBulge] = "LBulge";
  state_name[RBulge] = "RBulge";
  state_name[LRBulge] = "LRBulge";
  state_name[StemStart] = "StemStart";
  state_name[StemEnd] = "StemEnd";
  state_name[Multi] = "Multi";
  state_name[LMulti] = "LMulti";
  state_name[RMulti] = "RMulti";

  // create stacked dinucleotide PGroups
  const int A = alphabet().size();
  for (int i = 0; i < A*A; ++i)
    {
      sstring name;
      name << "stack" << opening_dinuc.index2word(i);
      stacked_dinuc[i] = pscore.new_alphabet_group (alphabet(), 2, name.c_str());
    }

  // set state types
  PFunc zero = 0;

  init_emit (Loop,    EmitXL, zero);
  init_emit (LBulge,  EmitXL, zero);
  init_emit (RBulge,  EmitXR, zero);
  init_emit (LRBulge, EmitXL, zero);

  state_type[StemStart] = Null;
  state_type[StemEnd]   = Null;

  for (int i = 0; i < A*A; ++i)
    init_emit (Stem+i, EmitXLR, zero);

  state_type[Multi]   = Bifurc;
  state_type[LMulti]  = Null;
  state_type[RMulti]  = Null;

  // set transitions
  // start state
  transition (Start, StemStart) = start_in_stem.y;
  transition (Start, StemEnd)   = start_in_stem.n;

  // loops
  transition (Loop, Loop) = loop_extend.y;
  transition (Loop, End)  = loop_extend.n;

  // bulges
  transition (LBulge,  LBulge)    = bulge_extend.y;
  transition (LBulge,  StemStart) = bulge_extend.n;

  transition (RBulge,  RBulge)    = bulge_extend.y;
  transition (RBulge,  StemStart) = bulge_extend.n;

  transition (LRBulge, LRBulge)   = bulge_extend.y;
  transition (LRBulge, RBulge)    = bulge_extend.n;

  // stems (including stem emissions)
  PFunc one = 1;
  vector<int> dinuc_i (2);
  vector<int> dinuc_j (2);
  for (dinuc_i[0] = 0; dinuc_i[0] < A; ++dinuc_i[0])
    for (dinuc_i[1] = 0; dinuc_i[1] < A; ++dinuc_i[1])
      {
	// calculate indices
	const int i = emit_xl_mul(EmitXLR) * dinuc_i[0] + emit_xr_mul(EmitXLR) * dinuc_i[1];
	const int dinuc_idx_i = opening_dinuc.intvec2index (dinuc_i);
	
	// set emission
	emit[Stem+i][i] = one;

	// set transitions
	transition (StemStart, Stem+i) = opening_dinuc[dinuc_idx_i];
	for (dinuc_j[0] = 0; dinuc_j[0] < A; ++dinuc_j[0])
	  for (dinuc_j[1] = 0; dinuc_j[1] < A; ++dinuc_j[1])
	    {
	      const int j = emit_xl_mul(EmitXLR) * dinuc_j[0] + emit_xr_mul(EmitXLR) * dinuc_j[1];
	      const int dinuc_idx_j = opening_dinuc.intvec2index (dinuc_j);
	      transition (Stem+j, Stem+i) = stem_extend.y * stacked_dinuc[dinuc_idx_j][dinuc_idx_i];
	    }
	transition (Stem+i, StemEnd) = stem_extend.n;

	// set state name
	state_name[Stem+i].clear();
	state_name[Stem+i] << "Stem" << CFG_alphabet.int2char (dinuc_i[0]) << CFG_alphabet.int2char (dinuc_i[1]);
      }

  // ends of stems
  transition (StemEnd, Loop)    = post_stem[0];
  transition (StemEnd, LBulge)  = post_stem[1] / 2;
  transition (StemEnd, RBulge)  = post_stem[1] / 2;
  transition (StemEnd, LRBulge) = post_stem[2];
  transition (StemEnd, Multi)   = post_stem[3];

  // multiloops
  bifurc[Multi] = Bifurcation (LMulti, RMulti);

  transition (LMulti, LBulge)    = multi_bulge_open.y;
  transition (LMulti, StemStart) = multi_bulge_open.n;

  transition (RMulti, Multi)     = multi_extend.y;
  transition (RMulti, StemStart) = multi_extend.n * multi_bulge_open.n * multi_bulge_open.n;
  transition (RMulti, LBulge)    = multi_extend.n * multi_bulge_open.y * multi_bulge_open.n;
  transition (RMulti, RBulge)    = multi_extend.n * multi_bulge_open.n * multi_bulge_open.y;
  transition (RMulti, LRBulge)   = multi_extend.n * multi_bulge_open.y * multi_bulge_open.y;

  // set unpaired nucleotide emissions
  for (int i = 0; i < A; ++i)
    emit[Loop][i] = emit[LBulge][i] = emit[RBulge][i] = emit[LRBulge][i] = PFunc(nuc[i]);

  // subtract off the null emission scores
  factor_in_null_model (null_emit);
  factor_in_null_extend (null_extend);
}

Dirichlet_prior Super_stem::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  // create the prior, assign idiotically simple pseudocounts to all PGroup's
  // (not quite Laplace, because k<1; but close)
  Dirichlet_prior prior (pscore);
  prior.assign_Laplace (start_in_stem, k);
  prior.assign_Laplace (loop_extend, k);
  prior.assign_Laplace (bulge_extend, k);
  prior.assign_Laplace (stem_extend, k);
  prior.assign_Laplace (multi_extend, k);
  prior.assign_Laplace (multi_bulge_open, k);
  prior.assign_Laplace (post_stem, k);
  prior.assign_Laplace (nuc, k);
  prior.assign_Laplace (opening_dinuc, k);
  for_const_contents (vector<Alphabet_group>, stacked_dinuc, pn)
    prior.assign_Laplace (*pn, k);
  // null model
  prior.assign_Laplace (null_emit, k);
  prior.assign_Laplace (null_extend, k);
  // return the prior
  return prior;
}
