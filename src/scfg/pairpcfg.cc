#include "scfg/pairpcfg.h"
#include "util/vector_output.h"

Pair_PCFG::Pair_PCFG (int states) :
  Pair_CFG<PFunc> (states, PFunc(0.0)),
  group_suffix (0)
{ }

Pair_PCFG::Pair_PCFG (int states, const vector<vector<sstring> >& group_suffix) :
  Pair_CFG<PFunc> (states, PFunc(0.0)),
  group_suffix (&group_suffix)
{ }

Pair_PCFG::Pair_PCFG (const Pair_PHMM& hmm)
  : Pair_CFG<PFunc> (hmm)
{
  // doesn't set group_suffix, because that's fucking private in Pair_PHMM. ho hum
}

Pair_PCFG::Pair_PCFG (const Pair_CFG_scores& cfg)
  : Pair_CFG<PFunc> (cfg, PFunc())
{
  // transition matrix
  for (int src = Start; src < cfg.states(); ++src)
    if (src != End)
      for (int dest = End; dest < cfg.states(); ++dest)
	if (dest != Start)
	  if (cfg.transition(src,dest) > -InfinityScore)
	    transition(src,dest) = PFunc (Score2Prob (cfg.transition(src,dest)));
  // emissions
  for (int s = 0; s < cfg.states(); ++s)
    for (int e = 0; e < (int) cfg.emit[s].size(); ++e)
      if (cfg.emit[s][e] > -InfinityScore)
	emit[s][e] = PFunc (Score2Prob (cfg.emit[s][e]));
}

void Pair_PCFG::set_group_suffix (const vector<vector<sstring> >& g_suffix)
{
  group_suffix = &g_suffix;
}

Pair_CFG_scores Pair_PCFG::eval_cfg_scores_gs (const PScores& var_scores) const
{
  ((Pair_PCFG*) this) -> set_group_suffix (var_scores.group_suffix);
  return eval_cfg_scores (var_scores);
}

Pair_CFG_scores Pair_PCFG::eval_cfg_scores (const PScores& var_scores) const
{
  if (CTAGGING(3,CFG_PFUNC))
    {
      CL << "Evaluating CFG scores from PFunc's:\n";
      show(CL);
    }
  if (CTAGGING(1,CFG_PSCORES))
    {
      CL << "Using the following PScores to evaluate CFG scores:\n";
      var_scores.show(CL);
    }
  Pair_CFG_scores cfg_scores (*this);
  cfg_scores.name = name;
  cfg_scores.state_name = state_name;
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	cfg_scores.transition(src,dest) = transition(src,dest).eval_sc(var_scores);
      for (int i = 0; i < (int) emit[src].size(); ++i)
	cfg_scores.emit[src][i] = emit[src][i].eval_sc(var_scores);
      cfg_scores.start[src] = start[src].eval_sc(var_scores);
      cfg_scores.end[src] = end[src].eval_sc(var_scores);
    }
  cfg_scores.start_to_end() = start_to_end().eval_sc(var_scores);
  if (CTAGGING(-1,CFG_SCORES))
    {
      CL << "Evaluated the following CFG scores from PFunc's:\n";
      cfg_scores.show (CL);
    }
  return cfg_scores;
}

void Pair_PCFG::inc_var_counts (PCounts& var_counts,
				const PScores& var_scores,
				const Pair_CFG_counts& cfg_counts,
				const Prob model_count) const
{
  assert_same_dimensions (cfg_counts);
  for (int src = 0; src < states(); src++)
    {
      for (int dest = 0; dest < states(); dest++)
	transition(src,dest).inc_var_counts (var_counts, var_scores, cfg_counts.transition(src,dest) * model_count);
      for (int i = 0; i < (int) emit[src].size(); ++i)
	emit[src][i].inc_var_counts (var_counts, var_scores, cfg_counts.emit[src][i] * model_count);
      start[src].inc_var_counts (var_counts, var_scores, cfg_counts.start[src] * model_count);
      end[src].inc_var_counts (var_counts, var_scores, cfg_counts.end[src] * model_count);
    }
  start_to_end().inc_var_counts (var_counts, var_scores, cfg_counts.start_to_end() * model_count);
}

void Pair_PCFG::factor_in_null_model (const Alphabet_group& null_emit)
{
  const int A = alphabet().size();
  for (int s = 0; s < states(); ++s)
    {
      const State_type t = state_type[s];
      if (!is_emit_type(t)) continue;
      const int A_xl = t & EmitXL ? A : 1;
      const int A_xr = t & EmitXR ? A : 1;
      const int A_yl = t & EmitYL ? A : 1;
      const int A_yr = t & EmitYR ? A : 1;
      for (int xl = 0; xl < A_xl; ++xl)
	for (int xr = 0; xr < A_xr; ++xr)
	  for (int yl = 0; yl < A_yl; ++yl)
	    for (int yr = 0; yr < A_yr; ++yr)
	      {
		const int idx = emit_xl_mul(t) * xl + emit_xr_mul(t) * xr + emit_yl_mul(t) * yl + emit_yr_mul(t) * yr;
		PFunc& e = emit[s][idx];
		if (!e.is_zero())
		  {
		    if (t & EmitXL) e /= null_emit[xl];
		    if (t & EmitXR) e /= null_emit[xr];
		    if (t & EmitYL) e /= null_emit[yl];
		    if (t & EmitYR) e /= null_emit[yr];
		  }
	      }
    }
}

void Pair_PCFG::factor_in_null_extend (const Boolean_group& null_extend)
{
  for (int d = 0; d < states(); ++d)
    {
      // factor null_extend.n out of Start transition, not End transition...
      // ...because Start transition only gets used once, but End may get used many times
      if (!transition(Start,d).is_zero())
	transition(Start,d) /= null_extend.n;
      // the rest involves null_extend.y, & so is only relevant for emit states
      const State_type t = state_type[d];
      if (!is_emit_type(t)) continue;
      PFunc f = 1;
      if (t & EmitXL) f *= null_extend.y;
      if (t & EmitXR) f *= null_extend.y;
      if (t & EmitYL) f *= null_extend.y;
      if (t & EmitYR) f *= null_extend.y;
      for (int s = Start; s < states(); ++s)
	if (s != End)
	  if (!transition(s,d).is_zero())
	    transition(s,d) /= f;
    }
  if (!transition(Start,End).is_zero())
    transition(Start,End) /= null_extend.n;
}

Pair_local_PCFG::Pair_local_PCFG (const Pair_PCFG& global,
				  PScope& pscope,
				  int mask_metascore_idx) :
  Pair_PCFG (global.states() + GlobalOffset),
  global (global),
  mask_metascore_idx (mask_metascore_idx)
{
  update();
}

void Pair_local_PCFG::update()
{
  // debugging output
  if (CTAGGING(-2,LOCAL_PCFG))
    {
      CL << "Building local model from the following global model:\n";
      global.show (CL);
    }
  // test that we've got the right number of states
  const int gs = global.states();
  if (states() != local_state(gs))
    THROWEXPR ("Global model has " << gs << " states, local model has " << states() << " != " << gs << " + " << GlobalOffset);
  // clear all transitions, emit vectors & metascore indices
  const PFunc zero (0.0);
  const PFunc one  (1.0);
  reset (zero);
  reset_meta();
  // set padding state_types
  init_emit (LFlankBif,  Null,     one);  // set to Bifurc in allow_lr_hits(), to avoid unnecessary Bifurc state
  init_emit (RFlankBif,  Null,     one);  // set to Bifurc in allow_lr_hits(), to avoid unnecessary Bifurc state
  init_emit (PadBif,     Null,     one);  // set to Bifurc in allow_multiple_hits(), to avoid unnecessary Bifurc state
  init_emit (PreFlank,   Null,     one);
  init_emit (PrePad,     Null,     one);
  init_emit (PreBif,     Null,     one);
  init_emit (XLFlank,    EmitXL,   one);
  init_emit (XRFlank,    EmitXR,   one);
  init_emit (XLRFlank,   EmitXLR,  zero);  // set to one in allow_lr_pad()
  init_emit (YLFlank,    EmitYL,   one);
  init_emit (YRFlank,    EmitYR,   one);
  init_emit (YLRFlank,   EmitYLR,  zero);  // set to one in allow_lr_pad()
  init_emit (XLPad,      EmitXL,   one);
  init_emit (XRPad,      EmitXR,   one);
  init_emit (XLRPad,     EmitXLR,  zero);  // set to one in allow_lr_pad()
  init_emit (YLPad,      EmitYL,   one);
  init_emit (YRPad,      EmitYR,   one);
  init_emit (YLRPad,     EmitYLR,  zero);  // set to one in allow_lr_pad()
  init_emit (LocalStart, Null,     one);

  // set transition "probabilities"
  transition (Start, PrePad) = 1;
  // do transitions for (XLFlank XRFlank XLRFlank)
  // likewise   (Y*Flank...),  (X*Pad...),  (Y*Pad...)
  for (int i = XLFlank; i <= XLPad; i += 6)
    {
      for (int j = i; j < i+6; ++j)  // forward transitions within block
	for (int k = j; k < i+6; ++k)
	  transition (j, k) = 1;
      for (int j = i; j < i+6; j += 3)  // loop back from XLRFlank to XLFlank, XRFlank
	{
	  transition (j+2, j)   = 1;
	  transition (j+2, j+1) = 1;
	}
    }
  for (int i = XLFlank; i <= YLRFlank; ++i)
    {
      transition (PreFlank, i) = 1;
      transition (i, End) = 1;
    }
  for (int i = XLPad; i <= YLRPad; ++i)
    {
      transition (PrePad, i) = 1;
      transition (i, PreBif) = 1;
      transition (i, LocalStart) = 1;
    }
  transition (PrePad, PreBif) = 1;
  transition (PrePad, LocalStart) = 1;

  // copy the global model
  transition (LocalStart, End) = global.transition (Start, End);
  for (int s = 0; s < gs; ++s)
    {
      int ls = local_state(s);
      // start & end transitions
      transition (LocalStart, ls) = global.transition (Start, s);
      transition (ls, End)        = global.transition (s, End);
      // other transitions
      for (int t = 0; t < gs; ++t)
	transition (local_state(t), ls) = global.transition(t,s);
      // emit profiles
      emit[ls] = global.emit[s];
      // metascore indices
      xlmeta_idx[ls] = global.xlmeta_idx[s];
      xrmeta_idx[ls] = global.xrmeta_idx[s];
      ylmeta_idx[ls] = global.ylmeta_idx[s];
      yrmeta_idx[ls] = global.yrmeta_idx[s];
    }
}

void Pair_local_PCFG::allow_multiple_hits()
{
  bifurc[PadBif] = Bifurcation (PrePad, PrePad);
  state_type[PadBif] = Bifurc;  // no need to call init_emit(), as it was a Null state before
  transition (PreBif, PadBif) = 1;
}

void Pair_local_PCFG::allow_lr_pad()
{
  bifurc[LFlankBif] = Bifurcation (PreFlank, PrePad);
  bifurc[RFlankBif] = Bifurcation (PrePad, PreFlank);

  transition (PreBif, LFlankBif) = 1;
  transition (PreBif, RFlankBif) = 1;
  
  PFunc one (1.0);
  init_emit (LFlankBif, Bifurc,  one);
  init_emit (RFlankBif, Bifurc,  one);
  init_emit (XLRFlank,  EmitXLR, one);
  init_emit (YLRFlank,  EmitYLR, one);
  init_emit (XLRPad,    EmitXLR, one);
  init_emit (YLRPad,    EmitYLR, one);
}

void Pair_local_PCFG::reset_mask (Named_profile& np) const
{
  if ((int) np.meta_sc.size() < mask_metascore_idx)
    CLOGERR << "Warning -- sequence '" << np.name << "' is missing metascores\n";
  while ((int) np.meta_sc.size() < mask_metascore_idx + 1) np.meta_sc.push_back (Metascore());
  np.meta_sc[mask_metascore_idx] = Metascore (np.size(), (Score) 0);
}

void Pair_local_PCFG::reset_masks (Sequence_database& db) const
{
  for_contents (Sequence_database, db, np) reset_mask (*np);
}

void Pair_local_PCFG::mask (Named_profile& np, const vector<Metaprob>& expected_metacounts) const
{
  const Metaprob& in_motif_prob = expected_metacounts[mask_metascore_idx];
  Metascore delta (in_motif_prob.size());
  
  for (int pos = 0; pos < (int) in_motif_prob.size(); ++pos)
    delta[pos] = Prob2Score (1.0 - in_motif_prob[pos]);
  
  if (CTAGGING(4,LOCALCFG_MASK)) CL << "Masking '" << np.name << "' with (" << delta << ")\n";
  np.add_metascores (mask_metascore_idx, delta);
}

void Pair_local_PCFG::get_mask_gff (GFF_list& gff_list, const vector<Metaprob>& expected_metacounts, Prob min_prob, const char* seqname, const char* seq, const char* source, const char* feature) const
{
  gff_list.acquire_mask (expected_metacounts[mask_metascore_idx], min_prob, seqname, seq, source, feature, 0);
}

Odds_PCFG::Odds_PCFG() : Pair_PCFG (0)
{ }

Odds_PCFG::Odds_PCFG (int states, Alphabet_group null_emit, Boolean_group null_extend) :
  Pair_PCFG (states),
  null_emit (null_emit),
  null_extend (null_extend)
{ }

Odds_PCFG::Odds_PCFG (int states, PScope& pscope) :
  Pair_PCFG (states),
  null_emit (pscope.new_alphabet_group (CFG_alphabet, 1, "nullEmit")),
  null_extend (pscope.new_boolean_group ("nullExtend"))
{ }

Trainable_PCFG::Trainable_PCFG (PScores& pscore, int states)
  : Odds_PCFG (states, pscore),
    first_pgroup_idx (pscore.pgroups()),
    pscore (pscore)
{ }

void Trainable_PCFG::init_my_pgroups()
{
  for (int g = first_pgroup_idx; g < pscore.pgroups(); ++g)
    my_pgroups.push_back (PGroup (g, pscore.group_size (g)));
}

Dirichlet_prior Trainable_PCFG::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  Dirichlet_prior prior (pscore);
  for_const_contents (vector<PGroup>, my_pgroups, pg)
    prior.assign_Laplace (*pg, k);
  return prior;
}
