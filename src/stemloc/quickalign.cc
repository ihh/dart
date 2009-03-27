#include "stemloc/quickalign.h"
#include "scfg/paircfg.h"
#include "stemloc/slkeywords.h"

Quick_align::Quick_align (PScores& ps, const Alphabet_group& ne)
  : Pair_PHMM (TotalStates, CFG_alphabet),
    pscore (&ps),
    is_global (pscore->new_boolean_group()),
    mat_pg (pscore->new_group (3, "match")),
    ins_pg (pscore->new_group (4, "insert")),
    gap_pg (pscore->new_group (3, "gap")),
    gapx_pg (pscore->new_group (2, "gapx")),
    null_emit (&ne),
    pair_nuc (CFG_alphabet.size())
{
  // set model name
  name = "Quick_align";

  // set state names
  state_name[LeftPadX] = "LeftPadX";
  state_name[LeftPadY] = "LeftPadY";
  state_name[Match] = "Match";
  state_name[InsX] = "InsX";
  state_name[InsY] = "InsY";
  state_name[Gap] = "Gap";
  state_name[GapX] = "GapX";
  state_name[GapY] = "GapY";
  state_name[RightPadX] = "RightPadX";
  state_name[RightPadY] = "RightPadY";

  // set is_global PVar names
  pscore->group_name[is_global.group_idx].clear();
  pscore->group_suffix[is_global.group_idx][0] = STEMLOC_LOCAL;
  pscore->group_suffix[is_global.group_idx][1] = STEMLOC_GLOBAL;

  // make global by default
  PVar local = is_global[0], global = is_global[1];
  ps[local] = -InfinityScore;
  ps[global] = 0;

  // add indel pgroups to mutable_pgroups
  mutable_pgroups.insert (mat_pg.group_idx);
  mutable_pgroups.insert (ins_pg.group_idx);
  mutable_pgroups.insert (gap_pg.group_idx);
  mutable_pgroups.insert (gapx_pg.group_idx);

  // create substitution PGroups & add to mutable_pgroups
  const int A = alphabet->size();
  for (int i = 0; i < A; ++i)
    {
      sstring name;
      name << "qsub" << null_emit->index2word(i);
      pair_nuc[i] = pscore->new_alphabet_group (*alphabet, 1, name.c_str());
      mutable_pgroups.insert (pair_nuc[i].group_idx);
    }

  // set state types
  PFunc one (1.);
  init_emit (LeftPadX, EmitX, one);
  init_emit (LeftPadY, EmitY, one);
  init_emit (Match, EmitXY, one);
  init_emit (InsX, EmitX, one);
  init_emit (InsY, EmitY, one);
  init_emit (Gap, EmitXY, one);
  init_emit (GapX, EmitX, one);
  init_emit (GapY, EmitY, one);
  init_emit (RightPadX, EmitX, one);
  init_emit (RightPadY, EmitY, one);

  // set transitions
  const double SmallProb = .00001;

  // local transitions
  //  transition (Start, Match) = local;    /* commented out because this is folded into global transition, below */
  transition (Start, LeftPadX) = local;
  transition (Start, LeftPadY) = local;

  transition (LeftPadX, LeftPadX) = local;
  transition (LeftPadX, LeftPadY) = local;
  transition (LeftPadX, Match) = local;
  transition (LeftPadX, End) = local;

  transition (LeftPadY, LeftPadY) = local;
  transition (LeftPadY, LeftPadX) = local * SmallProb;  // low Y->X probability encourages X-before-Y in padding regions, without disallowing Y-before-X in training alignments
  transition (LeftPadY, Match) = local;
  transition (LeftPadY, End) = local;

  transition (Match, RightPadX) = local;
  transition (Match, RightPadY) = local;  // NB outgoing transitions from Match sum to >1, since right (& left) pad supernormalised

  transition (RightPadX, RightPadX) = local;
  transition (RightPadX, RightPadY) = local;
  transition (RightPadX, End) = local;

  transition (RightPadY, RightPadY) = local;
  transition (RightPadY, RightPadX) = local * SmallProb;  // low Y->X probability encourages X-before-Y in padding regions, without disallowing Y-before-X in training alignments
  transition (RightPadY, End) = local;

  // global transitions
  transition (Start, Match) = local + global * mat_pg[0];
  transition (Start, InsX) = global * mat_pg[1] / 2;
  transition (Start, InsY) = global * mat_pg[1] / 2;
  transition (Start, Gap) = global * mat_pg[2];

  // remaining transitions
  transition (Start, End) = 1;

  transition (Match, Match) = mat_pg[0];
  transition (Match, InsX) = mat_pg[1] / 2;
  transition (Match, InsY) = mat_pg[1] / 2;
  transition (Match, Gap) = mat_pg[2];
  transition (Match, End) = 1;

  transition (InsX, Match) = ins_pg[0];
  transition (InsX, InsX) = ins_pg[1];
  transition (InsX, InsY) = ins_pg[2];
  transition (InsX, End) = ins_pg[3];

  transition (InsY, Match) = ins_pg[0];
  transition (InsY, InsY) = ins_pg[1];
  transition (InsY, InsX) = ins_pg[2];
  transition (InsY, End) = ins_pg[3];

  transition (Gap, Gap) = gap_pg[0];
  transition (Gap, GapX) = gap_pg[1] / 2;
  transition (Gap, GapY) = gap_pg[1] / 2;
  transition (Gap, Match) = gap_pg[2];

  transition (GapX, GapX) = gapx_pg[0];
  transition (GapX, Match) = gapx_pg[1];

  transition (GapY, GapY) = gapx_pg[0];
  transition (GapY, Match) = gapx_pg[1];

  // set emissions
  for (int i = 0; i < A; ++i)
    for (int j = 0; j < A; ++j)
      pair_emit[Match](i,j) = pair_nuc[i][j] / (*null_emit)[j];
}

set<int> Quick_align::non_pad_states() const
{
  set<int> states;
  states.insert (Match);
  states.insert (InsX);
  states.insert (InsY);
  states.insert (Gap);
  states.insert (GapX);
  states.insert (GapY);
  return states;
}

Dirichlet_prior Quick_align::default_prior() const
{
  const double k = .1;  // pseudocount for all Laplace priors
  // create the prior, assign idiotically simple pseudocounts to all PGroup's
  // (not quite Laplace, because k<1; but close)
  Dirichlet_prior prior (*pscore);
  prior.assign (mat_pg, Laplace_prior (mat_pg, k));
  prior.assign (ins_pg, Laplace_prior (ins_pg, k));
  prior.assign (gap_pg, Laplace_prior (gap_pg, k));
  prior.assign (gapx_pg, Laplace_prior (gapx_pg, k));
  for_const_contents (vector<Alphabet_group>, pair_nuc, pn)
    prior.assign (*pn, Laplace_prior (*pn, k));
  return prior;
}
