#include "handel/proghmm.h"

// Prog_HMM_base

Prog_HMM_base::Prog_HMM_base (int states, const Alphabet& alphabet)
  : Pair_PHMM (states, alphabet),
    pscore(),
    prior (pscore)
{ }

Prog_HMM_base::~Prog_HMM_base()
{ }

void Prog_HMM_base::assign_Laplace_prior()
{
  for (int g = 0; g < pscore.pgroups(); ++g)
    {
      const PGroup pg (g, pscore.group_size (g));
      prior.assign_Laplace (pg);
    }
}

// Prog_HMM

Prog_HMM::Prog_HMM (int indelTypes, int zones, const Alphabet& alphabet)
  : Prog_HMM_base (calcTotalStates (indelTypes, zones), alphabet),
    indelTypes (indelTypes),
    zones (zones),
    stop_pg (zones),
    indel_begin_pg (zones),
    indel_type_pg (zones),
    indel_extend_pg (zones, vector<Boolean_group> (indelTypes)),
    indel_swap_pg (zones, vector<Boolean_group> (indelTypes)),
    single_pg (zones),
    sub_pg (zones, vector<Alphabet_group> (alphabet.size())),
    start_zone_pg (pscore.new_group (zones + 1, "startZone")),
    next_zone_pg (zones)
{
  // create PGroups
  for (int z = 0; z < zones; ++z)
    {
      sstring stop_name, indel_begin_name, indel_type_name, indel_extend_name, indel_swap_name, single_name, sub_name, next_zone_name;
      stop_name << "zone" << z << "stop";
      indel_begin_name << "zone" << z << "indelBegin";
      indel_type_name << "zone" << z << "indelType";
      single_name << "zone" << z << "single";
      next_zone_name << "zone" << z << "next";

      stop_pg[z] = pscore.new_boolean_group (stop_name.c_str());
      indel_begin_pg[z] = pscore.new_boolean_group (indel_begin_name.c_str());
      indel_type_pg[z] = pscore.new_group (indelTypes, indel_type_name.c_str());
      single_pg[z] = pscore.new_alphabet_group (alphabet, single_name.c_str());
      next_zone_pg[z] = pscore.new_group (zones, next_zone_name.c_str());

      for (int i = 0; i < indelTypes; ++i)
	{
	  sstring extend_name, swap_name;
	  extend_name << "zone" << z << "indel" << i << "Extend";
	  swap_name << "zone" << z << "indel" << i << "Swap";
	  indel_extend_pg[z][i] = pscore.new_boolean_group (extend_name.c_str());
	  indel_swap_pg[z][i] = pscore.new_boolean_group (swap_name.c_str());
	}
      for (int c = 0; c < alphabet.size(); ++c)
	{
	  sstring s;
	  s << "zone" << z << "sub" << alphabet.int2char_uc(c);
	  sub_pg[z][c] = pscore.new_alphabet_group (alphabet, s.c_str());
	}
    }

  // set state types
  for (int z = 0; z < zones; ++z)
    {
      state_type[ZoneStart(z)] = Null;
      state_type[ZoneEnd(z)] = Null;
      init_emit (Match(z), EmitXY, PFunc());
      for (int i = 0; i < indelTypes; ++i)
	{
	  init_emit (Ins(z,i), EmitY, PFunc());
	  init_emit (Del(z,i), EmitX, PFunc());
	}
    }

  // set transitions
  transition (Start, End) = start_zone_pg[zones];

  for (int z = 0; z < zones; ++z)
    {
      // intra-zone & start/end transitions for zone z
      // from Start to ZoneStart(z)
      transition (Start, ZoneStart(z)) = start_zone_pg[z];

      // from ZoneStart(z) to Match(z)
      transition (ZoneStart(z), Match(z)) = indel_begin_pg[z].n;

      // from Match(z)
      transition (Match(z), Match(z))   = stop_pg[z].n * indel_begin_pg[z].n;
      transition (Match(z), ZoneEnd(z)) = stop_pg[z].y;
      
      // from ZoneEnd(z) to End
      transition (ZoneEnd(z), End) = next_zone_pg[z][z];

      // transitions involving Ins(z,i) and Del(z,i)
      for (int i = 0; i < indelTypes; ++i)
	{
	  const PFunc indel_begin  = indel_begin_pg[z].y * indel_type_pg[z][i] / 2;
	  const PFunc indel_mbegin = stop_pg[z].n * indel_begin;
	  const PFunc indel_extend = indel_extend_pg[z][i].y;
	  const PFunc indel_swap   = indel_extend_pg[z][i].n * indel_swap_pg[z][i].y;
	  const PFunc indel_match  = indel_extend_pg[z][i].n * indel_swap_pg[z][i].n * stop_pg[z].n;
	  const PFunc indel_end    = indel_extend_pg[z][i].n * indel_swap_pg[z][i].n * stop_pg[z].y;

	  // transitions to Ins(z,i) and Del(z,i)...
	  // ...from ZoneStart(z)
	  transition (ZoneStart(z), Ins(z,i)) = indel_begin;
	  transition (ZoneStart(z), Del(z,i)) = indel_begin;

	  // ...from Match(z)
	  transition (Match(z), Ins(z,i)) = indel_mbegin;
	  transition (Match(z), Del(z,i)) = indel_mbegin;

	  // transitions from Ins(z,i) and Del(z,i)
	  transition (Ins(z,i), Ins(z,i))   = indel_extend;
	  transition (Ins(z,i), Del(z,i))   = indel_swap;
	  transition (Ins(z,i), Match(z))   = indel_match;
	  transition (Ins(z,i), ZoneEnd(z)) = indel_end;

	  transition (Del(z,i), Del(z,i))   = indel_extend;
	  transition (Del(z,i), Ins(z,i))   = indel_swap;
	  transition (Del(z,i), Match(z))   = indel_match;
	  transition (Del(z,i), ZoneEnd(z)) = indel_end;
	}

      // inter-zone transitions from zone z-->z2
      for (int z2 = 0; z2 < zones; ++z2)
	if (z2 != z)
	  transition (ZoneEnd(z), ZoneStart(z2)) = next_zone_pg[z][z2];
    }
  
  // set emissions
  for (int z = 0; z < zones; ++z)
    for (int c = 0; c < alphabet.size(); ++c)
      {
	for (int i = 0; i < indelTypes; ++i)
	  {
	    single_emit[Ins(z,i)][c] = single_pg[z][c];
	    single_emit[Del(z,i)][c] = single_pg[z][c];
	  }
	for (int d = 0; d < alphabet.size(); ++d)
	  pair_emit[Match(z)](c,d) = (single_pg[z][c] * sub_pg[z][c][d] + single_pg[z][d] * sub_pg[z][d][c]) / 2;
      }

  // Set up the prior: tapered for start_zone_pg, next_zone_pg & indel_type_pg,
  //  increased weight on { sub_pg[z][c][c], indel_begin_pg[z][c].n } for lower z,
  //  increased weight on indel_extend_pg[z][i].y for higher i, flat otherwise.
  // This means that the lower zones are favoured & inherently prefer to have more matches & fewer indels,
  //  while lower indel types are also favoured & prefer to be shorter (but all these preferences are weak).
  assign_Laplace_prior();  // start with a flat prior for all PGroups by default
  const double zone_taper = .1;
  const double indel_taper = .1;
  const double indel_begin_weight_inc = .1;
  const double indel_extend_weight_inc = .1;
  const double match_weight_inc = .1;
  double indel_begin_weight = 0.;
  double match_weight = 0.;
  prior.assign_Laplace (start_zone_pg, 1., zone_taper);
  for (int z = zones - 1; z >= 0; --z)
    {
      prior.assign_Laplace (next_zone_pg[z], 1., zone_taper);
      prior.assign_Laplace (indel_type_pg[z], 1., indel_taper);

      const Multigroup indel_begin_mg (indel_begin_pg[z]);
      prior.assignment[indel_begin_mg].mixture.alpha[0][0][0] += indel_begin_weight;

      double indel_extend_weight = 0.;
      for (int i = 0; i < indelTypes; ++i)
	{
	  const Multigroup indel_extend_mg (indel_extend_pg[z][i]);
	  prior.assignment[indel_extend_mg].mixture.alpha[0][0][1] += indel_extend_weight;
	  indel_extend_weight += indel_extend_weight_inc;
	}

      for (int c = 0; c < alphabet.size(); ++c)
	{
	  const Multigroup sub_mg (sub_pg[z][c]);
	  prior.assignment[sub_mg].mixture.alpha[0][0][c] += match_weight;
	}

      indel_begin_weight += indel_begin_weight_inc;
      match_weight += match_weight_inc;
    }
}
