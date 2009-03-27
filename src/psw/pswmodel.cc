#include "psw/pswmodel.h"

PSW_vars::PSW_vars() { }
PSW_vars::PSW_vars (PScope& scope_ref, const Alphabet& alphabet_ref, int time_bins)
{
  init (scope_ref, alphabet_ref, time_bins);
}

void PSW_vars::init (PScope& scope_ref, const Alphabet& alphabet_ref, int tb)
{
  // reset
  scope = &scope_ref;
  alphabet = &alphabet_ref;
  time_bins = tb;
  a = b = c = d = e = f = g = h = vector<Boolean_group> (time_bins);
  subst = vector<vector<Alphabet_group> > (time_bins, vector<Alphabet_group> (alphabet->size()));
  // initialise
  aroot = scope->new_boolean_group ("SI");
  croot = scope->new_boolean_group ("II");
  emit = scope->new_alphabet_group (*alphabet, "emit");
  for (int t = 0; t < time_bins; ++t)
    {
      a[t] = new_boolean_group ("SI", t);
      b[t] = new_boolean_group ("SD", t);
      c[t] = new_boolean_group ("II", t);
      d[t] = new_boolean_group ("ID", t);
      e[t] = new_boolean_group ("MI", t);
      f[t] = new_boolean_group ("MD", t);
      g[t] = new_boolean_group ("DI", t);
      h[t] = new_boolean_group ("DD", t);
      for (int c = 0; c < alphabet->size(); ++c)
	{
	  sstring src_char ("subst");
	  src_char << alphabet->int2char_uc(c);
	  subst[t][c] = new_alphabet_group (src_char.c_str(), t);
	}
    }
}

Boolean_group PSW_vars::new_boolean_group (const char* prefix, int time)
{
  sstring tmp (prefix);
  tmp << (int) time;
  return scope->new_boolean_group (tmp.c_str());
}

Alphabet_group PSW_vars::new_alphabet_group (const char* prefix, int time)
{
  sstring tmp (prefix);
  tmp << (int) time;
  return scope->new_alphabet_group (*alphabet, tmp.c_str());
}

Zoned_PSW_vars::Zoned_PSW_vars (PScope& scope, const Alphabet& alphabet, int zones, int time_bins)
  : scope (scope), alphabet (alphabet), zones (zones), time_bins (time_bins),
    zone_vars (zones),
    zone_trans (zones + 1)
{
  vector<sstring> suffix (zones);
  for (int z = 0; z < zones; ++z)
    {
      zone_vars[z].init (scope, alphabet, time_bins);
      suffix[z] << (int) z;
    }
  zone_start = scope.new_group (zones + 1, "Z", suffix);
  for (int z = 0; z < zones; ++z)
    {
      sstring name;
      name << "Z" << z << "Z";
      zone_trans[z] = scope.new_group (zones + 1, name.c_str(), suffix);
    }
}

PSW_PHMM::PSW_PHMM (const Zoned_PSW_vars& vars, int time, bool allow_null_cycles)
  : Pair_PHMM (total_states(vars.zones), vars.alphabet),
    vars (vars), time (time), zones (vars.zones)
{
  const int A = alphabet->size();
  // set state types
  for (int z = 0; z < zones; ++z)
    {
      state_type[SS(z)] = Null;
      state_type[EE(z)] = Null;
      state_type[SI(z)] = EmitY;
      state_type[II(z)] = EmitY;
      state_type[IM(z)] = EmitXY;
      state_type[ID(z)] = EmitX;
    }
  // set zone emissions
  for (int z = 0; z < zones; ++z)
    {
      // initialise vectors/arrays
      single_emit[SI(z)] = vector<PFunc> (A);
      single_emit[II(z)] = vector<PFunc> (A);
      single_emit[ID(z)] = vector<PFunc> (A);
      pair_emit[IM(z)] = array2d<PFunc> (A, A);
      // fill vectors/arrays
      for (int c = 0; c < A; ++c)
	{
	  const PVar& emit_c = vars.zone_vars[z].emit[c];
	  single_emit[SI(z)][c] = single_emit[II(z)][c] = single_emit[ID(z)][c] = emit_c;
	  for (int d = 0; d < A; ++d)
	    pair_emit[IM(z)](c,d) = emit_c * vars.zone_vars[z].subst[time][c][d];
	}
    }
  // set intra-zone transitions
  for (int z = 0; z < zones; ++z)
    {
      // make local copies of relevant transducer vars
      const Boolean_group& a = vars.zone_vars[z].aroot;
      const Boolean_group& c = vars.zone_vars[z].croot;
      const Boolean_group& at = vars.zone_vars[z].a[time];
      const Boolean_group& bt = vars.zone_vars[z].b[time];
      const Boolean_group& ct = vars.zone_vars[z].c[time];
      const Boolean_group& dt = vars.zone_vars[z].d[time];
      const Boolean_group& et = vars.zone_vars[z].e[time];
      const Boolean_group& ft = vars.zone_vars[z].f[time];
      const Boolean_group& gt = vars.zone_vars[z].g[time];
      const Boolean_group& ht = vars.zone_vars[z].h[time];
      // transitions from SS
      if (allow_null_cycles)
	transition(SS(z),EE(z)) = at.n * a.n;
      transition(SS(z),SI(z)) = at.y;
      transition(SS(z),IM(z)) = at.n * a.y * bt.n;
      transition(SS(z),ID(z)) = at.n * a.y * bt.y;
      // transitions from SI
      transition(SI(z),EE(z)) = ct.n * a.n;
      transition(SI(z),II(z)) = ct.y;
      transition(SI(z),IM(z)) = ct.y * a.y * dt.n;
      transition(SI(z),ID(z)) = ct.y * a.y * dt.y;
      // transitions from II
      transition(II(z),EE(z)) = ct.n * c.n;
      transition(II(z),II(z)) = ct.y;
      transition(II(z),IM(z)) = ct.n * c.y * dt.n;
      transition(II(z),ID(z)) = ct.n * c.y * dt.y;
      // transitions from IM
      transition(IM(z),EE(z)) = et.n * c.n;
      transition(IM(z),II(z)) = et.y;
      transition(IM(z),IM(z)) = et.n * c.y * ft.n;
      transition(IM(z),ID(z)) = et.n * c.y * ft.y;
      // transitions from ID
      transition(ID(z),EE(z)) = gt.n * c.n;
      transition(ID(z),II(z)) = gt.y;
      transition(ID(z),IM(z)) = gt.n * c.y * ht.n;
      transition(ID(z),ID(z)) = gt.n * c.y * ht.y;
    }
  // set inter-zone transitions
  transition(Start,End) = vars.zone_start[zones];
  for (int zi = 0; zi < zones; ++zi)
    {
      transition(Start,SS(zi)) = vars.zone_start[zi];
      transition(EE(zi),End) = vars.zone_trans[zi][zones+1];
      for (int zj = 0; zj < zones; ++zj)
	transition(EE(zi),SS(zj)) = vars.zone_trans[zi][zj];
    }
}
