#include "scfg/cfgsexpr.h"
#include "scfg/cfgkeywords.h"
#include "seq/pkeywords.h"

// pseudocount for all Dirichlet priors
#define GRAMSET_DIRICHLET_PSEUDOCOUNT 0.1

void Gramset::clear()
{
  single_scfg_pscores = PScores();
  single_scfg_prior.clear();
  single_scfg_mutable_pgroups.clear();
  phmm_cfg_map.clear();
}

void PCFG_builder::init_null_pgroups (SExpr& model_sexpr, PScores& pscores, SymPVar& sym2pvar, Alphabet_group& null_emit, Boolean_group& null_extend, bool want_HMM)
{
  const Alphabet& alph = CFG_alphabet;

  if (SExpr* null_emit_sexpr = model_sexpr.find (CFG_NULL_EMIT, 1))
    {
      PGroup tmp_null_emit = init_pgroup (pscores, sym2pvar, null_emit_sexpr->value(), (set<int>*) 0, false, true, true);
      if (tmp_null_emit.group_size != alph.size())
	THROWEXPR ("In (" << *null_emit_sexpr << "):\nNull emit pgroup needs to be have " << alph.size() << " paramters (one for each alphabet token)");
      null_emit = Alphabet_group (tmp_null_emit.group_idx, alph);
    }
  else
    null_emit = pscores.new_alphabet_group (alph, CFG_NULL_EMIT);

  if (!want_HMM)
    {
      if (SExpr* null_extend_sexpr = model_sexpr.find (CFG_NULL_EXTEND, 1))
	{
	  PGroup tmp_null_extend = init_pgroup (pscores, sym2pvar, null_extend_sexpr->value(), (set<int>*) 0, false, true, true);
	  if (tmp_null_extend.group_size != 2)
	    THROWEXPR ("In (" << *null_extend_sexpr << "):\nNull extend pgroup needs to be boolean");
	  null_extend = Boolean_group (tmp_null_extend.group_idx);
	}
      else
	null_extend = pscores.new_boolean_group (CFG_NULL_EXTEND);
    }
}

Odds_PCFG PCFG_builder::init_grammar (SExpr& model_sexpr, PScores& pscores, set<int>& pad_states_ret, set<int>& mutable_pgroups_ret, bool want_single, bool want_HMM)
{
  // initialize PVars
  SymPVar sym2pvar;

  Alphabet_group null_emit;
  Boolean_group null_extend;
  init_null_pgroups (model_sexpr, pscores, sym2pvar, null_emit, null_extend, want_HMM);
  init_pgroups (pscores, sym2pvar, model_sexpr, CFG_CONST_BGROUP, (set<int>*) 0, false, true, true);
  init_pgroups (pscores, sym2pvar, model_sexpr, CFG_BGROUP, &mutable_pgroups_ret, false, true, true);

  // create type lookup
  map<sstring,State_type> type_lookup;
  type_lookup[sstring(CFG_TYPE_NULL)] = Null;
  type_lookup[sstring(CFG_TYPE_LEFT)] = EmitXL;
  type_lookup[sstring(CFG_TYPE_RIGHT)] = EmitXR;
  type_lookup[sstring(CFG_TYPE_PAIR)] = EmitXLR;

  type_lookup[sstring(CFG_TYPE_L)] = EmitXL;
  type_lookup[sstring(CFG_TYPE_R)] = EmitXR;
  type_lookup[sstring(CFG_TYPE_LR)] = EmitXLR;

  // initialize SCFG state descriptors
  SymIndex name2state ("State");
  name2state[sstring(CFG_START)] = Start;
  name2state[sstring(CFG_END)] = End;
  vector<State_type> state_type;
  vector<sstring> state_name;
  vector<SExpr*> sexpr_by_state;

  // loop through states
  const vector<SExpr*> state_sexprs = model_sexpr.find_all (CFG_STATE, 1);
  for_const_contents (vector<SExpr*>, state_sexprs, state_sexpr)
    {
      const int new_state_index = state_name.size();

      const sstring name = (**state_sexpr) (CFG_NAME).get_atom();
      if (name2state.find (name) != name2state.end())
	THROWEXPR ("In (" << model_sexpr << "):\nDuplicate state name " << name);

      if ((**state_sexpr).find (CFG_PAD))
	pad_states_ret.insert (new_state_index);

      State_type t = Null;
      SExpr& type_sexpr = (**state_sexpr).find_or_die (CFG_TYPE, 1);
      if (type_sexpr.find (CFG_BIFURCATION))
	t = Bifurc;
      else if (!type_sexpr.find (CFG_TYPE_NULL))
	{
	  SExpr *x_sexpr, *y_sexpr;
	  if ((x_sexpr = type_sexpr.find (CFG_X, 1)))
	    for_const_contents (list<SExpr>, x_sexpr->child, xtype_sexpr)
	      t = (State_type) (((int) t) | ((int) type_lookup[xtype_sexpr->get_atom()]));
	  if ((y_sexpr = type_sexpr.find (CFG_Y, 1)))
	    for_const_contents (list<SExpr>, y_sexpr->child, ytype_sexpr)
	      t = (State_type) (((int) t) | (((int) type_lookup[ytype_sexpr->get_atom()]) << ShiftYL));
	  if (!x_sexpr && !y_sexpr)
	    {
	      int lr_type = 0, match_type = 0;

	      if (want_HMM)
		lr_type = EmitXLR;
	      else
		lr_type =
		  type_sexpr.find(CFG_TYPE_LEFT)
		  ? EmitXL
		  : (type_sexpr.find(CFG_TYPE_RIGHT)
		     ? EmitXR
		     : (type_sexpr.find(CFG_TYPE_PAIR)
			? EmitXLR
			: -1));

	      if (want_single)
		match_type = (1 << ShiftXL);
	      else
		match_type = 
		  type_sexpr.find(CFG_TYPE_DELETE)
		  ? (1 << ShiftXL)
		  : (type_sexpr.find(CFG_TYPE_INSERT)
		     ? (1 << ShiftYL)
		     : (type_sexpr.find(CFG_TYPE_MATCH)
			? ((1 << ShiftXL) | (1 << ShiftYL))
			: -1));

	      sstring err;
	      if (lr_type < 0)
		err << "Is this state a l, r, lr or null type?\n";

	      if (match_type < 0)
		err << "Is this state a match, insert or delete type?\n";

	      if (err.size())
		THROWEXPR ("In (" << type_sexpr << "):" << err);

	      t = (State_type) (match_type * lr_type);
	    }
	}

      if (want_single && (t & EmitYLR))
	THROWEXPR ("In (" << type_sexpr << "):\nExpected a single-sequence model");

      if (want_HMM && (t & (EmitXR | EmitYR)))
	THROWEXPR ("In (" << type_sexpr << "):\nExpected a left-regular model, i.e. an HMM, not an SCFG");

      name2state[name] = new_state_index;
      state_name.push_back (name);
      state_type.push_back (t);
      sexpr_by_state.push_back (*state_sexpr);
    }

  // create the Odds_PCFG
  const int states = state_name.size();
  Odds_PCFG pcfg (states, null_emit, null_extend);
  if (SExpr* name_sexpr = model_sexpr.find (CFG_NAME))
    pcfg.name = name_sexpr->value().get_atom();
  else
    pcfg.name = sstring (want_HMM ? CFG_PAIR_HMM : (want_single ? CFG_SINGLE_SCFG : CFG_PAIR_SCFG));
  pcfg.state_name = state_name;
  pcfg.state_type = state_type;

  // initialize emissions & bifurcations
  const PFunc zero_pfunc (0.);
  for (int s = 0; s < states; ++s)
    {
      SExpr& ss = *sexpr_by_state[s];
      const State_type t = state_type[s];
      if (is_emit_type (t))
	{
	  pcfg.init_emit (s, t, zero_pfunc);
	  const vector<SExpr*> emit_sexprs = ss.find_all (CFG_EMIT, 1);
	  for_const_contents (vector<SExpr*>, emit_sexprs, emit_sexpr)
	    {
	      sstring chars;
	      if (t & EmitXLR)
		chars << token_list_to_string ((**emit_sexpr).find_or_die (CFG_X, 1), 1);
	      if (t & EmitYLR)
		chars << token_list_to_string ((**emit_sexpr).find_or_die (CFG_Y, 1), 1);

	      PFunc f = init_pfunc (sym2pvar, (**emit_sexpr).find_or_die(CFG_PROB,1), 1);

	      pcfg.emit_by_chars (s, chars.c_str()) = f;
	    }
	}
      else if (is_bifurc_type (t))
	{
	  SExpr& sb = ss.find_or_die(CFG_TYPE).find_or_die(CFG_BIFURCATION);
	  pcfg.init_bifurc (s,
			    name2state (sb(CFG_BIFURC_LEFT).get_atom(), sb),
			    name2state (sb(CFG_BIFURC_RIGHT).get_atom(), sb));
	}
    }

  // initialize transitions
  const vector<SExpr*> trans_sexprs = model_sexpr.find_all (CFG_TRANSITION, 1);
  for_const_contents (vector<SExpr*>, trans_sexprs, trans_sexpr)
    {
      PFunc f = init_pfunc (sym2pvar, (**trans_sexpr).find_or_die(CFG_PROB,1), 1);

      pcfg.transition (name2state ((**trans_sexpr)(CFG_FROM).get_atom(), **trans_sexpr),
		       name2state ((**trans_sexpr)(CFG_TO).get_atom(), **trans_sexpr))
	= f;
    }

  // return
  return pcfg;
}

Pair_PHMM PCFG_builder::init_hmm (SExpr& model_sexpr, PScores& pscores, Alphabet_group& null_emit_group_ret, set<int>& pad_states_ret, set<int>& mutable_pgroups_ret, bool want_single)
{
  Odds_PCFG pcfg = init_grammar (model_sexpr, pscores, pad_states_ret, mutable_pgroups_ret, want_single, true);
  Pair_HMM<PFunc> pfunc_hmm (pcfg.as_HMM());
  Pair_PHMM phmm (pfunc_hmm);
  null_emit_group_ret = pcfg.null_emit;
  return phmm;
}

void PCFG_builder::grammar2stream (ostream& out, const Pair_PCFG& pcfg, const PScores& pscores, const set<int>& mutable_pgroups, const char* tag, const char* indent, int null_emit_group_index, int null_extend_group_index, const set<int>* pad_states)
{
  const Alphabet& alph (CFG_alphabet);

  sstring ind (indent);
  ind << ' ';

  // model, name
  out << indent << '(' << SExpr_atom(tag) << '\n';
  out << ind << '(' << CFG_NAME << ' ' << SExpr_atom(pcfg.name) << ")\n";

  // params
  if (null_emit_group_index >= 0)
    {
      const vector<int> gvec (1, (int) null_emit_group_index);
      pscores2stream (out, pscores, CFG_NULL_EMIT, gvec, ind.c_str(), (PCounts*) 0, true);
    }

  if (null_extend_group_index >= 0)
    {
      const vector<int> gvec (1, (int) null_extend_group_index);
      pscores2stream (out, pscores, CFG_NULL_EXTEND, gvec, ind.c_str(), (PCounts*) 0, true);
    }

  vector<int> remaining_mutable_pgroups, remaining_const_pgroups;
  for (int g = 0; g < pscores.groups(); ++g)
    if (g != null_emit_group_index && g != null_extend_group_index)
      {
	if (mutable_pgroups.find (g) == mutable_pgroups.end())
	  remaining_const_pgroups.push_back (g);
	else
	  remaining_mutable_pgroups.push_back (g);
      }

  pscores2stream (out, pscores, CFG_CONST_BGROUP, remaining_const_pgroups, ind.c_str(), (PCounts*) 0, true);
  pscores2stream (out, pscores, CFG_BGROUP, remaining_mutable_pgroups, ind.c_str(), (PCounts*) 0, true);

  // states
  for (int s = 0; s < pcfg.states(); ++s)
    {
      // state, name
      out << ind << '(' << CFG_STATE << '\n';
      out << ind << " (" << CFG_NAME << ' ' << SExpr_atom(pcfg.state_name[s]) << ")\n";

      // padding?
      if (pad_states)
	if (pad_states->find(s) != pad_states->end())
	  out << ind << " (" << CFG_PAD << ")\n";

      // type
      const State_type t = pcfg.state_type[s];
      out << ind << " (" << CFG_TYPE;
      if (t == Bifurc)  // bifurcations
	out << " (" << CFG_BIFURCATION
	    << " (" << CFG_BIFURC_LEFT << ' ' << SExpr_atom(pcfg.state_name[pcfg.bifurc[s].l])
	    << ") (" << CFG_BIFURC_RIGHT << ' ' << SExpr_atom(pcfg.state_name[pcfg.bifurc[s].r]) << "))";

      else if (t == Null)
	out << " (" << CFG_TYPE_NULL << ")";

      else if (pcfg.is_emit_type (t))
	{
	  if (t & EmitXLR)
	    {
	      vector<sstring> s;
	      if (t & EmitXL)
		s.push_back (sstring (CFG_TYPE_L));
	      if (t & EmitXR)
		s.push_back (sstring (CFG_TYPE_R));
	      out << " (" << CFG_X << " " << s << ")";
	    }

	  if (t & EmitYLR)
	    {
	      vector<sstring> s;
	      if (t & EmitYL)
		s.push_back (sstring (CFG_TYPE_L));
	      if (t & EmitYR)
		s.push_back (sstring (CFG_TYPE_R));
	      out << " (" << CFG_Y << " " << s << ")";
	    }
	}

      else
	THROWEXPR ("Don't know how to encode state type " << t);

      out << ")\n";

      // emissions
      if (pcfg.is_emit_type (t))
	for (int e = 0; e < (int) pcfg.emit[s].size(); ++e)
	  {
	    const PFunc& f = pcfg.emit[s][e];
	    if (!f.is_null() && !f.is_zero())
	      {
		out << ind << " (" << CFG_EMIT;
		if (t & EmitXLR)
		  {
		    vector<char> s;
		    if (t & EmitXL)
		      s.push_back (alph.int2char (pcfg.emit_xl_idx (t, e)));
		    if (t & EmitXR)
		      s.push_back (alph.int2char (pcfg.emit_xr_idx (t, e)));
		    out << " (" << CFG_X << " " << s << ")";
		  }

		if (t & EmitYLR)
		  {
		    vector<char> s;
		    if (t & EmitYL)
		      s.push_back (alph.int2char (pcfg.emit_yl_idx (t, e)));
		    if (t & EmitYR)
		      s.push_back (alph.int2char (pcfg.emit_yr_idx (t, e)));
		    out << " (" << CFG_Y << " " << s << ")";
		  }

		out << " (" << CFG_PROB << ' ';
		pfunc2stream (out, pscores, f);
		out << "))\n";
	      }
	  }

      out << ind << ")  ; end " << CFG_STATE << ' ' << pcfg.state_name[s] << '\n';
    }

  // transitions
  const sstring start_string (CFG_START), end_string (CFG_END);
  for (int src = Start; src < pcfg.states(); ++src)
    if (src != End)
      for (int dest = End; dest < pcfg.states(); ++dest)
	if (dest != Start)
	  {
	    const PFunc& f = pcfg.transition (src, dest);
	    if (!f.is_null() && !f.is_zero())
	      {
		out << ind
		    << "(" << CFG_TRANSITION
		    << " (" << CFG_FROM << ' ' << SExpr_atom(src == Start ? start_string : pcfg.state_name[src])
		    << ") (" << CFG_TO << ' ' << SExpr_atom(dest == End ? end_string : pcfg.state_name[dest])
		    << ") (" << CFG_PROB << ' ';
		pfunc2stream (out, pscores, f);
		out << "))\n";
	      }
	  }

  // all done
  out << indent << ")  ;; end " << tag << ' ' << pcfg.name << '\n';
}


void PCFG_builder::hmm2stream (ostream& out, const Pair_PHMM& phmm, const PScores& pscores, const set<int>& mutable_pgroups, const char* tag, const char* indent, int null_emit_group_idx, const set<int>* pad_states)
{
  Pair_PCFG pcfg (phmm);
  grammar2stream (out, pcfg, pscores, mutable_pgroups, tag, indent, null_emit_group_idx, -1, pad_states);
}

void PCFG_builder::init_gramset (SExpr& gramset_sexpr, Gramset& gramset)
{
  gramset.clear();

  set<int> dummy_pad_states;
  gramset.single_scfg = init_grammar (gramset_sexpr.find_or_die (CFG_SINGLE_SCFG), gramset.single_scfg_pscores, dummy_pad_states, gramset.single_scfg_mutable_pgroups, true, false);
  gramset.single_scfg_prior.assign_Laplace (gramset.single_scfg_pscores, gramset.single_scfg_mutable_pgroups, GRAMSET_DIRICHLET_PSEUDOCOUNT);

  gramset.phmm_cfg_map.clear();
  const vector<SExpr*> phmm_cfg_sexprs = gramset_sexpr.find_all (CFG_PARAM_SET, 1);
  if (phmm_cfg_sexprs.empty())
    THROWEXPR (CFG_GRAMSET << " file needs at least one " << CFG_PARAM_SET);
  for_const_contents (vector<SExpr*>, phmm_cfg_sexprs, phmm_cfg_sexpr)
    {
      PHMM_CFG phmm_cfg;

      phmm_cfg.hmm = init_hmm ((**phmm_cfg_sexpr).find_or_die (CFG_PAIR_HMM), phmm_cfg.hmm_pscores, phmm_cfg.hmm_null_emit, phmm_cfg.hmm_pad_states, phmm_cfg.hmm_mutable_pgroups, false);
      phmm_cfg.hmm_prior.assign_Laplace (phmm_cfg.hmm_pscores, phmm_cfg.hmm_mutable_pgroups, GRAMSET_DIRICHLET_PSEUDOCOUNT);

      phmm_cfg.cfg = init_grammar ((**phmm_cfg_sexpr).find_or_die (CFG_PAIR_SCFG), phmm_cfg.cfg_pscores, dummy_pad_states, phmm_cfg.cfg_mutable_pgroups, false, false);
      phmm_cfg.cfg_prior.assign_Laplace (phmm_cfg.cfg_pscores, phmm_cfg.cfg_mutable_pgroups, GRAMSET_DIRICHLET_PSEUDOCOUNT);

      const sstring& name = (**phmm_cfg_sexpr) (CFG_NAME).get_atom();
      gramset.phmm_cfg_map[name] = phmm_cfg;
    }
}

void PCFG_builder::gramset2stream (ostream& out, const Gramset& gramset, const char* tag)
{
  out << "(" << tag << '\n';
  grammar2stream (out, gramset.single_scfg, gramset.single_scfg_pscores, gramset.single_scfg_mutable_pgroups, CFG_SINGLE_SCFG, " ", gramset.single_scfg.null_emit.group_idx, gramset.single_scfg.null_extend.group_idx);
  for_const_contents (PHMM_CFG_map, gramset.phmm_cfg_map, name_phmm_cfg)
    {
      out << " (" << CFG_PARAM_SET
	  << "\n  (" << CFG_NAME << ' ' << SExpr_atom(name_phmm_cfg->first) << ")\n";
      hmm2stream (out, name_phmm_cfg->second.hmm, name_phmm_cfg->second.hmm_pscores, name_phmm_cfg->second.hmm_mutable_pgroups, CFG_PAIR_HMM, "  ", name_phmm_cfg->second.hmm_null_emit.group_idx, &name_phmm_cfg->second.hmm_pad_states);
      grammar2stream (out, name_phmm_cfg->second.cfg, name_phmm_cfg->second.cfg_pscores, name_phmm_cfg->second.cfg_mutable_pgroups, CFG_PAIR_SCFG, "  ", name_phmm_cfg->second.cfg.null_emit.group_idx, name_phmm_cfg->second.cfg.null_extend.group_idx);
      out << " )  ; end " << CFG_PARAM_SET << ' ' << name_phmm_cfg->first << '\n';
    }
  out << ")\n";
}
