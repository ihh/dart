#include "seq/psexpr.h"
#include "seq/pkeywords.h"
#include "seq/alignment.h"
#include "util/vector_output.h"

using std::fixed;

PFunc_builder::SymIndex::SymIndex (const char* st)
  : map<sstring,int>(),
    symbol_type (st)
{ }

int PFunc_builder::SymIndex::operator() (const sstring& key, const SExpr& containing_sexpr) const
{
  map<sstring,int>::const_iterator iter = find (key);
  if (iter == end())
    THROWEXPR ("In (" << containing_sexpr << "):\n" << symbol_type << ' ' << key << " not declared");
  return iter->second;
}

void PFunc_builder::alphabet2stream (ostream& out, const Alphabet& alph)
{
  // print alphabet
  out << ";; Alphabet " << alph.name << "\n;;\n";
  out << '(' << PK_ALPHABET << '\n';
  out << " (" << PK_NAME << ' ' << alph.name << ")\n";
  out << " (" << PK_ALPHABET_TOKEN << " (";
  for (int i = 0; i < alph.size(); ++i)
    {
      if (i > 0)
	out << ' ';
      out << alph.int2char (i);
    }
  out << "))\n";
  if (alph.has_complement())
    {
      out << " (" << PK_ALPHABET_COMPLEMENT << " (";
      for (int i = 0; i < alph.size(); ++i)
	{
	  if (i > 0)
	    out << ' ';
	  out << alph.int2char (alph.complement (i));
	}
      out << "))\n";
    }
  const sstring deg (alph.degenerate_chars());
  const char wildcard_char = alph.score2char (alph.flat_score_map (0));
  for (int i = 0; i < (int) deg.size(); ++i)
    if (deg[i] != wildcard_char)
      {
	// assume all degenerate possibilities are equiprobable
	out << " (" << PK_ALPHABET_EXTEND << " (" << PK_TO << ' ' << deg[i] << ')';
	const Symbol_weight_map& swm = alph.char2weight (deg[i]);
	for_const_contents (Symbol_weight_map, swm, sw)
	  out << " (" << PK_FROM << ' ' << alph.int2char (sw->first) << ')';
	out << ")\n";
      }
  out << " (" << PK_ALPHABET_WILDCARD << ' ' << wildcard_char << ")\n";
  out << ")  ;; end alphabet " << alph.name << "\n\n";
}


void PFunc_builder::init_alphabet (Alphabet& alph, SExpr& sexpr)
{
  sstring name, token;

  if (sexpr.find (PK_NAME))
    name = sexpr (PK_NAME).get_atom();

  token = token_list_to_string (sexpr (PK_ALPHABET_TOKEN));
  for_const_contents (sstring, token, c)
    if (Alignment::is_gap_char(*c))
      CLOGERR << "Warning: alphabet character '" << *c << "' is also a gap character and will be treated as a gap!\n";

  alph.reset (name.c_str(), token.size());
  alph.init_chars (token.c_str(), "");

  if (sexpr.find (PK_ALPHABET_COMPLEMENT))
    {
      const vector<sstring> dummy_class_label (1, sstring());
      const sstring comp = token_list_to_string (sexpr (PK_ALPHABET_COMPLEMENT), alph, dummy_class_label);
      alph.reset (name.c_str(), token.size());
      alph.init_chars (token.c_str(), comp.c_str());
    }

  typedef map<char,sstring> DegCharMap;
  sstring deg;
  DegCharMap deg2char;

  bool got_wildcard = false;
  char wildcard = '*';
  if (sexpr.find (PK_ALPHABET_WILDCARD))
    {
      const sstring& tok = sexpr(PK_ALPHABET_WILDCARD).get_atom();
      assert_valid_token (tok, sexpr);
      deg.push_back (tok[0]);
      wildcard = tok[0];
    }

  vector<SExpr*> extend_sexpr = sexpr.find_all (PK_ALPHABET_EXTEND, 1);
  for_const_contents (vector<SExpr*>, extend_sexpr, ex)
    {
      const sstring& to_char = (**ex) (PK_TO).get_atom();
      assert_valid_token (to_char, sexpr);

      vector<SExpr*> from_sexpr = (**ex).find_all (PK_FROM, 1);
      sstring from_chars;
      for_const_contents (vector<SExpr*>, from_sexpr, f)
	{
	  const sstring& from_char = (**f).value().get_atom();
	  assert_valid_token (from_char, alph, sexpr);
	  from_chars.push_back (from_char[0]);
	}

      deg.push_back (to_char[0]);
      deg2char[to_char[0]] = from_chars;
    }
  if (deg.size())
    {
      alph.init_degen_chars (deg.c_str());
      if (got_wildcard)
	alph.set_wildcard (wildcard);
      for_const_contents (DegCharMap, deg2char, dc)
	alph.set_degen_char (dc->first, dc->second.c_str());
      alph.init_degen_complement();
    }
}

int PFunc_builder::get_state_class (const sstring& atom, const vector<sstring>& class_alph)
{
  int c;
  for (c = 0; c < (int) class_alph.size(); ++c)
    if (class_alph[c] == atom)
      break;
  if (c == (int) class_alph.size())
    THROWEXPR ("Label '" << atom << "' is not a valid state class (options are " << class_alph << ")");
  return c;
}

void PFunc_builder::assert_valid_token (const sstring& atom, const SExpr& sexpr)
{
  if (atom.size() != 1)
    THROWEXPR ("In (" << sexpr << "):\nBad alphabet token '" << atom << "' (alphabet tokens must be single characters)");
}

void PFunc_builder::assert_valid_token (const sstring& atom, const Alphabet& alph, const SExpr& sexpr)
{
  assert_valid_token (atom, sexpr);
  if (!alph.contains_strict (atom[0]))
    THROWEXPR ("In (" << sexpr << "):\nAlphabet '" << alph.name << "' does not contain token '" << atom << "'");
}

sstring PFunc_builder::token_list_to_string (SExpr& token_list, int offset)
{
  const vector<sstring> toks = token_list.atoms_to_strings (offset);
  for_const_contents (vector<sstring>, toks, tok)
    assert_valid_token (*tok, token_list);
  return sstring::join (toks, "");
}

sstring PFunc_builder::token_list_to_string (SExpr& token_list, const Alphabet& alph, const vector<sstring>& class_alph, int offset)
{
  vector<sstring> toks = token_list.atoms_to_strings (offset);
  sstring s, state_class;
  if (class_alph.size() > 1)
    {
      state_class = toks.back();
      toks.pop_back();
      get_state_class (state_class, class_alph);  // assert state class is valid
    }
  for_const_contents (vector<sstring>, toks, tok)
    {
      assert_valid_token (*tok, alph, token_list);
      s << *tok;
    }
  s << state_class;
  return s;
}

void PFunc_builder::pscores2stream (ostream& out, const PScores& pscores, const PCounts* pcounts, bool use_bitscores)
{
  set<int> mutable_pgroups;
  for (int g = 0; g < pscores.groups(); ++g)
    mutable_pgroups.insert (g);

  pscores2stream (out, pscores, mutable_pgroups, pcounts, use_bitscores);
}

void PFunc_builder::pscores2stream (ostream& out, const PScores& pscores, const set<int>& mutable_pgroups, const PCounts* pcounts, bool use_bitscores)
{
  vector<vector<int> > pgroup_subset (4);
  for (int g = 0; g < pscores.groups(); ++g)
    {
      const int is_mutable = mutable_pgroups.find (g) == mutable_pgroups.end() ? 0 : 1;
      const int is_rate = pscores.group_size(g) > 1 ? 0 : 1;
      const int subset_index = is_mutable + 2*is_rate;
      pgroup_subset[subset_index].push_back (g);
    }
  for (int is_mutable = 0; is_mutable < 2; ++is_mutable)
    for (int is_rate = 0; is_rate < 2; ++is_rate)
      {
	const int subset_index = is_mutable + 2*is_rate;
	if (pgroup_subset[subset_index].size())
	  {
	    out << '\n';
	    const char* tag = is_rate
	      ? (is_mutable ? PK_RATE : PK_CONST_RATE)
	      : (is_mutable ? PK_PGROUP : PK_CONST_PGROUP);
	    pscores2stream (out, pscores, tag, pgroup_subset[subset_index], " ", pcounts, use_bitscores);
	  }
      }
}

sstring PFunc_builder::score_sexpr (Score sc, bool use_bitscores)
{
  sstring s;
  if (use_bitscores)
    {
      if (sc <= -InfinityScore)
	s << PK_INFINITE;
      else if (sc == 0)
	s << 0;
      else
	{
	  s.precision(3);
	  s << fixed << -Score2Bits (sc);
	}
    }
  else
    s << Score2Prob (sc);
  return s;
}

void PFunc_builder::pscores2stream (ostream& out, const PScores& pscores, const char* tag, const vector<int>& pgroups_to_show, const char* indent, const PCounts* pcounts, bool use_bitscores)
{
  if (pgroups_to_show.size())
    {
      out << indent << '(' << tag << '\n';
      for_const_contents (vector<int>, pgroups_to_show, g)
	{
	  out << indent << " (";
	  const int sz = pscores.group_size(*g);
	  if (sz == 1)
	    {
	      out << pscores.group_suffix[*g][0] << ' ' << score_sexpr (pscores[PVar(*g,0)], use_bitscores);
	      // commented out due to clutter...
#if 0
	      if (pcounts)
		{
		  print_count (out, (*pcounts)[PVar(*g,0)]);
		  print_time (out, pcounts->wait[*g]);
		}
#endif
	    }
	  else
	    for (int v = 0; v < sz; ++v)
	      {
		if (v > 0)
		  out << '\n' << indent << "  ";
		out << '(' << pscores.group_suffix[*g][v] << ' ' << score_sexpr (pscores[PVar(*g,v)], use_bitscores);

		// commented out due to clutter...
#if 0
		if (pcounts)
		  out << " (" << PK_COUNT << ' ' << (*pcounts)[PVar(*g,v)] << ')';
#endif

		out << ')';
	      }
	  out << ")\n";
	}
      out << indent << " )  ;; end " << tag << '\n';
    }
}

void PFunc_builder::pcounts2stream (ostream& out, const PCounts& pcounts, const char* tag, const PCounts* baseline_pcounts,
				   bool interpret_single_element_pgroups_as_rate_variables, bool print_zero_counts)
{
  out << " (" << tag << "\n";
  for (int g = 0; g < pcounts.groups(); ++g)
    {
      const int sz = pcounts.group_size(g);
      if (sz == 1 && interpret_single_element_pgroups_as_rate_variables)
	{
	  double pcount = pcounts[PVar(g,0)];
	  double ptime = pcounts.wait[g];
	  if (baseline_pcounts)
	    {
	      pcount -= (*baseline_pcounts)[PVar(g,0)];
	      ptime -= baseline_pcounts->wait[g];
	    }
	  if (pcount > 0. || ptime > 0.)
	    out << "  (" << pcounts.group_suffix[g][0]
		<< ' ' << pcount
		<< ' ' << ptime
		<< ")\n";
	}
      else
	for (int v = 0; v < sz; ++v)
	  {
	    double pcount = pcounts[PVar(g,v)];
	    if (baseline_pcounts)
	      pcount -= (*baseline_pcounts)[PVar(g,v)];
	    if (pcount > 0. || print_zero_counts)
	      out << "  (" << pcounts.group_suffix[g][v]
		  << ' ' << pcount
		  << ")\n";
	  }
    }
  out << "  )  ;; end " << tag << "\n";
}

void PFunc_builder::pfunc2stream (ostream& out, const PScores& pscores, const PFunc& pfunc)
{
  pfunc.show (out, &pscores.group_suffix, true);
}


void PFunc_builder::init_pgroups (PScores& pscores, SymPVar& sym2pvar, SExpr& grammar_sexpr, const char* tag, set<int>* mutable_pgroups, bool force_rate, bool disallow_rate, bool use_bitscores)
{
  const vector<SExpr*> all_sexpr = grammar_sexpr.find_all (tag, 1);
  for_const_contents (vector<SExpr*>, all_sexpr, sexpr)
    {
      const vector<SExpr*> all_child_sexpr = (*sexpr)->values();
      for_const_contents (vector<SExpr*>, all_child_sexpr, c_sexpr)
	init_pgroup (pscores, sym2pvar, **c_sexpr, mutable_pgroups, force_rate, disallow_rate, use_bitscores);
    }
}

PGroup PFunc_builder::init_pgroup (PScores& pscores, SymPVar& sym2pvar, SExpr& pgroup_sexpr, set<int>* mutable_pgroups, bool force_rate, bool disallow_rate, bool use_bitscores)
{
  // get all param names & values
  vector<sstring> suffix;
  vector<Score> score;
  // check if first two child S-expressions are atoms, not lists; if so, this is a rate declaration
  const bool looks_like_rate = pgroup_sexpr.child.size() >= 2 && pgroup_sexpr[0].is_atom() && pgroup_sexpr[1].is_atom();
  if (disallow_rate && looks_like_rate)
    THROWEXPR ("In " << pgroup_sexpr << "\nSingle-element probability groups are not allowed\n");
  if (force_rate && !looks_like_rate)
    THROWEXPR ("In " << pgroup_sexpr << "\nThis does not look like a rate parameter declaration\n");
  if (force_rate || looks_like_rate)
    {
      suffix.push_back (pgroup_sexpr[0].get_atom());
      score.push_back (Prob2Score (pgroup_sexpr[1].get_atom().to_nonneg_double_strict()));
    }
  else  // treat it as a group of probability parameters (KNOWN BUG: internal representation means that single-element PGroups are still treated as rates)
    for_contents (list<SExpr>, pgroup_sexpr.child, child_iter)
    {
      suffix.push_back (child_iter->tag());
      const sstring dbl_string = child_iter->value().get_atom();
      Score sc;
      if (dbl_string == sstring (PK_INFINITE))
	{
	  if (!use_bitscores)
	    THROWEXPR ("Illegal infinite value in assignment " << pgroup_sexpr);
	  sc = -InfinityScore;
	}
      else if (use_bitscores)
	sc =  -Bits2Score(dbl_string.to_double());
      else
	sc = Prob2Score (dbl_string.to_nonneg_double_strict());
      score.push_back (sc);
    }

  // populate PScores object
  PGroup g = pscores.new_group (suffix.size(), 0, suffix);
  for (int v = 0; v < (int) score.size(); ++v)
    {
      if (sym2pvar.find (suffix[v]) != sym2pvar.end())
	THROWEXPR ("Multiple definitions for parameter " << suffix[v]);
      pscores[g[v]] = score[v];
      sym2pvar[suffix[v]] = g[v];
    }

  // add to mutable_pgroups
  if (mutable_pgroups)
    mutable_pgroups->insert (g.group_idx);

  // return
  return g;
}

const sstring plus_string("+"), times_string("*"), div_string("/"), const_string (PK_CONST), pow_string ("^");
PFunc PFunc_builder::init_pfunc (const SymPVar& sym2pvar, SExpr& pfunc_sexpr, int offset)
{
  if (CTAGGING(-2,INIT_PFUNC))
    CL << "Parsing PFunc S-expression: " << pfunc_sexpr << '\n';

  PFunc f;
  if (pfunc_sexpr.is_atom())
    {
      // handle constants
      SymPVar::const_iterator sym2pvar_iter = sym2pvar.find (pfunc_sexpr.get_atom());
      if (sym2pvar_iter != sym2pvar.end())
	f = sym2pvar_iter->second;
      else
	f = pfunc_sexpr.get_atom().to_nonneg_double_strict ("Undefined parameter or malformed expression");
    }
  else
    {
      // handle exponentiation
      if (pfunc_sexpr.child.size() == 3
	  && pfunc_sexpr[1].is_atom()
	  && pfunc_sexpr[1].atom == pow_string
	  && pfunc_sexpr[2].is_atom())
	{
	  const PFunc h = init_pfunc (sym2pvar, pfunc_sexpr[0], 0);
	  const double exponent = pfunc_sexpr[2].atom.to_double_strict ("Exponent is not a numeric constant");
	  f = PFunc (h, exponent);
	}
      else
	{
	  // process other types of expression
	  PFunc g;
	  bool div = false;
	  SExprIter child_iter = pfunc_sexpr.child.begin();
	  while (offset-- > 0)
	    ++child_iter;
	  for (; child_iter != pfunc_sexpr.child.end(); ++child_iter)
	    {
	      if (child_iter->is_atom() && child_iter->atom == plus_string)
		{
		  // +
		  if (div)
		    THROWEXPR ("Bad binary arithmetic operator sequence (/+) in S-expression " << pfunc_sexpr);
		  if (f.is_null())
		    f = g;
		  else if (!g.is_null())
		    f += g;
		  g = PFunc();
		}
	      else if (child_iter->is_atom() && child_iter->atom == times_string)
		{
		  // *
		  if (div)
		    THROWEXPR ("Bad binary arithmetic operator sequence (/*) in S-expression " << pfunc_sexpr);
		  // ignore multiplication signs
		}
	      else if (child_iter->is_atom() && child_iter->atom == div_string)
		{
		  // /
		  if (div)
		    THROWEXPR ("Bad binary arithmetic operator sequence (//) in S-expression " << pfunc_sexpr);
		  div = true;
		}
	      else if (child_iter->is_atom() && child_iter->atom == const_string)
		{
		  // #
		  if (++child_iter == pfunc_sexpr.child.end())
		    THROWEXPR ("Isolated '" << const_string << "' in S-expression " << pfunc_sexpr);

		  // inert PVar
		  PFunc h = init_pfunc (sym2pvar, *child_iter, 0);
		  h.make_inert();
	      
		  // default: multiply
		  g = g.is_null() ? h : (div ? g/h : g*h);
		  div = false;
		}
	      else if (child_iter->is_atom() && child_iter->atom == pow_string)
		{
		  THROWEXPR ("In " << pfunc_sexpr << ":\nExponents must appear in the form (F ^ C) where F is a function and C is a constant\nFurthermore, if F is a nontrivial function, it must itself be enclosed by parentheses");
		}
	      else
		{
		  // default: multiply
		  const PFunc h = init_pfunc (sym2pvar, *child_iter, 0);
		  g = g.is_null() ? h : (div ? g/h : g*h);
		  div = false;
		}
	    }
	  if (f.is_null())
	    f = g;
	  else if (!g.is_null())
	    f += g;
	}
    }

  if (CTAGGING(-2,INIT_PFUNC))
    {
      CL << "Returning PFunc: ";
      f.show (CL);
      CL << '\n';
    }

  return f;
}


void PFunc_builder::print_count (ostream& out, double count)
{
  // commented out due to clutter...
#if 0
  out << " (" << PK_COUNT << ' ' << count << ')';
#endif
}

void PFunc_builder::print_time (ostream& out, double t)
{
  // commented out due to clutter...
#if 0
  out << " (" << PK_TIME << ' ' << t << ')';
#endif
}
