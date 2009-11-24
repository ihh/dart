#include <map>
#include "ecfg/ecfgsexpr.h"
#include "util/vector_output.h"
#include "util/svisitor.h"
#include "seq/pkeywords.h"

// dummy class label for hidden alphabets
#define DUMMY_CLASS_LABEL "dummy_class_label"

int ECFG_builder::token_list_to_state (SExpr& token_list, const Alphabet& alph, int word_len, const vector<sstring>& class_alph)
{
  vector<sstring> toks = token_list.atoms_to_strings();
  const int expected_word_len = word_len + (class_alph.size() > 1 ? 1 : 0);
  if ((int) toks.size() != expected_word_len)
    THROWEXPR ("In (" << token_list << "):\nList contains " << toks.size() << " tokens; expected " << expected_word_len);
  int c = 0;
  if (class_alph.size() > 1)
    {
      c = get_state_class (toks.back(), class_alph);
      toks.pop_back();
    }
  int mul = 1, state = 0;
  for_const_contents (vector<sstring>, toks, tok)
    {
      assert_valid_token (*tok, alph, token_list);
      state += mul * alph.char2int_strict ((*tok)[0]);
      mul *= alph.size();
    }
  state += c * mul;
  return state;
}

const char* ECFG_builder::policy2string (Update_policy type)
{
  switch (type)
    {
    case Rev:
      return EG_POLICY_REV;

    case Irrev:
      return EG_POLICY_IRREV;

    case Rind:
      return EG_POLICY_RIND;

    case Parametric:
      return EG_PARAMETRIC;

    default:
      break;
    }

  return 0;
}

ECFG_enum::Update_policy ECFG_builder::string2policy (const sstring& policy_string)
{
  for (int i = 0; i < TotalPolicies; ++i)
    if (policy_string == sstring (policy2string ((Update_policy) i)))
      return (Update_policy) i;
  THROWEXPR ("Bad update policy '" << policy_string << "'");
  return TotalPolicies;  // unreachable
}

void ECFG_builder::init_chain_classes (sstring& class_row, vector<sstring>& class_alph, const int terms, SExpr& chain_sexpr)
{
  SExpr* class_sexpr = chain_sexpr.find (EG_CHAIN_CLASS);
  if (class_sexpr)
    {
      class_alph = (*class_sexpr) (EG_TRANSFORM_LABEL).atoms_to_strings();
      class_row = (*class_sexpr) (EG_TRANSFORM_ROW).get_atom();
      for_contents (vector<sstring>, class_alph, cl)
	{
	  if (cl->size() == 1)
	    *cl = sstring (terms, (*cl)[0]);
	  else if ((int) cl->size() != terms)
	    THROWEXPR ("In '" << *class_sexpr << "':\nChain has " << terms << " terminals, but class label '" << *cl << "' has " << cl->size() << " characters");
	}
    }
  else
    class_alph.push_back (sstring (DUMMY_CLASS_LABEL));  // add a dummy class label, so there's at least one "hidden" class
}

void ECFG_builder::init_chain (ECFG_matrix_set& ems, SymIndex& term2chain, const SymPVar& sym2pvar, SExpr& chain_sexpr, double tres)
{
  list<SExpr>& state_list = chain_sexpr (EG_CHAIN_TERMINAL).child;

  sstring class_row;
  vector<sstring> class_alph;
  init_chain_classes (class_row, class_alph, state_list.size(), chain_sexpr);

  const int chain_index = ems.add_matrix (state_list.size(),
					  string2policy (chain_sexpr (EG_CHAIN_POLICY).get_atom()),
					  class_alph.size(),
					  tres);

  ECFG_chain& chain = ems.chain[chain_index];
  chain.class_labels = class_alph;
  chain.class_row = class_row;

  const Alphabet& alph = ems.alphabet;
  chain.state.clear();
  for_contents (list<SExpr>, state_list, state)
    {
      if (!state->is_atom())
	THROWEXPR ("In (" << chain_sexpr << ")\nBad state name");
      const sstring& state_name = state->atom;
      if (term2chain.find (state_name) != term2chain.end())
	THROWEXPR ("In (" << chain_sexpr << ")\nDuplicate chain state symbol (" << state_name << ")\nSExpr: " << chain_sexpr);
      chain.state.push_back (state_name);
      term2chain[state_name] = chain_index;
    }
  if (chain.state.size() == 0)
    THROWEXPR ("In (" << chain_sexpr << "): terminal list empty");

  vector<SExpr*> pi_list = chain_sexpr.find_all (EG_CHAIN_INITIAL, 1);
  for_contents (vector<SExpr*>, pi_list, pi_sexpr)
    {
      const int state = token_list_to_state ((**pi_sexpr)(EG_CHAIN_STATE), alph, chain.word_len, class_alph); 
      SExpr& prob_sexpr = (*pi_sexpr)->find_or_die (EG_PROB, 1);

      // Default behavior for duplicated initial probs is to add them together
      if (chain.is_parametric)
	{
	  PFunc& pi_pfunc = chain.matrix_funcs->pi[state];
	  PFunc new_pfunc = init_pfunc (sym2pvar, prob_sexpr, 1);
	  if (pi_pfunc.is_null())
	    pi_pfunc = new_pfunc;
	  else
	    pi_pfunc += new_pfunc;
	}
      else
	chain.matrix->pi[state] += prob_sexpr.value().get_atom().to_nonneg_double_strict();
    }

  chain.matrix->init_update_flags (false);  // clear all update flags
  vector<SExpr*> mutate_list = chain_sexpr.find_all (EG_CHAIN_MUTATE, 1);
  for_contents (vector<SExpr*>, mutate_list, mutate_sexpr)
    {
      const int from_state = token_list_to_state ((**mutate_sexpr) (EG_FROM), alph, chain.word_len, class_alph);
      const int to_state = token_list_to_state ((**mutate_sexpr) (EG_TO), alph, chain.word_len, class_alph);
      SExpr& rate_sexpr = (*mutate_sexpr)->find_or_die (EG_RATE, 1);

      // Default behavior for duplicated rates is to add them together
      if (chain.is_parametric)
	{
	  PFunc& rate_pfunc = chain.matrix_funcs->X[0] (from_state, to_state);
	  PFunc new_pfunc = init_pfunc (sym2pvar, rate_sexpr, 1);

	  if (rate_pfunc.is_null())
	    rate_pfunc = new_pfunc;
	  else
	    rate_pfunc += new_pfunc;
	}
      else
	{
	  const double rate = rate_sexpr.value().get_atom().to_nonneg_double_strict();
	  chain.matrix->X[0] (from_state, to_state) += rate;
	  chain.matrix->X[0] (from_state, from_state) -= rate;
	  chain.matrix->X_update_flag[0] (from_state, to_state) = true;
	}
    }

  // initialise alphabet, update eigenvectors, etc.
  Alphabet dummy_alphabet ("", chain.matrix->A);
  chain.matrix->init_alphabet (dummy_alphabet);
  chain.matrix->update();
}

void ECFG_builder::init_hybrid_chain (ECFG_matrix_set& ems, SymIndex& term2chain, const SymPVar& sym2pvar, SExpr& hybrid_chain_sexpr)
{
  // Note that all submodels in a hybrid chain MUST have the same update-policy (& hence the same basis-projection subclass of EM_matrix_base)
  // See dart/src/hsm/em_matrix_base.h for a discussion of why this is so.

  list<SExpr>& state_list = hybrid_chain_sexpr (EG_CHAIN_TERMINAL).child;

  sstring class_row;
  vector<sstring> class_alph;
  init_chain_classes (class_row, class_alph, state_list.size(), hybrid_chain_sexpr);

  const int chain_index = ems.add_matrix (state_list.size(), Hybrid, class_alph.size());

  ECFG_chain& chain = ems.chain[chain_index];
  chain.class_labels = class_alph;
  chain.class_row = class_row;

  chain.state.clear();
  for_contents (list<SExpr>, state_list, state)
    {
      if (!state->is_atom())
	THROWEXPR ("In (" << hybrid_chain_sexpr << ")\nBad state name");
      const sstring& state_name = state->atom;
      if (term2chain.find (state_name) != term2chain.end())
	THROWEXPR ("In (" << hybrid_chain_sexpr << ")\nDuplicate chain state symbol (" << state_name << ")\nSExpr: " << hybrid_chain_sexpr);
      chain.state.push_back (state_name);
      term2chain[state_name] = chain_index;
    }
  if (chain.state.size() == 0)
    THROWEXPR ("In (" << hybrid_chain_sexpr << "): terminal list empty");

  chain.gs_tag = hybrid_chain_sexpr (EG_TRANSFORM_ROW).get_atom();

  const vector<SExpr*> components = hybrid_chain_sexpr.find_or_die (EG_HYBRID_COMPONENTS).values();
  if (components.size() == 0)
    THROWEXPR ("In (" << hybrid_chain_sexpr << "): component chain list empty");

  vector<int> component_chain_indices;
  for_const_contents (vector<SExpr*>, components, component_sexpr)
    {
      const vector<sstring> component_terminals = (**component_sexpr) (EG_CHAIN_TERMINAL).atoms_to_strings();
      if (component_terminals.size() == 0)
	THROWEXPR ("In (" << **component_sexpr << "): component terminal list empty");

      int component_chain_index = -1;
      ECFG_chain* cpt_chain = 0;
      for (int pos = 0; pos < (int) component_terminals.size(); ++pos)
	{
	  if (term2chain.find (component_terminals[pos]) == term2chain.end())
	    THROWEXPR ("In (" << **component_sexpr << "): unknown terminal");

	  const int term_chain_index = term2chain[component_terminals[pos]];
	  if (pos > 0 && term_chain_index != component_chain_index)
	    THROWEXPR ("In (" << **component_sexpr << "): can't mix & match terminals from different chains in one component chain");

	  component_chain_index = term_chain_index;
	  cpt_chain = &ems.chain[component_chain_index];

	  if (cpt_chain->state.size() != component_terminals.size())
	    THROWEXPR ("In (" << **component_sexpr << "): number of terminals in component chains doesn't match number of terminals in hybrid chain");

	  if (cpt_chain->state[pos] != component_terminals[pos])
	    THROWEXPR ("In (" << **component_sexpr << "): order of terminals in component chain expression doesn't match order in original chain declaration; this is not allowed");
	}

      // check that hidden class labels match
      if (cpt_chain->classes > 1 && chain.classes == 1)
	THROWEXPR ("In (" << **component_sexpr << "): component chain has hidden classes (" <<  cpt_chain->class_labels << "), whereas hybrid chain has no hidden classes");

      if (cpt_chain->classes == 1 && chain.classes > 1)
	THROWEXPR ("In (" << **component_sexpr << "): component chain has no hidden classes, whereas hybrid chain has hidden classes (" <<  chain.class_labels << ")");

      if (cpt_chain->classes != chain.classes || cpt_chain->class_labels != chain.class_labels)
	THROWEXPR ("In (" << **component_sexpr << "): component chain has hidden classes (" <<  cpt_chain->class_labels << "), whereas hybrid chain has different hidden classes (" <<  chain.class_labels << ")");

      // all component chains must be parametric
      // two reasons for this:
      //
      // (1) the weird mix of inheritance and templating in EM_matrix_base,
      //     along with the way the legacy code was hacked to introduce hybrid-chains,
      //     means that all rate matrices for a single Column_matrix must have the same type
      //
      // (2) without pseudocounts, component chains can often have zero start[] counts,
      //     which leads to 'nan' values and numerical Newmat crashes.
      //
      if (cpt_chain->type != Parametric)
	THROWEXPR ("In (" << hybrid_chain_sexpr << "): all component chains must have " << EG_PARAMETRIC << ' ' << EG_CHAIN_POLICY);

      const sstring gs_value = (**component_sexpr) (EG_TRANSFORM_LABEL).get_atom();
      if (chain.gs_tag_value_chain_index.find (gs_value) != chain.gs_tag_value_chain_index.end())
	  THROWEXPR ("In (" << hybrid_chain_sexpr << "): duplicate component chain label " << gs_value);

      chain.gs_values.push_back (gs_value);

      chain.gs_tag_value_chain_index[gs_value] = component_chain_index;
    }

  // if gs_tag begins with a special character, check that there are exactly two values, corresponding to true & false
  if (chain.gs_tag[0] == ECFG_IMPLICIT_GS_TAG_EQUALS || chain.gs_tag[0] == ECFG_IMPLICIT_GS_TAG_ANCESTOR)
    {
      bool implicit_gs_vals_ok = false;
      if (chain.gs_values.size() == 2)
	if ((chain.gs_values[0] == sstring (ECFG_IMPLICIT_GS_VALUE_TRUE) && chain.gs_values[1] == sstring (ECFG_IMPLICIT_GS_VALUE_FALSE))
	    || (chain.gs_values[1] == sstring (ECFG_IMPLICIT_GS_VALUE_TRUE) && chain.gs_values[0] == sstring (ECFG_IMPLICIT_GS_VALUE_FALSE)))
	  implicit_gs_vals_ok = true;
      if (!implicit_gs_vals_ok)
	THROWEXPR ("In (" << hybrid_chain_sexpr << "): hybrid chains with row labels beginning with the special character '"
		   << chain.gs_tag[0] << "' must have exactly two labels, '" << ECFG_IMPLICIT_GS_VALUE_TRUE << "' and '" << ECFG_IMPLICIT_GS_VALUE_FALSE << "'");
    }
}

void ECFG_builder::init_chain_given_alphabet (EM_matrix_base& hsm, const Alphabet& alph, SExpr& chain_sexpr, int required_pseudoterms)
{
  ECFG_matrix_set ems (alph);
  SymIndex term2chain;
  SymPVar sym2pvar;

  init_chain (ems, term2chain, sym2pvar, chain_sexpr.find_or_die (EG_CHAIN));

  const ECFG_chain& chain = ems.chain[0];
  if (chain.word_len != required_pseudoterms)
    THROWEXPR ("Chain has " << chain.word_len << " pseudoterminals where " << required_pseudoterms << " were expected");

  const EM_matrix_base& read_hsm = *chain.matrix;

  hsm.init_matrix (read_hsm.C, read_hsm.A);
  hsm.assign_matrix_params (read_hsm);

  hsm.init_alphabet (((EM_matrix_base&) read_hsm).alphabet());
  hsm.update();
}

void ECFG_builder::init_chain_and_alphabet (Alphabet& alph, EM_matrix_base& hsm, SExpr& alph_chain_sexpr)
{
  Alphabet base_alph ("", 1);
  init_alphabet (base_alph, alph_chain_sexpr.find_or_die (PK_ALPHABET));

  ECFG_matrix_set ems (base_alph);
  SymIndex term2chain;
  SymPVar sym2pvar;
  init_chain (ems, term2chain, sym2pvar, alph_chain_sexpr.find_or_die (EG_CHAIN));

  const ECFG_chain& chain = ems.chain[0];
  const EM_matrix_base& read_hsm = *chain.matrix;

  alph.init_hidden (base_alph, chain.class_labels);

  hsm.init_matrix (read_hsm.C, read_hsm.A);
  hsm.assign_matrix_params (read_hsm);
  hsm.init_alphabet (alph);
  hsm.update();
}

ECFG_builder::ECFG_symbol_sequence::ECFG_symbol_sequence (const ECFG_matrix_set& ems, const SymIndex& nonterm2state, const SymIndex& term2chain, SExpr& symseq_sexpr)
{
  sexpr = &symseq_sexpr;
  chain = 0;
  chain_index = -1;
  int pos = 0;
  for_contents (list<SExpr>, symseq_sexpr.child, lhs)
    {
      if (!lhs->is_atom())
	THROWEXPR ("In (" << symseq_sexpr << ")\nSymbols on LHS of rules must be atoms");
      sstring lhs_name = lhs->atom;
      bool comp = false;
      if (lhs_name[0] == ECFG_complement_character)
	{
	  comp = true;
	  lhs_name.erase (lhs_name.begin());
	}
      if (term2chain.find (lhs_name) != term2chain.end())
	{
	  // symbol is an evolving terminal on an ECFG_chain, so find the chain and token
	  const int c = term2chain.find(lhs_name)->second;  // chain
	  const vector<sstring>& state = ems.chain[c].state;
	  const int tok = find (state.begin(), state.end(), lhs_name) - state.begin();

	  if (chain)  // test for consistency with previous terminals
	    {
	      if (chain_index != c)
		THROWEXPR ("In (" << symseq_sexpr << ")\nChain mismatch");
	      if (term_pos.find (tok) != term_pos.end())
		THROWEXPR ("In (" << symseq_sexpr << ")\nDuplicate token");
	    }
	  else
	    {
	      // this is the first terminal in the sequence; initialise vars
	      chain = &ems.chain[c];
	      chain_index = c;
	    }

	  // record this symbol
	  is_term.push_back (true);
	  sym.push_back (tok);
	  term_pos[tok] = pos;
	  pos_comped[pos] = comp;
	}

      else if (!comp)
	{
	  bool postemit = false;
	  bool found = nonterm2state.find (lhs_name) != nonterm2state.end();
	  if (!found)
	    {
	      const char c = lhs_name.back();
	      if (c == ECFG_post_emit_character || c == ECFG_deprecated_post_emit_character)
		{
		  postemit = true;
		  lhs_name.chop();
		  found = nonterm2state.find (lhs_name) != nonterm2state.end();
		}
	    }
	  if (found)
	    {
	      // symbol is a nonterminal state; record it
	      const int state = nonterm2state.find (lhs_name)->second;
	      is_term.push_back (false);
	      sym.push_back (state);
	      pos_nonterm[pos] = state;
	      pos_postemit[pos] = postemit;
	      pos_name[pos] = lhs_name;
	    }
	  else
	    THROWEXPR ("In (" << symseq_sexpr << "):\nSymbol (" << lhs_name << ") unknown");
	}

      // increment position
      ++pos;
    }
}

ECFG_builder::ECFG_rule_block::ECFG_rule_block (const ECFG_matrix_set& ems,
						const SymIndex& nonterm2state,
						const SymIndex& term2chain,
						bool parametric_transitions,
						const SymPVar& sym2pvar,
						SExpr& transform_sexpr)
  : nonterm2state (nonterm2state),
    term2chain (term2chain),
    parametric_transitions (parametric_transitions),
    sym2pvar (sym2pvar),
    sexpr (&transform_sexpr)
{
  // parse lhs and rhs
  lhs = ECFG_symbol_sequence (ems, nonterm2state, term2chain, transform_sexpr (EG_FROM));
  rhs = ECFG_symbol_sequence (ems, nonterm2state, term2chain, transform_sexpr (EG_TO));
}

ECFG_builder::SymIndex ECFG_builder::init_nonterm2state (const SymIndex& term2chain, SExpr& grammar_sexpr)
{
  SymIndex nonterm2state;
  int n = 0;

  const vector<SExpr*> nonterm = grammar_sexpr.find_all (EG_NONTERMINAL, 1);
  for_const_contents (vector<SExpr*>, nonterm, nt)
    {
      const sstring nt_name = (**nt) (EG_NAME).get_atom();
      if (nonterm2state.find (nt_name) == nonterm2state.end())
	nonterm2state[nt_name] = n++;
    }

  const vector<SExpr*> trans = grammar_sexpr.find_all (EG_TRANSFORM, 1);
  for_const_contents (vector<SExpr*>, trans, t)
    {
      SExpr& lhs = (**t) (EG_FROM);
      for_const_contents (list<SExpr>, lhs.child, lhs_sym)
	{
	  sstring lsym = lhs_sym->atom;
	  if (lsym[0] == ECFG_complement_character)
	    lsym.erase (lsym.begin());
	  else if (lsym.back() == ECFG_post_emit_character || lsym.back() == ECFG_deprecated_post_emit_character)
	    lsym.erase (--lsym.end());
	  if (term2chain.find (lsym) != term2chain.end())
	    continue;
	  if (nonterm2state.find (lsym) == nonterm2state.end())
	    nonterm2state[lsym] = n++;
	}
    }

  return nonterm2state;
}

void ECFG_builder::ECFG_rule_block::parse (ECFG_scores& ecfg)
{
  // log
  if (CTAGGING(-1,ECFG_PARSE_RULE))
    CL << "Parsing rule block: " << *sexpr << "\n";

  // check there's exactly one nonterminal on LHS
  if (lhs.pos_nonterm.size() != 1)
    THROWEXPR ("In (" << *sexpr << ")\nLHS of '" << EG_TRANSFORM << "' rule must have exactly one nonterminal");

  // figure out what LHS nonterminal is
  const int lhs_nonterm_pos = lhs.pos_nonterm.begin()->first;
  const int lhs_nonterm = lhs.pos_nonterm[lhs_nonterm_pos];
  ECFG_state_info& info = ecfg.state_info[lhs_nonterm];
  const sstring& lhs_name = lhs.pos_name[lhs_nonterm_pos];

  // if LHS nonterminal name ends in post-emit character, check it has been initialised as an emit state, and if not, flag it as such
  const bool from_postemit = lhs.pos_postemit[lhs_nonterm_pos];
  if (from_postemit)
    {
      const sstring post_emit_name (ECFG_post_emit_character);
      if (info.name.size() && info.name != post_emit_name)
	{
	  if (info.total_size() == 0)
	    THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: LHS nonterminal ends in \""
		       << ECFG_post_emit_character << "\" or \"" << ECFG_deprecated_post_emit_character
		       << "\", but is not an emit state; call me easily confused, but I'm bailing");
	}
      else
	{
	  info.name = post_emit_name;   // flag as an uninitialized emit state
	  CTAG(1,ECFG_PARSE_RULE) << "Flagging " << lhs_name << " as an emit state\n";
	}
    }

  // figure out what kind of state this is: emit, bifurcation, null
  if (rhs.chain)  // emission?
    {
      if (lhs.chain_index >= 0 && lhs.chain_index != rhs.chain_index)
	THROWEXPR ("In (" << *sexpr << ")\nChain mismatch on LHS and RHS of '" << EG_TRANSFORM << "' rule");
      if ((int) rhs.term_pos.size() != rhs.chain->word_len)
	THROWEXPR ("In (" << *sexpr << ")\nChain terminals unaccounted for in '" << EG_TRANSFORM << "' rule");
      if ((int) lhs.term_pos.size() == rhs.chain->word_len)
	THROWEXPR ("In (" << *sexpr << ")\nAll chain terminals appear on LHS of '" << EG_TRANSFORM << "' rule");
      if (rhs.pos_nonterm.size() != 1)
	THROWEXPR ("In (" << *sexpr << ")\nFor emissions, must have exactly one nonterminal on RHS of '" << EG_TRANSFORM << "' rule");

      // check RHS nonterm is post-emit version of LHS nonterm
      const int rhs_nonterm_pos = rhs.pos_nonterm.begin()->first;
      if (lhs_nonterm != rhs.pos_nonterm[rhs_nonterm_pos]
	  || lhs.pos_postemit[lhs_nonterm_pos]
	  || !rhs.pos_postemit[rhs_nonterm_pos])
	THROWEXPR ("In (" << *sexpr << ")\nEmit rules must be of the form ... X ... -> ... X" << ECFG_post_emit_character << " ... where X is a nonterminal");
      // check context-sensitive terminals appear in the same order on LHS & RHS
      for (int lpos = 0; lpos < lhs_nonterm_pos; ++lpos)
	if (!rhs.is_term[lpos] || rhs.sym[lpos] != lhs.sym[lpos] || rhs.pos_comped[lpos] != lhs.pos_comped[lpos])
	  THROWEXPR ("In (" << *sexpr << ")\nLHS/RHS mismatch at position " << lpos);
      for (int lpos = (int) lhs.sym.size() - 1; lpos > lhs_nonterm_pos; --lpos)
	{
	  const int rpos = rhs.sym.size() - lhs.sym.size() + lpos;
	  if (!rhs.is_term[rpos] || rhs.sym[rpos] != lhs.sym[lpos] || rhs.pos_comped[rpos] != lhs.pos_comped[lpos])
	    THROWEXPR ("In (" << *sexpr << ")\nLHS/RHS mismatch at position " << lpos << " of LHS, " << rpos << " of RHS");
	}
      // Check state not already initialised
      if (info.name.size() && info.name != sstring (ECFG_post_emit_character))
	THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: describes emit state, but state already initialised");

      // log
      if (CTAGGING(-1,ECFG_PARSE_RULE))
	CL << "Initializing " << lhs_name << " as an emit state\n";

      // Initialise the emit state
      const int r_context = lhs.sym.size() - lhs_nonterm_pos - 1;
      info = ECFG_state_info (lhs_nonterm_pos, rhs_nonterm_pos - lhs_nonterm_pos,
			      rhs.sym.size() - rhs_nonterm_pos - 1 - r_context, r_context);
      info.name = lhs_name;
      info.matrix = rhs.chain_index;
      CTAG(1,ECFG_PARSE_RULE) << "Initializing " << lhs_name << " as an emit state\n";

      // initialise emission
      // since term_pos is sorted by increasing chain token, 'mul' will be correctly set
      int mul = 1;
      for_const_contents (ECFG_symbol_sequence::IntIntMap, rhs.term_pos, tp)
	{
	  const int pos = tp->second;
	  int chain_pos = pos > rhs_nonterm_pos ? (pos - 1) : pos;
	  info.mul[chain_pos] = mul;
	  info.comp[chain_pos] = rhs.pos_comped[pos];
	  mul *= ecfg.matrix_set.alphabet.size();
	}

      // initialise indel params
      init_gaps (info, sym2pvar, sexpr);

      // initialise the annotation
      typedef map<sstring,sstring> String_string_map;
      String_string_map deterministic_annot;
      const vector<SExpr*> all_annot_sexpr = sexpr->find_all (EG_TRANSFORM_ANNOTATE, 1);
      for_const_contents (vector<SExpr*>, all_annot_sexpr, annot_sexpr)
	{
	  const sstring& row = (**annot_sexpr) (EG_TRANSFORM_ROW).get_atom();
	  if ((*annot_sexpr)->find (EG_ANNOTATE_EMIT) != 0) {
	    // handle probabilistic annotations
	    const vector<SExpr*> all_emit_sexpr = (*annot_sexpr)->find_all (EG_ANNOTATE_EMIT, 1);
	    for_const_contents (vector<SExpr*>, all_emit_sexpr, emit_sexpr) {
	      const sstring multi_label = sstring::join ((**emit_sexpr) (EG_TRANSFORM_LABEL).atoms_to_strings(), "");
	      PFunc pfunc = init_pfunc (sym2pvar, (**emit_sexpr) (EG_PROB), 1);
	      info.annot[row][multi_label] = pfunc;
	    }

	  } else {
	    // handle old-format deterministic annotations
	    const sstring& col = (**annot_sexpr) (EG_TRANSFORM_COLUMN).get_atom();
	    const sstring& label = (**annot_sexpr) (EG_TRANSFORM_LABEL).get_atom();

	    // check chain terminal
	    if (term2chain.find (col) == term2chain.end())
	      THROWEXPR ("In (" << **annot_sexpr << "):\nChain terminal (" << col << ") unknown");

	    const int chain_idx = term2chain.find(col)->second;
	    if (chain_idx != rhs.chain_index)
	      THROWEXPR ("In (" << **annot_sexpr << "):\nChain terminal (" << col << ") is not on the same chain as the terminals in this rule");

	    // get column corresponding to named chain terminal
	    const vector<sstring>& state = ecfg.matrix_set.chain[rhs.chain_index].state;
	    const int term_idx = find (state.begin(), state.end(), col) - state.begin();
	    int col_idx = rhs.term_pos[term_idx];
	    if (col_idx > rhs_nonterm_pos)
	      --col_idx;

	    // check column is emitted rather than context
	    if (col_idx < info.l_context || col_idx >= rhs.chain->word_len - info.r_context)
	      THROWEXPR ("In (" << **annot_sexpr << "):\nCannot annotate chain terminal (" << col << ") as it is contextual, rather than emitted, in this rule");

	    // check label is a single char
	    if (label.size() != 1)
	      THROWEXPR ("In (" << **annot_sexpr << "):\nLabel (" << label << ") must be a single character");

	    // annotate
	    if (deterministic_annot.find (row) == deterministic_annot.end())
	      deterministic_annot[row] = sstring (info.emit_size(), ECFG_annotation_wildcard);
	    deterministic_annot[row][col_idx - info.l_context] = label[0];
	  }
	}

      // convert old-format deterministic annotations into new probabilistic-annotation data structure
      if (deterministic_annot.size()) {
	for_const_contents (String_string_map, deterministic_annot, ss) {
	  if (info.annot.find (ss->first) != info.annot.end())
	    THROWEXPR ("Cannot mix deterministic and probabilistic annotations (tag: " << ss->first << ")");
	  info.annot[ss->first][ss->second] = PFunc(1.);
	}
      }

      // min_len, max_len, infix, prefix, suffix
      if (SExpr* minlen_sexpr = sexpr->find (EG_TRANSFORM_MINLEN, 1))
	info.min_len = minlen_sexpr->value().get_atom().to_int();

      if (sexpr->find (EG_TRANSFORM_INFIX, 1))
	info.infix = true;
      else if (SExpr* maxlen_sexpr = sexpr->find (EG_TRANSFORM_MAXLEN, 1))
	info.max_len = maxlen_sexpr->value().get_atom().to_int();

      if (sexpr->find (EG_TRANSFORM_PREFIX, 1))
	info.prefix = true;

      if (sexpr->find (EG_TRANSFORM_SUFFIX, 1))
	info.suffix = true;
    }

  else  // bifurcation/null?
    {
      if (lhs.chain)
	THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: has terminals on LHS but not on RHS");
      if (sexpr->find (EG_TRANSFORM_ANNOTATE))
	THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: is not an emission, but contains '" << EG_TRANSFORM_ANNOTATE << "' clause");
      switch (rhs.pos_nonterm.size())
	{
	case 0:
	case 1:  // null
	  if (!from_postemit)
	    {
	      // Check state not already initialised as emit or bifurc state
	      if (info.name.size() && (info.total_size() || info.bifurc))
		THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: describes null state, but state already initialised as emission or bifurcation");
	      // set state type to null
	      info = ECFG_null_state_info();
	      info.name = lhs_name;
	      CTAG(1,ECFG_PARSE_RULE) << "Initializing " << lhs_name << " as a null state\n";
	    }
	  break;
	case 2:  // bifurcation
	  // Check state not already initialised
	  if (info.name.size())
	    THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: describes bifurcation state, but state already initialised");
	  // Initialise bifurc state
	  info = ECFG_bifurc_state_info (rhs.sym[0], rhs.sym[1]);
	  info.name = lhs_name;
	  CTAG(1,ECFG_PARSE_RULE) << "Initializing " << lhs_name << " as a bifurcation state\n";
	  break;
	default:  // error
	  THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: has more than two nonterminals on RHS");
	  break;
	}
    }

  // Set transition scores
  if (!rhs.chain)
    {
      typedef map<int,bool> IntBoolMap;
      for_const_contents (IntBoolMap, rhs.pos_postemit, pd)
	if (pd->second)
	  THROWEXPR ("In (" << *sexpr << ")\nBad '" << EG_TRANSFORM << "' rule: not an emission, but RHS contains a nonterminal ending with post-emit characters \""
		     << ECFG_post_emit_character << "\" or [deprecated] \"" << ECFG_deprecated_post_emit_character << "\"");

      // get probability, or PFunc
      Prob prob = 1.;
      PFunc pfunc (prob);
      SExpr* prob_sexpr = sexpr->find (EG_PROB);
      sstring dest_state_text, trans_prob_text;  // for logging
      if (prob_sexpr)
	{
	  if (parametric_transitions)
	    {
	      pfunc = init_pfunc (sym2pvar, *prob_sexpr, 1);
	      pfunc.show (trans_prob_text);
	      trans_prob_text << " (PFunc)";
	    }
	  else
	    {
	      prob = prob_sexpr->value().get_atom().to_nonneg_double_strict();
	      trans_prob_text << prob << " (const)";
	    }
	}

      // figure out dest state
      switch (rhs.pos_nonterm.size())
	{
	case 0:  // transition to end state
	  dest_state_text = "End";
	  if (parametric_transitions)
	    ecfg.trans_funcs.transition (lhs_nonterm, End) = pfunc;
	  else
	    ecfg.transition (lhs_nonterm, End) = Prob2Score (prob);
	  break;
	  
	case 1:  // transition to null state
	  {
	    const int rhs_nonterm = rhs.pos_nonterm.begin()->second;
	    dest_state_text = rhs.pos_name.begin()->second;
	    if (parametric_transitions)
	      ecfg.trans_funcs.transition (lhs_nonterm, rhs_nonterm) = pfunc;
	    else
	      ecfg.transition (lhs_nonterm, rhs_nonterm) = Prob2Score (prob);
	    break;
	  }

	case 2:  // bifurcation
	  break;

	default:  // something else; complain (although earlier tests should catch syntax errors)
	  THROWEXPR ("In (" << *sexpr << ")\nBad RHS block (this error should be unreachable)");
	  break;
	}

      // log, if loggable
      if (CTAGGING(-1,ECFG_PARSE_RULE) && dest_state_text.size())
	CL << "Adding transition from " << lhs_name << " to " << dest_state_text << " with probability " << trans_prob_text << "\n";
    }
}

void ECFG_builder::init_gaps (ECFG_state_info& info, const SymPVar& sym2pvar, SExpr* sexpr)
{
  SExpr* gap_sexpr;
  if (sexpr->find (EG_TRANSFORM_NO_GAPS, 1))
    info.gaps_ok = info.wild_gaps = info.indels = false;
  else if (sexpr->find (EG_TRANSFORM_STRICT_GAPS, 1))
    {
      info.gaps_ok = true;
      info.wild_gaps = false;
      info.indels = false;
    }
  else if ((gap_sexpr = sexpr->find (EG_TRANSFORM_GAP_MODEL, 1)))
    {
      info.gaps_ok = true;
      info.wild_gaps = sexpr->find (EG_TRANSFORM_STRICT_GAPS, 1) ? false : true;
      info.indels = true;
      info.has_parametric_indels = gap_sexpr->find (EG_PARAMETRIC) != 0;

      if (info.has_parametric_indels)
	{
	  info.link_extend_func = init_pfunc (sym2pvar, gap_sexpr->find_or_die (EG_TRANSFORM_EXTEND_PROB), 1);
	  info.link_end_func = init_pfunc (sym2pvar, gap_sexpr->find_or_die (EG_TRANSFORM_END_PROB), 1);
	  info.ins_rate_func = init_pfunc (sym2pvar, gap_sexpr->find_or_die (EG_TRANSFORM_INS_RATE), 1);
	  info.del_rate_func = init_pfunc (sym2pvar, gap_sexpr->find_or_die (EG_TRANSFORM_DEL_RATE), 1);
	}
      else
	{
	  info.link_extend = (*gap_sexpr) (EG_TRANSFORM_EXTEND_PROB).get_atom().to_nonneg_double_strict();
	  info.link_end = 1. - info.link_extend;
	  info.ins_rate = (*gap_sexpr) (EG_TRANSFORM_INS_RATE).get_atom().to_nonneg_double_strict();
	  info.del_rate = (*gap_sexpr) (EG_TRANSFORM_DEL_RATE).get_atom().to_nonneg_double_strict();
	}
    }
  else if (sexpr->find (EG_TRANSFORM_IGNORE_GAPS, 1))
    {
      info.gaps_ok = true;
      info.wild_gaps = true;
      info.indels = false;
    }
}

ECFG_scores* ECFG_builder::init_ecfg (const Alphabet& alph, SExpr& grammar_sexpr, double tres)
{
  // create & run the validator
  SExpr_validator ecfg_validator
    ("Grammar->('"EG_GRAMMAR" GrammarProperty*);"
     "GrammarProperty->Name|('"EG_META" Wild)|('"EG_TRANSIENT_META" Wild)|('"EG_UPDATE_RULES" Atom)|('"EG_UPDATE_RATES" Atom)|ParametricFlag|('"EG_PSEUDOCOUNTS" Count*)|('"PK_RATE" ParameterValue*)|('"PK_PGROUP" ParameterGroup*)|('"PK_CONST_RATE" ParameterValue*)|('"PK_CONST_PGROUP" ParameterGroup*)|('"EG_PARAMS" ParameterGroup*)|('"EG_CONST" ParameterGroup*)|('"EG_NONTERMINAL" NontermProperty*)|QualifiedSummationDirective|('"EG_HYBRID_CHAIN" HybridChainProperty*)|('"EG_CHAIN" ChainProperty*)|('"EG_TRANSFORM" RuleProperty*);"
     "Name->('"EG_NAME" Atom);"
     "ParametricFlag->('"EG_PARAMETRIC" End);"
     "MinimumLength->('"EG_TRANSFORM_MINLEN" Atom);"
     "MaximumLength->('"EG_TRANSFORM_MAXLEN" Atom);"
     "SummationDirective->('"EG_TRANSFORM_SUM_FROM" End);"
     "QualifiedSummationDirective->('"EG_TRANSFORM_SUM_FROM" Atom);"
     "PrefixConstraint->('"EG_TRANSFORM_PREFIX" End);"
     "SuffixConstraint->('"EG_TRANSFORM_SUFFIX" End);"
     "InfixConstraint->('"EG_TRANSFORM_INFIX" End);"
     "SourceStateList->('"EG_FROM" (Atom*));"
     "DestinationStateList->('"EG_TO" (Atom*));"
     "ProbabilityExpression->('"EG_PROB" Function);"
     "Function->Wild;"
     "ParameterValue->(Atom Atom);"
     "ParameterValuePair->(Atom ValuePair);"
     "ValuePair->Atom Atom;"
     "ParameterGroup->(ParameterValue*);"
     "Count->ParameterValue|ParameterValuePair;"
     "NontermProperty->Name|MinimumLength|MaximumLength|SummationDirective|PrefixConstraint|SuffixConstraint|InfixConstraint;"
     "RuleProperty->SourceStateList|DestinationStateList|ProbabilityExpression|('"EG_TRANSFORM_ANNOTATE" AnnotationProperty*)|MinimumLength|MaximumLength|InfixConstraint|PrefixConstraint|SuffixConstraint|SummationDirective|('"EG_TRANSFORM_NO_GAPS" End)|('"EG_TRANSFORM_STRICT_GAPS" End)|('"EG_TRANSFORM_IGNORE_GAPS" End)|('"EG_TRANSFORM_GAP_MODEL" GapModelProperty*);"
     "AnnotationProperty->AnnotationRow|AnnotationColumn|AnnotationLabel|ProbabilisticAnnotation;"
     "AnnotationRow->('"EG_TRANSFORM_ROW" Atom);"
     "AnnotationColumn->('"EG_TRANSFORM_COLUMN" Atom);"
     "AnnotationLabel->('"EG_TRANSFORM_LABEL" Atom);"
     "AnnotationLabelList->('"EG_TRANSFORM_LABEL" (Atom*));"
     "ProbabilisticAnnotation->('"EG_ANNOTATE_EMIT" ProbabilisticAnnotationProperty*);"
     "ProbabilisticAnnotationProperty->AnnotationLabelList|ProbabilityExpression;"
     "GapModelProperty->ParametricFlag|('"EG_TRANSFORM_STRICT_GAPS" End)|('"EG_TRANSFORM_EXTEND_PROB" Function)|('"EG_TRANSFORM_END_PROB" Function)|('"EG_TRANSFORM_INS_RATE" Function)|('"EG_TRANSFORM_DEL_RATE" Function);"
     "HybridChainProperty->TerminalList|AnnotationRow|('"EG_HYBRID_COMPONENTS" HybridChainComponent*);"
     "TerminalList->('"EG_TERMINAL" (Atom*));"
     "HybridChainComponent->(HybridChainComponentProperty HybridChainComponentProperty);"
     "HybridChainComponentProperty->AnnotationLabel|TerminalList;"
     "ChainProperty->TerminalList|('"EG_CHAIN_CLASS" HiddenClassDescription)|InitialStateProbability|MutationRate|('"EG_CHAIN_POLICY" UpdatePolicy);"
     "HiddenClassDescription->AnnotationRow|AnnotationLabelList;"
     "UpdatePolicy->'"EG_POLICY_REV"|'"EG_POLICY_IRREV"|'"EG_POLICY_RIND"|'"EG_PARAMETRIC";"
     "InitialStateProbability->('"EG_CHAIN_INITIAL" InitialStateProbabilityProperty*);"
     "InitialStateProbabilityProperty->StateList|ProbabilityExpression;"
     "StateList->('"EG_CHAIN_STATE" (Atom*));"
     "MutationRate->('"EG_CHAIN_MUTATE" MutationRateProperty*);"
     "MutationRateProperty->SourceStateList|DestinationStateList|RateExpression;"
     "RateExpression->('"EG_RATE" Function);"
);
  ecfg_validator.parse(grammar_sexpr);

  // pointer to ECFG_scores
  ECFG_scores* ecfg = 0;

  // create symbol lookup tables
  SymIndex term2chain, nonterm2state;
  SymPVar sym2pvar;

  // do two passes: the first is to count states (i.e. nonterminals)
  // christ, this is horribly inefficient... a lot of initialization is done twice, unnecessarily. EUCH.
  int states = 0;
  for (int pass = 0; pass < 2; ++pass)
    {
      // initialise ECFG
      if (ecfg)
	delete ecfg;
      ecfg = new ECFG_scores (alph, states);
      for (int s = 0; s < states; ++s)
	ecfg->state_info[s].name.clear();  // clear state names

      // initialise name
      ecfg->name = grammar_sexpr (EG_NAME).get_atom();

      // initialise meta
      const vector<SExpr*> all_meta_sexpr = grammar_sexpr.find_all (EG_META, 1);
      for_const_contents (vector<SExpr*>, all_meta_sexpr, meta_sexpr)
	ecfg->meta.push_back (**meta_sexpr);

      // initialise PScores
      sym2pvar.clear();

      // first old-style (params, const) ... should probably issue a warning here
      init_pgroups (ecfg->pscores, sym2pvar, grammar_sexpr, EG_PARAMS, &ecfg->mutable_pgroups, false, false);
      init_pgroups (ecfg->pscores, sym2pvar, grammar_sexpr, EG_CONST, (set<int>*) 0, false, false);

      // now the new style (pgroup, rate, const-pgroup, const-rate)
      init_pgroups (ecfg->pscores, sym2pvar, grammar_sexpr, PK_PGROUP, &ecfg->mutable_pgroups, false, true);
      init_pgroups (ecfg->pscores, sym2pvar, grammar_sexpr, PK_RATE, &ecfg->mutable_pgroups, true, false);
      init_pgroups (ecfg->pscores, sym2pvar, grammar_sexpr, PK_CONST_PGROUP, (set<int>*) 0, false, true);
      init_pgroups (ecfg->pscores, sym2pvar, grammar_sexpr, PK_CONST_RATE, (set<int>*) 0, true, false);

      // initialise PCounts
      const vector<SExpr*> all_pc_sexpr = grammar_sexpr.find_all (EG_PSEUDOCOUNTS, 1);
      if (all_pc_sexpr.size())
	ecfg->pcounts = PCounts (ecfg->pscores);   // zero the pseudocounts initially
      for_const_contents (vector<SExpr*>, all_pc_sexpr, pseudocounts_sexpr)
	{
	  const vector<SExpr*> all_pseudocounts_sexpr = (*pseudocounts_sexpr)->values();
	  for_const_contents (vector<SExpr*>, all_pseudocounts_sexpr, pcount_sexpr)
	    {
	      const sstring& pvar_name = (*pcount_sexpr)->tag();
	      if (sym2pvar.find (pvar_name) != sym2pvar.end())
		{
		  const PVar& pv = sym2pvar[pvar_name];
		  if ((**pcount_sexpr).child.size() < 2)
		    THROWEXPR("In (" << **pcount_sexpr << "):\nPseudocounts for rate and probability parameters need an event count");
		  ecfg->pcounts[pv] = (**pcount_sexpr)[1].get_atom().to_nonneg_double_strict();
		  if (ecfg->pscores.group_size (pv.group_idx) == 1) {  // rate parameter?
		    if ((**pcount_sexpr).child.size() < 3)
		      THROWEXPR("In (" << **pcount_sexpr << "):\nPseudocounts for rate parameters comprise both an event count and a time period");
		    ecfg->pcounts.wait[pv.group_idx] = (**pcount_sexpr)[2].get_atom().to_nonneg_double_strict();
		  }
		}
	      else
		CLOGERR << "Warning: a pseudocount is specified for parameter " << pvar_name << " but there is no corresponding parameter definition\n";
	    }
	}

      // initialise chains
      term2chain.clear();
      const vector<SExpr*> all_chain_sexpr = grammar_sexpr.find_all (EG_CHAIN, 1);
      for_const_contents (vector<SExpr*>, all_chain_sexpr, chain_sexpr)
	init_chain (ecfg->matrix_set, term2chain, sym2pvar, **chain_sexpr, tres);

      // initialise hybrid chains
      const vector<SExpr*> all_hybrid_chain_sexpr = grammar_sexpr.find_all (EG_HYBRID_CHAIN, 1);
      for_const_contents (vector<SExpr*>, all_hybrid_chain_sexpr, hybrid_chain_sexpr)
	init_hybrid_chain (ecfg->matrix_set, term2chain, sym2pvar, **hybrid_chain_sexpr);

      // initialise nonterminals
      nonterm2state = init_nonterm2state (term2chain, grammar_sexpr);

      // count nonterminals
      states = nonterm2state.size();
      if (states == 0)
	THROWEXPR ("In (" << grammar_sexpr << "):\nGrammar doesn't have any nonterminals");
    }

  // parametric?
  ecfg->has_parametric_transitions = grammar_sexpr.find (EG_PARAMETRIC) != 0;

  // parse "transform" blocks
  const vector<SExpr*> all_transform_sexpr = grammar_sexpr.find_all (EG_TRANSFORM, 1);
  for_const_contents (vector<SExpr*>, all_transform_sexpr, transform_sexpr)
    {
      ECFG_rule_block rule_block (ecfg->matrix_set, nonterm2state, term2chain,
				  ecfg->has_parametric_transitions, sym2pvar, **transform_sexpr);
      rule_block.parse (*ecfg);
    }

  // look for isolated "sum-from" modifiers; flag any sum states as such
  const vector<SExpr*> all_sum_state_sexpr = grammar_sexpr.find_all (EG_TRANSFORM_SUM_FROM, 1);
  for_const_contents (vector<SExpr*>, all_sum_state_sexpr, sum_state_sexpr)
    if ((*sum_state_sexpr)->has_value())
      {
	const sstring& sum_state_name = (*sum_state_sexpr)->value().get_atom();
	if (nonterm2state.find (sum_state_name) != nonterm2state.end())
	  ecfg->state_info[nonterm2state[sum_state_name]].sum_state = true;
      }

  // look for "nonterminal" blocks containing "sum-from", "minlen", "maxlen", "infix", "prefix", "suffix" modifiers
  const vector<SExpr*> nonterm_blocks = grammar_sexpr.find_all (EG_NONTERMINAL, 1);
  for_const_contents (vector<SExpr*>, nonterm_blocks, nt_ptr)
    {
      SExpr& nt_sexpr = **nt_ptr;
      const sstring nt_name = nt_sexpr(EG_NAME).get_atom();
      const int state = nonterm2state[nt_name];
      ECFG_state_info& info = ecfg->state_info[state];

      if (nt_sexpr.find (EG_TRANSFORM_SUM_FROM))
	info.sum_state = true;

      if (nt_sexpr.find (EG_TRANSFORM_MINLEN))
	info.min_len = nt_sexpr(EG_TRANSFORM_MINLEN).get_atom().to_int();

      if (nt_sexpr.find (EG_TRANSFORM_INFIX))
	info.infix = true;
      else if (nt_sexpr.find (EG_TRANSFORM_MAXLEN))
	info.max_len = nt_sexpr(EG_TRANSFORM_MAXLEN).get_atom().to_int();

      if (nt_sexpr.find (EG_TRANSFORM_PREFIX))
	info.prefix = true;

      if (nt_sexpr.find (EG_TRANSFORM_SUFFIX))
	info.suffix = true;

      // initialise indel params
      init_gaps (info, sym2pvar, &nt_sexpr);

      // initialise GFFs
      const vector<SExpr*> all_gff_sexpr = nt_sexpr.find_all (EG_GFF, 1);
      for_const_contents (vector<SExpr*>, all_gff_sexpr, gff_sexpr)
	init_gff (ecfg, info, *gff_sexpr);
    }

  // add start transition to first state (slightly inefficient, this)
  if (ecfg->has_parametric_transitions)
    ecfg->trans_funcs.transition (Start, 0) = 1.;
  else
    ecfg->transition (Start, 0) = 0;

  // flag GFF annotations
  const vector<SExpr*> all_gff_sexpr = grammar_sexpr.find_all (EG_GFF, 1);
  for_const_contents (vector<SExpr*>, all_gff_sexpr, gff_sexpr)
    {
      const sstring& nonterm_name = (**gff_sexpr) (EG_GFF_NONTERM).get_atom();
      if (nonterm2state.find (nonterm_name) == nonterm2state.end())
	THROWEXPR ("In " << *gff_sexpr << ": nonterminal " << nonterm_name << " not found");
      ECFG_state_info& info = ecfg->state_info[nonterm2state[nonterm_name]];

      init_gff (ecfg, info, *gff_sexpr);
    }

  // initialise update_rates and update_rules
  SExpr* update_rates_sexpr = grammar_sexpr.find (EG_UPDATE_RATES, 1);
  if (update_rates_sexpr && update_rates_sexpr->has_value())
    ecfg->update_rates = update_rates_sexpr->value().get_atom().to_int();
  else
    ecfg->update_rates = true;

  SExpr* update_rules_sexpr = grammar_sexpr.find (EG_UPDATE_RULES, 1);
  if (update_rules_sexpr && update_rules_sexpr->has_value())
    ecfg->update_rules = update_rules_sexpr->value().get_atom().to_int();
  else
    ecfg->update_rules = true;

  // evaluate PFunc's, in case any subsequent routine tries to use ECFG_scores
  ecfg->eval_funcs();

  // return
  return ecfg;
}

void ECFG_builder::init_gff (ECFG_scores* ecfg, ECFG_state_info& info, SExpr* gff_sexpr)
{
  GFF gff;

  if (gff_sexpr->find (EG_GFF_SOURCE, 1))
    gff.source = (*gff_sexpr) (EG_GFF_SOURCE).get_atom();
  else
    gff.source = ecfg->name;  // default value for "source" field is grammar name

  if (gff_sexpr->find (EG_GFF_TYPE, 1))
    gff.feature = (*gff_sexpr) (EG_GFF_TYPE).get_atom();
  else
    gff.feature = info.name;  // default value for "type" field is nonterminal name

  if (gff_sexpr->find (EG_GFF_STRAND, 1))
    gff.strand = GFF::string2strand ((*gff_sexpr) (EG_GFF_STRAND).get_atom());

  if (gff_sexpr->find (EG_GFF_FRAME, 1))
    gff.frame = GFF::string2frame ((*gff_sexpr) (EG_GFF_FRAME).get_atom());

  info.gff.push_back (gff);
}

void ECFG_builder::init_grammars (Alphabet& alph, vector<ECFG_scores*>& ecfgs, SExpr& grammars_sexpr, const Tree_alignment_database* align_db, double tres)
{
  // expand includes
  SExpr_file_operations file_ops;
  file_ops.preorder_visit (grammars_sexpr);

  // initialise alphabet
  init_alphabet (alph, grammars_sexpr.find_or_die (PK_ALPHABET));

  // expand remaining macros
  expand_macros (grammars_sexpr, alph, align_db);

  // init all grammars
  const vector<SExpr*> all_grammar_sexpr = grammars_sexpr.find_all (EG_GRAMMAR);
  for_const_contents (vector<SExpr*>, all_grammar_sexpr, ecfg_sexpr)
    ecfgs.push_back (init_ecfg (alph, **ecfg_sexpr, tres));
}

void ECFG_builder::load_xgram_alphabet_and_grammars (const sstring& filename, Alphabet& alph, vector<ECFG_scores*>& ecfgs, const Tree_alignment_database* align_db, int max_subseq_len, double tres)
{
  // open grammars file & read it
  ifstream ecfg_sexpr_file (filename.c_str());
  if (!ecfg_sexpr_file)
    THROWEXPR ("Grammar/alphabet file not found: '" << filename << "'");
  sstring ecfg_sexpr_string, s;
  while (ecfg_sexpr_file && !ecfg_sexpr_file.eof())
    {
      s.getline (ecfg_sexpr_file);
      ecfg_sexpr_string << s;
    }
  ecfg_sexpr_file.close();

  // create SExpr
  SExpr ecfg_sexpr (ecfg_sexpr_string.begin(), ecfg_sexpr_string.end());

  // delegate
  load_xgram_alphabet_and_grammars (ecfg_sexpr, alph, ecfgs, align_db, max_subseq_len, tres);
}

void ECFG_builder::load_xgram_alphabet_and_grammars (SExpr& ecfg_sexpr, Alphabet& alph, vector<ECFG_scores*>& ecfgs, const Tree_alignment_database* align_db, int max_subseq_len, double tres)
{
  // initialise ECFGs & alphabet
  init_grammars (alph, ecfgs, ecfg_sexpr, align_db, tres);

  // set infix lengths
  for_contents (vector<ECFG_scores*>, ecfgs, ecfg)
    (*ecfg)->set_infix_len (max_subseq_len);
}

void ECFG_builder::expand_macros (SExpr& grammars_sexpr, const Alphabet& alph, const Tree_alignment_database* align_db)
{
  SExpr_macros macros;

  vector<sstring> tokens;
  for (int sym = 0; sym < alph.size(); ++sym)
    {
      sstring tok;
      tok << alph.int2char (sym);
      tokens.push_back (tok);
    }
  macros.foreach[sstring (EG_FOREACH_TOKEN)] = tokens;

  sstring tokens_str;
  tokens_str << alph.size();
  macros.replace[sstring (EG_TOKENS)] = tokens_str;

  if (align_db)
    {
      if (align_db->size() > 0)
	{
	  if (align_db->size() > 1)
	    CTAG(5,ECFG) << "Warning: more than one alignment in database. Using tree from first alignment for tree-related macros.\n";

	  const Alignment& align = align_db->tree_align[0]->align;
	  const PHYLIP_tree& tree = align_db->tree_align[0]->tree;

	  sstring nodes_str, leaves_str, branches_str, ancestors_str, columns_str;

	  nodes_str << tree.nodes();
	  leaves_str << tree.leaves();
	  ancestors_str << tree.internals();
	  branches_str << tree.branches();
	  columns_str << align.columns();

	  macros.replace[sstring (EG_NODES)] = nodes_str;
	  macros.replace[sstring (EG_BRANCHES)] = branches_str;
	  macros.replace[sstring (EG_LEAVES)] = leaves_str;
	  macros.replace[sstring (EG_ANCESTORS)] = ancestors_str;
	  macros.replace[sstring (EG_COLUMNS)] = columns_str;

	  vector<sstring> nodes, leaves, branches, ancestors;
	  for (int n = 0; n < tree.nodes(); ++n)
	    if (tree.node_name[n].size())
	      {
		const sstring& nn = tree.node_name[n];

		nodes.push_back (nn);
		if (n != tree.root)
		  branches.push_back (nn);

		if (tree.is_leaf (n))
		  leaves.push_back (nn);
		else
		  ancestors.push_back (nn);
	      }

	  macros.foreach[sstring (EG_FOREACH_NODE)] = nodes;
	  macros.foreach[sstring (EG_FOREACH_BRANCH)] = branches;
	  macros.foreach[sstring (EG_FOREACH_LEAF)] = leaves;
	  macros.foreach[sstring (EG_FOREACH_ANCESTOR)] = ancestors;
	}
      else
	CTAG(5,ECFG) << "Warning: since I have no alignments from which to extract trees or other data, tree- and alignment-related macros will be ignored.\n";
    }

  macros.preorder_visit (grammars_sexpr);

  SExpr_list_operations list_ops;
  list_ops.postorder_visit (grammars_sexpr);

  if (CTAGGING(4,ECFG ECFG_MACROS))
    CL << "Grammar S-expression following macro expansion:\n" << grammars_sexpr << '\n';
}

void ECFG_builder::ecfg2stream (ostream& out, const Alphabet& alph, const ECFG_scores& ecfg, const ECFG_counts* counts)
{
  // open grammar block
  out << ";; Grammar " << ecfg.name << "\n;;\n";
  out << '(' << EG_GRAMMAR << '\n';

  // print name
  out << " (" << EG_NAME << ' ' << ecfg.name << ")\n";

  // print meta-info
  if (ecfg.meta.size())
    {
      out << " (" << EG_META;
      for_const_contents (list<SExpr>, ecfg.meta, meta)
	out << "\n  (" << *meta << ")";
      out << ")\n";
    }

  if (ecfg.transient_meta.size())
    {
      out << " (" << EG_TRANSIENT_META;
      for_const_contents (list<SExpr>, ecfg.transient_meta, meta)
	out << "\n  (" << *meta << ")";
      out << ")\n";
    }

  // parametric?
  if (ecfg.has_parametric_transitions)
    out << " (" << EG_PARAMETRIC << ")\n";

  // print update_rules and update_rates
  out << " (" << EG_UPDATE_RATES << ' ' << (ecfg.update_rates ? 1 : 0) << ")\n";
  out << " (" << EG_UPDATE_RULES << ' ' << (ecfg.update_rules ? 1 : 0) << ")\n";

  // Print nonterminal declarations; buffer up transformation rules
  out << "\n ;; Nonterminal declarations and modifiers\n\n";
  sstring trans_block;
  trans_block << "\n ;; Transformation rules for grammar symbols\n";

  // if the start state has transitions to multiple states, then add a dummy null state & print start rules
  // otherwise, sort the states beginning with the first accessible nonterminal from the start state
  vector<int> state_order (ecfg.states()), starts;
  for (int s = 0; s < ecfg.states(); ++s)
    {
      const bool is_start =
	ecfg.has_parametric_transitions
	? !ecfg.trans_funcs.transition (Start, s).is_null()
	: ecfg.transition (Start, s) > -InfinityScore;

      state_order[s] = s;
      if (is_start)
	starts.push_back (s);
    }

  const int n_starts = starts.size();
  const bool start_to_end =
    ecfg.has_parametric_transitions
    ? !ecfg.trans_funcs.transition (Start, End).is_null()
    : ecfg.transition (Start, End) > -InfinityScore;

  if (n_starts == 0 && !start_to_end)
    {
      THROWEXPR ("Start state is isolated, very strange grammar.... refusing to print it... meltdown...");
    }
  else if (n_starts == 1 && !start_to_end)
    swap (state_order[0], state_order[starts[0]]);
  else
    {
      // right now we assume there is no other state named ECFG_default_start_nonterminal.
      // this is pretty fragile...
      const sstring start_name (ECFG_default_start_nonterminal);  // could do a bit more work here, to ensure name is unique
      out << " (" << EG_NONTERMINAL << " (" << EG_NAME << ' ' << start_name << "))\n";
      trans_block << "\n ;; Nonterminal " << start_name << "\n ;;\n";
      if (start_to_end)
	starts.push_back (End);
      for_const_contents (vector<int>, starts, s_iter)
	{
	  trans_block << " ("
		      << EG_TRANSFORM << " ("
		      << EG_FROM << " (" << start_name << ")) ("
		      << EG_TO << " ("
		      << (*s_iter == End ? "" : ecfg.state_info[*s_iter].name.c_str())
		      << ")) (" << EG_PROB << ' ';
	  if (ecfg.has_parametric_transitions)
	    pfunc2stream (trans_block, ecfg.pscores, ecfg.trans_funcs.transition (Start, *s_iter));
	  else
	    trans_block << Score2Prob (ecfg.transition (Start, *s_iter));
	  trans_block << "))\n";
	}
    }

  // OK, Start state is done; now print other states
  for_const_contents (vector<int>, state_order, s_iter)
    {
      const int s = *s_iter;
      const ECFG_state_info& info = ecfg.state_info[s];

      out << " (" << EG_NONTERMINAL << " (" << EG_NAME << ' ' << info.name << ')';
      if (info.min_len > info.emit_size())
	out << " (" << EG_TRANSFORM_MINLEN << ' ' << info.min_len << ')';
      if (info.infix)
	out << " (" << EG_TRANSFORM_INFIX << ')';
      else if (info.max_len >= 0)
	out << " (" << EG_TRANSFORM_MAXLEN << ' ' << info.max_len << ')';
      if (info.prefix)
	out << " (" << EG_TRANSFORM_PREFIX << ')';
      if (info.suffix)
	out << " (" << EG_TRANSFORM_SUFFIX << ')';
      if (info.sum_state)
	out << " (" << EG_TRANSFORM_SUM_FROM << ')';
      for_const_contents (list<GFF>, info.gff, gff)
	{
	  out << " (" << EG_GFF;
	  if (gff->source != ecfg.name)
	    out << " (" << EG_GFF_SOURCE << ' ' << gff->source << ')';
	  if (gff->feature != info.name)
	    out << " (" << EG_GFF_TYPE << ' ' << gff->feature << ')';
	  if (gff->strand != GFF_enum::NoStrand)
	    out << " (" << EG_GFF_STRAND << ' ' << (gff->strand > 0 ? '+' : '-') << ')';
	  if (gff->frame != GFF_enum::NoFrame)
	    out << " (" << EG_GFF_FRAME << ' ' << gff->frame << ')';
	  out << ')';
	}
      out << ")\n";

      trans_block << "\n ;; Nonterminal " << info.name << "\n ;;\n";
      if (info.bifurc)
	{
	  trans_block << " ("
		      << EG_TRANSFORM << " ("
		      << EG_FROM << " (" << info.name << ")) ("
		      << EG_TO << " (" << ecfg.state_info[info.ldest].name << ' ' << ecfg.state_info[info.rdest].name
		      << ")))\n";
	}
      else  // emit or null
	{
	  if (info.emit_size())
	    {
	      trans_block << " (" << EG_TRANSFORM << " (" << EG_FROM << " (";
	      const ECFG_chain& chain = ecfg.matrix_set.chain[info.matrix];
	      // figure out pos2term map
	      vector<int> pos2term (chain.word_len, -1);
	      for (int mul = 1, i = 0; i < chain.word_len; mul *= alph.size(), ++i)
		for (int k = 0; k < (int) info.mul.size(); ++k)
		  if (info.mul[k] == mul)
		    pos2term[k] = i;
	      // print LHS of emit rule
	      sstring sep;
	      for (int i = 0; i < info.l_context; ++i)
		{
		  trans_block << sep << chain.state[pos2term[i]];
		  sep = " ";
		}
	      trans_block << sep << info.name;
	      for (int i = info.r_context - 1; i >= 0; --i)
		trans_block << ' ' << chain.state[pos2term[chain.word_len - 1 - i]];
	      trans_block << ")) (" << EG_TO << " (";
	      // print RHS of emit rule
	      sep = "";
	      for (int i = 0; i < info.l_context + info.l_emit; ++i)
		{
		  trans_block << sep << chain.state[pos2term[i]];
		  sep = " ";
		}
	      trans_block << sep << info.name << ECFG_post_emit_character;
	      for (int i = info.r_emit + info.r_context - 1; i >= 0; --i)
		trans_block << ' ' << chain.state[pos2term[chain.word_len - 1 - i]];
	      trans_block << "))";

	      // gap model
	      if (!info.gaps_ok)
		trans_block << " (" << EG_TRANSFORM_NO_GAPS << ")";
	      else if (info.indels)
		{
		  if (!info.wild_gaps)
		    trans_block << " (" << EG_TRANSFORM_STRICT_GAPS << ")";
		  trans_block << "\n  ("
		      << EG_TRANSFORM_GAP_MODEL << " (";
		  if (info.has_parametric_indels)
		    {
		      trans_block << EG_PARAMETRIC << ") (" << EG_TRANSFORM_EXTEND_PROB << ' ';
		      pfunc2stream (trans_block, ecfg.pscores, info.link_extend_func);
		      if (counts)
			print_count (trans_block, counts->link_extend_count[s]);

		      trans_block << ") (" << EG_TRANSFORM_END_PROB << ' ';
		      pfunc2stream (trans_block, ecfg.pscores, info.link_end_func);
		      if (counts)
			print_count (trans_block, counts->link_end_count[s]);

		      trans_block << ") (" << EG_TRANSFORM_INS_RATE << ' ';
		      pfunc2stream (trans_block, ecfg.pscores, info.ins_rate_func);
		      if (counts)
			{
			  print_count (trans_block, counts->ins_count[s]);
			  print_time (trans_block, counts->ins_wait[s]);
			}

		      trans_block << ") (" << EG_TRANSFORM_DEL_RATE << ' ';
		      pfunc2stream (trans_block, ecfg.pscores, info.del_rate_func);
		      if (counts)
			{
			  print_count (trans_block, counts->del_count[s]);
			  print_time (trans_block, counts->del_wait[s]);
			}

		      trans_block << "))";
		    }
		  else
		    {
		      trans_block << EG_TRANSFORM_EXTEND_PROB << ' ' << info.link_extend;
		      if (counts)
			print_count (trans_block, counts->link_extend_count[s]);

		      trans_block << ") (" << EG_TRANSFORM_END_PROB << ' ' << info.link_end;
		      if (counts)
			print_count (trans_block, counts->link_end_count[s]);

		      trans_block << ") (" << EG_TRANSFORM_INS_RATE << ' ' << info.ins_rate;
		      if (counts)
			{
			  print_count (trans_block, counts->ins_count[s]);
			  print_time (trans_block, counts->ins_wait[s]);
			}

		      trans_block << ") (" << EG_TRANSFORM_DEL_RATE << ' ' << info.del_rate;
		      if (counts)
			{
			  print_count (trans_block, counts->del_count[s]);
			  print_time (trans_block, counts->del_wait[s]);
			}

		      trans_block << "))";
		    }
		}
	      else if (info.wild_gaps)
		trans_block << " (" << EG_TRANSFORM_IGNORE_GAPS << ")";
	      else   // !info.wild_gaps && !info.indels
		trans_block << " (" << EG_TRANSFORM_STRICT_GAPS << ")";

	      // annotation
	      for_const_contents (ECFG_state_info::ECFG_state_annotation, info.annot, row_annot)
		{
		  const sstring& row = row_annot->first;
		  if (row_annot->second.size() == 1 && row_annot->second.begin()->second.is_one()) {
		    // deterministic
		    const sstring& label = row_annot->second.begin()->first;
		    for (int i = 0; i < (int) label.size(); ++i)
		      if (label[i] != ECFG_annotation_wildcard)
			trans_block << "\n  ("
				    << EG_TRANSFORM_ANNOTATE << " ("
				    << EG_TRANSFORM_ROW << ' ' << row << ") ("
				    << EG_TRANSFORM_COLUMN << ' ' << chain.state[pos2term[info.l_context + i]] << ") ("
				    << EG_TRANSFORM_LABEL << ' ' << label[i] << "))";
		  } else {
		    // probabilistic
		    trans_block << "\n  ("
				<< EG_TRANSFORM_ANNOTATE << " ("
				<< EG_TRANSFORM_ROW << ' ' << row << ")";
		    for_const_contents (ECFG_state_info::String_prob_dist, row_annot->second, label_prob) {
		      const sstring& label = label_prob->first;
		      const PFunc& prob = label_prob->second;
		      trans_block << "\n   (" << EG_ANNOTATE_EMIT << " (" << EG_TRANSFORM_LABEL << " (";
		      for (int i = 0; i < (int) label.size(); ++i) {
			if (i > 0)
			  trans_block << ' ';
			trans_block << label[i];
		      }
		      trans_block << ")) (" << EG_PROB << ' ';
		      pfunc2stream (trans_block, ecfg.pscores, prob);
		      trans_block << "))";
		    }
		    trans_block << ")";
		  }
		}

	      // end of emit rule
	      trans_block << ")\n";
	    }

	  // print transitions
	  for (int d = End; d < ecfg.states(); ++d)
	    if (d != Start)
	      {
		const bool print_transition =
		  ecfg.has_parametric_transitions
		  ? !ecfg.trans_funcs.transition (s, d).is_null()
		  : ecfg.transition (s, d) > -InfinityScore;

		if (print_transition)
		  {
		    trans_block << " (" << EG_TRANSFORM << " (" << EG_FROM << " (" << info.name;
		    if (info.emit_size())  // add post-emit character
		      trans_block << ECFG_post_emit_character;
		    trans_block << ")) (" << EG_TO << " (";
		    if (d != End)
		      trans_block << ecfg.state_info[d].name;

		    trans_block << ")) (" << EG_PROB << ' ';
		    if (ecfg.has_parametric_transitions)
		      pfunc2stream (trans_block, ecfg.pscores, ecfg.trans_funcs.transition (s, d));
		    else
		      trans_block << Score2Prob (ecfg.transition (s, d));
		    if (counts)
		      print_count (trans_block, counts->transition (s, d));
		    trans_block << "))\n";
		  }
	      }
	}
    }

  // print transformation rule block
  out << trans_block;

  // print PScores & PCounts
  if (ecfg.pscores.groups())
    {
      out << "\n ;; Probability and rate parameters\n";
      pscores2stream (out, ecfg.pscores, ecfg.mutable_pgroups, counts ? &counts->var_counts : (const PCounts*) 0);
    }
  if (ecfg.pcounts.groups())
    {
      out << "\n ;; Pseudocounts for probability and rate parameters\n\n";
      pcounts2stream (out, ecfg.pcounts, EG_PSEUDOCOUNTS, (const PCounts*) 0, true, false);
    }
  if (counts && ecfg.update_rates)
    {
      out << "\n ;; Expected counts and wait times for probability and rate parameters\n\n";
      pcounts2stream (out, counts->var_counts, EG_EXPECTED_COUNTS, ecfg.pcounts.groups() ? &ecfg.pcounts : (const PCounts*) 0, true, true);

      out << "\n ;; Expected counts and wait time for Markov chain substitution models\n\n";
      out << " (" << EG_CHAIN_COUNTS << '\n';
      chain_counts2stream (out, alph, ecfg, *counts);
      out << " )  ;; end " << EG_CHAIN_COUNTS << '\n';
    }

  // print chains
  out << "\n ;; Markov chain substitution models\n";
  for (int n_chain = 0; n_chain < (int) ecfg.matrix_set.chain.size(); ++n_chain)
    {
      const ECFG_chain& chain = ecfg.matrix_set.chain[n_chain];
      if (chain.type == Hybrid)
	{
	  // hybrid-chain
	  out << "\n (" << EG_HYBRID_CHAIN << '\n';
	  out << "  (" << EG_CHAIN_TERMINAL << " (" << chain.state << "))\n";
	  if (chain.classes > 1)
	    {
	      out << "  (" << EG_CHAIN_CLASS << "\n";
	      out << "   (" << EG_TRANSFORM_ROW << ' ' << chain.class_row << ")\n";
	      out << "   (" << EG_TRANSFORM_LABEL << " (" << chain.class_labels << ")))\n";
	    }
	  out << "  (" << EG_TRANSFORM_ROW << ' ' << chain.gs_tag << ")\n";

	  out << "  (" << EG_HYBRID_COMPONENTS << '\n';
	  for_const_contents (vector<sstring>, chain.gs_values, gs_value)
	    out << "   ((" << EG_TRANSFORM_LABEL << ' ' << *gs_value
		<< ") (" << EG_CHAIN_TERMINAL << " (" << ecfg.matrix_set.chain[((map<sstring,int>&) chain.gs_tag_value_chain_index)[*gs_value]].state  // cast away const
		<< ")))\n";

	  out << "  )\n";
	  out << " )  ;; end hybrid chain " << chain.state << "\n";
	}
      else
	{
	  // regular chain
	  const Update_statistics* stats = counts ? &counts->stats[n_chain] : (const Update_statistics*) 0;

	  if (chain.matrix->C > 1)
	    THROWEXPR ("Don't know how to display EM_matrix objects with hidden states");
	  out << "\n (" << EG_CHAIN << '\n';
	  out << "  (" << EG_CHAIN_POLICY << ' ' << chain.matrix->update_policy() << ")\n";
	  out << "  (" << EG_CHAIN_TERMINAL << " (" << chain.state << "))\n";
	  if (chain.classes > 1)
	    {
	      out << "  (" << EG_CHAIN_CLASS << "\n";
	      out << "   (" << EG_TRANSFORM_ROW << ' ' << chain.class_row << ")\n";
	      out << "   (" << EG_TRANSFORM_LABEL << " (" << chain.class_labels << ")))\n";
	    }

	  out << "\n  ;; initial probability distribution\n";
	  for (int s = 0; s < chain.matrix->m(); ++s)
	    {
	      const bool print_pi =
		chain.is_parametric
		? !chain.matrix_funcs->pi[s].is_null()
		: chain.matrix->pi[s] > 0.;
	  
	      if (print_pi)
		{
		  out << "  (" << EG_CHAIN_INITIAL << " (" << EG_CHAIN_STATE << " ";
		  print_state (out, s, chain.word_len, alph, chain.class_labels);
		  out << ") (" << EG_PROB << ' ';

		  if (chain.is_parametric)
		    pfunc2stream (out, ecfg.pscores, chain.matrix_funcs->pi[s]);
		  else
		    out << chain.matrix->pi[s];
		  out << ')';

		  if (stats)
		    print_count (out, stats->s[s]);

		  out << ")\n";
		}
	    }

	  out << "\n  ;; mutation rates\n";
	  for (int s = 0; s < chain.matrix->m(); ++s)
	    {
	      // commented out due to clutter...
#if 0
	      if (stats)
		{
		  out << "  (" << EG_WAIT << " (" << EG_CHAIN_STATE << ' ';
		  print_state (out, s, chain.word_len, alph, chain.class_labels);
		  out << ')';
		  print_time (out, stats->w[s]);
		  out << ")\n";
		}
#endif

	      for (int d = 0; d < chain.matrix->m(); ++d)
		{
		  const bool print_rate =
		    chain.is_parametric
		    ? !chain.matrix_funcs->X[0](s,d).is_null()
		    : (s != d && chain.matrix->X[0](s,d) > 0.);

		  if (print_rate)
		    {
		      out << "  (" << EG_CHAIN_MUTATE << " (" << EG_FROM << " ";
		      print_state (out, s, chain.word_len, alph, chain.class_labels);
		      out << ") (" << EG_TO << " ";
		      print_state (out, d, chain.word_len, alph, chain.class_labels);
		      out << ") (" << EG_RATE << ' ';

		      if (chain.is_parametric)
			pfunc2stream (out, ecfg.pscores, chain.matrix_funcs->X[0](s,d));
		      else
			out << chain.matrix->X[0](s,d);

		      out << ')';

		      if (stats && s != d)
			print_count (out, stats->u(s,d));

		      out << ")\n";
		    }
		}
	    }
	  out << " )  ;; end chain " << chain.state << "\n";
	}
    }

  // close grammar block
  out << "\n)  ;; end grammar " << ecfg.name << "\n\n";
}

void ECFG_builder::print_state (ostream& out, int state, int wordlen, const Alphabet& alph, const vector<sstring>& class_alph)
{
  out << '(';
  int mul, pos;
  for (mul = 1, pos = 0; pos < wordlen; mul *= alph.size(), ++pos)
    {
      if (pos > 0)
	out << ' ';
      out << alph.int2char ((state / mul) % alph.size());
    }
  if (class_alph.size() > 1)
    out << ' ' << class_alph[(state / mul) % class_alph.size()];
  out << ')';
}

void ECFG_builder::grammars2stream (ostream& out, const Alphabet& alph, const vector<ECFG_scores*>& grammars)
{
  const vector<ECFG_trainer*> trainer (grammars.size(), (ECFG_trainer*) 0);
  grammars2stream (out, alph, grammars, &trainer[0]);
}

void ECFG_builder::grammars2stream (ostream& out, const Alphabet& alph, const vector<ECFG_scores*>& grammars, const ECFG_trainer* const* trainer)
{
  for (int n = 0; n < (int) grammars.size(); ++n)
    ecfg2stream (out, alph, *grammars[n], (trainer ? trainer[n] : 0) ? &trainer[n]->counts : (const ECFG_counts*) 0);
  alphabet2stream (out, alph);
}

void ECFG_builder::chain_counts2stream (ostream& out, const Alphabet& alph, const ECFG_scores& ecfg, const ECFG_counts& counts)
{
  for (int n_chain = 0; n_chain < (int) ecfg.matrix_set.chain.size(); ++n_chain)
    {
      const ECFG_chain& chain = ecfg.matrix_set.chain[n_chain];
      if (chain.type != Hybrid)
	if (chain.matrix->C == 1)  // don't barf here if C>1; leave it to chain display code
	  {
	    const Update_statistics& stats = counts.stats[n_chain];

	    out << "  ((" << EG_CHAIN_TERMINAL << " (" << chain.state << "))\n";

	    for (int s = 0; s < chain.matrix->m(); ++s)
	      {
		out << "   (" << EG_CHAIN_INITIAL << " (" << EG_CHAIN_STATE << " ";
		print_state (out, s, chain.word_len, alph, chain.class_labels);
		out << ") (" << PK_COUNT << ' ' << stats.s[s] << "))\n";
	      }

	    for (int s = 0; s < chain.matrix->m(); ++s)
	      {
		out << "   (" << EG_WAIT << " (" << EG_CHAIN_STATE << ' ';
		print_state (out, s, chain.word_len, alph, chain.class_labels);
		out << ") (" << PK_TIME << ' ' << stats.w[s] << "))\n";

		for (int d = 0; d < chain.matrix->m(); ++d)
		  if (s != d)
		    {
		      out << "   (" << EG_CHAIN_MUTATE << " (" << EG_FROM << " ";
		      print_state (out, s, chain.word_len, alph, chain.class_labels);
		      out << ") (" << EG_TO << " ";
		      print_state (out, d, chain.word_len, alph, chain.class_labels);
		      out << ") (" << PK_COUNT << ' ' << stats.u(s,d) << "))\n";
		    }
	      }

	    out << "  )  ;; end chain " << chain.state << "\n";
	  }
    }
}
