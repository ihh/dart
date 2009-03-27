#ifndef ECFG_SEXPR
#define ECFG_SEXPR

#include "util/sexpr.h"
#include "ecfg/ecfg.h"
#include "ecfg/ecfgdp.h"
#include "ecfg/ecfgkeywords.h"

#include "seq/psexpr.h"

// class to build an ECFG (ECFG_scores) from an S-expression (SExpr)
struct ECFG_builder : ECFG_enum, PFunc_builder
{
  // output methods
  static void ecfg2stream (ostream& out, const Alphabet& alph, const ECFG_scores& ecfg, const ECFG_counts* counts = 0);
  static void chain_counts2stream (ostream& out, const Alphabet& alph, const ECFG_scores& ecfg, const ECFG_counts& counts);

  static void grammars2stream (ostream& out, const Alphabet& alph, const vector<ECFG_scores*>& grammars);
  static void grammars2stream (ostream& out, const Alphabet& alph, const vector<ECFG_scores*>& grammars, const ECFG_trainer* const* trainer);

  // top-level init method
  // initialise an Alphabet and a bunch of grammars
  static void init_grammars (Alphabet& alph, vector<ECFG_scores*>& ecfgs, SExpr& grammars_sexpr, const Tree_alignment_database* align_db = 0, double tres = DEFAULT_TIMEPOINT_RES);
  static void load_xgram_alphabet_and_grammars (const sstring& filename, Alphabet& alph, vector<ECFG_scores*>& ecfgs, const Tree_alignment_database* align_db = 0, int max_subseq_len = -1, double tres = DEFAULT_TIMEPOINT_RES);
  static void expand_macros (SExpr& grammars_sexpr, const Alphabet& alph, const Tree_alignment_database* align_db = 0);

  // initialise a single ECFG_scores from an EG_GRAMMAR expression
  static ECFG_scores* init_ecfg (const Alphabet& alph, SExpr& pgroup_sexpr, double tres = DEFAULT_TIMEPOINT_RES);

  // initialise an ECFG_chain in an ECFG_matrix_set (ems) from an EG_CHAIN expression (chain_sexpr),
  // updating a symbol-to-chain-index map (term2chain)
  static void init_chain (ECFG_matrix_set& ems, SymIndex& term2chain, const SymPVar& sym2pvar, SExpr& chain_sexpr, double tres = DEFAULT_TIMEPOINT_RES);

  // initialize hidden class labels for a hybrid chain
  static void init_chain_classes (sstring& class_row, vector<sstring>& class_alph, const int n_terminals, SExpr& chain_sexpr);

  // initialize a hybrid chain
  static void init_hybrid_chain (ECFG_matrix_set& ems, SymIndex& term2chain, const SymPVar& sym2pvar, SExpr& hybrid_chain_sexpr);

  // initialise a chain and an alphabet together. This routine is used by Handel
  static void init_chain_and_alphabet (Alphabet& alph, EM_matrix_base& hsm, SExpr& alph_chain_sexpr);

  // initialise a chain given an alphabet. Used by evoldoer
  static void init_chain_given_alphabet (EM_matrix_base& hsm, const Alphabet& alph, SExpr& chain_sexpr, int required_pseudoterms = 1);

  // initialise nonterminal symbols
  static SymIndex init_nonterm2state (const SymIndex& term2chain, SExpr& grammar_sexpr);

  // initialise gap model for a state
  // NB the constructor for ECFG_state_info defaults to "gaps-ok"
  static void init_gaps (ECFG_state_info& info, const SymPVar& sym2pvar, SExpr* sexpr);

  // initialise a GFF block
  static void init_gff (ECFG_scores* ecfg, ECFG_state_info& info, SExpr* gff_sexpr);

  // structure describing a symbol sequence (LHS or RHS of a rule)
  struct ECFG_symbol_sequence
  {
    // typedefs
    typedef map<int,int> IntIntMap;

    // data
    vector<bool> is_term;  // position-to-terminal-flag map
    vector<int> sym;   // position-to-symbol map (both terminals and nonterminals)
    IntIntMap term_pos;  // symbol-to-position map for terminals
    IntIntMap pos_nonterm;  // position-to-symbol map for nonterminals
    map<int,bool> pos_postemit;  // position-to-postemit-flag map for nonterminals (true if nonterminal is post-emit, e.g. "X'")
    map<int,sstring> pos_name;  // position-to-name map for nonterminals
    map<int,bool> pos_comped;  // position-to-complemented-flag map for terminals (true if term is "~S" rather than "S")
    int chain_index;
    const ECFG_chain* chain;
    SExpr* sexpr;

    // constructors
    // default constructor
    ECFG_symbol_sequence() { }
    // initialise an ECFG_symbol_sequence from a symbol sequence SExpr (symseq_sexpr)
    ECFG_symbol_sequence (const ECFG_matrix_set& ems, const SymIndex& nonterm2state, const SymIndex& term2chain, SExpr& symseq_sexpr);
  };

  // structure describing a "(from ...)" block
  struct ECFG_rule_block
  {
    // data
    ECFG_symbol_sequence lhs, rhs;
    const SymIndex& nonterm2state;
    const SymIndex& term2chain;
    bool parametric_transitions;
    const SymPVar& sym2pvar;
    SExpr* sexpr;

    // constructor
    // initialise an ECFG_rule_block from a from SExpr (from_sexpr)
    ECFG_rule_block (const ECFG_matrix_set& ems,
		     const SymIndex& nonterm2state,
		     const SymIndex& term2chain,
		     bool parametric_transitions,
		     const SymPVar& sym2pvar,
		     SExpr& transform_sexpr);

    // parse method
    void parse (ECFG_scores& ecfg);
  };

  // parser helper methods
  static int token_list_to_state (SExpr& token_list, const Alphabet& alph, int word_len, const vector<sstring>& class_alph);

  static const char* policy2string (Update_policy type);
  static Update_policy string2policy (const sstring& policy_string);

  static void print_state (ostream& out, int state, int wordlen, const Alphabet& alph, const vector<sstring>& class_alph);
};

#endif /* ECFG_SEXPR */

