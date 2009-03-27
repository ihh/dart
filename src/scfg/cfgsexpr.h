#ifndef CFG_SEXPR_INCLUDED
#define CFG_SEXPR_INCLUDED

#include "util/sexpr.h"
#include "hmm/pairphmm.h"
#include "scfg/pairpcfg.h"
#include "scfg/cfgkeywords.h"
#include "seq/psexpr.h"

// Stemloc grammar set typedefs
struct PHMM_CFG
{
  // data
  Pair_PHMM hmm;  // pair HMM used for pre-aligning
  Odds_PCFG cfg;  // pair SCFG used for aligning-and-folding
  Alphabet_group hmm_null_emit;  // null model for HMM (we don't have an Odds_PHMM class to hold this, ugh)
  set<int> hmm_pad_states;  // non-pad states in the HMM (again, no class to hold this; formerly held in Quick_align)
  PScores hmm_pscores, cfg_pscores;
  set<int> hmm_mutable_pgroups, cfg_mutable_pgroups;
  Dirichlet_prior hmm_prior, cfg_prior;
  // constructor
  PHMM_CFG() : hmm(0,CFG_alphabet) { }
};
typedef map<sstring,PHMM_CFG> PHMM_CFG_map;

// Stemloc grammar set
struct Gramset
{
  Odds_PCFG single_scfg;   // SCFG used for pre-folding
  PScores single_scfg_pscores;
  Dirichlet_prior single_scfg_prior;
  set<int> single_scfg_mutable_pgroups;

  PHMM_CFG_map phmm_cfg_map;  // pair HMMs for pre-aligning and corresponding pair SCFGs for pre-folding

  void clear();  // reset everything
};

// class to build an RNA Pair SCFG (Pair_PCFG), or a Gramset, from an S-expression (SExpr)
struct PCFG_builder : Pair_CFG_state_type_enum, Grammar_state_enum, PFunc_builder
{
  // SExprIter
  typedef SExpr::SExprIter SExprIter;

  // method to read the null_emit and null_extend PGroup's
  static void init_null_pgroups (SExpr& model_sexpr, PScores& pscores, SymPVar& sym2pvar, Alphabet_group& null_emit, Boolean_group& null_extend, bool want_HMM = false);

  // method to read an SCFG
  static Odds_PCFG init_grammar (SExpr& model_sexpr, PScores& pscores, set<int>& pad_states_ret, set<int>& mutable_pgroups_ret, bool want_single = false, bool want_HMM = false);

  // method to read an HMM
  static Pair_PHMM init_hmm (SExpr& model_sexpr, PScores& pscores, Alphabet_group& null_emit_group_ret, set<int>& pad_states_ret, set<int>& mutable_pgroups_ret, bool want_single = false);

  // method to write an SCFG
  static void grammar2stream (ostream& out, const Pair_PCFG& pcfg, const PScores& pscores, const set<int>& mutable_pgroups, const char* tag = CFG_PAIR_SCFG, const char* indent = "", int null_emit_group_index = -1, int null_extend_group_index = -1, const set<int>* pad_states = 0);

  // method to write an HMM
  static void hmm2stream (ostream& out, const Pair_PHMM& phmm, const PScores& pscores, const set<int>& mutable_pgroups, const char* tag = CFG_PAIR_HMM, const char* indent = "", int null_emit_group_index = -1, const set<int>* pad_states = 0);

  // method to read a Gramset
  static void init_gramset (SExpr& grammar_set_sexpr, Gramset& gramset);

  // method to write a Gramset
  static void gramset2stream (ostream& out, const Gramset& gramset, const char* tag = CFG_GRAMSET);

  // write a short program to save sldefaults.h as gramset files

  // next, convert Stemloc to use this class:
  // -- use string name to select param_set, instead of integer index
  // -- use I/O methods from this class
  // AND/OR write another short program to convert Telegraph format grammars into the new format

  // finally, test! (loading, saving & training)
};


#endif /* CFG_SEXPR_INCLUDED */
