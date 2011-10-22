#ifndef PSEXPR_INCLUDED
#define PSEXPR_INCLUDED

#include "util/sexpr.h"

#include "seq/pvar.h"
#include "seq/pfunc.h"

// class to build PScores, PFunc and Alphabet objects from an S-expression (SExpr)
struct PFunc_builder
{
  // SExprIter
  typedef SExpr::SExprIter SExprIter;

  // symbol-to-index map; throws an exception if operator() called for a nonexistent key
  struct SymIndex : map<sstring,int>
  {
    sstring symbol_type;  // initialized by constructor, used for error message only
    SymIndex (const char* st = "Symbol");
    int operator() (const sstring& key, const SExpr& containing_sexpr) const;
  };

  // symbol-to-PVar map
  typedef map<sstring,PVar> SymPVar;

  // output
  // PScores, PFunc, PCounts
  static void pscores2stream (ostream& out, const PScores& pscores, const PCounts* pcounts = 0, bool use_bitscores = false);
  static void pscores2stream (ostream& out, const PScores& pscores, const char* tag, const vector<int>& pgroups_to_show, const char* indent = "", const PCounts* pcounts = 0, bool use_bitscores = false);
  static void pscores2stream (ostream& out, const PScores& pscores, const set<int>& mutable_pgroups, const PCounts* pcounts = 0, bool use_bitscores = false);  // uses mutable_pgroups to figure out const & pgroup/rate tags

  static sstring mutable_pscores2string (const PScores& pscores, const set<int>& mutable_pgroups);

  static void pcounts2stream (ostream& out, const PCounts& pcounts, const char* tag, const PCounts* baseline_pcounts = 0,
			      bool interpret_single_element_pgroups_as_rate_variables = true,
			      bool print_zero_counts = false);

  static void pfunc2stream (ostream& out, const PScores& pscores, const PFunc& pfunc);

  static sstring score_sexpr (Score sc, bool use_bitscores);
  static sstring score_sexpr (FScore sc, bool use_bitscores);

  // Alphabet
  static void alphabet2stream (ostream& out, const Alphabet& alph);

  // input
  // PScores and PFuncs
  static PFunc init_pfunc (const SymPVar& sym2pvar, SExpr& pfunc_sexpr, int offset = 0);
  static PGroup init_pgroup (PScores& pscores, SymPVar& sym2pvar, SExpr& pgroup_sexpr, set<int>* mutable_pgroups = 0, bool force_rate = false, bool disallow_rate = false, bool use_bitscores = false);
  // init_pgroups: mutable_pgroups can be used to record indices of all PGroup's with this tag. Typically this is used to distinguish mutable from constant pgroups (e.g. in xrate)
  static void init_pgroups (PScores& pscores, SymPVar& sym2pvar, SExpr& grammar_sexpr, const char* tag, set<int>* mutable_pgroups = 0, bool force_rate = false, bool disallow_rate = false, bool use_bitscores = false);
  // wrapper method for const/mutable probs/rates
  static void init_pgroups_and_rates (PScores& pscores, SymPVar& sym2pvar, SExpr& sexpr, set<int>* mutable_pgroups = 0, bool use_bitscores = false);

  // initialise an Alphabet (alph) from an EG_ALPHABET expression (alphabet_sexpr)
  static void init_alphabet (Alphabet& alph, SExpr& alphabet_sexpr, bool allow_multi_char_tokens = false);

  // helpers
  static int get_state_class (const sstring& atom, const vector<sstring>& class_alph);

  static void assert_valid_token (const sstring& atom, const SExpr& sexpr);
  static void assert_valid_token (const sstring& atom, const Alphabet& alph, const SExpr& sexpr);

  static sstring token_list_to_string (SExpr& token_list, int offset = 0);
  static sstring token_list_to_string (SExpr& token_list, const Alphabet& alph, const vector<sstring>& class_alph, int offset = 0);

  static void print_count (ostream& out, double count);
  static void print_time (ostream& out, double t);
};


#endif /* PSEXPR_INCLUDED */
