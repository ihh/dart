#ifndef ONTO_SEXPR_INCLUDED
#define ONTO_SEXPR_INCLUDED

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include <libguile.h>
#else
#error Terminatrix requires guile! Please install guile from http://www.gnu.org/software/guile/ and re-run the configure script.
#endif

#include "ecfg/ecfgsexpr.h"
#include "seq/psexpr.h"
#include "util/svisitor.h"


// Terminatrix class
struct Terminatrix
{
  // primary object interface to Scheme
  SExpr_Scheme_evaluator& scheme;

  // SExpr's
  SExpr init_sexpr, model_sexpr, tree_db_sexpr, knowledge_sexpr;

  // symbol tables
  PFunc_builder::SymPVar sym2pvar;
  PFunc_builder::SymIndex term2chain;

  // params
  PScores pscores;
  PCounts pcounts, var_counts;
  set<int> mutable_pgroups;

  // alphabet namespace
  // (alph_list is a bit wasteful; we already store every Alphabet in alph_dict. Could just store a list of Alphabet names instead...)
  const Alphabet dummy_alph;  // the ECFG_matrix_set needs a default Alphabet. Keep this around to pacify it
  list<Alphabet> alph_list;
  Alphabet_dictionary alph_dict;

  // continuous-time Markov chain
  ECFG_matrix_set matrix_set;

  // results
  Update_statistics stats;

  // status variables
  bool got_counts;  // true if model has been trained

  // constructor
  Terminatrix (SExpr_Scheme_evaluator& scheme);  // the SExpr_Scheme_evaluator is assumed to have been initialized

  // helpers
  void eval_funcs();  // populates the ECFG_chain (owned by the ECFG_matrix_set) with numeric values, by evaluating functions

  // accessors
  ECFG_chain& chain();
  EM_matrix_base& rate_matrix();
};

// Terminatrix I/O adapter
struct Terminatrix_builder : ECFG_builder
{
  // load
  static void init_terminatrix (Terminatrix& terminatrix, SExpr& terminatrix_sexpr);
  // save
  static void terminatrix2stream (ostream& out, Terminatrix& terminatrix);

  // input helpers
  static void init_terminatrix_params (Terminatrix& terminatrix, SExpr& terminatrix_params_sexpr);
  static void init_terminatrix_model (Terminatrix& terminatrix, SExpr& terminatrix_model_sexpr);
  static void init_terminatrix_member_sexpr (SExpr& member_sexpr, SExpr& parent_sexpr, const char* tag);

  // output helpers
  static void terminatrix_params2stream (ostream& out, Terminatrix& terminatrix);
  static void terminatrix_model2stream (ostream& out, Terminatrix& terminatrix);
  static void terminatrix_member_sexpr2stream (ostream& out, SExpr& member_sexpr);
};

// Terminatrix family visitor
struct Terminatrix_family_visitor
{
  // global info
  Terminatrix& terminatrix;

  // info on the current family
  int current_family_index;
  SCM current_name_scm, current_newick_scm;
  sstring current_name;
  PHYLIP_tree *current_tree;

  // constructor
  Terminatrix_family_visitor (Terminatrix& term)
    : terminatrix(term),
      current_family_index(-1),
      current_name_scm(SCM_BOOL_F),
      current_newick_scm(SCM_BOOL_F),
      current_tree(NULL)
  { }

  // virtual methods
  // destructor
  virtual ~Terminatrix_family_visitor() { }

  // map-reduce
  // the result of the map operation is not explicitly represented; it is assumed to be described by the current state of the subclass
  virtual scm_t_bits zero() { return SCM_UNPACK (SCM_BOOL_F); }  // guaranteed to be called before any families are visited. "What result will you get if there are no trees?"
  virtual void map_current() { }  // guaranteed to be called exactly once for every family, right before reduce(). "Prepare your internal state to reflect the current family"
  virtual scm_t_bits reduce (scm_t_bits previous) { return previous; }  // guaranteed to be called after map_current(). "Combine your internal state with the results so far"
  virtual SCM finalize (scm_t_bits result) { return SCM_PACK (result); }  // guaranteed to be called after all families visited. "Convert the results to a SCM object"

  // map-reduce method
  SCM map_reduce();
};

struct Terminatrix_concatenator
{
  scm_t_bits reduce (scm_t_bits previous) { SCM cons = scm_cons (current_mapped_scm(), SCM_PACK(previous)); return SCM_UNPACK(cons); }
  scm_t_bits zero() { return SCM_UNPACK (scm_list_n (SCM_UNDEFINED)); }  // creates an empty list
  virtual SCM current_mapped_scm() = 0;  // guaranteed to be called after map_current()
};

struct Terminatrix_EM_visitor : Terminatrix_family_visitor
{
  // global info
  Update_statistics stats;  // uses these in preference to Terminatrix's; copy across if needed
  // info on the current family
  Column_matrix current_colmat;
  // methods
  void map_current() { initialize_current_colmat(); }
  void initialize_current_colmat();
};

struct Terminatrix_log_evidence : Terminatrix_EM_visitor, Terminatrix_concatenator
{
  SCM current_mapped_scm()
  {
    Terminatrix_EM_visitor::current_colmat.fill_up (Terminatrix_family_visitor::terminatrix.rate_matrix(), *(Terminatrix_family_visitor::current_tree));
    return scm_from_double (Terminatrix_EM_visitor::current_colmat.total_log_likelihood());
  }
};

#endif /* ONTO_SEXPR_INCLUDED */
