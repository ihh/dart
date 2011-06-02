#ifndef ONTO_SEXPR_INCLUDED
#define ONTO_SEXPR_INCLUDED

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include <libguile.h>
#else
#error Termx requires guile! Please install guile from http://www.gnu.org/software/guile/ and re-run the configure script.
#endif

#include "ecfg/ecfgsexpr.h"
#include "ecfg/guile-ecfg.h"
#include "seq/psexpr.h"
#include "util/svisitor.h"
#include "guile/assmap.h"
#include "ontology/onto_keywords.h"


// Termx class
// this C++/Scheme hybrid, and its helper classes, propagate knowledge from an ontological knowledge base
// over a particular kind of graphical model: a continuous-time Markov chain over a phylogenetic tree.
// the class is documented below with reference to
//  "families" (the trees)
//  "genes" (the things at the nodes of the tree)
//  "terms" (the members of the ontology, that can be associated with the genes)
// however these are somewhat arbitrary labels for the values that the random variables can take
// (e.g. it works perfectly well as a substitution model for nucleotides, codons, amino acids, other discrete characters)
// With a view to future parallelization, the core map-reduce operation common to all the usual
// variations on the sum-product algorithm (evidence, posteriors, EM, etc) has been abstracted out.
struct Termx
{
  // Scheme functions that are used to create & initialize the knowledge-base
  SCM init_scm,  // called when the Termx is initialized. Return result is discarded
    model_scm,     // use to generate the model. Must return a quoted S-expression of the form (model (alphabet ...)+ (chain ...))
    tree_db_scm,   // use to generate the phylogenetic tree database. Must return a Scheme association list mapping names (i.e. strings) to trees (i.e. newick smobs)
    knowledge_scm; // must return a function that will be called like (knowledge familyName geneName termTuple) and must return #t (gene-term association allowed) or #f (gene-term association disallowed)

  SCM knowledge_func_scm;  // return value of knowledge_scm

  // symbol tables (only really used by Termx_builder class; should probably move them there)
  PFunc_builder::SymPVar sym2pvar;  // numerical parameters
  PFunc_builder::SymIndex term2chain;  // term variables in the Markov chain

  // numerical parameters of the model
  PScores pscores;
  PCounts pcounts, var_counts;
  set<int> mutable_pgroups;

  // alphabet namespace
  // in our nomenclature, corrupted by ripping this code out of the guts of a sequence analysis package,
  // an "Alphabet" class models a named, unstructured set of ontology terms (represented as string "tokens").
  const Alphabet dummy_alph;  // the ECFG_matrix_set needs a default Alphabet. Keep this around to pacify it
  list<Alphabet> alph_list;  // TODO: alph_list is a bit wasteful; we already store every Alphabet in alph_dict. Could just store a list of Alphabet names here instead...
  Alphabet_dictionary alph_dict;  // mapping from alphabet name to the corresponding Alphabet

  // continuous-time Markov chain
  // for icky legacy reasons, this is wrapped in a class from the ECFG (Evolutionary Context-Free Grammar) library.
  ECFG_matrix_set matrix_set;

  // summary statistics of EM computations
  Update_statistics stats;

  // status variables, used to throttle output
  bool got_counts;  // true if var_counts have been populated

  // constructors
  // These assume Guile has been initialized, and the Newick smob added.
  Termx();
  Termx (SCM scm);
  Termx (const SExpr& sexpr);

  // destructor
  ~Termx();
  void unprotect (SCM scm);  // destructor helper

  // helpers
  void eval_funcs();  // populates the ECFG_chain with numeric values, by evaluating PFunc's (parsed algebraic functions); also resets stats

  // accessors
  ECFG_chain& chain();
  EM_matrix_base& rate_matrix();
  const Alphabet& default_alphabet();
};

// Termx family visitor
// This class and its subclasses form a sort of bastard map-reduce framework.
// The input is a list of "families" (named phylogenetic trees).
// For efficiency, the base class (Termx_family_visitor) does reduction in place, i.e. a left-fold.
// Default behavior of this class is to return #f, discarding values in the reduction.
// To get closer to a true map-reduce, override current_mapped_scm() and reduce_scm() methods in derived class (Termx_concatenator).
// The idealized, "pure" map-reduce functions would be as follows
//     map : (string current_name, PHYLIP_tree current_tree) -> SCM
//  reduce : (SCM, SCM) -> SCM
struct Termx_family_visitor
{
  // global info
  Termx& termx;

  // info on the current family
  int current_family_index;
  SCM current_name_scm, current_newick_scm;
  sstring current_name;
  PHYLIP_tree *current_tree;

  // constructor
  Termx_family_visitor (Termx& term)
    : termx(term),
      current_family_index(-1),
      current_name_scm(SCM_BOOL_F),
      current_newick_scm(SCM_BOOL_F),
      current_tree(NULL)
  { }

  // virtual methods
  // destructor
  virtual ~Termx_family_visitor() { }

  // map-reduce virtual methods
  // the order of calls is as follows:
  //  - zero();
  //  - for each family:
  //     - init_current();
  //     - reduce().
  //  - finalize().
  // the result of the map operation is not explicitly represented by this class; for that, use the derived class Termx_concatenator.
  // in this class, the intermediate results of the fold are represented as the integral type, scm_t_bits, rather than the guile internal type SCM,
  // to allow derived classes to pass other types instead of a Scheme object.
  // practically, this just means a few calls to SCM_PACK and SCM_UNPACK.
  // the derived class Termx_concatenator overrides this behavior too, passing everything as SCM.
  virtual scm_t_bits zero() { return SCM_UNPACK (SCM_BOOL_F); }  // guaranteed to be called before any families are visited. "What result will you get if there are zero families?"
  virtual void init_current() { }  // guaranteed to be called exactly once for every family, right before reduce(). "Prepare your internal state to reflect the current family"
  virtual scm_t_bits reduce (scm_t_bits previous) { return previous; }  // guaranteed to be called after init_current(). "Combine your internal state with the results so far"
  virtual SCM finalize (scm_t_bits result) { return SCM_PACK (result); }  // guaranteed to be called after all families visited. "Convert the results to a SCM object"

  // map-reduce methods
  SCM map_reduce_scm();
  SExpr map_reduce_sexpr() { return scm_to_sexpr (map_reduce_scm()); }

  // Termx convenience wrappers
  EM_matrix_base& rate_matrix() { return termx.rate_matrix(); }
  ECFG_chain& chain() { return termx.chain(); }
  int number_of_states() { return termx.rate_matrix().number_of_states(); }
};

// Termx_concatenator
// A virtual class returning a cons of the map results.
// This class also commits the map-reduction to passing SCM objects.
// Overrides the zero() and reduce(previous) methods.
// Introduces new virtual methods current_mapped_scm(), reduce_scm(), zero_scm(), finalize_scm
struct Termx_concatenator : virtual Termx_family_visitor
{
  Termx_concatenator (Termx& term) : Termx_family_visitor(term) { }
  scm_t_bits reduce (scm_t_bits previous) {  // delegate to reduce_scm(). intended final
    SCM previous_scm = SCM_PACK(previous);
    SCM reduced_scm = reduce_scm (current_mapped_scm(), previous_scm);
    return SCM_UNPACK(reduced_scm);
  }
  scm_t_bits zero() {  // delegate to zero_scm(). intended final
    SCM z_scm = zero_scm();
    return SCM_UNPACK (z_scm);
  }
  virtual SCM finalize (scm_t_bits result) {
    SCM final_scm = Termx_family_visitor::finalize (result);
    return finalize_scm (final_scm);
  }
  virtual SCM current_mapped_scm() = 0;  // guaranteed to be called after init_current()
  virtual SCM reduce_scm (SCM a, SCM b) { return scm_cons (a, b); }  // does a cons
  virtual SCM zero_scm() { return scm_list_n (SCM_UNDEFINED); }  // creates an empty list
  virtual SCM finalize_scm (SCM result) { return result; }  // does nothing
};

// Termx_keyed_concatenator
// A virtual class similar to Termx_concatenator,
// but returning a list of the form ((family1-id family1-results) (family2-id family2-results) ...)
// rather than a list of the form (family1-results family2-results ...)
struct Termx_keyed_concatenator : virtual Termx_concatenator
{
  Termx_keyed_concatenator (Termx& term) : Termx_family_visitor(term), Termx_concatenator(term) { }
  scm_t_bits reduce (scm_t_bits previous)  // delegate to reduce_scm(). intended final
  {
    SCM previous_scm = SCM_PACK(previous);
    SCM current_value_scm = current_mapped_scm();
    SCM current_scm = scm_cons (current_name_scm, current_value_scm);
    SCM reduced_scm = reduce_scm (current_scm, previous_scm);
    return SCM_UNPACK(reduced_scm);
  }
};

// Termx_EM_visitor
// A class that initializes a Column_matrix (DP matrix for a phylogenetic continuous-time Markov chain)
// by querying the Termx "knowledge function".
struct Termx_EM_visitor : virtual Termx_family_visitor
{
  // global info
  Update_statistics stats;  // used by this class in preference to Termx's; copy across if needed
  PCounts var_counts;  // initialized to pseudocounts; used by this class in preference to Termx's; copy across if needed
  // info on the current family
  Column_matrix current_colmat;
  // methods
  Termx_EM_visitor (Termx& term)
    : Termx_family_visitor(term),
      stats(term.rate_matrix().number_of_states()),
      var_counts(term.pscores)
  {
    // evaluate all initial probability & mutation rate functions
    term.eval_funcs();
  }
  void init_current() { initialize_current_colmat(); }  // intended final
  void initialize_current_colmat();  // called by init_current()
  SCM node_name_scm (int node) {
      const sstring& node_name = current_tree->node_name[node];
      return scm_from_locale_string (node_name.c_str());
  }
  bool knowledge_func (int node, int state) {
    SCM knowledge_result_scm = scm_call_3 (termx.knowledge_func_scm, current_name_scm, node_name_scm (node), state_tuple_scm (state));
    if (!scm_boolean_p (knowledge_result_scm))
      THROWEXPR ("(" TERMX_KNOWLEDGE_SCM " " << current_name << " " << current_tree->node_name[state] << " (" << state_tuple(state) << ")) should return #t or #f, but it didn't");
    return knowledge_result_scm == SCM_BOOL_T;
  }
  vector<sstring> state_tuple (int state) {
    return termx.chain().get_symbol_tokens (state, termx.alph_dict, termx.default_alphabet());
  }
  SCM state_tuple_scm (int state) {
    const vector<sstring> tuple = state_tuple(state);
    return vector_to_scm (tuple);
  }
  // Column_matrix wrappers
  void fill_up() {
    current_colmat.fill_up (Termx_family_visitor::rate_matrix(),
			    *(Termx_family_visitor::current_tree));
  }
  void fill_down() {
    current_colmat.fill_down (Termx_family_visitor::rate_matrix(),
			      *(Termx_family_visitor::current_tree),
			      stats);
  }
  Loge current_log_evidence() {
    return current_colmat.total_log_likelihood();
  }
  Prob node_post_prob (int node, int state) {
    return current_colmat.node_post_prob (node, state,
					  *(Termx_family_visitor::current_tree),
					  Termx_family_visitor::rate_matrix());
  }
  void accumulate_chain_counts() {
    stats.transform (Termx_family_visitor::rate_matrix(), false);
    Termx_family_visitor::chain().inc_var_counts (stats, var_counts, termx.pscores);
  }
  void accumulate_pseudocounts() {
    var_counts += termx.pcounts;
  }
  static SCM pscores_to_scm (const PScores& pscores, const char* tag = TERMX_PARAMS) {
    sstring pscores_str;
    pscores_str << '(' << tag << ' ';
    PFunc_builder::pscores2stream (pscores_str, pscores);
    pscores_str << ')';
    SExpr pscores_sexpr (pscores_str.begin(), pscores_str.end());
    CTAG(1,TERMX) << "In pscores_to_scm: " << pscores_sexpr.to_string() << '\n';
    return sexpr_to_scm (&pscores_sexpr);
  }
  static PScores pscores_from_scm (SCM scm) {
    SExpr sexpr = scm_to_sexpr (scm);
    PScores pscores;
    PFunc_builder::SymPVar sym2pvar;
    PFunc_builder::init_pgroups_and_rates (pscores, sym2pvar, sexpr(TERMX_PARAMS));
    return pscores;
  }
  static SCM pcounts_to_scm (const PCounts& pcounts, const char* tag = TERMX_PARAM_COUNTS) {
    sstring pcounts_str;
    PFunc_builder::pcounts2stream (pcounts_str, pcounts, tag);
    SExpr pcounts_sexpr (pcounts_str.begin(), pcounts_str.end());
    CTAG(1,TERMX) << "In pcounts_to_scm: " << pcounts_sexpr.to_string() << '\n';
    return sexpr_to_scm (&pcounts_sexpr);
  }
};

// Termx_log_evidence
// The current_mapped_scm function runs the Column_matrix implementation of Felsenstein's algorithm.
struct Termx_log_evidence : Termx_keyed_concatenator, Termx_EM_visitor
{
  Termx_log_evidence (Termx& term)
    : Termx_family_visitor(term),
      Termx_concatenator(term),
      Termx_keyed_concatenator(term),
      Termx_EM_visitor(term)
  { }
  SCM current_mapped_scm()
  {
    Termx_EM_visitor::fill_up();
    return scm_list_1 (scm_from_double (Termx_EM_visitor::current_log_evidence()));
  }
  SCM finalize_scm (SCM result)
  {
    return scm_list_2 (string_to_scm (TERMX_TERMX),
		       scm_cons (string_to_scm (TERMX_LOG_EVIDENCE),
				 result));
  }
};

// Termx_prediction
// The current_mapped_scm function runs the Column_matrix implementations of the sum-product message-passing algorithms
// (Felsenstein's pruning algorithm, then Elston & Stewart's peeling algorithm, or however you want to think about it),
// and returns a table of posterior probabilities for the possible states at the various nodes.
struct Termx_prediction : Termx_keyed_concatenator, Termx_EM_visitor
{
  Termx_prediction (Termx& term)
    : Termx_family_visitor(term),
      Termx_concatenator(term),
      Termx_keyed_concatenator(term),
      Termx_EM_visitor(term)
  { }
  SCM current_mapped_scm()
  {
    Termx_EM_visitor::fill_up();
    Termx_EM_visitor::fill_down();
    SCM by_node = scm_list_n (SCM_UNDEFINED);  // empty list
    for (Phylogeny::Node node = 0; node < current_tree->nodes(); ++node)
      if (current_tree->node_name[node].size())
	{
	  SCM by_state = scm_list_n (SCM_UNDEFINED);  // empty list
	  for (int state = 0; state < Termx_family_visitor::number_of_states(); ++state)
	    if (knowledge_func (node, state))
	      by_state = scm_cons (scm_list_2 (state_tuple_scm (state),
					       scm_from_double (Termx_EM_visitor::node_post_prob (node, state))),
				   by_state);
	  by_node = scm_cons (scm_list_2 (string_to_scm (current_tree->node_name[node]),
					  by_state),
			      by_node);
	}
    return by_node;
  }
  SCM finalize_scm (SCM result)
  {
    return scm_list_2 (string_to_scm (TERMX_TERMX),
		       scm_cons (string_to_scm (TERMX_POSTERIOR),
				 result));
  }
};

// Termx_learning_step
struct Termx_learning_step : Termx_EM_visitor
{
  // data
  PScores next_pscores;
  Loge total_log_ev;
  // constructor
  Termx_learning_step (Termx& term)
    : Termx_family_visitor(term),
      Termx_EM_visitor(term),
      next_pscores (term.pscores),
      total_log_ev (0.)
  { }
  // overrides
  scm_t_bits reduce (scm_t_bits previous) {  // intended final
    Termx_EM_visitor::fill_up();
    Termx_EM_visitor::fill_down();
    NatsPMulAcc (total_log_ev, current_log_evidence());
    return Termx_family_visitor::reduce (previous);
  }
  SCM finalize (scm_t_bits result) {  // intended final
    Termx_EM_visitor::accumulate_chain_counts();
    Termx_EM_visitor::accumulate_pseudocounts();
    var_counts.optimise (next_pscores, Termx_family_visitor::termx.mutable_pgroups);
    /*
    return scm_list_2 (string_to_scm (TERMX_TRAINING_STEP),
		       scm_list_2 (string_to_scm (TERMX_LOG_EVIDENCE),
				   scm_from_double (total_log_ev)));
    */
    return scm_list_5 (string_to_scm (TERMX_TRAINING_STEP),
		       scm_list_2 (string_to_scm (TERMX_LOG_EVIDENCE),
				   scm_from_double (total_log_ev)),
		       pscores_to_scm (Termx_family_visitor::termx.pscores, TERMX_PARAMS),
		       pcounts_to_scm (Termx_EM_visitor::var_counts, TERMX_PARAM_COUNTS),
		       pscores_to_scm (next_pscores, TERMX_NEXT_PARAMS));
  }
};

// Termx_trainer
struct Termx_trainer
{
  // data
  Loge max_total_log_ev;
  PScores argmax_total_log_ev;
  SCM log_scm;
  // constructor (everything happens in here)
  Termx_trainer (Termx& term, int max_learning_steps)
  {
    log_scm = scm_list_n (SCM_UNDEFINED);  // empty list
    for (int n_step = 0; n_step < max_learning_steps; ++n_step)
      {
	Termx_learning_step step (term);
	SCM step_log_scm = step.map_reduce_scm();  // must call map_reduce_scm() to ensure step.next_pscores is set correctly
	if (n_step == 0 || step.total_log_ev > max_total_log_ev)
	  {
	    max_total_log_ev = step.total_log_ev;
	    argmax_total_log_ev = step.next_pscores;
	  }
	else
	  break;
	CTAG(6,TERMX) << "EM iteration #" << (n_step + 1) << ": log-evidence " << max_total_log_ev << '\n';
	log_scm = scm_cons (step_log_scm, log_scm);
	term.pscores = step.next_pscores;
      }
    term.pscores = argmax_total_log_ev;
    scm_gc_protect_object (log_scm);
  }
  ~Termx_trainer() {
    scm_gc_unprotect_object (log_scm);
  }
  // accessors
  SCM final_scm() {
    SCM pscores_scm = Termx_EM_visitor::pscores_to_scm (argmax_total_log_ev);
    // TODO: optionally include training log
    return scm_cons (string_to_scm (TERMX_TERMX),
		     pscores_scm);
  }
  SExpr final_sexpr() {
    SCM scm = final_scm();
    return scm_to_sexpr (scm);
  }
};

// Termx I/O adapter
struct Termx_builder : ECFG_builder
{
  // load
  static void init_termx (Termx& termx, SExpr& termx_sexpr);
  static void init_termx (Termx& termx, SCM termx_scm);
  static void init_termx (Termx& termx, const Ass_map& ass_map);
  // save
  static void termx2stream (ostream& out, Termx& termx);

  // input helpers
  static void init_termx_params (Termx& termx, SExpr& termx_params_sexpr);
  static void init_termx_model (Termx& termx, SExpr& termx_model_sexpr);
  static void init_termx_member_scm (SCM& member_scm, const Ass_map& parent_ass_map, const char* tag);

  // output helpers
  static void termx_params2stream (ostream& out, Termx& termx);
  static void termx_model2stream (ostream& out, Termx& termx);
  static void termx_member_scm2stream (ostream& out, SCM& member_scm, const char* tag);
};

// guile methods
SCM termx_evidence (SCM termx_scm);
SCM termx_prediction (SCM termx_scm);
SCM termx_learn (SCM max_steps_scm, SCM termx_scm);

void init_termx_primitives (void);

#endif /* ONTO_SEXPR_INCLUDED */
