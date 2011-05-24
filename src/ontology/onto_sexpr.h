#ifndef ONTO_SEXPR_INCLUDED
#define ONTO_SEXPR_INCLUDED

#if defined(GUILE_INCLUDED) && GUILE_INCLUDED
#include <libguile.h>
#else
#error Terminatrix requires guile! Please install guile from http://www.gnu.org/software/guile/ and re-run the configure script.
#endif

#include "ecfg/ecfgsexpr.h"
#include "ecfg/guile-ecfg.h"
#include "seq/psexpr.h"
#include "util/svisitor.h"
#include "guile/assmap.h"
#include "ontology/onto_keywords.h"


// Terminatrix class
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
struct Terminatrix
{
  // Scheme functions that are used to create & initialize the knowledge-base
  SCM init_scm,  // called when the Terminatrix is initialized. Return result is discarded
    model_scm,     // use to generate the model. Must return a quoted S-expression of the form (model (alphabet ...)+ (chain ...))
    tree_db_scm,   // use to generate the phylogenetic tree database. Must return a Scheme association list mapping names (i.e. strings) to trees (i.e. newick smobs)
    knowledge_scm; // must return a function that will be called like (knowledge familyName geneName termTuple) and must return #t (gene-term association allowed) or #f (gene-term association disallowed)

  SCM knowledge_func_scm;  // return value of knowledge_scm

  // symbol tables (only really used by Terminatrix_builder class; should probably move them there)
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
  bool got_counts;  // true if pcounts have been populated

  // constructor
  // Assumes Guile has been initialized, and the Newick smob added.
  Terminatrix();

  // helpers
  void eval_funcs();  // populates the ECFG_chain (owned by the ECFG_matrix_set) with numeric values, by evaluating PFunc's (parsed algebraic functions)

  // accessors
  ECFG_chain& chain();
  EM_matrix_base& rate_matrix();
  const Alphabet& default_alphabet();
};

// Terminatrix family visitor
// This class and its subclasses form a sort of bastard map-reduce framework.
// The input is a list of "families" (named phylogenetic trees).
// For efficiency, the base class (Terminatrix_family_visitor) does reduction in place, i.e. a left-fold.
// Default behavior of this class is to return #f, discarding values in the reduction.
// To get closer to a true map-reduce, override current_mapped_scm() and reduce_scm() methods in derived class (Terminatrix_concatenator).
// The idealized, "pure" map-reduce functions would be as follows
//     map : (string current_name, PHYLIP_tree current_tree) -> SCM
//  reduce : (SCM, SCM) -> SCM
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

  // map-reduce virtual methods
  // the order of calls is as follows:
  //  - zero();
  //  - for each family:
  //     - init_current();
  //     - reduce().
  //  - finalize().
  // the result of the map operation is not explicitly represented by this class; for that, use the derived class Terminatrix_concatenator.
  // in this class, the intermediate results of the fold are represented as the integral type, scm_t_bits, rather than the guile internal type SCM,
  // to allow derived classes to pass other types instead of a Scheme object.
  // practically, this just means a few calls to SCM_PACK and SCM_UNPACK.
  // the derived class Terminatrix_concatenator overrides this behavior too, passing everything as SCM.
  virtual scm_t_bits zero() { return SCM_UNPACK (SCM_BOOL_F); }  // guaranteed to be called before any families are visited. "What result will you get if there are zero families?"
  virtual void init_current() { }  // guaranteed to be called exactly once for every family, right before reduce(). "Prepare your internal state to reflect the current family"
  virtual scm_t_bits reduce (scm_t_bits previous) { return previous; }  // guaranteed to be called after init_current(). "Combine your internal state with the results so far"
  virtual SCM finalize (scm_t_bits result) { return SCM_PACK (result); }  // guaranteed to be called after all families visited. "Convert the results to a SCM object"

  // map-reduce methods
  SCM map_reduce_scm();
  SExpr map_reduce_sexpr() { return scm_to_sexpr (map_reduce_scm()); }

  // Terminatrix convenience wrappers
  EM_matrix_base& rate_matrix() { return terminatrix.rate_matrix(); }
  ECFG_chain& chain() { return terminatrix.chain(); }
  int number_of_states() { return terminatrix.rate_matrix().number_of_states(); }
};

// Terminatrix_concatenator
// A virtual class returning a cons of the map results.
// This class also commits the map-reduction to passing SCM objects.
// Overrides the zero() and reduce(previous) methods.
// Introduces new virtual methods current_mapped_scm(), reduce_scm(), zero_scm(), finalize_scm
struct Terminatrix_concatenator : virtual Terminatrix_family_visitor
{
  Terminatrix_concatenator (Terminatrix& term) : Terminatrix_family_visitor(term) { }
  scm_t_bits reduce (scm_t_bits previous) { SCM reduction = reduce_scm (current_mapped_scm(), SCM_PACK(previous)); return SCM_UNPACK(reduction); }  // delegate to reduce_scm(). intended final
  scm_t_bits zero() { return SCM_UNPACK (zero_scm()); }  // delegate to zero_scm(). intended final
  virtual SCM finalize (scm_t_bits result) { return finalize_scm (Terminatrix_family_visitor::finalize (result)); }
  virtual SCM current_mapped_scm() = 0;  // guaranteed to be called after init_current()
  virtual SCM reduce_scm (SCM a, SCM b) { return scm_cons (a, b); }  // does a cons
  virtual SCM zero_scm() { return scm_list_n (SCM_UNDEFINED); }  // creates an empty list
  virtual SCM finalize_scm (SCM result) { return result; }  // does nothing
};

// Terminatrix_keyed_concatenator
// A virtual class similar to Terminatrix_concatenator,
// but returning a list of the form ((family1-id family1-results) (family2-id family2-results) ...)
// rather than a list of the form (family1-results family2-results ...)
struct Terminatrix_keyed_concatenator : virtual Terminatrix_concatenator
{
  Terminatrix_keyed_concatenator (Terminatrix& term) : Terminatrix_family_visitor(term), Terminatrix_concatenator(term) { }
  scm_t_bits reduce (scm_t_bits previous)  // delegate to reduce_scm(). intended final
  {
    SCM previous_scm = SCM_PACK(previous);
    SCM current_scm = scm_list_2 (current_name_scm, current_mapped_scm());
    SCM reduced_scm = reduce_scm (current_scm, previous_scm);
    return SCM_UNPACK(reduced_scm);
  }
};

// Terminatrix_EM_visitor
// A class that initializes a Column_matrix (DP matrix for a phylogenetic continuous-time Markov chain)
// by querying the Terminatrix "knowledge function".
struct Terminatrix_EM_visitor : virtual Terminatrix_family_visitor
{
  // global info
  Update_statistics stats;  // used by this class in preference to Terminatrix's; copy across if needed
  PCounts var_counts;  // initialized to pseudocounts; used by this class in preference to Terminatrix's; copy across if needed
  // info on the current family
  Column_matrix current_colmat;
  // methods
  Terminatrix_EM_visitor (Terminatrix& term) : Terminatrix_family_visitor(term) { clear(); }
  void init_current() { initialize_current_colmat(); }  // intended final
  void initialize_current_colmat();  // called by init_current()
  SCM node_name_scm (int node) {
      const sstring& node_name = current_tree->node_name[node];
      return scm_from_locale_string (node_name.c_str());
  }
  bool knowledge_func (int node, int state) {
    SCM knowledge_result_scm = scm_call_3 (terminatrix.knowledge_func_scm, current_name_scm, node_name_scm (node), state_tuple_scm (state));
    if (!scm_boolean_p (knowledge_result_scm))
      THROWEXPR ("(" TERMINATRIX_KNOWLEDGE_SCM " " << current_name << " " << current_tree->node_name[state] << " (" << state_tuple(state) << ")) should return #t or #f, but it didn't");
    return knowledge_result_scm == SCM_BOOL_T;
  }
  vector<sstring> state_tuple (int state) {
    return terminatrix.chain().get_symbol_tokens (state, terminatrix.alph_dict, terminatrix.default_alphabet());
  }
  SCM state_tuple_scm (int state) {
    const vector<sstring> tuple = state_tuple(state);
    return vector_to_scm (tuple);
  }
  // Column_matrix wrappers
  void fill_up() {
    current_colmat.fill_up (Terminatrix_family_visitor::rate_matrix(),
			    *(Terminatrix_family_visitor::current_tree));
  }
  void fill_down() {
    current_colmat.fill_down (Terminatrix_family_visitor::rate_matrix(),
			      *(Terminatrix_family_visitor::current_tree),
			      stats);
  }
  Loge total_log_likelihood() {
    return current_colmat.total_log_likelihood();
  }
  Prob node_post_prob (int node, int state) {
    return current_colmat.node_post_prob (node, state,
					  *(Terminatrix_family_visitor::current_tree),
					  Terminatrix_family_visitor::rate_matrix());
  }
  void inc_var_counts() {
    Terminatrix_family_visitor::chain().inc_var_counts (stats, var_counts, terminatrix.pscores);
  }
  void clear() { stats.clear(); var_counts = terminatrix.pcounts; }
};

// Terminatrix_log_evidence
// The current_mapped_scm function runs the Column_matrix implementation of Felsenstein's algorithm.
struct Terminatrix_log_evidence : Terminatrix_keyed_concatenator, Terminatrix_EM_visitor
{
  Terminatrix_log_evidence (Terminatrix& term)
    : Terminatrix_family_visitor(term),
      Terminatrix_concatenator(term),
      Terminatrix_keyed_concatenator(term),
      Terminatrix_EM_visitor(term)
  { }
  SCM current_mapped_scm()
  {
    Terminatrix_EM_visitor::fill_up();
    return scm_from_double (Terminatrix_EM_visitor::total_log_likelihood());
  }
  SCM finalize_scm (SCM result)
  {
    return scm_list_2 (string_to_scm (TERMINATRIX_TERMINATRIX),
		       scm_cons (string_to_scm (TERMINATRIX_LOG_EVIDENCE),
				 result));
  }
};

// Terminatrix_prediction
// The current_mapped_scm function runs the Column_matrix implementations of the sum-product message-passing algorithms
// (Felsenstein's pruning algorithm, then Elston & Stewart's peeling algorithm, or however you want to think about it),
// and returns a table of posterior probabilities for the possible states at the various nodes.
struct Terminatrix_prediction : Terminatrix_keyed_concatenator, Terminatrix_EM_visitor
{
  Terminatrix_prediction (Terminatrix& term)
    : Terminatrix_family_visitor(term),
      Terminatrix_concatenator(term),
      Terminatrix_keyed_concatenator(term),
      Terminatrix_EM_visitor(term)
  { }
  SCM current_mapped_scm()
  {
    Terminatrix_EM_visitor::fill_up();
    Terminatrix_EM_visitor::fill_down();
    SCM by_node = scm_list_n (SCM_UNDEFINED);  // empty list
    for (Phylogeny::Node node = 0; node < current_tree->nodes(); ++node)
      if (current_tree->node_name[node].size())
	{
	  SCM by_state = scm_list_n (SCM_UNDEFINED);  // empty list
	  for (int state = 0; state < Terminatrix_family_visitor::number_of_states(); ++state)
	    if (knowledge_func (node, state))
	      by_state = scm_cons (scm_list_2 (state_tuple_scm (state),
					       scm_from_double (Terminatrix_EM_visitor::node_post_prob (node, state))),
				   by_state);
	  by_node = scm_cons (scm_list_2 (string_to_scm (current_tree->node_name[node]),
					  by_state),
			      by_node);
	}
    return by_node;
  }
  SCM finalize_scm (SCM result)
  {
    return scm_list_2 (string_to_scm (TERMINATRIX_TERMINATRIX),
		       scm_cons (string_to_scm (TERMINATRIX_POSTERIOR),
				 result));
  }
};

// Terminatrix I/O adapter
struct Terminatrix_builder : ECFG_builder
{
  // load
  static void init_terminatrix (Terminatrix& terminatrix, SExpr& terminatrix_sexpr);
  static void init_terminatrix (Terminatrix& terminatrix, SCM terminatrix_scm);
  static void init_terminatrix (Terminatrix& terminatrix, const Ass_map& ass_map);
  // save
  static void terminatrix2stream (ostream& out, Terminatrix& terminatrix);

  // input helpers
  static void init_terminatrix_params (Terminatrix& terminatrix, SExpr& terminatrix_params_sexpr);
  static void init_terminatrix_model (Terminatrix& terminatrix, SExpr& terminatrix_model_sexpr);
  static void init_terminatrix_member_scm (SCM& member_scm, const Ass_map& parent_ass_map, const char* tag);

  // output helpers
  static void terminatrix_params2stream (ostream& out, Terminatrix& terminatrix);
  static void terminatrix_model2stream (ostream& out, Terminatrix& terminatrix);
  static void terminatrix_member_scm2stream (ostream& out, SCM& member_scm, const char* tag);
};

// guile methods
SCM terminatrix_evidence (SCM terminatrix_scm);
SCM terminatrix_prediction (SCM terminatrix_scm);
void init_terminatrix_primitives (void);

#endif /* ONTO_SEXPR_INCLUDED */
