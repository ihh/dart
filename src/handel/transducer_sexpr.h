#ifndef TRANSDUCER_SEXPR_INCLUDED
#define TRANSDUCER_SEXPR_INCLUDED

#include "handel/transducer.h"
#include "seq/psexpr.h"

// name of automatically created dummy singleton transducers
#define Singleton_transducer_name   "*SINGLETON"

// This class started as an S-Expression input adapter for Transducer<PFunc>,
// as a clean test of the Handel MCMC code, which was proving hard to debug.
//
struct Transducer_SExpr_file : TSpaceEnum, Transducer_state_type_enum, Grammar_state_enum, PFunc_builder
{
  // typedefs
  typedef Pair_transducer_funcs Pair_trans_funcs;  // pairwise transducer
  typedef map<sstring,Alphabet_group> SymPGroup;  // name lookup for multidimensional token-indexed array labels (Alphabet_group's)
  typedef map<sstring,Pair_trans_funcs> SymTrans;  // name lookup for pairwise transducers
  typedef map<int,Score_profile> NodeProfileMap;  // map from node indices to observed sequence bit-profiles (Score_profile's)
  typedef map<int,vector<int> > NodePathMap;  // map from node indices to state paths responsible for the emission of that node
  typedef map<sstring,int> StateIndexMap;  // name lookup for states

  // data
  // alphabet
  vector<sstring> alphabet;
  map<sstring,int> sym2tok;  // map of alphabet symbols to integer tokens
  int alphabet_size;

  // parameters
  PScores pscores;
  set<PVar> defined_pvars;  // list of PVars whose actual values are specified
  SymPVar sym2pvar;  // map of transition labels to PVar's
  SymPGroup sym2pgroup;  // map of emit labels to Alphabet_group's

  // transducers
  SymTrans trans_dict;  // dictionary of transducers

  // tree
  ETree etree;
  vector<sstring> branch_name;
  map<int,sstring> tape_name;  // node names default to "$0" (for node -1), "$1" (for node 0), "$2" (for node 1), etc.

  // branch->transducer mapping (complete; indexed by branch)
  vector<sstring> branch_trans_name;
  vector<Pair_trans_funcs> branch_trans;

  // tape->sequence mapping (partial)
  NodeProfileMap prof;

  // branch->statepath mappings (partial)
  // path: state paths that will be preserved
  // old_path: state paths that will be resampled
  NodePathMap path, old_path;

  // rank->branch mapping (partial): specifies order of partial compositions for Redelings-Suchard MCMC kernel
  typedef vector<vector<int> > RedSuchSchedule;
  RedSuchSchedule proposal_branches;

  // branch->banding coefficient mapping
  map<int,double> band_coeff;

  // subtrees induced by branch->statepath mapping
  // these variables are initialized by setup_clique()
  map<int,int> uncons_clique, cons_clique;  // map from nodes to clique indices. each node can be a member of up to two cliques: constrained & unconstrained
  set<int> hinge_nodes;  // "hinges" are nodes that are in members of both types of clique, i.e. on the perimeter (so they need to be pruned down to profiles)
  vector<bool> branch_cons;  // OTOH, each *branch* must either be constrained or unconstrained. branch_cons[N] is true if branch to node N is constrained
  vector<set<int> > clique;  // clique[C] = set of nodes in clique #C
  int free_clique;  // only one clique can correspond to an unconstrained subtree. this is the index of that clique. if negative, then entire tree is constrained

  // composite paths for constrained subtrees
  // used only by show_tree output method
  typedef map<int,int> Branch_transition_map;
  struct Composite_path_step
  {
    vector<int> emit_nodes;
    Branch_transition_map branch_trans;
  };
  struct Composite_path
  {
    vector<Composite_path_step> steps;
    vector<Pair_transducer_scores> tmp_branch_sc;  // with singleton at root
  };
  typedef map<int,Composite_path> Clique_path;
  Clique_path clique_path;

  // tapes to reconstruct
  set<int> nodes_to_reconstruct;

  // composite transducer name
  sstring composite_name;

  // methods
  // constructors
  Transducer_SExpr_file() : etree(0) { }  // override default ETree constructor, which automatically creates one node
  Transducer_SExpr_file (SExpr& transducer_file_sexpr);
  Transducer_SExpr_file (const Transducer_SExpr_file& transducer_sexpr_file);

  // the following 'explicit' form of the constructor does the following:
  //  - puts all PVar's into defined_pvars[]
  //  - calls setup_cliques(), autoname_tree(), setup_composite_name()
  //  - initializes sym2tok, alphabet_size, trans_dict, branch_trans_name, composite_name
  Transducer_SExpr_file (const vector<sstring>& alphabet,
			 const PScores& pscores,
			 const ETree& etree,
			 const vector<Pair_trans_funcs>& branch_trans,
			 const NodeProfileMap& prof,
			 const NodePathMap& path,
			 const NodePathMap& old_path);

  // initializer (helper for explicit & copy constructors)
  void init();  // DO NOT CALL DIRECTLY; for use by constructors only

  // builder methods
  // method to read PVar defs
  void read_pvars (SExpr& transducer_file_sexpr, const char* keyword, bool is_bitscore);
  void set_pvar (const PVar& pvar, Score sc);

  // method to resolve S-expressions to PVar's
  PVar resolve_pvar (SExpr& var_sexpr);

  // method to resolve alphabet tokens
  int resolve_token (const sstring& tok);

  // recursive tree building method
  void build_node (SExpr& node, int parent, const map<sstring,sstring>& transducer_renaming_map);

  // method to rename tapes, branches & transducers
  void rebuild_tree_names (const map<int,sstring>& tape_name, const vector<sstring>& branch_name, const vector<sstring>& branch_trans_name);

  // methods to auto-name anonymous tree nodes & branches
  void autoname_tree();
  static sstring auto_tape_name (int node);
  sstring auto_branch_name (int node);

  // method to auto-build transducer dictionary
  void autobuild_trans_dict();

  // static method to munge PScores::group_suffix
  static void munge_group_suffix (PScores& pscores, PGroup& pg, const sstring& emit_label);

  // method to set up clique info
  void setup_cliques();

  // method to set composite name
  void setup_composite_name();

  // sequence pointer accessor for DP
  vector<Score_profile*> make_seq_vec();

  // accessor for vector of Pair_transducer_scores
  vector<Pair_transducer_scores> make_branch_sc();

  // vector of observed nodes
  vector<int> observed_nodes();

  // free clique set
  const set<int>& free_clique_set() const { return clique[free_clique]; }

  // EHMM builder
  inline EHMM_transducer_funcs build_composite()
  {
    EHMM_transducer_funcs ehmm_funcs (etree, branch_trans, &pscores, true);
    ehmm_funcs.tape_name = tape_name;
    return ehmm_funcs;
  }

  // method to peel a single constrained clique and return a profile of the hinge node (or other node, if prof_node is specified)
  // also computes the clique's score and adds to cumulative_score
  Score_profile peel_clique (int n_clique, const vector<Pair_transducer_scores>& branch_sc, Loge& clique_loglike, int prof_node = -1, bool normalize_prof_node = true);

  // method to peel off all constrained cliques and return a new Transducer_SExpr_file object.
  // the ETree of the new Transducer_SExpr_file is the subtree corresponding to the nodes in free_clique_set().
  // also returns the total log-likelihood of all the peeled cliques in peeling_loglike.
  Transducer_SExpr_file peel_constrained (Loge& peeling_loglike, bool normalize_prof_nodes = true);

  // method to forward-simulate, i.e. generate sequences and alignment paths
  void simulate();

  // methods to compute likelihood of a given path
  Loge get_path_loglike (const NodePathMap& node_path_map) const;
  NodePathMap get_path_map (const EHMM_transducer_funcs& ehmm, const vector<int>& composite_ehmm_trace) const;

  // method to automatically generate a proposal schedule
  void autopropose();

  // methods to get band diameter using banding coefficients
  // formula is: band_diameter = effective_tree_banding_coefficient * sqrt(rms_sequence_length)
  // where effective_tree_banding_coefficient = max_{leaf nodes L} sum_{nodes N from root to L} banding_coefficient[N]
  bool has_all_banding_coefficients() const;  // true only if all (unconstrained) branches have banding coefficients
  long get_band_diameter() const;  // calculated using above formula

  // output methods
  // output method for alphabet, parameters and transducers
  void show_defs (ostream& out);

  // output method for tree
  void show_tree (ostream& out, int base_indent = 0, bool show_composite_paths = false);

  // output methods for composed EHMM
  void show_composite (EHMM_transducer_funcs& ehmm, ostream& out);
  void show_composite_name_format (ostream& out);

  // output method for a Score_profile
  void show_prof (const Score_profile& prof_sc, ostream& out, int base_indent = 0);
  void show_ssm (const Symbol_score_map& ssm, ostream& out);

  // output methods for state paths
  template<class T>
  void show_path (const Transducer<T>& ehmm_sc, const vector<int>& path, const char* path_tag, ostream& out);

  template<class T>
  void show_path_types (const Transducer<T>& ehmm_sc, const vector<int>& path, ostream& out);

  // output methods for state paths with score & decomposition
  void show_path_with_tree (const EHMM_transducer_funcs& ehmm, const vector<int>& path, Loge path_ll, ostream& out, int base_indent = 0);
  void show_path_with_tree (const NodePathMap& path, Loge path_ll, ostream& out, int base_indent = 0);

  // output method for subtrees
  void show_clique (int n_clique, bool constrained, ostream& out, bool show_composite_paths = false);

  // output methods for bit scores
  sstring score_sexpr (Score sc);
  sstring score_sexpr (Loge ll);

  // adapter methods
  // method to return the Stockholm alignment implied by path & prof
  Stockade stockade();
};


// base class for singleton transducers (insert-only)
struct Singleton_transducer_base
{
  inline int insert_state() const { return 0; }
  inline int wait_state() const { return 1; }
};

// dummy singleton transducer whose path weight PFunc's are all 1
struct Singleton_transducer_funcs : Pair_transducer_funcs, Singleton_transducer_base
{
  Singleton_transducer_funcs (int alphabet_size = 1);
};

// score version of the singleton transducer; Score's are all 0
struct Singleton_transducer_scores : Pair_transducer_scores, Singleton_transducer_base
{
  Singleton_transducer_scores (int alphabet_size = 1);
};


// template method defs (thanks C++, a delight to work with you, as always)
template<class T>
void Transducer_SExpr_file::show_path (const Transducer<T>& ehmm_sc, const vector<int>& path, const char* path_tag, ostream& out)
{
  out << '(' << path_tag << " (";
  for (int pos = 0; pos < (int) path.size(); ++pos)
    out << (pos == 0 ? "" : " ") << ehmm_sc.sexpr_state_name (path[pos]);
  out << "))";
}

template<class T>
void Transducer_SExpr_file::show_path_types (const Transducer<T>& ehmm_sc, const vector<int>& path, ostream& out)
{
  out << '(' << TSEXPR_TYPE << " (";
  for (int pos = 0; pos < (int) path.size(); ++pos)
    {
      sstring state_type_sexpr;
      ehmm_sc.get_state_type_sexpr (path[pos], state_type_sexpr);
      out << (pos == 0 ? "" : " ") << state_type_sexpr;
    }
  out << "))";
}


#endif /* TRANSDUCER_SEXPR_INCLUDED */
