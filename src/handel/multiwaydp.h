#ifndef MULTIWAYDP_INCLUDED
#define MULTIWAYDP_INCLUDED

#include "handel/transducer_sexpr.h"

// structure describing a parameterization of the pruning/peeling matrix
struct Transducer_tree_emission
{
  // data
  const vector<const Symbol_score_map*>& node_scores;
  const vector<int>& branch_trans_state;
  ENode inserter;
  const vector<ENode>& deleters;
  const vector<ENode>& matchers;
  sstring state_name;

  // constructor
  Transducer_tree_emission (const vector<const Symbol_score_map*>& node_scores,
			    const vector<int>& branch_trans_state,
			    ENode inserter,
			    const vector<ENode>& deleters,
			    const vector<ENode>& matchers);
};

// derived Transducer_tree_emission class that owns some data
struct Transducer_tree_emission_auto
{
  // data
  vector<ENode> my_deleters, my_matchers;

  // constructor
  // initializes a TEmission object from branch_trans_state[], then calls summarize_emission() to get inserter, deleters[] & matchers[]
  Transducer_tree_emission_auto (const vector<const Symbol_score_map*>& node_scores,
				 const vector<int>& branch_trans_state);

};

// Score-based Felsenstein pruning & peeling matrix for transducers
struct Transducer_peeler
{
  // data
  int alphabet_size;

  // Let N be an ENode tree index and I,J,K be alphabet symbols.
  // Let T[N] be all observations in subtree rooted at node N.
  // Let Parent[N] be the parent of N.
  // F_matrix[N][I] = Prob2Score P({observations in T[N]} | state at N is I)
  // E_matrix[N][J] = Prob2Score P({observations in T[N]} | state at Parent[N] is J)
  // G_matrix[N][K] = Prob2Score P({observations NOT in T[N]} AND state at N is K)
  vector<vector<Score> > F_matrix, E_matrix, G_matrix;
  Score pruning_sc;  // Prob2Score P(all observations)

  // tree-related stuff
  ETree etree;
  ESubtrees subtrees;

  // tree->transducer mapping
  vector<Pair_transducer_scores> branch_transducer;
  vector<Pair_transducer_counts*> branch_trans_counts;

  // init methods
  void init (int alphabet_size,
	     const ETree& etree,
	     const vector<Pair_transducer_scores>& branch_transducer);

  void init (int alphabet_size,
	     const ETree& etree,
	     const vector<Pair_transducer_scores>& branch_transducer,
	     vector<Pair_transducer_counts>& branch_trans_counts);

  // prune method
  void prune (const Transducer_tree_emission& emit);

  // pruning matrix stochastic traceback
  vector<int> sample_trace (const Transducer_tree_emission& emit);

  // peel method; call after prune
  void peel (const Transducer_tree_emission& emit, double weight = 1.);  // weight is for counts

  // peeling matrix accessor
  Symbol_score_map post_sc (ENode node, bool normalize = true) const;
};

// Transducer_DP_base
struct Transducer_DP_base : Transducer_state_type_enum
{
  // typedefs
  typedef multi_array<Score> Cell_multi_array;

  // data that must be initialized by caller
  EHMM_transducer_scores* trans_sc;  // transducer
  vector<Score_profile*> seq;  // sequences. some of these pointers may be null, indicating unobserved (i.e. null) nodes
  bool want_counts;  // false by default; set to true if you want to do peeling & accumulate EM counts

  // internal data
  // final score
  Score end_sc;

  // alphabet size
  int alphabet_size;

  // indices of observed sequences, i.e. non-null entries in seq[] vector (set up by alloc())
  vector<int> observed_seqs;

  // annealing "temperature", set to 1.0 by default (probably due for deprecation)
  double kT;

  // work structures for Felsenstein pruning (set up by alloc())
  Transducer_peeler peeler;  // main DP matrix for pruning & peeling
  vector<const Symbol_score_map*> node_scores;
  vector<ENode> inserter;  //  indexed by state
  vector<vector<ENode> > deleters, matchers;  //  indexed by state
  vector<Pair_transducer_counts> branch_trans_counts;  // branch transducer counts; indexed by output node, only allocated if want_counts is true
  Transducer_counts ehmm_transition_counts;  // EHMM transition counts

  // incoming & outgoing states (set up by alloc())
  vector<vector<int> > incoming;  // incoming states by dest state
  vector<vector<int> > outgoing;  // outgoing states by src state

  // constructor
  Transducer_DP_base();

  // virtual destructor (does nothing)
  virtual ~Transducer_DP_base() { }

  // alloc method: sets up colmat, observed_seqs
  // requires caller to have set trans_sc, seq & tree (if tree is being used)
  virtual void alloc();

  // helper method to calculate emit scores by Felsenstein pruning
  Score calc_cell_emit (int state, const vector<int>& seq_coords, Prob weight = 1.);

  // accessors
  inline int states() const { return trans_sc->states(); }
  inline int sequences() const { return observed_seqs.size(); }
};

// Transducer_DP_matrix
struct Transducer_DP_matrix : Transducer_DP_base
{
  // data
  vector<int> dim;  // dim[i] = 1 + observed_seqs[i].size()
  vector<int> seqzero;  // vector of zeroes with same # of elements as dim[]

  // DP matrix
  vector<Cell_multi_array> cell_sc;  // indexed as cell_sc[state][seqcoords]

  // internal vars
  array2d<int> seq_delta;  // indexed as seq_delta[state][observed_seq_index]
  unsigned int n_types;  // number of possible state types = 1 << trans_sc->nodes()

  // Constructor: does nothing.
  // To initialise, set member variables (tree, trans_sc, seq) then call alloc()
  Transducer_DP_matrix() { }

  // alloc method
  // initialises observed_seq[], dim[] and cell
  void alloc();

  // accessors & helpers
  double cells() const;

  // display method
  void show (ostream& out);
};

// virtual interface to Transducer_forward_matrix
struct Transducer_forward_matrix_interface
{
  virtual void fill() = 0;
  virtual vector<int> sample_traceback() = 0;
  virtual ~Transducer_forward_matrix_interface() { }
};

// Transducer_forward_matrix
// NB assumes no null states!
struct Transducer_forward_matrix : Transducer_DP_matrix, Transducer_forward_matrix_interface
{
  // constructor
  Transducer_forward_matrix();

  // fill method
  void fill();

  // sample state path
  vector<int> sample_traceback();
};

// Transducer_Viterbi_matrix
// can handle null states just fine.
struct Transducer_Viterbi_matrix : Transducer_DP_matrix
{
  // constructor
  Transducer_Viterbi_matrix();

  // fill method
  void fill();

  // get Viterbi state path
  vector<int> traceback();
};

// Transducer_backward_matrix
// cell_sc[STATE][COORDS] includes score of emission from STATE at COORDS
typedef vector<multi_array<Prob> > Multiseq_prob_matrix;
struct Transducer_backward_matrix : Transducer_DP_matrix
{
  // data
  Transducer_forward_matrix fwd;  // the forward matrix
  Multiseq_prob_matrix cell_match_by_type;  // matrix of postprobs for optimal accuracy
  vector<unsigned int> compressed_state_type;  // mapping from trans_sc state type to a "compressed" state type (bits corresponding to unobserved seqs are removed)
  unsigned int n_compressed_types;  // range of compressed type
  bool fwd_initialized;  // flag indicating whether forward matrix has already been alloc'd and filled

  // optimal accuracy data
  vector<double> reward;  // mapping from [compressed] state type to reward
  multi_array<double> optacc;

  // constructor
  Transducer_backward_matrix (bool want_counts = false);  // set want_counts to collect EM-style Baum-Welch & peeling counts

  // derived alloc method
  // unless fwd_initialized flag has been set, this method copies tree, trans_sc and seq to fwd & calls fwd.alloc() before superclass alloc()
  void alloc();

  // method to intialize fwd from pre-filled Forward matrix
  // This method should be called before alloc().
  // NB: this leaves "private but exposed" parts of fwd uninitialized, while destroying the original Transducer_forward_matrix object,
  // so it should only be used for Forward-Backward AFTER sampling forward tracebacks from the original object.
  void init_fwd (Transducer_forward_matrix& existing_fwd);

  // fill method
  void fill();

  // method to increment expected PVar counts
  // (requires that caller set want_counts during construction or pre-alloc initialization)
  // "elim_sc", if non-null, must be a downcast pointer to the same object as "trans_sc".
  // If it's left null, then method won't do null state correction.
  void inc_var_counts (const EHMM_transducer_funcs& ehmm_funcs,
		       PCounts& var_counts,
		       const PScores& var_scores,
		       const Eliminated_EHMM_transducer_scores* elim_sc = 0,  // set this to do null state correction
		       Prob weight = 1.);

  // optimal accuracy methods
  void init_sumpairs_reward();  // set up reward using sum-of-pairs scoring
  void alloc_optacc();
  void fill_optacc();
  Alignment_path optacc_traceback();

  // optimal accuracy display
  void show_optacc (ostream& out);
};


#endif /* MULTIWAYDP_INCLUDED */
