#ifndef EM_MATRIX_BASE
#define EM_MATRIX_BASE

#include "newmat/newmat.h"
#include "tree/tree_alignment.h"
#include "tree/irrev_diag_matrix_factory.h"
#include "util/piper.h"
#include "util/vector_output.h"
#include "seq/pfunc.h"

// default resolution for timepoint cache
#define DEFAULT_TIMEPOINT_RES .01

// EQUIVALENT EIGENVALUE HACK: IH, 4/20/2005
// Numerical resolution for eigenvalues.
// If any two eigenvalues are closer together than this, they're forced to be equal.
// This should catch problems with "approximately equal" eigenvalues that were being missed.
#define EM_matrix_eigenvalue_tolerance 0.00001

// The basic parameters of an EM_matrix_base.
// This is a basic substitution model over a finite alphabet, with a hidden dynamic "site class" variable.
// (Holmes & Rubin, JMB, 2002. "An Expectation Maximisation Algorithm for Training Hidden Substitution Models".)
// The alphabet has A residues, while the site class variable can take C values.
// Thus, there are C*A states in the model. (NB different sites along the sequence are assumed to be independent.)
// If C=1, then the model reduces to a general reversible finite-alphabet substitution model.
// For example, DNA has A=4,C=1: the state space is {a,c,g,t}.
// DNA with two site classes has A=4,C=2: the state space is {a1,a2,c1,c2,g1,g2,t1,t2}.
// Instantaneous substitution events can change the site class or the residue, but not both.
// Thus, in the above two-class DNA model, events a1->a2 or a1->c1 are allowed, but a1->c2 is not
// (though a1 could mutate indirectly to c2 via multiple events, e.g. a1->c1->c2 or a1->a2->c2).
// This means that the full rate matrix is sparse; in this class, we store only the rates that can be nonzero,
// i.e. the intra-class rates (a1->c1) and the inter-class rates (a1->a2).
//
template <class T>
struct EM_matrix_params_template
{
  // data
  int C;  // site classes
  int A;  // alphabet size
  vector<T> pi;  // equilibrium frequencies
  vector<array2d<T> > X;  // intra-class rate matrices; outer vector indexed by class
  vector<array2d<T> > Y;  // inter-class rate matrices; outer vector indexed by residue

  // methods
  inline void init_matrix_template (int new_C, int new_A, T default_val)
  {
    C = new_C;
    A = new_A;
  
    pi = vector<T> (A*C, default_val);
    X = vector<array2d<T> > (C, array2d<T> (A, A, default_val));
    Y = vector<array2d<T> > (A, array2d<T> (C, C, default_val));
  }
};

struct EM_matrix_params : EM_matrix_params_template<double>
{
  // initialise method
  void init_matrix_params (int new_C, int new_A)
  { init_matrix_template (new_C, new_A, 0.); }

  // assign method
  void assign_matrix_params (const EM_matrix_params& params)
  { *this = params; }
};

// EM_matrix_base is a base class containing common methods for the various ways of
// estimating EM_matrix_params using the rate matrix Expectation Maximisation algorithm.
// (Holmes & Rubin, JMB, 2002. "An Expectation Maximisation Algorithm for Training Hidden Substitution Models".)
struct EM_matrix_base : Irrev_diagonalised_matrix_factory, EM_matrix_params, Piper
{
  // data
  // primary data (NB we also inherit data from EM_matrix_params)
  const Tree_alignment_database* align_db;  // alignment database
  double timepoint_res;  // resolution of time for Timepoint_cache

  // secondary data
  Matrix R;  // full rate matrix
  vector<Loge>  log_pi;  // logs of eqm freqs

  // by-timepoint matrix cache
  struct Timepoint_data
  {
    array2d<Loge> M;  // M(t)_{ij} = \log (exp(R*t)_{ij})
    array2d<Complex> J;  // J(k,l)(t)
    Timepoint_data (int m = 0);
  };
  typedef map<int,Timepoint_data> Timepoint_cache;
  Timepoint_cache timepoint_cache;

  // EM params
  double em_min_inc;  // min fractional log-likelihood improvement
  int em_max_iter;  // max iterations

  // Update flags: these can be used to keep some matrix elements fixed (e.g. at zero).
  // Initially, these are all set to true, but may later be cleared.
  // This mechanism is mostly used by input routines, which clear the update flag on encountering zero rates
  // (this was introduced to prevent zero rates from becoming finite due to numerical imprecision).
  // Note that diagonal elements are always updated by default, so the on-diagonal flags are ignored.
  // Note also that all these flags are ignored by the RIND_EM_matrix derived class.
  vector<array2d<short> > X_update_flag;  // flags indicating whether to update elements of X
  vector<array2d<short> > Y_update_flag;  // flags indicating whether to update elements of Y

  // Alphabet
  Alphabet hidden_alphabet;

  // methods
  // constructor
  EM_matrix_base (int C,
		  int A,
		  int max_fork = 1,
		  const Tree_alignment_database* align_db = 0,
		  double timepoint_res = DEFAULT_TIMEPOINT_RES);

  // virtual destructor
  virtual ~EM_matrix_base() { }

  // descriptor (for S-expression phylogrammar format; see source file "ecfg/ecfgsexpr.h")
  virtual const char* update_policy() const { return "uninitialised-update-policy"; }

  // initialisers
  void init_matrix (int C, int A, const Tree_alignment_database* align_db = 0, double timepoint_res = DEFAULT_TIMEPOINT_RES);
  void init_update_flags (bool default_flag);
  void init_cache (const Tree_alignment_database* align_db = 0, double new_timepoint_res = DEFAULT_TIMEPOINT_RES);
  void add_branch_lengths_to_cache (const PHYLIP_tree& tree);
  virtual void init_alphabet (const Alphabet& base_alphabet);

  // assignment method: initialise a matrix with C=1 from a Substitution_matrix_factory; calls init_matrix & update
  void assign (Substitution_matrix_factory& submat_factory);

  // accessors
  const Timepoint_data& timepoint_data (double t) const;

  // trivial methods
  inline int ca (int c, int a) const { return c * A + a; }  // index of state for class c, residue a
  inline int m() const { return C*A; }  // total number of states
  inline int state_class (int i) const { return i / A; }
  inline int state_residue (int i) const { return i % A; }

  // inherited from Diagonalised_matrix_factory
  const Alphabet& alphabet();

  // method to discretize branch lengths
  inline int discrete_time (double time) const
  {
    int tp = (int) (time / timepoint_res + .5);

    // IH, 4/20/2005
    // Explicitly guard against nonzero time values being
    // rounded down to zero.
    // Proposed by Nick Goldman.
    if (tp == 0 && time > 0)
      tp = 1;

    return tp;
  }

  // update methods for secondary data
  void update_R();  // updates R from (X,Y)
  void update_log_pi();  // updates log_pi[] and irrev_prior[] from pi[]
  virtual void diagonalize() = 0;  // delegate to subclass-specific diagonalization method
  void sanitize_eigenvalues();  // ensures no eigenvalues are positive
  void update_timepoint_cache();  // updates timepoint_cache from (S,mu,U); call after update_S()
  void update();  // calls all update methods, including diagonalize()

  // struct for EM update statistics
  struct Update_statistics
  {
    int states;
    Loge log_likelihood;  // log-likelihood
    vector<double> s;  // expected state start count
    vector<double> w;  // expected state wait time
    array2d<double> u;  // expected transition usage
    array2d<Complex> DJU;  // expected transition usage in eigenvector space
    array2d<Complex> DJU_rev;  // version of DJU with all branches flipped
    // methods
	// constructors
    Update_statistics (int states = 0);
    void clear (double pseud_init = 0., double pseud_mutate = 0., double pseud_wait = 0.);  // resets all counts
    void clear_DJU();  // clears DJU only
    void transform (const EM_matrix_base& hsm, bool symmetrise);  // transforms from DJU-space into u-space
    void check_waits_transitions (const EM_matrix_base& hsm);	// check for negative or zero waits, transitions
    friend ostream& operator<< (ostream& o, const Update_statistics& stats);
    void send (ostream& out) const;
    void receive (istream& in);
  };

  // wrapper for transform_waits_transitions
  // subclass-specific interpretation of "symmetrise" flag: reversible heeds it, irreversible ignores it
  virtual void transform_symmetrised_waits_transitions (Update_statistics& stats, bool symmetrise) const = 0;

  // back-transformation of wait times & transition usages from eigenspace
  void transform_waits_transitions (Update_statistics& stats, bool symmetrise) const;

  // interface to phylo_em.h
  // returns a matrix whose on-diagonal elements correspond to wait times, and off-diagonals are usage counts.
  array2d<double> clean_phylo_EM (int src, int dest, double T) const;

  // up/down algorithm for a single column
  struct Column_matrix
  {
    // The "signposts" for the Column_matrix.
    // These indicate cliques, gapped nodes, allowed states, etc.
    vector<int>          gapped;    // flag saying whether each node is gapped
    vector<vector<int> > allowed;   // list of states allowed at each node
    vector<int>          root;      // pointer to the clique root for each node
    vector<int>          leaf;      // flag saying whether each node is a leaf
    vector<int>          clique;    // list of clique root nodes

    // Parameters
    int nodes;
    int states;

    // U, D, L tables
    vector<vector<Loge> > U;  // U[node][state] = log P(subtree rooted at node 'node' | node 'node is in state 'state')
    vector<vector<Loge> > D;  // D[node][state] = log P(everything *except* subtree rooted at node 'node' | parent of 'node' is in state 'state')
    vector<Loge>          L;  // L[node] = log P(clique rooted at 'node')
    Loge tll;  // total log-likelihood

    // inferred class label string
    sstring class_label;
    inline bool estimate_class_labels() const { return class_label.size() != 0; }

    // allocator
    void alloc (int _nodes, int _states, bool alloc_class_labels = FALSE);

    // inline method to clear DP arrays
    inline void clear()
    {
      for (int n = 0; n < nodes; ++n)
	{
	  L[n] = -InfinityLoge;
	  for (int s = 0; s < states; ++s)
	    U[n][s] = D[n][s] = -InfinityLoge;
	}
      tll = -InfinityLoge;
    }

    // four alternative initialise methods
    void initialise (const PHYLIP_tree& tree, const vector<const Symbol_score_map*>& node_scores);  // null pointers signify gaps
    void initialise (const PHYLIP_tree& tree, const Alphabet& alphabet, const char* node_chars);
    void initialise (const Tree_alignment& tree_align, int column, const vector<int>& seq_coords, const Symbol_score_map* wildcard);

    // method to set up "signpost" arrays (gapped[], allowed[], root[], leaf[], clique[])
    // this method is called by the other initialise methods
    void init_signposts (const PHYLIP_tree& tree);

    // up/down algorithm
    // up...
    void fill_up (const EM_matrix_base& hsm,
		  const PHYLIP_tree& tree,
		  int column = -1);

    void fill_up (const EM_matrix_base& hsm,
		  const PHYLIP_tree& tree,
		  const sstring& position_descriptor);

    void fill_up (const vector<const EM_matrix_base*>& hsm,
		  const PHYLIP_tree& tree,
		  const sstring& position_descriptor);

    // down...
    void fill_down (const EM_matrix_base& hsm,
		    const PHYLIP_tree& tree,
		    Update_statistics& stats,
		    int column = -1,
		    double weight = 1.);

    void fill_down (const EM_matrix_base& hsm,
		    const PHYLIP_tree& tree,
		    Update_statistics& stats,
		    const sstring& position_descriptor,
		    double weight = 1.);

    void fill_down (const vector<const EM_matrix_base*>& hsm_vec,
		    const PHYLIP_tree& tree,
		    const vector<Update_statistics*>* stats_vec,
		    const sstring& position_descriptor,
		    double weight = 1.);

    // total log-likelihood accessor
    inline Loge total_log_likelihood() const { return tll; }

    // posterior probability accessors
    // (somewhat clumsily implemented)
    inline Prob node_post_prob (int node, int state, const PHYLIP_tree& tree,
				const EM_matrix_base& hsm) const
    {
      if (node == tree.root)
	return Nats2Prob (NatsPMul3 (hsm.log_pi[state], U[node][state], -L[node]));
      Prob p = 0.;
      for (int parent_state = 0; parent_state < hsm.m(); ++parent_state)
	p += branch_post_prob (node, parent_state, state, tree, hsm);
      return p;
    }

    inline Prob branch_post_prob (int node,
				  int parent_state,
				  int node_state,
				  const PHYLIP_tree& tree,
				  const EM_matrix_base& hsm) const
    {
      const int parent = tree.parent[node];
      const double branch_length = tree.branch_length (parent, node);
      return Nats2Prob (NatsPMul (NatsPMul (U[node][node_state], D[node][parent_state]),
				  NatsPMul (hsm.timepoint_data (branch_length).M (parent_state, node_state), -L[root[node]])));
    }

    // method to set total_log_likelihood(), called by fill_up()
    void calculate_total_log_likelihood();  // sums over cliques
  };

  // Column_matrix for each column of alignment
  struct Alignment_matrix
  {
    const EM_matrix_base& hsm;
    const Tree_alignment& tree_align;
    const PHYLIP_tree& tree;
    const Alignment& align;
    const int states;
    Column_matrix colmat;
    Loge total_loglike;  // summed over columns
    Symbol_score_map wildcard_ssm;
    // constructor
    Alignment_matrix (const EM_matrix_base& hsm, const Tree_alignment& ta, bool alloc_class_labels = FALSE);
    // up/down wrappers
    void fill_up();
    void fill_up_down (Update_statistics& stats);
    // up/down helpers
    void fill_initialised_colmat (int column);
    void check_finite() const;  // prints log messages
    inline Loge total_log_likelihood() const { return total_loglike; }
    // assignment operator, to keep vector<> happy
    Alignment_matrix& operator= (const Alignment_matrix& aln);
    // class label output
    vector<sstring> get_class_labels (Update_statistics& stats);  // calls fill_up, fill_down
    void display_class_labels (ostream& out);  // calls fill_up, fill_down
  };

  // up/down algorithm
  void up_down (Update_statistics& stats, bool likelihood_only, bool symmetrise, bool infer_class_labels = FALSE) const;
  Update_statistics get_stats_unequilibrated (bool symmetrise, bool infer_class_labels = FALSE) const;
  Loge log_likelihood() const;

  // class inference
  void display_classes (ostream& out) const;
  
  // EM algorithm(s)
  void randomise (double prior_dev = .1, double intra_min = .05, double intra_dev = .01, double inter_min = .005, double inter_dev = .001);
  void init_flat() { randomise (0., 1./(double)m(), 0., 1./(double)m(), 0.); }
  Loge iterate_quick_EM (bool intra = TRUE, bool inter = TRUE, int forgive = 20, bool infer_class_labels = FALSE);  // returns best log-likelihood
  virtual Update_statistics single_quick_EM (bool intra = TRUE, bool inter = TRUE, bool infer_class_labels = FALSE) = 0;  // does partial M-step only
  virtual void quick_M (const Update_statistics& stats, bool intra = TRUE, bool inter = TRUE) = 0;  // the M-step

  void normalise_substitution_rate();

  // I/O methods
  void write (ostream& out) const;
  void read (istream& in);
  void guess_alphabet();  // guesses an alphabet based on value of C

  // file format help method
  static bool hsm_help (Opts_list* ol);
};

typedef EM_matrix_base::Update_statistics Update_statistics;
typedef EM_matrix_base::Column_matrix Column_matrix;

#endif /* EM_MATRIX_BASE */
