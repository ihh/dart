#ifndef TRANSMAT_INCLUDED
#define TRANSMAT_INCLUDED

#include "seq/alignment.h"
#include "seq/local.h"
#include "seq/pfunc.h"

// Grammar state enumeration. A trivial class purely for inheriting from
struct Grammar_state_enum
{
  enum { Start = -1, End = -2, UndefinedState = -3 };    // indices of special states
  static inline int hmm_state_enum_idx (int state) { return state + 2; }   // "effective" state (i.e. leaving room for Start & End)
};

// transition matrix template
template <class T, class T_matrix = array2d<T> >
struct Transition_matrix : Grammar_state_enum, Stream_saver
{
  typedef T entry_type;

  typedef typename T_matrix::Row_view    T_row;
  typedef typename T_matrix::Column_view T_column;

  T_matrix _matrix;

  T_row    start;   // start[STATE] == transition (Start, STATE)
  T_column end;     //   end[STATE] == transition (STATE, End)

  static inline int tm_index (int state) { return state + 2; }  // backward compatibility wrapper
  static inline int tm_inverse (int state) { return state - 2; }  // backward compatibility wrapper

  // transition accessors

  inline int tm_states() const { return _matrix.xsize() - 2; }

  inline T& transition (int source, int dest)
  {
#ifdef DART_DEBUG
    if (source == End || dest == Start)
      DART_DEBUG_ERROR ("HMMs can't have *->Start or End->* transitions");
#endif /* DART_DEBUG */
    return _matrix (tm_index(source), tm_index(dest));
  }

  inline const T& transition (int source, int dest) const
  {
#ifdef DART_DEBUG
    if (source == End || dest == Start)
      DART_DEBUG_ERROR ("HMMs can't have *->Start or End->* transitions");
#endif /* DART_DEBUG */
    return _matrix (tm_index(source), tm_index(dest));
  }

  T& start_to_end() { return _matrix (tm_index(Start), tm_index(End)); }
  const T& start_to_end() const { return _matrix (tm_index(Start), tm_index(End)); }

  // constructors

  Transition_matrix() :
    _matrix (tm_index(0), tm_index(0)),
    start (_matrix, tm_index(Start), tm_index(0), tm_index(0)),
    end (_matrix, tm_index(End), tm_index(0), tm_index(0))
  { }

  Transition_matrix (int states) :
    _matrix (tm_index(states), tm_index(states)),
    start (_matrix, tm_index(Start), tm_index(0), tm_index(states)),
    end (_matrix, tm_index(End), tm_index(0), tm_index(states))
  { }
  
  Transition_matrix (int states, T t) :
    _matrix (tm_index(states), tm_index(states), t),
    start (_matrix, tm_index(Start), tm_index(0), tm_index(states)),
    end (_matrix, tm_index(End), tm_index(0), tm_index(states))
  { }

  // copy constructor must make its own start & end views onto transition matrix...

  Transition_matrix (const Transition_matrix<T,T_matrix>& tm) :
    _matrix (tm._matrix),
    start (_matrix, tm_index(Start), tm_index(0), tm_index(tm.tm_states())),
    end (_matrix, tm_index(End), tm_index(0), tm_index(tm.tm_states()))
  { }

  // virtual destructor
  virtual ~Transition_matrix<T,T_matrix>() { }

  // assignment operator must rebuild start & end views
  Transition_matrix<T,T_matrix>& operator= (const Transition_matrix<T,T_matrix>& tm)
  {
    assign_transition_matrix (tm);
    return *this;
  }

  void assign_transition_matrix (const Transition_matrix<T,T_matrix>& tm)
  {
    _matrix = tm._matrix;
    start = T_row (_matrix, tm_index(Start), tm_index(0), tm_index(tm_states()));
    end = T_column (_matrix, tm_index(End), tm_index(0), tm_index(tm_states()));
  }

  // equality test operator
  bool tm_equals (const Transition_matrix<T,T_matrix>& tm) const { return _matrix == tm._matrix; }

  // reset method
  void reset_transitions (const T& t) { template_for_contents (T_matrix, _matrix, tm) *tm = t; }

  // output method
  void show_transitions (ostream& o) const;

  // the show_transitions() method needs a bunch of virtuals
  // a templated base class might be cleaner, but i can't be arsed to code that right now...
  // (probably a better way to handle this would be with template specializations, c.f. Transducer_score_manipulator<> in "src/handel/transducer.h")
  virtual const char* element_descriptor() const { return "element"; }
  virtual int         element_width() const { return 1; }
  virtual void        show_element (const T& element, ostream& o) const { o << '?'; }

  // methods to render the transition matrix as a graphviz dot-format file
  virtual void print_dotfile (ostream& out) const;
  virtual void print_dotfile_nodes (ostream& out) const;
  virtual void print_dotfile_edges (ostream& out) const;
  virtual void print_dotfile_node (ostream& out, int state) const;
  virtual void print_dotfile_edge (ostream& out, int src, int dest) const;

  virtual sstring dotfile_graph_name() const;
  virtual sstring dotfile_node_id (int state) const;
  virtual sstring dotfile_node_label (int state) const;
  virtual sstring dotfile_edge_weight (int src, int dest) const;

  virtual map<sstring,sstring> dotfile_node_attrs (int state) const;
  virtual map<sstring,sstring> dotfile_edge_attrs (int src, int dest) const;

  // helper methods to determine if an element is null
  // ugly technique this... implementing a static method for every conceivable template type.... euuugh
  // (probably a better way to handle this would be with template specializations, c.f. Transducer_score_manipulator<> in "src/handel/transducer.h")
  inline static bool is_non_null (const sstring& s) { return s.size() > 0; }
  inline static bool is_non_null (const Score& sc) { return sc > -InfinityScore; }
  inline static bool is_non_null (const Prob& pr) { return pr > 0; }  // DODGY: potential for type clash between Prob & Loge
  inline static bool is_non_null (const PFunc& f) { return !f.is_null() && !f.is_zero(); }
};

typedef Transition_matrix<Score> Transition_scores;
typedef Transition_matrix <Score, array2d <Score, array2d_sparse_vector<Score> > > Sparse_transition_scores;

typedef Transition_matrix<Prob> Transition_counts;
typedef Transition_matrix <Prob, array2d <Prob, array2d_sparse_vector<Prob> > > Sparse_transition_counts;

typedef Transition_counts Transition_probs;

typedef Transition_matrix<PFunc> Transition_funcs;
typedef Transition_matrix <PFunc, array2d <PFunc, array2d_sparse_vector<PFunc> > > Sparse_transition_funcs;

template <class Transition_probs_base>
struct Concrete_transition_probs_template : Transition_probs_base
{
  // constructors
  // see also Transition_methods::score2prob
  Concrete_transition_probs_template (int size, Prob t) : Transition_probs_base (size, t) { }
  Concrete_transition_probs_template (const Transition_probs& m) : Transition_probs_base (m) { }

  // show methods
  const char* element_descriptor() const { return "probability"; }
  int         element_width() const { return 10; }
  void        show_element (const Prob& element, ostream& o) const
  {
    o.width(10);
    o << element;
  }
};

typedef Concrete_transition_probs_template<Transition_probs> Concrete_transition_probs;
typedef Concrete_transition_probs_template<Transition_counts> Concrete_transition_counts;

template <class Transition_scores_base>
struct Concrete_transition_scores_template : Transition_scores_base
{
  // constructors
  // see also Transition_methods::prob2score
  Concrete_transition_scores_template (int size, Score t) : Transition_scores_base (size, t) { }
  Concrete_transition_scores_template (const Transition_scores& m) : Transition_scores_base (m) { }

  // show methods
  const char* element_descriptor() const { return "score"; }
  int         element_width() const { return 10; }
  void        show_element (const Score& element, ostream& o) const
  {
    o.width(10);
    o << element;
  }
};

typedef Concrete_transition_scores_template<Transition_scores> Concrete_transition_scores;
typedef Concrete_transition_scores_template<Sparse_transition_scores> Concrete_sparse_transition_scores;

// static helper methods for Transition_scores
struct Transition_methods : Grammar_state_enum
{
  static vector<int> make_global_path (const vector<int>& local_path);  // adds Start & End states to beginning & end of path
  static vector<int> make_local_path (const vector<int>& global_path);  // strips Start & End states
  static vector<int> sample_local_path (const Transition_scores& ts, int start_state, const set<int>& end_states);
  static vector<int> sample_global_path (const Transition_scores& ts);
  static vector<int> consensus_local_path (const Transition_scores& ts, int start_state, const set<int>& end_states);
  static vector<int> consensus_global_path (const Transition_scores& ts);
  static Score path_transition_score (const Transition_scores& ts, const vector<int>& state_path);
  static void add_transition_counts_from_path (Transition_counts& tc, const vector<int>& state_path);

  // sort methods
  // for each state, find all the states that have transitions into it:
  template<class TScoreMatrix>
  static vector<vector<int> > incoming_states (const TScoreMatrix& ts);

  // for each state, find all the states in (selection) that it has transitions into, keeping the order specified by (selection):
  template<class TScoreMatrix>
  static vector<vector<int> > selected_outgoing_states (const TScoreMatrix& ts, const vector<int>& selection);

  // for each state, find all the states in (selection) that have transitions into it, keeping the order specified by (selection):
  template<class TScoreMatrix>
  static vector<vector<int> > selected_incoming_states (const TScoreMatrix& ts, const vector<int>& selection);

  // do a topological sort of null states (sources first), throwing an exception if a cycle is detected
  template<class TScoreMatrix>
  static vector<int> topological_sort (const TScoreMatrix& ts, const vector<int>& null_states, bool suppress_cycle_warning = false);

  // eliminate null states
  // returns a matrix with the same dimensions, but transitions between "emit" (i.e. non-null) states
  // are increased to effectively incorporate sum over all paths via null states.
  //
  // if zero2null is true, then transition probs to null states are set to zero;
  //  otherwise, they are set to the sum-over-null-paths probability.
  //
  // if s is non-null, then the sum-over-null-paths probabilities from null to emit states are stored in the null->emit transitions of *s.
  //
  // See (Holmes, 2003) or comments in C++ file for algebra of state elimination and restoral.
  static Concrete_transition_probs eliminate (const Transition_probs& tp, const vector<int>& null_states, bool zero2null = true, Concrete_transition_probs* s = 0);

  // sample path through null states
  // given an apparent path through emit states in a "null-states-eliminated" matrix,
  // samples a hidden path (possibly involving null states) through the original matrix
  // if choose_ML_path==true, the maximum-likelihood path is returned instead
  static vector<int> sample_eliminated (const Transition_probs& tp_orig,
					const Transition_probs& tp_elim,
					const vector<int>& null_states,
					const vector<int>& elim_path,
					bool choose_ML_path = false);

  // get transition counts involving null states
  // given a set of counts of apparent transitions between emit states in a "null-states-eliminated" matrix,
  // returns posterior expectations of counts for the original matrix (including transitions to/from null states)
  // TODO: implement this method!
  static Concrete_transition_counts count_eliminated (const Transition_probs& tp_orig,
						      const Transition_probs& tp_elim,
						      const Transition_probs& s_mx,
						      const vector<int>& null_states,
						      const vector<int>& emit_states,
						      const Transition_counts& elim_counts);

  // Score<-->Prob conversion
  static Concrete_transition_probs score2prob (const Transition_scores& ts);
  static Concrete_transition_scores prob2score (const Transition_probs& tp);
};

// template method code
#include "inline_transmat.h"

#endif
