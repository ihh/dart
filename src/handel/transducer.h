#ifndef TRANSDUCER_INCLUDED
#define TRANSDUCER_INCLUDED

#include <iomanip>
using namespace std;

#include "handel/ehmm.h"
#include "handel/transducer_keywords.h"
#include "seq/pkeywords.h"

#include "util/multi_array.h"
#include "seq/psexpr.h"
#include "hmm/transmat.h"
#include "hmm/pairhmm.h"
#include "hsm/em_matrix_base.h"

// default transducer state names
#define Default_name_for_transducer_start_state Grammar_start_state_name
#define Default_name_for_transducer_end_state   Grammar_end_state_name
#define Undefined_transducer_state_name         Grammar_undef_state_name

// max number of tries for feeding a parent sequence through a transducer
// (some kind of DP algorithm would be more robust..)
#define TRANSDUCER_SAMPLE_MAX_TRIES 1000

// graphviz dotfile strings for emissions & states
#define Transducer_emit_bgcolor   0xff0000
#define Transducer_null_bgcolor   0x7f5f5f

#define Transducer_match_bgcolor  0xff00ff
#define Transducer_insert_bgcolor 0x00ffff
#define Transducer_delete_bgcolor 0xffff00

#define Transducer_wait_bgcolor   0xa0a0a0

#define Transducer_undef_bgcolor  0xffffff

#define Transducer_start_bgcolor  0xa0a0a0
#define Transducer_end_bgcolor    0xa0a0a0

#define Transducer_white_bgcolor  0xffffff

#define Transducer_fontname       "Courier"

// Multiple-sequence, single-terminal transducer state type enum
// Nodes (i.e. sequences) are numbered as per "ehmm.h":
//  -1 (root), 0 (subroot), 1...
struct Transducer_state_type_enum
{
  typedef unsigned int State_type;
  enum { TransducerStartType = -1, TransducerEndType = -2, TransducerUndefinedType = -3, TransducerWaitType = 0, TransducerDeleteType = 1, TransducerInsertType = 2, TransducerMatchType = 3 };
  static int type_node_emit (State_type t, int node)
  { return (t >> (node + 1)) & 1; }
};
typedef Transducer_state_type_enum::State_type Transducer_state_type;

// Transducer score/PFunc abstraction classes
// logic for PFunc multiplication is different from Score multiplication
// the following functions will be statically bound
//  to the right type in the different template classes
template<class T> struct Transducer_score_manipulator { };

template<> struct Transducer_score_manipulator<Score>
{
  static inline void set_to_zero (Score& x) { x = -InfinityScore; }
  static inline void set_to_one (Score& x) { x = 0; }
  static inline void pmul_acc (Score& x, const Score& y)
  { ScorePMulAcc (x, y); }
  static inline void psum_acc (Score& x, const Score& y)
  { ScorePSumAcc (x, y); }
  static inline bool is_nonzero (const Score& x)  { return x > -InfinityScore; }
  static inline bool is_one (const Score& x)  { return x == 0; }

  static inline sstring pval_string (const Score& x)
  { sstring s; s << "2^(1/1000 * "; ShowScore(x,s); s << ")"; return s; }

  static inline sstring score_label_string (const Score& x, const PScores* dummy_pscores)
  { return sstring(); }

  // ignore eval_sc flag: Score containers don't need to evaluate PFunc's in order to show Score values
  static inline sstring score_bitvalue_string (const Score& x, const PScores* dummy_pscores, bool eval_sc)
  { sstring s; if (x > -InfinityScore) s << -Score2Bits(x); else s << PK_INFINITE; return s; }
};

template<> struct Transducer_score_manipulator<PFunc>
{
  static inline void set_to_zero (PFunc& x) { x = 0.; }
  static inline void set_to_one (PFunc& x) { x = 1.; }
  static inline void pmul_acc (PFunc& x, const PFunc& y) { if (y.is_null() || y.is_zero()) set_to_zero(x); else if (x.is_one()) x = y; else if (!x.is_null() && !x.is_zero() && !y.is_one()) x *= y; }
  static inline void psum_acc (PFunc& x, const PFunc& y) { if (x.is_zero()) x = y; else if (!y.is_null() && !y.is_zero()) x += y; }
  static inline bool is_nonzero (const PFunc& x)  { return !x.is_null() && !x.is_zero(); }
  static inline bool is_one (const PFunc& x)  { return x.is_one(); }

  static inline sstring pval_string (const PFunc& x)
  { sstring s; x.show(s); return s; }

  static inline sstring score_label_string (const PFunc& x, const PScores* pscores)
  { sstring s; PFunc_builder::pfunc2stream (s, *pscores, x); return s; }

  static inline sstring score_bitvalue_string (const PFunc& x, const PScores* pscores, bool eval_sc)
  { sstring s; if (eval_sc) { const Score sc = x.eval_sc (*pscores); if (sc > -InfinityScore) s << -Score2Bits(sc); else s << PK_INFINITE; } return s; }
};

template<> struct Transducer_score_manipulator<Prob>
{
  static inline void set_to_zero (Prob& x) { x = 0.; }
  static inline void set_to_one (Prob& x) { x = 1.; }
  static inline void pmul_acc (Prob& x, const Prob& y)
  { x *= y; }
  static inline void psum_acc (Prob& x, const Prob& y)
  { x += y; }
  static inline bool is_nonzero (const Prob& x)  { return x > 0.; }
  static inline bool is_one (const Prob& x)  { return x == 1.; }

  static inline sstring pval_string (const Prob& x)
  { sstring s; s << x; return s; }

  static inline sstring score_label_string (const Prob& x, const PScores* dummy_pscores)
  { return sstring(); }

  // ignore eval_sc flag: Prob containers don't need to evaluate PFunc's in order to show Prob values
  static inline sstring score_bitvalue_string (const Prob& x, const PScores* dummy_pscores, bool eval_sc)
  { sstring s; if (x > 0.) s << -Prob2Bits(x); else s << PK_INFINITE; return s; }
};

// Multiple-sequence, single-terminal transducer
// (Single-terminal transducers can be associated with
//  homogeneous substitution processes on Felsenstein trees)
template<class T>
struct Transducer
  : Transition_matrix<T>,
    Transducer_state_type_enum,
    Transducer_score_manipulator<T>
{
  // typedefs
  typedef T entry_type;

  // state types
  vector<State_type> state_type;

  // compact state names indicating branch states & tree structure
  vector<sstring> state_name;
  sstring start_name, end_name;
  map<int,sstring> tape_name;  // names of tapes, indexed by node; can be left blank

  // state descriptors for "nodes" in graphviz dotfiles
  map<int,sstring> node_html;  // HTML TABLE's
  map<int,sstring> node_record;  // native dotfile record format

  // edge descriptors
  Transition_matrix<sstring> edge_label;

  // wrapper to access state name; safe for Start & End states
  sstring get_state_name (int state) const
  {
    if (state >= 0 && state < states())
      return state_name[state];
    if (state == Grammar_state_enum::Start)
      return start_name;
    if (state == Grammar_state_enum::End)
      return end_name;
    return Undefined_transducer_state_name;
  }

  // wrapper to access tape name
  sstring get_tape_name (int n) const
  {
    sstring s;
    if (tape_name.find (n) == tape_name.end())
      s << n + 1;
    else
      s = ((map<int,sstring>&) tape_name)[n];  // cast away const
    return s;
  }

  // wrapper to get state name --> index map
  map<sstring,int> state_index() const
  {
    map<sstring,int> si;
    si[start_name] = Grammar_state_enum::Start;
    si[end_name] = Grammar_state_enum::End;
    for (int s = 0; s < (int) state_name.size(); ++s)
      si[state_name[s]] = s;
    return si;
  }

  // wrapper to get token --> tape index map
  map<sstring,int> tape_index() const
  {
    map<sstring,int> ti;
    for (map<int,sstring>::const_iterator tn_ptr = tape_name.begin(); tn_ptr != tape_name.end(); ++tn_ptr)
      ti[tn_ptr->second] = tn_ptr->first;
    return ti;
  }

  // number of nodes
  inline int nodes() const
  {
    State_type all = 0;
    for_const_contents (vector<State_type>, state_type, t)
      all |= *t;
    int n = 0;
    while (all) { ++n; all >>= 1; }
    return n;
  };

  // max node index is nodes()-2, since min is -1
  inline int max_node_index() const { return nodes() - 2; }

  // transducer is "joint" if there are no emissions to the root node (-1),
  //  and "conditional" otherwise
  inline bool is_joint() const
  {
    for_const_contents (vector<State_type>, state_type, t)
      if (type_node_emit (*t, -1))  // test for emissions to node -1
	return false;
    return true;
  }

  inline bool is_conditional() const
  {
    return !is_joint();
  }

  // constructors
  Transducer (int states)
    : Transition_matrix<T> (states),
      state_type (states, 0),
      state_name (states),
      start_name (Default_name_for_transducer_start_state),
      end_name (Default_name_for_transducer_end_state),
      edge_label (states)
  { }

  Transducer (int states, T t)
    : Transition_matrix<T> (states, t),
      state_type (states, 0),
      state_name (states),
      start_name (Default_name_for_transducer_start_state),
      end_name (Default_name_for_transducer_end_state),
      edge_label (states)
  { }


  // templated copy constructors
  Transducer<T> (const Transducer<T>& trans)
    : Transition_matrix <T> (trans),
      state_type (trans.state_type),
      state_name (trans.state_name),
      start_name (trans.start_name),
      end_name (trans.end_name), 
      tape_name (trans.tape_name), 
      node_html (trans.node_html), 
      node_record (trans.node_record),
      edge_label (trans.edge_label)
  { }

  template<class S>
  Transducer<T> (const Transducer<S>& trans, const T& val)
    : Transition_matrix<T> (trans.states(), val),
      state_type (trans.state_type),
      state_name (trans.state_name),
      start_name (trans.start_name),
      end_name (trans.end_name), 
      tape_name (trans.tape_name), 
      node_html (trans.node_html), 
      node_record (trans.node_record),
      edge_label (trans.edge_label)
  { }

  // resize method
  void resize (int states, T t)
  {
    Transition_matrix<T> transmat (states, t);
    ((Transition_matrix<T>*) this)->assign_transition_matrix (transmat);

    state_type = vector<State_type> (states, 0);
    state_name = vector<sstring> (states);

    Transition_matrix<sstring> new_edge_label (states);
    edge_label.assign_transition_matrix (new_edge_label);
  }

  // accessors
  inline int states() const
  { return state_type.size(); }   // excludes start & end states

  // dx helper
  inline int node_emit (int state, int node) const
  { return type_node_emit (state_type[state], node); }

  // helper to convert state path to alignment path
  Alignment_path state2align_path (const vector<int>& state_path, bool include_subroot = false) const
  {
    const int rows = include_subroot ? nodes() : (nodes() - 1);
    Alignment_path a (rows);
    vector<bool> col_data (rows);
    for_const_contents (vector<int>, state_path, s)
      if (*s >= 0)
	{
	  for (int n = 0; n < rows; ++n)
	    col_data[n] = node_emit (*s, include_subroot ? (n-1) : n);
	  a.append_column (col_data);
	}
    return a;
  }

  // helper to get lists of emit & null states for a given set of observed nodes
  // NB Start & End are excluded from the null states list
  void get_emit_and_null_states (const vector<int>& observed_nodes, vector<int>& emit_states, vector<int>& null_states) const
  {
    for (int s = 0; s < states(); ++s)
      {
	const Transducer_state_type_enum::State_type type = state_type[s];
	bool is_null = true;
	for_const_contents (vector<int>, observed_nodes, obs)
	  if (Transducer_state_type_enum::type_node_emit (type, *obs))
	    {
	      is_null = false;
	      break;
	    }
	if (is_null)
	  null_states.push_back (s);
	else
	  emit_states.push_back (s);
      }
  }

  // show method
  void show (ostream& o) const
  {
    o << "Transducer state types:\n";
    for (int s = 0; s < states(); ++s)
      {
	o << "State " << s;
	if (state_name[s].size())
	  o << '[' << state_name[s] << ']';

	unsigned int t = state_type[s];
	if (t)
	  {
	    o << " emits to";
	    for (int seq = 0; t != 0; t >>= 1, ++seq)
	      if (t & 1)
		o << ' ' << seq - 1;
	  }
	else
	  o << " is null";
	o << "\n";
      }
    Transition_matrix<T>::show_transitions(o);
  }

  // show_sexpr helpers
  virtual const char* transducer_keyword() const { return TSEXPR_TRANSDUCER; }

  virtual void get_state_type_sexpr (int state, sstring& sts) const
  {
    int syms = 0;
    sts << '(';
    if (state >= 0)
      {
	for (int n = -1; n < nodes(); ++n)
	  if (type_node_emit (state_type[state], n))
	    sts << (syms++ ? " " : "")
		<< get_tape_name(n);
      }
    sts << ')';
  }

  virtual void get_emit_label_sexpr (int state, sstring& ls) const { }
  virtual sstring transition_label_sexpr (const T& t, const PScores* ps) const { return ((Transducer_score_manipulator<T>*) this)->score_label_string (t, ps); }
  virtual bool has_emit_labels() const { return false; }

  // show_sexpr method
  // if pscores is non-null, it is used to display transition labels
  // if eval_sc is true, then PFunc's will be evaluated and bitscores displayed
  void show_sexpr (ostream& out, const sstring& name, const PScores* pscores = 0, bool eval_sc = true)
  {
    out << '(' << transducer_keyword() << '\n';
    if (name.size())
      out << " (" << TSEXPR_NAME << ' ' << name << ")\n";
    out << '\n';

    for (int state = Grammar_state_enum::End; state < states(); ++state)
      {
	out << " (" << TSEXPR_STATE
	    << " (" << TSEXPR_NAME << ' ' << sexpr_state_name(state) << ')';
	if (state >= 0)
	  {
	    sstring state_type_sexpr, label_sexpr;
	    get_state_type_sexpr (state, state_type_sexpr);
	    get_emit_label_sexpr (state, label_sexpr);
	    out << " (" << TSEXPR_TYPE << ' ' << state_type_sexpr << ')';
	    if (label_sexpr.size())
	      out << " (" << TSEXPR_LABEL << ' ' << label_sexpr << ')';
	  }
	else
	  out << " (" << TSEXPR_TYPE << ' '
	      << (state == Grammar_state_enum::Start ? TSEXPR_START : TSEXPR_END)
	      << ')';
	out << ")\n";
      }

    out << '\n';
    for (int src = Grammar_state_enum::Start; src < states(); ++src)
      {
	bool src_visible = false;
	for (int dest = Grammar_state_enum::End; dest < states(); ++dest)
	  if (dest != Grammar_state_enum::Start)
	    {
	      T& t = Transition_matrix<T>::transition (src, dest);
	      if (((Transducer_score_manipulator<T>*) this)->is_nonzero (t))
		{
		  src_visible = true;
		  out << " (" << TSEXPR_TRANSITION
		      << " (" << TSEXPR_FROM << ' ' << sexpr_state_name(src) << ')'
		      << " (" << TSEXPR_TO << ' ' << sexpr_state_name(dest) << ')';
		  if (pscores != 0 && !((Transducer_score_manipulator<T>*) this)->is_one(t))
		    {
		      const sstring bvs = ((Transducer_score_manipulator<T>*) this)->score_bitvalue_string (t, pscores, eval_sc);
		      if (bvs.size())
			out << " (" << TSEXPR_BITVALUE << ' ' << bvs << ')';
		      const sstring tls = transition_label_sexpr (t, pscores);
		      if (tls.size())
			out << " (" << TSEXPR_LABEL << ' ' << tls << ')';
		    }
		  out << ")\n";
		}
	    }
	if (src_visible)
	  out << '\n';
      }

    out << ")\n";
  }

  // helper method to get state name in an S-expression format. EHMM_transducer overrides this with a disgustingly-implemented hack
  virtual sstring sexpr_state_name (int state) const { return get_state_name (state); }

  // virtual Transition_matrix methods....
  // historically these have been implemented separately for each subclass,
  //  but I don't think we really need them,
  //  so I'm putting braindead code in here
  // (probably a better way to handle this would be with template specializations, c.f. Transducer_score_manipulator<>)
  const char* element_descriptor() const { return "entries"; }
  int         element_width() const { return 10; }

  // dotfile stuff
  sstring dotfile_node_label (int state) const
  {
    sstring label;
    label << "<<TABLE CELLSPACING=\"0\" CELLBORDER=\"0\"><TR><TD>";
#if 0
    label << "<TABLE BORDER=\"0\" CELLBORDER=\"0\">";
    const int n_nodes = nodes();
    for (int n = is_joint() ? 0 : -1; n < n_nodes - 1; ++n)
      {
	const bool emits = state >= 0 ? type_node_emit (state_type[state], n) : false;
	label.fill ('0');
	label << "<TR><TD BGCOLOR=\"#"
	      << setbase(16) << setw(6)
	      << (emits ? Transducer_emit_bgcolor : Transducer_null_bgcolor)
	      << setbase(10)
	      << "\">"
	      << (emits ? '*' : '-')
	      << "</TD></TR>";
      }
    label << "</TABLE></TD><TD>";
#endif
    if (node_html.find (state) != node_html.end())
      label << node_html.find(state)->second;
    else
      label << Transition_matrix<T>::dotfile_node_label (state);  // call super
    label << "</TD></TR></TABLE>";
    label << ">";
    return label;
  }
};

// Pair_transducer
template<class T>
struct Pair_transducer : Transducer<T>, TSpaceEnum
{
  // typedefs
  typedef T entry_type;
  typedef Transducer<T> transducer_type;

  // state emit labels (matrices for match states, vectors for insert/delete states, ignored for wait states)
  vector<sstring> emit_label;
  int alphabet_size;

  // pair_emit: emission scores/functions
  // for X-absorption Y-emission (match states), entry is pair_emit(X,Y)
  // for X-absorptions (delete states), entry is pair_emit(X,0)
  // for Y-emissions (insert states), entry is pair_emit(Y,0)
  vector<array2d<T> >  pair_emit;

  // accessors for pair_emit
  inline T& match_val (int state, int xsym, int ysym) { return static_match_val (pair_emit[state], xsym, ysym); }
  inline T& insert_val (int state, int ysym) { return static_insert_val (pair_emit[state], ysym); }
  inline T& delete_val (int state, int xsym) { return static_delete_val (pair_emit[state], xsym); }

  static inline T& static_match_val (array2d<T>& pair_emit_array, int xsym, int ysym) { return pair_emit_array (xsym, ysym); }
  static inline T& static_insert_val (array2d<T>& pair_emit_array, int ysym) { return pair_emit_array (ysym, 0); }
  static inline T& static_delete_val (array2d<T>& pair_emit_array, int xsym) { return pair_emit_array (xsym, 0); }

  // allocator for pair_emit
  void alloc_pair_emit (T t)
  {
    pair_emit = vector<array2d<T> > (Transducer<T>::states());
    for (int s = 0; s < Transducer<T>::states(); ++s)
      if (Transducer<T>::state_type[s] == Transducer_state_type_enum::TransducerInsertType || Transducer<T>::state_type[s] == Transducer_state_type_enum::TransducerDeleteType)
	pair_emit[s] = array2d<T> (alphabet_size, 1, t);
      else if (Transducer<T>::state_type[s] == Transducer_state_type_enum::TransducerMatchType)
	pair_emit[s] = array2d<T> (alphabet_size, alphabet_size, t);
  }

  // constructors
  Pair_transducer (int states = 0) : Transducer<T> (states), emit_label (states), alphabet_size (1) { }
  Pair_transducer (int states, T t) : Transducer<T> (states, t), emit_label (states), alphabet_size (1) { }

  // templated copy constructors
  Pair_transducer<T> (const Pair_transducer<T>& pair_trans)
    : Transducer<T> (pair_trans),
      emit_label (pair_trans.emit_label),
      alphabet_size (pair_trans.alphabet_size),
      pair_emit (pair_trans.pair_emit)
  { }

  template<class S>
  Pair_transducer<T> (const Pair_transducer<S>& pair_trans, const T& val)
    : Transducer<T> (pair_trans, val),
      emit_label (pair_trans.emit_label),
      alphabet_size (pair_trans.alphabet_size)
  {
    alloc_pair_emit (val);
  }

  // TState to Transducer_state_type map
  inline static Transducer_state_type tstate2type (TState tstate)
  {
    switch (tstate)
      {
      case TransWait:
	return (Transducer_state_type) Transducer_state_type_enum::TransducerWaitType;
      case TransDelete:
	return (Transducer_state_type) Transducer_state_type_enum::TransducerDeleteType;
      case TransInsert:
	return (Transducer_state_type) Transducer_state_type_enum::TransducerInsertType;
      case TransMatch:
	return (Transducer_state_type) Transducer_state_type_enum::TransducerMatchType;
      default:
	THROWEXPR ("Can't map TState "
		   << tstate << " to pair transducer state type");
	break;
      }
    return (Transducer_state_type) 0;
  }

  // Transducer_state_type to TState map
  inline static TState type2tstate (Transducer_state_type type)
  {
    switch (type)
      {
      case 0:
	return (TState) TransWait;
      case 1:
	return (TState) TransDelete;
      case 2:
	return (TState) TransInsert;
      case 3:
	return (TState) TransMatch;
      default:
	THROWEXPR ("Can't map pair transducer state type "
		   << type << " to TState");
	break;
      }
    return (TState) TransUndef;
  }

  // display method for emit labels
  void show_pair_emit (ostream& out) const
  {
    out << "Emit labels:\n";
    for (int s = 0; s < transducer_type::states(); ++s)
      switch (transducer_type::state_type[s])
	{
	case Transducer_state_type_enum::TransducerInsertType:
	  out << "Insert state " << transducer_type::state_name[s] << " [" << s << "]\n";
	  for (int i = 0; i < alphabet_size; ++i)
	    out << ' ' << ((Pair_transducer<T>&)*this).insert_val(s,i);  // cast away const
	  out << '\n';
	  break;
	case Transducer_state_type_enum::TransducerDeleteType:
	  out << "Delete state " << transducer_type::state_name[s] << " [" << s << "]\n";
	  for (int i = 0; i < alphabet_size; ++i)
	    out << ' ' << ((Pair_transducer<T>&)*this).delete_val(s,i);  // cast away const
	  out << '\n';
	  break;
	case Transducer_state_type_enum::TransducerMatchType:
	  out << "Match state " << transducer_type::state_name[s] << " [" << s << "]\n";
	  for (int i = 0; i < alphabet_size; ++i)
	    {
	      for (int j = 0; j < alphabet_size; ++j)
		out << ' ' << ((Pair_transducer<T>&)*this).match_val(s,i,j);  // cast away const
	      out << '\n';
	    }
	  break;
	default:
	  break;
	}
  }

  // show_sexpr helpers
  void get_state_type_sexpr (int state, sstring& sts) const
  {
    switch (state)
      {
      case Grammar_state_enum::Start: sts << TSEXPR_START; break;
      case Grammar_state_enum::End: sts << TSEXPR_END; break;
      default:
	switch (Transducer<T>::state_type[state])
	  {
	  case Transducer_state_type_enum::TransducerWaitType: sts << TSEXPR_WAIT; break;
	  case Transducer_state_type_enum::TransducerInsertType: sts << TSEXPR_INSERT; break;
	  case Transducer_state_type_enum::TransducerDeleteType: sts << TSEXPR_DELETE; break;
	  case Transducer_state_type_enum::TransducerMatchType: sts << TSEXPR_MATCH; break;
	  default: THROWEXPR ("Unknown state type"); break;
	  }
	break;
      }
  }

  bool has_emit_labels() const
  {
    for (int s = 0; s < Transducer<T>::states(); ++s)
      if (emit_label[s].size())
	return true;
    return false;
  }

  void get_emit_label_sexpr (int state, sstring& ls) const
  {
    ls = emit_label[state];
  }
};

// EHMM_transducer
template<class T, class PairTrans = Pair_transducer<T> >
struct EHMM_transducer : TSpaceEnum, Transducer<T>
{
  // typedefs
  typedef T entry_type;
  typedef Transducer<T> transducer_type;
  typedef PairTrans pair_transducer_type;
  typedef vector<vector<int> > Branch_trans_states;

  // typedefs for build()
  typedef map<sstring,vector<int> > States_by_EHMM_type;

  // tree and branch transducers; alphabet size
  ETree etree;
  vector<pair_transducer_type> branch_transducer;
  int alphabet_size;

  // lookup of transducer states by state index and ETree-node
  Branch_trans_states branch_trans_states;

  // lookup of EHMM emission/absorption profiles by state index
  vector<TEmission> emission;
  vector<TTerm> absorption;

  // lookup of EHMM state types by state index
  map<int,EState> state2estate;

  // override virtual display method
  void show_element (const Score& sc, ostream& out) const { ShowScore(sc,out); }
  void show_element (const PFunc& f, ostream& out) const { f.show (out); }

  // helper objects for build()

  // context structure for states
  struct State_context
  {
    // references
    const int s;
    const ETree& etree;
    const vector<pair_transducer_type>& branch_transducer;
    const Branch_trans_states& branch_trans_states;
    const vector<TEmission>& emission;
    const EState& estate;

    // helpers for subclass
    sstring get_state_name (int n) const
    {
      if (s == Grammar_state_enum::Start)
	return branch_transducer[n].start_name;
      if (s == Grammar_state_enum::End)
	return branch_transducer[n].end_name;
      return branch_transducer[n].get_state_name (branch_trans_states[s][n]);
    }

    TTerm get_tterm (int n) const { return s < 0 ? TTermNull : emission[s].tterm[n]; }

    // constructors
    State_context (const int s,
		   const ETree& etree,
		   const vector<pair_transducer_type>& branch_transducer,
		   const Branch_trans_states& branch_trans_states,
		   const vector<TEmission>& emission,
		   const EState& estate)
      : s (s),
	etree (etree),
	branch_transducer (branch_transducer),
	branch_trans_states (branch_trans_states),
	emission (emission),
	estate (estate)
    { }

    State_context (const State_context& context)
      : s (context.s),
	etree (context.etree),
	branch_transducer (context.branch_transducer),
	branch_trans_states (context.branch_trans_states),
	emission (context.emission),
	estate (context.estate)
    { }
  };

  // State_label
  struct State_label : State_context
  {
    // data
    vector<sstring> buffer;  // partial name buffer
    vector<ENode> children;  // children of current node
    int n;  // running node index

    // accessors
    sstring& current() { return buffer[n]; }
    sstring n_name() const { return State_context::get_state_name (n); }
    TTerm tterm() const { return State_context::get_tterm (n); }

    // n_emits
    unsigned long n_emits() const
      {
	if (!n_changed())
	  return false;

	const pair_transducer_type& branch_trans = State_context::branch_transducer[n];
	if (State_context::s == Grammar_state_enum::Start
	    || State_context::s == Grammar_state_enum::End)
	  return false;
	
	const int branch_state = State_context::branch_trans_states[State_context::s][n];
	if (branch_state == Grammar_state_enum::Start
	    || branch_state == Grammar_state_enum::End)
	  return false;

	const Transducer_state_type branch_type = branch_trans.state_type[branch_state];
	return
	  branch_type == Transducer_state_type_enum::TransducerInsertType
	  || branch_type == Transducer_state_type_enum::TransducerMatchType;
      }

    // n_absorbs
    unsigned long n_absorbs() const
      {
	if (!n_changed())
	  return false;

	const pair_transducer_type& branch_trans = State_context::branch_transducer[n];
	if (State_context::s == Grammar_state_enum::Start
	    || State_context::s == Grammar_state_enum::End)
	  return false;
	
	const int branch_state = State_context::branch_trans_states[State_context::s][n];
	if (branch_state == Grammar_state_enum::Start
	    || branch_state == Grammar_state_enum::End)
	  return false;

	const Transducer_state_type branch_type = branch_trans.state_type[branch_state];
	return
	  branch_type == Transducer_state_type_enum::TransducerMatchType
	  || branch_type == Transducer_state_type_enum::TransducerDeleteType;
      }

    // n_changed
    bool n_changed() const
    {
      const vector<ENode> chang = State_context::estate.collapsedChanges();
      return find (chang.begin(), chang.end(), (ENode) n) != chang.end();
    }

    // emit_color
    unsigned long emit_color() const
    {
      return n_emits() ? Transducer_emit_bgcolor : Transducer_null_bgcolor;
    }

    // emit_color
    unsigned long absorb_color() const
    {
      return n_absorbs() ? Transducer_emit_bgcolor : Transducer_null_bgcolor;
    }

    unsigned long fade_color (unsigned long int col)
    {
      unsigned long int newcol = 0;
      for (int b = 0; b < 3; ++b)
	{
	  const int cpt = (col >> (8*b)) & 0xff;
	  const int undef = (Transducer_undef_bgcolor >> (8*b)) & 0xff;
	  const int newcpt = (cpt + undef) >> 1;
	  newcol |= (newcpt & 0xff) << (8*b);
	}
      return newcol;
    }

    // branch_color
    unsigned long branch_color()
    {
      const pair_transducer_type& branch_trans = State_context::branch_transducer[n];

      if (State_context::s == Grammar_state_enum::Start)
	return Transducer_start_bgcolor;
      if (State_context::s == Grammar_state_enum::End)
	return Transducer_end_bgcolor;

      const int branch_state = State_context::branch_trans_states[State_context::s][n];
      if (branch_state == Grammar_state_enum::Start)
	return Transducer_start_bgcolor;

      unsigned long col = 0;
      switch (branch_trans.state_type[branch_state])
	{
	case Transducer_state_type_enum::TransducerMatchType:
	  col = Transducer_match_bgcolor;
	  break;
	case Transducer_state_type_enum::TransducerInsertType:
	  col = Transducer_insert_bgcolor;
	  break;
	case Transducer_state_type_enum::TransducerDeleteType:
	  col = Transducer_delete_bgcolor;
	  break;
	case Transducer_state_type_enum::TransducerWaitType:
	  col = Transducer_wait_bgcolor;
	  break;
	default:
	  col = Transducer_undef_bgcolor;
	  break;
	}

      // fade unchanged nodes in collapsed EHMM state
      if (!n_changed())
	col = fade_color (col);

      // return
      return col;
    }  // end branch_color

    // constructor
    State_label (const State_context& context)
      : State_context (context),
	buffer (context.etree.nodes())
    {
      if (!State_context::etree.nodes())
	THROWEXPR ("Empty ETree in EHMM_transducer::State_label");
    }

    // virtual destructor, assignment
    virtual ~State_label() { }
    virtual State_label& operator= (const State_label& l) { return *this; }

    // virtual build
    virtual void build_node() { }

    // label
    const sstring& label()
    {
      for (n = State_context::etree.nodes() - 1; n >= 0; --n)
	{
	  children = State_context::etree.children (n);
	  build_node();
	}
      return buffer[0];
    }
  };

  // various node labels
  struct Node_name : State_label
  {
    Node_name (const State_context& context) : State_label (context) { }
    void build_node()
    {
      State_label::current() << State_label::n_name();
      if (State_label::children.size())
	{
	  State_label::current() << '(';
	  for (unsigned int i = 0; i < State_label::children.size(); ++i)
	    State_label::current() << (i > 0 ? "," : "") << State_label::buffer[State_label::children[i]];   // ACHTUNG, HACK! comma is stripped out in Transducer_SExpr.... ugh
	  State_label::current() << ')';
	}
    }
  };

  // Node records
  struct Node_rec : State_label
  {
    Node_rec (const State_context& context) : State_label (context) { }
    void build_node()
    {
      State_label::current() << "{" << State_label::n_name();
      for (unsigned int i = 0; i < State_label::children.size(); ++i)
	State_label::current() << (i ? "|" : "{") << State_label::buffer[State_label::children[i]];
      State_label::current() << "}";
    }
  };

  // HTML table for dotfiles
  struct Node_table : State_label
  {
    Node_table (const State_context& context) : State_label (context) { }
    void build_node()
    {
      State_label::current().fill ('0');
      const bool absorbs = State_label::n_absorbs();
      const bool emits = State_label::n_emits();

      State_label::current()
	<< "<TABLE BORDER=\"0\" CELLSPACING=\"0\" CELLBORDER=\"0\"><TR>";

      // absorb at root?
      if (State_label::n == 0)
	State_label::current()
	  << "<TD BORDER=\"0\" BGCOLOR=\"#"
	  << setbase(16) << setw(6)
	  << State_label::absorb_color()
	  << setbase(10)
	  << "\">"
	  << (absorbs ? '*' : '-')
	  << "</TD>";

      State_label::current()
	<< "<TD BORDER=\""
	<< (State_label::n_changed() ? 1 : 0)
	<< "\" BGCOLOR=\"#"
	<< setbase(16) << setw(6)
	<< State_label::branch_color()
	<< setbase(10)
	<< "\">"
	<< State_label::n_name()
	<< "</TD>"

	<< "<TD BORDER=\"0\" BGCOLOR=\"#"
	<< setbase(16) << setw(6)
	<< State_label::emit_color()
	<< setbase(10)
	<< "\">"
	<< (emits ? '*' : '-')
	<< "</TD>";

      if (State_label::children.size())
	{
	  State_label::current()
	    << "<TD><TABLE BORDER=\"0\" CELLSPACING=\"0\" CELLBORDER=\"0\">";

	  for_const_contents (vector<ENode>, State_label::children, c)
	    State_label::current()
	    << "<TR><TD>"
	    << State_label::buffer[*c]
	    << "</TD></TR>";

	  State_label::current()
	    << "</TABLE></TD>";
	}
      State_label::current()
	<< "</TR></TABLE>";
    }
  };

  // context structure for transitions
  struct Transition_context
  {
    // data
    State_context src_context, dest_context;
    const CollapsedETrans& etrans;

    // constructors
    Transition_context (const State_context& src_context,
			const State_context& dest_context,
			const CollapsedETrans& etrans)
      : src_context (src_context),
	dest_context (dest_context),
	etrans (etrans)
    { }

    Transition_context (const Transition_context& context)
      : src_context (context.src_context),
	dest_context (context.dest_context),
	etrans (context.etrans)
    { }
  };

  // Transition formulae in LaTeX
  struct Transition_formula : Transition_context
  {
    // data
    vector<sstring> terms_to_multiply;
    vector<sstring> terms_to_sum;

    // constructor
    Transition_formula (const Transition_context& context)
      : Transition_context (context)
    { }

    // label
    sstring label() { return sstring::join (terms_to_multiply, " "); }

    // builder methods
    sstring trans_term (int node, int src, int dest)
    {
      Transition_context::src_context.n = node;
      Transition_context::dest_context.n = node;

      sstring term;
      term << "t^{" << node
	   << "}(" << Transition_context::src_context.n_name()
	   << "," << Transition_context::dest_context.n_name()
	   << ")";

      return term;
    }

    void multiply (int node, int src, int dest)
    {
      terms_to_multiply.push_back (trans_term (node, src, dest));
    }

    void sum (int node, int start, int mid, int end)
    {
      sstring term;
      term << trans_term (node, start, mid) << ' ' << trans_term (node, mid, end);
      terms_to_sum.push_back (term);
    }

    void clear_sum()
    {
      if (terms_to_sum.size() > 0)
	{
	  sstring term;
	  if (terms_to_sum.size() > 1)
	    term << '(';
	  term << sstring::join (terms_to_sum, "+");
	  if (terms_to_sum.size() > 1)
	    term << ')';
	  terms_to_multiply.push_back (term);
	  terms_to_sum.clear();
	}
    }

  };

  // helper class for eliminating states (by simply dropping them) from an EHMM_transducer object
  struct Trimmed_EHMM_transducer : EHMM_transducer<T,PairTrans>
  {
    // constructors
    Trimmed_EHMM_transducer (const EHMM_transducer<T,PairTrans>& ehmm_trans,
			     const vector<int>& states_to_eliminate)
      : EHMM_transducer<T,PairTrans> (ehmm_trans.states() - (int) states_to_eliminate.size())
    {
      // figure out states to keep
      vector<int> states_to_keep;
      const set<int> elim_set (states_to_eliminate.begin(), states_to_eliminate.end());
      for (int s = 0; s < ehmm_trans.states(); ++s)
	if (elim_set.find (s) == elim_set.end())
	  states_to_keep.push_back (s);
      if ((int) states_to_keep.size() != EHMM_transducer<T,PairTrans>::states())
	THROWEXPR ("Oops, miscalculated number of states to keep");

      // copy over elements from the original EHMM, one by one
      EHMM_transducer<T,PairTrans>::start_to_end() = ehmm_trans.start_to_end();
      for (int i = Grammar_state_enum::End; i <= Grammar_state_enum::Start; ++i)
	EHMM_transducer<T,PairTrans>::state2estate[i] = ((EHMM_transducer<T,PairTrans>&)ehmm_trans).state2estate[i];
      for (int i = 0; i < (int) states_to_keep.size(); ++i)
	{
	  const int ei = states_to_keep[i];

	  EHMM_transducer<T,PairTrans>::state_type[i] = ehmm_trans.state_type[ei];
	  EHMM_transducer<T,PairTrans>::state_name[i] = ehmm_trans.state_name[ei];
	  EHMM_transducer<T,PairTrans>::node_record[i] = ((map<int,sstring>&) ehmm_trans.node_record)[ei];  // cast away const
	  EHMM_transducer<T,PairTrans>::node_html[i] = ((map<int,sstring>&) ehmm_trans.node_html)[ei];  // cast away const

	  EHMM_transducer<T,PairTrans>::branch_trans_states[i] = ehmm_trans.branch_trans_states[ei];
	  EHMM_transducer<T,PairTrans>::emission[i] = ehmm_trans.emission[ei];
	  EHMM_transducer<T,PairTrans>::absorption[i] = ehmm_trans.absorption[ei];
	  state2estate[i] = ((EHMM_transducer<T,PairTrans>&)ehmm_trans).state2estate[ei];

	  EHMM_transducer<T,PairTrans>::start[i] = ehmm_trans.start[ei];
	  EHMM_transducer<T,PairTrans>::end[i] = ehmm_trans.end[ei];
	  for (int j = 0; j < (int) states_to_keep.size(); ++j)
	    {
	      const int ej = states_to_keep[j];
	      EHMM_transducer<T,PairTrans>::transition (i, j) = ehmm_trans.transition (ei, ej);
	    }
	}

      EHMM_transducer<T,PairTrans>::start_name = ehmm_trans.start_name;
      EHMM_transducer<T,PairTrans>::end_name = ehmm_trans.end_name;

      EHMM_transducer<T,PairTrans>::etree = ehmm_trans.etree;
      EHMM_transducer<T,PairTrans>::branch_transducer = ehmm_trans.branch_transducer;
      EHMM_transducer<T,PairTrans>::alphabet_size = ehmm_trans.alphabet_size;

      EHMM_transducer<T,PairTrans>::tape_name = ehmm_trans.tape_name;

      // print some log messages
      if (states_to_eliminate.size()
	  && CTAGGING(5,TRANSDUCER TRANSDUCER_ELIM TRANSDUCER_ELIM_STATES))
	{
	  CL << "Silently dropping the following states:\n";
	  for_const_contents (vector<int>, states_to_eliminate, s)
	    CL << ' ' << ehmm_trans.get_state_name(*s) << '\n';
	}
    }
  };

  // build() --- build method for transducer composition
  // builds a multi-sequence transducer from pair transducers,
  // using the conditional[*] collapsed EHMM.
  // entry "branch_transducer[n]" is the conditionally-normalised[*]
  // pair transducer that transforms node "etree.parent[n]" into node "n".
  // [*] if branch_transducer[0] is jointly normalised then joint EHMM,
  // instead of conditional EHMM, will be built.
  void build()
  {
    const bool is_cond = branch_transducer[0].is_conditional();
    if (CTAGGING(3,TRANSDUCER))
      CL << "Building " << (is_cond ? "conditional" : "joint") << " EHMM for ETree " << etree.etree2string(true) << " with " << etree.nodes() << " nodes\n";

    if (CTAGGING(-1,TRANSDUCER))
      CL << "ETree parent[] vector: (" << etree.parent << ")\n";

    if (!etree.valid())
      THROWEXPR ("Invalid ETree in EHMM_transducer::build");

    if (CTAGGING(2,TRANSDUCER BRANCH_TRANSDUCERS))
      {
	CL << "Displaying branch transducers for ETree " << etree.etree2string(true) << ":\n";
	for (int n = 0; n < (int) branch_transducer.size(); ++n)
	  {
	    CL << "Branch from node " << etree.parent[n] << " to node " << n << ":\n";
	    branch_transducer[n].show (CL);
	  }
      }

    // make EStateList and CollapsedEMatrix
    const TSpace tspace (1);  // single-character transducer state space
    CollapsedESpaceSubsets subsets (tspace, etree);

    const EStateList estatelist =
      is_cond
      ? subsets.conditional (tspace)
      : subsets.joint (tspace);
    const CollapsedEMatrix ematrix (estatelist, tspace, etree);

    // States_by_EHMM_type is a map from estate2string() string identifiers,
    // to vectors of actual composite-transducer states having that type.
    // We make this map now.
    States_by_EHMM_type states_by_ehmm_type;
    for_const_contents (EStateList, estatelist, estate)
      {
	// check for TransStart and TransEnd EHMM state types; handle these separately
	if (estate->tstate[0] == TransStart || estate->tstate[0] == TransEnd)
	  {
	    // it's not all-TransStart or all-TransEnd if we find a dissenter, comrades
	    bool found_dissenter = false;
	    for (int n = 1; n < etree.nodes(); ++n)
	      if (estate->tstate[n] != estate->tstate[0])
		{
		  found_dissenter = true;
		  break;
		}
	    if (!found_dissenter)
	      {
		// No dissenters: all TransStart or TransEnd.
		// Initialize the list of states of this type to be Start or End, respectively.
		const int next_state_index =
		  estate->tstate[0] == TransStart
		  ? Grammar_state_enum::Start
		  : Grammar_state_enum::End;
		states_by_ehmm_type[estate->estate2string()]
		  = vector<int> (1, next_state_index);
		state2estate[next_state_index] = *estate;
		continue;
	      }
	  }

	// proceed to find all matching EHMM states
	const vector<vector<int> > estate_branch_trans_states
	  = get_branch_trans_states (branch_transducer, estate->tstate);  // call helper method

	vector<int> state_indices;  // make list of states
	const TEmission temission = estate->emission (tspace);
	const TTerm tabsorption = estate->absorption (tspace);
	for_const_contents (vector<vector<int> >,
			    estate_branch_trans_states,
			    ebts)
	  {
	    const int next_state_index = branch_trans_states.size();
	    state_indices.push_back (next_state_index);
	    branch_trans_states.push_back (*ebts);
	    emission.push_back (temission);
	    absorption.push_back (tabsorption);
	    state2estate[next_state_index] = *estate;
	  }

	// record
	states_by_ehmm_type[estate->estate2string()] = state_indices;
      }

    // log states_by_ehmm_type and branch_trans_states
    if (CTAGGING(2,TRANSDUCER))
      {
	for_const_contents (States_by_EHMM_type, states_by_ehmm_type, type_states)
	  CL << "States of EHMM type " << type_states->first << ": " << type_states->second << "\n";
	for (int s = 0; s < (int) branch_trans_states.size(); ++s)
	  {
	    CL << "State " << s << ": branch transducer states (" << branch_trans_states[s] << "), emission";
	    for_const_contents (vector<TTerm>, emission[s].tterm, tt)
	      if (*tt == TSpaceEnum::TTermNull)
		CL << " -";
	      else
		CL << ' ' << *tt;
	    CL << "\n";
	  }
      }

    // now we have a list of all states, so resize
    entry_type zero_entry;
    ((Transducer_score_manipulator<T>*) this)->set_to_zero (zero_entry);
    transducer_type::resize (branch_trans_states.size(), zero_entry);

    // make state types, names, dotfile labels
    for (int s = Grammar_state_enum::End; s < transducer_type::states(); ++s)
      {
	// state context for names & dotfile labels
	State_context state_context (s, etree, branch_transducer, branch_trans_states, emission, state2estate[s]);

	// name & type
	const sstring nn = Node_name (state_context).label();
	Transducer_state_type t = 0;

	// extra stuff for non-Start-or-End states
	if (s >= 0)
	  {
	    if (absorption[s] != TTermNull)
	      t |= 1;
	    for (int n = 0; n < etree.nodes(); ++n)
	      if (emission[s].tterm[n] != TTermNull)
		t |= 1 << (n + 1);

	    transducer_type::state_type[s] = t;
	    transducer_type::state_name[s] = nn;
	  }

	// graphviz dotfile labels
	transducer_type::node_record[s] = Node_rec (state_context).label();
	transducer_type::node_html[s] = Node_table (state_context).label();

	// log message
	CTAG(-1,TRANSDUCER) << "State " << s << "(" << nn << ") has type " << t << "\n";
      }

    // start, end names
    State_context start_state_context (Grammar_state_enum::Start, etree, branch_transducer, branch_trans_states, emission, state2estate[Grammar_state_enum::Start]);
    transducer_type::start_name = Node_name (start_state_context).label();

    State_context end_state_context (Grammar_state_enum::End, etree, branch_transducer, branch_trans_states, emission, state2estate[Grammar_state_enum::End]);
    transducer_type::end_name = Node_name (end_state_context).label();

    // for each branch transducer, make a list of wait states
    vector<vector<int> > wait_states (etree.nodes());
    for (int n = 0; n < etree.nodes(); ++n)
      {
	for (int s = 0; s < branch_transducer[n].states(); ++s)
	  if (branch_transducer[n].state_type[s] == 0)
	    wait_states[n].push_back (s);
      }

    // make transitions
    for_const_contents (vector<CollapsedETrans>, ematrix.etrans, etrans)
      {
	// loop through all (src,dest) pairs with the correct types
	const sstring src_type_string = etrans->src().estate2string();
	const sstring dest_type_string = etrans->dest().estate2string();

	if (CTAGGING(-1,TRANSDUCER))
	  CL << "Considering all transitions from type '" << src_type_string << "' to type '" << dest_type_string << "'\n";

	// loop over src composite state
	for_const_contents (vector<int>, states_by_ehmm_type[src_type_string], src)
	  {
	    // create src_context
	    State_context src_context (*src, etree, branch_transducer, branch_trans_states, emission, state2estate[*src]);
	    // loop over dest composite state
	    for_const_contents (vector<int>,
				states_by_ehmm_type[dest_type_string],
				dest)
	      {
		// log
		if (CTAGGING(-2,TRANSDUCER))
		  CL << "Considering transition from " << transducer_type::get_state_name(*src) << " to " << transducer_type::get_state_name(*dest) << "\n";

		// create dest_context
		State_context dest_context (*dest, etree, branch_transducer, branch_trans_states, emission, state2estate[*dest]);

		// create trans_context
		// commented out because currently unused
		// Transition_context trans_context (src_context, dest_context, *etrans);

		// initialise trans_prob
		entry_type trans_prob;
		((Transducer_score_manipulator<T>*) this)->set_to_one (trans_prob);

		// The amazing tree-winding loop!
		// Loop over nodes, multiplying branch probabilities to get transition probabilities.
		for (int n = 0; n < etree.nodes(); ++n)
		  {
		    // retrieve the EHMM representation of the path: the list of state types.
		    const TPath& tpath = etrans->tpath[n];

		    // get start & endpoints of the path, in actual branch-transducer state space
		    const int path_start = *src == Grammar_state_enum::Start ? *src : branch_trans_states[*src][n];
		    const int path_end = *dest == Grammar_state_enum::End ? *dest : branch_trans_states[*dest][n];

		    // print contemplative log message
		    CTAG(-3,TRANSDUCER) << "Node " << n << ": considering transition from " << path_start << " to " << path_end << "\n";

		    // initialize branch_prob
		    entry_type branch_prob;
		    ((Transducer_score_manipulator<T>*) this)->set_to_one (branch_prob);

		    // what kind of transition path is this?
		    switch (tpath.tstate.size())
		      {
		      case 1:
			// empty path: no transition allowed on this branch
			if (path_start == path_end)
			  CL << "(no change in state, so no contribution to branch prob)\n";
			else
			  ((Transducer_score_manipulator<T>*) this)->set_to_zero (branch_prob);
			break;

		      case 2:
			// direct transition
			branch_prob = branch_transducer[n].transition (path_start, path_end);
			CL << "Branch prob is " << ((Transducer_score_manipulator<T>*) this)->pval_string(branch_prob) << " (direct)\n";
			break;

		      case 3:
			{
			  // indirect transition: sum over intermediate null states
			  // The intermediate state must be a wait state.
			  if (tpath.tstate[1] != TransWait)
			    THROWEXPR ("In EHMM: two-step transition via non-wait state");

			  // clear branch_prob
			  ((Transducer_score_manipulator<T>*) this)->set_to_zero (branch_prob);

			  // loop over intermediate states, summing branch_prob
			  for_const_contents (vector<int>, wait_states[n], path_mid)
			  {
			    const entry_type& start_to_mid = branch_transducer[n].transition (path_start, *path_mid);
			    const entry_type& mid_to_end = branch_transducer[n].transition (*path_mid, path_end);

			    if (((Transition_matrix<T>*) this)->is_non_null (start_to_mid) && ((Transition_matrix<T>*) this)->is_non_null (mid_to_end))
			      {
				entry_type via_prob (start_to_mid);
				((Transducer_score_manipulator<T>*) this)->pmul_acc (via_prob, mid_to_end);
				((Transducer_score_manipulator<T>*) this)->psum_acc (branch_prob, via_prob);

				CL << "Added " << ((Transducer_score_manipulator<T>*) this)->pval_string(via_prob) << " to branch prob (transition via wait state " << *path_mid << ")\n";
			      }
			  }
			  
			  CL << "Branch prob is " << ((Transducer_score_manipulator<T>*) this)->pval_string(branch_prob) << " (via wait states)\n";
			  break;
			}

		      default:
			// don't know this type of path
			THROWEXPR ("Bad path type in EHMM transition");
			break;
		      }

		    // multiply transition prob by prob for branch
		    ((Transducer_score_manipulator<T>*) this)->pmul_acc (trans_prob, branch_prob);
		  }

		// store composite transducer transition probability
		if (((Transducer_score_manipulator<T>*) this)->is_nonzero (trans_prob))
		  transducer_type::transition (*src, *dest) = trans_prob;

		if (CTAGGING(-2,TRANSDUCER))
		  CL << "Transition from " << transducer_type::get_state_name(*src) << " to " << transducer_type::get_state_name(*dest) << " has prob " << ((Transducer_score_manipulator<T>*) this)->pval_string(trans_prob) << "\n";
	      }
	  }
      }

    // drop inaccessible states
    Trimmed_EHMM_transducer trimmed_ehmm (*this, inaccessible_states());
    swap (*this, (EHMM_transducer<T,PairTrans>&) trimmed_ehmm);
  }

  // try to find inaccessible states
  vector<int> inaccessible_states()
  {
    // figure out which states are inaccessible
    // first look for paths from Start
    vector<int> start_accessible (Transducer<T>::states(), (int) 0);
    set<int> visited;
    stack<int> src_states;
    src_states.push (Grammar_state_enum::Start);
    while (src_states.size())
      {
	const int src = src_states.top();
	src_states.pop();

	if (visited.find(src) == visited.end())
	  {
	    visited.insert (src);
	    for (int dest = 0; dest < Transducer<T>::states(); ++dest)
	      if (((Transducer_score_manipulator<T>*) this)->is_nonzero (Transducer<T>::transition(src,dest)))
		{
		  start_accessible[dest] = 1;
		  src_states.push (dest);
		}
	  }
      }

    // now look for paths from End
    vector<int> end_accessible (Transducer<T>::states(), (int) 0);
    visited.clear();
    stack<int> dest_states;
    dest_states.push (Grammar_state_enum::End);
    while (dest_states.size())
      {
	const int dest = dest_states.top();
	dest_states.pop();

	if (visited.find(dest) == visited.end())
	  {
	    visited.insert (dest);
	    for (int src = 0; src < Transducer<T>::states(); ++src)
	      if (((Transducer_score_manipulator<T>*) this)->is_nonzero (Transducer<T>::transition(src,dest)))
		{
		  end_accessible[src] = 1;
		  dest_states.push (src);
		}
	  }
      }

    // make vector of inaccessible states
    vector<int> inacc_states;
    for (int s = 0; s < Transducer<T>::states(); ++s)
      if (!start_accessible[s] || !end_accessible[s])
	inacc_states.push_back (s);

    // return
    return inacc_states;
  }

  // helper method to figure out whether a particular transition involves a change of state at a particular node
  // return codes:
  //   0 for no transition (branch transducer stays in the same state)
  //   1 for single transition
  //   2 for double transition via wait state
  //  -1 if composite transition is illegal
  int count_branch_transitions (int src, int dest, int node) const
  {
    const ESubtrees esubtrees (etree);
    const TSpace dummy_tspace (1);

    EHMM_transducer* mutable_this = (EHMM_transducer*) this;  // cast away const
    const CollapsedETrans etrans (mutable_this->state2estate[src], mutable_this->state2estate[dest], dummy_tspace, etree, esubtrees);
    if (!etrans.tpath.size())  // tpath empty => invalid transition
      return -1;

    // retrieve the EHMM representation of the path: the list of state types.
    const TPath& tpath = etrans.tpath[node];

    // get start & endpoints of the path, in actual branch-transducer state space
    const int path_start = src == Grammar_state_enum::Start ? src : branch_trans_states[src][node];
    const int path_end = dest == Grammar_state_enum::End ? dest : branch_trans_states[dest][node];

    // print contemplative log message
    CTAG(-3,TRANSDUCER) << "Checking for branch transition at node " << node
			<< " from " << Transducer<T>::get_state_name(src)
			<< " to " << Transducer<T>::get_state_name(dest)
			<< "\n";

    // what kind of transition path is this?
    int trans_at_node = -1;
    switch (tpath.tstate.size())
      {
      case 1:
	// empty path: no transition allowed on this branch
	if (path_start == path_end)
	  trans_at_node = 0;
	else
	  trans_at_node = -1;
	break;

      case 2:
	// direct transition
	trans_at_node = 1;
	break;

      case 3:
	{
	  // indirect transition: sum over intermediate null states
	  // The intermediate state must be a wait state.
	  if (tpath.tstate[1] != TransWait)
	    trans_at_node = -1;  // two-step transition via non-wait state
	  else
	    trans_at_node = 2;
	  break;
	}

      default:
	// don't know this type of path
	trans_at_node =  -1;  // bad path type
	break;
      }

    // return
    return trans_at_node;
  }


  // constructors: these call build()
  EHMM_transducer (const ETree& etree,
		   const vector<pair_transducer_type>& branch_transducer)
    : Transducer<T> (0),
      etree (etree),
      branch_transducer (branch_transducer),
      alphabet_size (branch_transducer.front().alphabet_size)
  {
    build();
  }

  EHMM_transducer (const ETree& etree,
		   const pair_transducer_type& pair_trans)
    : Transducer<T> (0),
      etree (etree),
      branch_transducer (etree.nodes(), pair_trans),
      alphabet_size (pair_trans.alphabet_size)
  {
    build();
  }

  EHMM_transducer (int states)
    : Transducer<T> (states),
      alphabet_size (0),
      branch_trans_states (states),
      emission (states, TEmission (0)),
      absorption (states)
  { }

  // helper: returns all combinations of branch transducer states
  // matching a particular set of state types
  static vector<vector<int> > get_branch_trans_states
  (const vector<pair_transducer_type>& branch_transducer,
   const vector<TState>& tstate)
  {
    if (CTAGGING(0,TRANSDUCER))
      {
	CL << "get_branch_trans_states: tstate = ";
	for_const_contents (vector<TState>, tstate, ts)
	  CL << TSpaceEnum::tstate2char (*ts);
	CL << "\n";
      }

    vector<vector<int> > valid_transducer_states_by_node ((int) tstate.size());
    vector<int> count_begin (tstate.size()), count_end (count_begin);
    for (int n = 0; n < (int) tstate.size(); ++n)
      {
	if (tstate[n] == TransStart)
	  {
	    valid_transducer_states_by_node[n].push_back
	      (Grammar_state_enum::Start);
	    count_end[n] = 1;
	  }
	else if (tstate[n] == TransEnd)
	  {
	    valid_transducer_states_by_node[n].push_back
	      (Grammar_state_enum::End);
	    count_end[n] = 1;
	  }
	else
	  for (int s = 0; s < (int) branch_transducer[n].states(); ++s)
	    if (branch_transducer[n].state_type[s] == (Transducer_state_type_enum::State_type) TSpaceEnum::xy_emit (tstate[n]))
	      {
		valid_transducer_states_by_node[n].push_back (s);
		++count_end[n];
	      }

	if (CTAGGING(0,TRANSDUCER))
	    CL << "get_branch_trans_states, node n=" << n << ": tstate[n]=" << TSpaceEnum::tstate2char(tstate[n]) << ", valid_transducer_states_by_node[n]=(" << valid_transducer_states_by_node[n] << ")\n";
      }

    bool matching_states_exist = true;
    for_const_contents (vector<vector<int> >, valid_transducer_states_by_node, vts)
      if (vts->size() == 0)
	{
	  matching_states_exist = false;
	  break;
	}

    vector<vector<int> > branch_trans_states;
    if (matching_states_exist)
      {
	vector<int> transducer_state ((int) tstate.size());
	for (Counter counter (count_begin, count_begin, count_end);
	     counter < count_end; ++counter)
	  {
	    for (int n = 0; n < (int) tstate.size(); ++n)
	      transducer_state[n] = valid_transducer_states_by_node[n][counter[n]];
	    branch_trans_states.push_back (transducer_state);
	  }
      }

    return branch_trans_states;
  }

  // helper to find emission subtrees
  void summarize_emission (int state,
			   ENode& inserter,
			   vector<ENode>& deleters,
			   vector<ENode>& absorbers) const
  {
    emission[state].summarize (etree, inserter, deleters, absorbers);
  }


  // show_sexpr helpers
  const char* transducer_keyword() const { return TSEXPR_COMPOSITE; }

  void get_state_type_sexpr (int state, sstring& sts) const
  {
    switch (state)
      {
      case Grammar_state_enum::Start:
	sts << TSEXPR_START;
	break;
      case Grammar_state_enum::End:
	sts << TSEXPR_END;
	break;

      default:
	switch (Transducer<T>::state_type[state])
	  {
	  case Transducer_state_type_enum::TransducerWaitType:
	    sts << TSEXPR_WAIT;
	    break;
	  default:
	    Transducer<T>::get_state_type_sexpr (state, sts);
	    break;
	  }
	break;
      }
  }

  void get_emit_label_sexpr (int state, sstring& label_sexpr) const
  {
    if (Transducer<T>::state_type[state])
      {
	if (has_emit_labels())
	  {
	    // hacky: check for emit_root == -1
	    int emit_root = emission[state].emitter();
	    if (emit_root == 0)
	      {
		const Transducer_state_type root_type = branch_transducer[0].state_type[branch_trans_states[state][0]];
		if (root_type == Transducer_state_type_enum::TransducerMatchType
		    || root_type == Transducer_state_type_enum::TransducerDeleteType)  // absorb anything from the root?
		  emit_root = -1;
	      }
	    label_sexpr = emit_label (state, emit_root);
	  }
      }
  }

  sstring transition_label_sexpr (const T& t, const PScores* ps) const
  {
    sstring ls;
    sstring sls = ((Transducer_score_manipulator<T>*) this)->score_label_string (t, ps);
    if (sls.size())
      ls << '(' << sls << ')';
    return ls;
  }

  sstring emit_label (int state, ENode node) const
  {
    bool absorbs = false, emits = false;
    sstring suffix, label;
    bool need_paren = false;   // true if we need parentheses in the TSEXPR_SUM expression
    if (node < 0)
      {
	suffix << emit_label (state, 0);
	emits = true;
      }
    else
      {
	const pair_transducer_type& trans = branch_transducer[node];
	const int trans_state_index = branch_trans_states[state][node];
	const sstring& trans_state_emit_label = trans.emit_label[trans_state_index];
	const Transducer_state_type trans_state_type = trans.state_type[trans_state_index];

	absorbs = trans_state_type & Transducer_state_type_enum::TransducerDeleteType;
	emits = trans_state_type & Transducer_state_type_enum::TransducerInsertType;

	if (trans_state_emit_label.size() && (absorbs || emits))
	  {
	    suffix << '(' << trans_state_emit_label;
	    if (absorbs)
	      suffix << ' ' << Transducer<T>::get_tape_name (etree.parent[node]);
	    if (emits)
	      suffix << ' ' << Transducer<T>::get_tape_name (node);
	    suffix << ')';
	  }

	const vector<ENode> kids = etree.children (node);
	for_const_contents (vector<ENode>, kids, c)
	  {
	    const int cs = branch_trans_states[state][*c];
	    if (cs >= 0)
	      if (branch_transducer[*c].state_type[cs] != Transducer_state_type_enum::TransducerWaitType)
		{
		  const sstring child_label = emit_label (state, *c);
		  if (child_label.size())
		    {
		      if (suffix.size())
			{
			  suffix << " * ";
			  need_paren = true;
			}
		      suffix << child_label;
		    }
		}
	  }
      }
    if (suffix.size())
      {
	if (emits)
	  label << '(' << TSEXPR_SUM
		<< ' ' << Transducer<T>::get_tape_name (node) << ' '
		<< (need_paren ? "(" : "")
		<< suffix
		<< (need_paren ? ")" : "")
		<< ')';
	else
	  ((basic_string<char>&)label).swap (suffix);
      }
    return label;
  }

  bool has_emit_labels() const
  {
    for (int n = 0; n < (int) branch_transducer.size(); ++n)
      for (int s = 0; s < branch_transducer[n].states(); ++s)
	if (branch_transducer[n].state_type[s] != 0 && branch_transducer[n].emit_label[s].size() > 0)
	  return true;
    return false;
  }

  // extremely hacky overriding helper method that strips out commas from EHMM state names.... euch
  sstring sexpr_state_name (int state) const
  {
    sstring s;
    s << '(' << Transducer<T>::get_state_name (state) << ')';
    for (int i = 0; i < (int) s.size(); ++i)
      if (s[i] == ',')
	s[i] = ' ';
    return s;
  }

  // state path conversion
  vector<vector<int> > branch_paths (const vector<int>& composite_path) const
  {
    vector<vector<int> > branch_paths (etree.nodes());

    if (composite_path.size())
      for (int n = etree.nodes() - 1; n >= 0; --n)
	{
	  vector<int> branch_path;
	  for (int pos = 0; pos < (int) composite_path.size(); ++pos)
	    if (composite_path[pos] >= 0)
	      {
		const int src = pos == 0 ? Grammar_state_enum::Start : composite_path[pos-1];
		const int dest = composite_path[pos];
		const int n_branch_trans = count_branch_transitions (src, dest, n);
		if (n_branch_trans < 0)
		  THROWEXPR ("Illegal transition from " << sexpr_state_name(src) << " to " << sexpr_state_name(dest));
		if (n_branch_trans > 0
		    && branch_transducer[n].state_type[branch_trans_states[composite_path[pos]][n]]
		    != Transducer_state_type_enum::TransducerWaitType)      // exclude wait states
		  branch_path.push_back (branch_trans_states[composite_path[pos]][n]);
	      }
	  swap (branch_paths[n], branch_path);
	}

    return branch_paths;
  }

  vector<int> composite_path (const vector<vector<int> >& branch_paths) const
  {
    // invert branch_trans_states
    map<vector<int>,int> composite_state;
    for (int s = 0; s < Transducer<T>::states(); ++s)
      composite_state[branch_trans_states[s]] = s;
    composite_state[vector<int> (etree.nodes(), Grammar_state_enum::Start)] = Grammar_state_enum::Start;

    // create composite path
    vector<int> path_cursor (etree.nodes(), 0);
    vector<int> branch_state (etree.nodes(), (int) Grammar_state_enum::Start);
    vector<int> composite_path;
    while (1)
      {
	// record current composite state
	vector<sstring> branch_state_name;
	for (int n = 0; n < (int) branch_state.size(); ++n)
	  branch_state_name.push_back (branch_transducer[n].get_state_name (branch_state[n]));
	if (composite_state.find (branch_state) == composite_state.end())
	  THROWEXPR ("Can't find composite state (" << branch_state_name << ')');

	composite_path.push_back (composite_state[branch_state]);
	CTAG(2,TRANSDUCER TRANSDUCER_COMPOSITE_PATH) << "Added state (" << branch_state_name << ") to path\n";

	// find next branch transducer to advance: highest numbered insert state
	ENode inserter = -1;
	bool found_ins = false;
	for (int n = etree.nodes() - 1; n >= 0; --n)
	  {
	    const pair_transducer_type& btn = branch_transducer[n];
	    int next_state;
	    if (path_cursor[n] < (int) branch_paths[n].size())
	      {
		next_state = branch_paths[n][path_cursor[n]];
		if (btn.state_type[next_state] == Transducer_state_type_enum::TransducerInsertType)
		  {
		    inserter = n;
		    found_ins = true;
		    break;
		  }
	      }
	    else
	      next_state = Grammar_state_enum::End;
	    // no insert, so check that node is in a wait state
	    const int last_state = branch_state[n];
	    if (last_state < 0 ? true : btn.state_type[last_state] != Transducer_state_type_enum::TransducerWaitType)
	      {
		// place node in an appropriate wait state
		vector<int> wait_states, allowed_wait_states;
		for (int s = 0; s < btn.states(); ++s)
		  if (btn.state_type[s] == Transducer_state_type_enum::TransducerWaitType)
		    {
		      wait_states.push_back (s);
		      if (((Transducer_score_manipulator<T>*) this)->is_nonzero (btn.transition (last_state, s)) && ((Transducer_score_manipulator<T>*) this)->is_nonzero (btn.transition (s, next_state)))
			allowed_wait_states.push_back (s);
		    }
		if (allowed_wait_states.size())
		  branch_state[n] = allowed_wait_states[0];
		else if (wait_states.size())
		  {
		    CLOGERR << "Can't find a wait state between " << btn.state_name[last_state] << " & " << btn.state_name[next_state] << "; using " << btn.state_name[wait_states[0]] << "\n";
		    branch_state[n] = wait_states[0];
		  }
		else
		  THROWEXPR ("Can't find any wait states in transducer; freaking out");
	      }
	  }
	// if we didn't find any inserters, just head to the End state
	if (!found_ins)
	  break;

	// advance branches, propagating emissions down the tree
	stack<int> branches_to_advance;
	branches_to_advance.push (inserter);
	while (!branches_to_advance.empty())
	  {
	    const int n = branches_to_advance.top();
	    if (path_cursor[n] >= (int) branch_paths[n].size())
	      THROWEXPR ("Specified state path for branch " << n << " is inconsistent (too short?)");

	    const int s = branch_paths[n][path_cursor[n]];
	    const int t = branch_transducer[n].state_type[s];

	    branch_state[n] = s;
	    branches_to_advance.pop();
	    ++path_cursor[n];

	    // is there an emission at this node?
	    if (t == Transducer_state_type_enum::TransducerInsertType || t == Transducer_state_type_enum::TransducerMatchType)
	      {
		// put kids on "advance branch" stack
		const vector<ENode> kids = etree.children (n);
		for_const_reverse_contents (vector<ENode>, kids, c)
		  branches_to_advance.push (*c);
	      }
	  }
      }

    // finish
    composite_path.push_back (Grammar_state_enum::End);
    return composite_path;
  }

};

struct Pair_transducer_scores : Pair_transducer<Score>
{
  // constructors
  Pair_transducer_scores (int states)
    : Pair_transducer<Score> (states, -InfinityScore)
  { }

  // helpers
  // calculate the score for a branch
  // (hopefully with efficiency like handel/tkfdata.cc, i.e. O(L),
  //  rather than hmm/pairhmm.cc, i.e. O(L^2))
  Score pairwise_path_score (const Pairwise_path& path);

  // calculate the effective transition score between two non-wait states (summed over intermediate wait states)
  Score effective_trans_score (int src, int dest) const;

  // sample a child path, conditional on a parent path
  // max_tries is the maximum number of attempts before emitting an empty child sequence (for transducers which can get stuck)
  void sample (const Digitized_biosequence& parent_dsq, vector<int>& parent_child_path, Digitized_biosequence& child_dsq, int max_tries = TRANSDUCER_SAMPLE_MAX_TRIES);

  // convert to Pair_HMM_scores
  Pair_HMM_scores pair_hmm (const Alphabet& alph) const;

  // string stuff
  void show_element (const Score& element, ostream& o) const
  { ShowScore (element, o); }
};

// typedefs
typedef Transducer<Score> Transducer_scores;
typedef Transducer<PFunc> Transducer_funcs;
typedef Transducer<Prob> Transducer_counts;

// methods for working with transducers
struct Transducer_methods : Grammar_state_enum, Transducer_state_type_enum
{
  // populate corresponding Transducer<Score> object
  static void populate_transducer_scores (const Transducer<PFunc>& trans_func, Transducer<Score>& trans_sc, const PScores& pscore)
  {
    trans_sc.state_type = trans_func.state_type;
    trans_sc.state_name = trans_func.state_name;
    trans_sc.start_name = trans_func.start_name;
    trans_sc.end_name = trans_func.end_name;
    trans_sc.node_html = trans_func.node_html;
    trans_sc.node_record = trans_func.node_record;
    trans_sc.edge_label = trans_func.edge_label;
    trans_sc.tape_name = trans_func.tape_name;

    // transitions
    for (int src = Start; src < trans_func.states(); ++src)
      for (int dest = End; dest < trans_func.states(); ++dest)
	if (dest != Start)
	  {
	    const PFunc f = trans_func.transition (src, dest);
	    trans_sc.transition (src, dest) = f.is_null() ? -InfinityScore : f.eval_sc (pscore);
	  }
  }

  // graphviz dotfiles (PFunc)
  // virtual method returning dotfile edge attributes
  template<class T>
  static map<sstring,sstring> static_dotfile_edge_attrs (const Transducer<T>& trans_funcs, int src, int dest, const PScores* pscores = 0, bool no_group_prefix = false)
  {
    // there used to be more here; currently this function just converts the transition score to a probabilistic weight
    // we probably should just override the edge_node_label to achieve this, but this way allows for further decoration of edges

    map<sstring,sstring> attr;

    // edge label
    const PFunc& f = trans_funcs.transition (src, dest);
    if (!f.is_one())
      {
	sstring& s = attr[sstring("label")];
	s << '"';
	if (pscores)
	  f.show (s, &pscores->group_suffix, no_group_prefix);
	else
	  f.show (s);
	s << '"';
      }

    // font
    attr[sstring("fontname")] = Transducer_fontname;

    // return
    return attr;
  }

  // virtual method returning dotfile node attributes
  template<class T>
  static map<sstring,sstring> static_dotfile_node_attrs (const Transducer<T>& trans, int state)
  {
    map<sstring,sstring> attr;

    // label
    attr[sstring("label")] = trans.dotfile_node_label(state);

    // font
    attr[sstring("fontname")] = Transducer_fontname;

    // return
    return attr;
  }

  // graphviz dotfiles (Score)
  // virtual method returning dotfile edge attributes
  static map<sstring,sstring> dotfile_edge_attrs (const Transducer<Score>& trans_sc, int src, int dest)
  {
    // there used to be more here; currently this function just converts the transition score to a probabilistic weight
    // we probably should just override the edge_node_label to achieve this, but this way allows for further decoration of edges

    map<sstring,sstring> attr;

    // edge label
    const Score sc = trans_sc.transition (src, dest);
    if (sc != 0)
      {
	const Prob p = Score2Prob (sc);
	attr[sstring("label")] << '"' << p << '"';
      }

    // font
    attr[sstring("fontname")] = Transducer_fontname;

    // return
    return attr;
  }

  // methods for updating PCounts via PFunc's
  // emissions and transitions are handled separately:
  //  emit counts are associated with emit labels of branch transducers
  //  transition counts are associated with composite EHMM transitions
  static void inc_transition_counts (const Transducer<Prob>& transition_counts,
				     const Transducer<PFunc>& transition_funcs,
				     PCounts& var_counts,
				     const PScores& var_scores,
				     const Prob weight = 1.);

  static void inc_emit_label_counts (const Pair_transducer<Prob>& transducer_counts,
				     const Pair_transducer<PFunc>& transducer_funcs,
				     PCounts& var_counts,
				     const PScores& var_scores,
				     const Prob weight = 1.);

};

struct Pair_transducer_funcs : Pair_transducer<PFunc>
{
  // constructors
  Pair_transducer_funcs (int states = 0)
    : Pair_transducer<PFunc> (states, PFunc())
  { }

  // the following form of the constructor "converts" a Pair_transducer_scores into a fully parametric Pair_transducer_funcs
  // by automatically creating & assigning appropriate scalar & alphabet-subscripted variables in pscores
  Pair_transducer_funcs (const Pair_transducer_scores& pair_trans_sc, PScores& pscores, const vector<sstring>& alphabet, const char* pvar_prefix = "");

  // eval
  Pair_transducer_scores eval_sc (const PScores& var_sc) const
  {
    Pair_transducer_scores pair_sc (states());
    Transducer_methods::populate_transducer_scores (*this, pair_sc, var_sc);

    // emissions
    pair_sc.emit_label = emit_label;
    pair_sc.alphabet_size = alphabet_size;
    pair_sc.state_type = state_type;
    pair_sc.alloc_pair_emit ((Score) 0);

    for (int state = 0; state < states(); ++state)
      for (int x = 0; x < pair_emit[state].xsize(); ++x)
	for (int y = 0; y < pair_emit[state].ysize(); ++y)
	  {
	    const PFunc f = pair_emit[state](x,y);
	    pair_sc.pair_emit[state](x,y) = f.is_null() ? -InfinityScore : f.eval_sc (var_sc);
	  }

    return pair_sc;
  }

  // string stuff
  map<sstring,sstring> dotfile_edge_attrs (int src, int dest) const
  { return Transducer_methods::static_dotfile_edge_attrs (*this, src, dest); }

  map<sstring,sstring> dotfile_node_attrs (int state) const
  { return Transducer_methods::static_dotfile_node_attrs (*this, state); }

  void show_element (const PFunc& element, ostream& o) const
  { element.show (o); }
};

struct Pair_transducer_counts : Pair_transducer<Prob>
{
  // constructors
  Pair_transducer_counts (int states = 0)
    : Pair_transducer<Prob> (states, 0.)
  { }

  // templated copy constructor
  template<class T>
  Pair_transducer_counts (const Pair_transducer<T>& pair_trans)
    : Pair_transducer<Prob> (pair_trans, 0.)
  { }

  // string stuff
  map<sstring,sstring> dotfile_edge_attrs (int src, int dest) const
  { return map<sstring,sstring>(); }

  void show_element (const Score& element, ostream& o) const
  { ShowScore (element, o); }
};


struct EHMM_transducer_scores : EHMM_transducer<Score,Pair_transducer_scores>
{
  // constructors
  EHMM_transducer_scores (const ETree& etree,
			  const vector<pair_transducer_type>& branch_transducers)
    : EHMM_transducer<Score,Pair_transducer_scores> (etree, branch_transducers)
  { }
  
  EHMM_transducer_scores (const ETree& etree,
		   const pair_transducer_type& branch_transducer)
    : EHMM_transducer<Score,Pair_transducer_scores> (etree, branch_transducer)
  { }

  EHMM_transducer_scores (int states = 0)
    : EHMM_transducer<Score,Pair_transducer_scores> (states)
  { }

  // string stuff
  map<sstring,sstring> dotfile_edge_attrs (int src, int dest) const
  { return Transducer_methods::dotfile_edge_attrs (*this, src, dest); }

  map<sstring,sstring> dotfile_node_attrs (int state) const
  { return Transducer_methods::static_dotfile_node_attrs (*this, state); }

  void show_element (const Score& element, ostream& o) const
  { ShowScore (element, o); }
};

// class for eliminating null states (by summing them out) from an EHMM_transducer_scores object
struct Eliminated_EHMM_transducer_scores : EHMM_transducer_scores
{
  // data
  const EHMM_transducer_scores* ehmm_trans_sc;  // original
  const vector<int>* states_to_eliminate;
  vector<int> states_to_keep;
  Concrete_transition_probs loopy_probs, loop_exit, elim_probs;
  Concrete_transition_scores elim_scores;

  // constructors
  Eliminated_EHMM_transducer_scores (const EHMM_transducer_scores& ehmm_trans_sc,
				     const vector<int>& states_to_eliminate);

  // transducer name
  const char* transducer_keyword() const { return TSEXPR_ELIMINATED; }

  // convert a state path back into original EHMM state indices, randomlyly adding eliminated states
  // if choose_ML_path==true, the ML path is taken instead
  vector<int> sample_eliminated (const vector<int>& path,
				 bool choose_ML_path = false);

  // convert a transition counts matrix back to counts of transitions in the original EHMM
  Transducer_counts count_eliminated (const Transducer_counts& transition_counts) const;
};

struct EHMM_transducer_funcs : EHMM_transducer<PFunc,Pair_transducer_funcs>
{
  // data
  const PScores* pscores;
  bool no_group_prefix;

  // constructors
  EHMM_transducer_funcs (int states = 0)
    : EHMM_transducer<PFunc,Pair_transducer_funcs> (states),
      pscores (0),
      no_group_prefix (false)
  { }

  EHMM_transducer_funcs (const ETree& etree,
			 const vector<pair_transducer_type>& branch_transducers,
			 PScores* pscores = 0,
			 bool no_group_prefix = false)
    : EHMM_transducer<PFunc,Pair_transducer_funcs> (etree, branch_transducers),
      pscores (pscores),
      no_group_prefix (no_group_prefix)
  { }
  
  EHMM_transducer_funcs (const ETree& etree,
			 const pair_transducer_type& branch_transducer)
    : EHMM_transducer<PFunc,Pair_transducer_funcs> (etree, branch_transducer),
      pscores (0),
      no_group_prefix (false)
  { }

  // eval
  EHMM_transducer_scores eval_sc (const PScores& var_sc) const
  {
    vector<Pair_transducer_scores> branch_sc;
    for_const_contents (vector<pair_transducer_type>, branch_transducer, btrans)
      branch_sc.push_back (btrans->eval_sc (var_sc));
    EHMM_transducer_scores ehmm_sc (etree, branch_sc);
    Transducer_methods::populate_transducer_scores (*this, ehmm_sc, var_sc);
    ehmm_sc.branch_trans_states = branch_trans_states;
    ehmm_sc.emission = emission;
    ehmm_sc.absorption = absorption;
    ehmm_sc.state2estate = state2estate;
    return ehmm_sc;
  }

  // string stuff
  map<sstring,sstring> dotfile_edge_attrs (int src, int dest) const
  { return Transducer_methods::static_dotfile_edge_attrs (*this, src, dest, pscores, no_group_prefix); }

  map<sstring,sstring> dotfile_node_attrs (int state) const
  { return Transducer_methods::static_dotfile_node_attrs (*this, state); }

  void show_element (const PFunc& element, ostream& o) const
  { element.show (o); }
};

#endif /* TRANSDUCER_INCLUDED */
