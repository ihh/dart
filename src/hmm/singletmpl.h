#ifndef SINGLETMPL_INCLUDED
#define SINGLETMPL_INCLUDED

#include "seq/alignment.h"
#include "seq/local.h"
#include "hmm/transmat.h"

// import Grammar_state_enum
typedef Grammar_state_enum HMM_state_enum;

// default names
#define Single_HMM_default_name         "Singlet_HMM"
#define Single_HMM_default_state_prefix "Singlet_HMM_state_"

// State type enumeration for single HMMs.
// Each state in a single HMM can have one of two types: emit or null.
// This base class carries a list of state types,
// plus an optional alphabet (mandatory for Telegraph i/o)
class Single_state_typing
{
protected:
  const Alphabet* _alphabet;

public:
  enum State_type { Null = 0, Emit = 1 };
  vector<State_type> state_type;

  // constructors

  Single_state_typing() : _alphabet(0), state_type() { }
  Single_state_typing (int states, const Alphabet& alph) : _alphabet (&alph), state_type (states, (State_type) Emit) { }
  Single_state_typing (const Single_state_typing& base) : _alphabet (base._alphabet), state_type (base.state_type) { }

  // assignment operator

  Single_state_typing& operator= (const Single_state_typing& base)
  {
    _alphabet = base._alphabet;
    state_type = base.state_type;
    return *this;
  }

  // accessors
  const Alphabet& alphabet() const { if (_alphabet==0) THROWEXPR("Alphabet undefined"); return *_alphabet; }
  inline int states() const { return state_type.size(); }   // excludes start & end states

  // reset method
  void reset_state_type() { state_type = vector<State_type> (states(), (State_type) Emit); }
};

// A template class for a transition matrix, a single-HMM state typing and a set of emit profiles
// Applications of this template include arrays of HMM scores, Baum-Welch counts and PFunc's (symbolic differentiable functions)
template<class T>
class Single_HMM : public Single_state_typing, public Transition_matrix<T>
{
public:
  sstring name;
  vector<sstring> state_name;

  vector<vector<T> > emit;

  // constructors

  Single_HMM() :
    Single_state_typing(),
    Transition_matrix<T>(),
    name (Single_HMM_default_name),
    state_name(),
    emit()
  { }

  Single_HMM (int states, const Alphabet& alphabet) :
    Single_state_typing (states, alphabet),
    Transition_matrix<T> (states),
    name (Single_HMM_default_name),
    state_name (states, sstring (Single_HMM_default_state_prefix)),
    emit (states, vector<T> (alphabet.size()))
  {
    for (int s = 0; s < states; ++s)
      state_name[s] << (s + 1);
  }
  
  Single_HMM (int states, const Alphabet& alphabet, T t) :
    Single_state_typing (states, alphabet),
    Transition_matrix<T> (states, t),
    name (Single_HMM_default_name),
    state_name (states, sstring (Single_HMM_default_state_prefix)),
    emit (states, vector<T> (alphabet.size(), t))
  {
    for (int s = 0; s < states; ++s)
      state_name[s] << (s + 1);
  }

  // copy constructors

  Single_HMM (const Single_HMM<T>& hmm) :
    Single_state_typing (hmm),
    Transition_matrix<T> (hmm),
    name (hmm.name),
    state_name (hmm.state_name),
    emit (hmm.emit)
  { }

  template<class S>
  Single_HMM (const Single_HMM<S>& hmm) :
    Single_state_typing (hmm.states(), hmm.alphabet()),
    Transition_matrix<T> (hmm.states()),
    name (hmm.name),
    state_name (hmm.state_name),
    emit (hmm.states(), vector<T> (hmm.alphabet().size()))
  {
    state_type = hmm.state_type;
  }

  template<class S>
  Single_HMM (const Single_HMM<S>& hmm, T t) :
    Single_state_typing (hmm.states(), hmm.alphabet()),
    Transition_matrix<T> (hmm.states(), t),
    name (hmm.name),
    state_name (hmm.state_name),
    emit (hmm.states(), vector<T> (hmm.alphabet().size(), t))
  {
    state_type = hmm.state_type;
  }

  // assignment operator

  Single_HMM<T>& operator= (const Single_HMM<T>& hmm)
  {
    ((Single_state_typing&) *this) = hmm;
    ((Transition_matrix<T>&) *this) = hmm;
    name = hmm.name;
    state_name = hmm.state_name;
    emit = hmm.emit;
    return *this;
  }

  // reset method

  void reset (const T& t)
  {
    reset_transitions (t);
    emit = vector<vector<T> > (states(), vector<T> (alphabet().size(), t));
  }

  // debugging output method
  void show (ostream& o) const;

  // read/write methods
  void write (ostream& out) const;
  void read (istream& in);

protected:

  // dimensional comparison

  template<class S>
  bool same_dimensions (const Single_HMM<S>& s) const
    {
      if (state_type != s.state_type) return 0;
      for (int i = 0; i < states(); ++i)
	if (emit[i].size() != s.emit[i].size()) return 0;
      return 1;
    }

  template<class S>
  void assert_same_dimensions (const Single_HMM<S>& s) const
  {
    if (!same_dimensions(s)) THROW Standard_exception ("Single HMM size mismatch");
  }

};

template<class T>
struct Single_meta_HMM : Single_HMM<T>
{
  // the metascore vector is a hack that allows certain states to be flagged as incurring site- and sequence-specific scores.
  // the original reason for this is to allow masking, but it also allows for Genewise-style splice site prediction,
  // among other things.
  vector<vector<int> > metascore_idx;   // metascore_idx[state] = index of metascore vector used by state <state>

  // helpers

  bool uses_metascores() const { for_const_contents (vector<vector<int> >, metascore_idx, m) if (m->size()) return 1; return 0; }
  void assert_no_metascores() const { if (uses_metascores()) THROW Standard_exception ("Can't handle metascores"); }
  int  max_metascore_idx() const
  {
    int max_i = -1;
    for_const_contents (vector<vector<int> >, metascore_idx, i)
      if (i->size())
	{
	  const int state_max = *(max_element (i->begin(), i->end()));
	  if (state_max > max_i) max_i = state_max;
	}
    return max_i;
  }

  // reset method
  void reset_meta() { metascore_idx = vector<vector<int> > (this->states(), vector<int>()); }

  // constructors

  Single_meta_HMM() :
    Single_HMM<T>()
  { }
  
  Single_meta_HMM (int states, const Alphabet& alphabet) :
    Single_HMM<T> (states, alphabet),
    metascore_idx (states, vector<int>())
  { }
  
  Single_meta_HMM (int states, const Alphabet& alphabet, T t) :
    Single_HMM<T> (states, alphabet, t),
    metascore_idx (states, vector<int>())
  { }

  template<class S>
  Single_meta_HMM (const Single_meta_HMM<S>& hmm) :
    Single_HMM<T> (hmm),
    metascore_idx (hmm.metascore_idx)
  { }

  template<class S>
  Single_meta_HMM (const Single_meta_HMM<S>& hmm, T t) :
    Single_HMM<T> (hmm),
    metascore_idx (hmm.metascore_idx)
  { }
};


// Single_HMM_scores contains transition and emission scores
// It is used as a template by other Single_HMM classes for dimensions, topology etc
//
// NB for backward-compatibility reasons, paths through Single_HMM's are all assumed to be global
// & *implicitly* have "ghost" Start states at the beginning and End states at the end.
//
// This may change, but for now, the wrapper calls to Transition_methods have to call
// make_global_path() and make_local_path() to add & remove these states appropriately.
//
struct Single_HMM_scores : Single_meta_HMM<Score>
{
  // constructor
  Single_HMM_scores();
  Single_HMM_scores (int states, const Alphabet& alphabet);

  // templated copy constructor

  template<class S>
  Single_HMM_scores (const Single_meta_HMM<S>& hmm) :
    Single_meta_HMM<Score> (hmm, -InfinityScore)
  { }

  // assignment
  Single_HMM_scores& operator= (const Single_HMM_scores& s);

  // emit methods
  vector<int> sample_state_path() const;  // wrapper for Transition_methods
  void sample_sequence (const vector<int>& state_path, Digitized_biosequence& dsq) const;       // sample sequence conditional on state path

  vector<int> consensus_state_path() const;  // wrapper for Transition_methods
  void consensus_sequence (const vector<int>& state_path, Digitized_biosequence& dsq) const;
  void make_profile (const vector<int>& state_path, Score_profile& prof_sc) const;

  // sparseness methods
  vector<int>           null_states_unsorted() const;   // null states
  vector<int>           null_states() const;            // null states, sorted topologically (croaks on null cycles)
  vector<int>           emit_states() const;            // non-null states
  vector<vector<int> >  incoming_states() const;        // for each state, find all the states that have transitions into it
  vector<vector<int> >  selected_outgoing_states (const vector<int>& selection) const;   // for each state, find all the states in (selection) that it has transitions into, keeping the order specified by (selection)

  // path score calculation
  Score path_transition_score (const vector<int>& state_path) const;  // wrapper for Transition_methods
  Score path_emit_score (const vector<int>& state_path, const Meta_profile& np) const;  // uses prof_sc's
  Score path_score (const vector<int>& state_path, const Meta_profile& np) const;  // uses prof_sc's
  Score path_emit_score_dsq (const vector<int>& state_path, const Named_profile& np) const;  // uses dsq's
  Score path_score_dsq (const vector<int>& state_path, const Named_profile& np) const;  // uses dsq's

  // misc methods
  void scale_all_scores (double beta);     // multiplies all scores by beta (useful for simulated annealing)
  
  const char* element_descriptor() const { return "scores"; }
  int  element_width() const { return 10; }
  void show_element (const Score& element, ostream& o) const { ShowScore (element, o); }
};

// Single_HMM_derivatives contains derivatives (w.r.t. some variable) of transition probabilities
//
struct Single_HMM_derivatives : Single_HMM<double>
{
  sstring wrt;  // derivative is "with respect to" (wrt) the parameter described in this string
  Single_HMM_derivatives (int states, const Alphabet& alphabet);
  Single_HMM_derivatives (int states, const char* x, const Alphabet& alphabet);   // x is a string describing the parameter w.r.t. which the derivatives are taken
  const char* element_descriptor() const { return wrt.c_str(); }
  int  element_width() const { return 10; }
  void show_element (const double& element, ostream& o) const { o << element; }
};

// Single_HMM_mask structure is used to restrict Single_HMM_scores parameter updates (e.g. during simple expectation maximization).
// The classes in "em.h" provide a more flexible way of achieving this.
//
struct Single_HMM_mask : Single_HMM<int>
{
  Single_HMM_mask (const Single_HMM_scores& hmm, bool init_flag = 1);
  void set_incoming_transitions (int state);   // NB sets start transition
  void clear_incoming_transitions (int state); // NB clears start transition
  void set_outgoing_transitions (int state);   // NB sets end transition
  void clear_outgoing_transitions (int state); // NB clears end transition
  void set_emit (int state);
  void clear_emit (int state);
  void set_all_start();
  void clear_all_start();
  void set_all_end();
  void clear_all_end();
  void set_all_transitions();    // NB calls set_all_start() & set_all_end()
  void clear_all_transitions();  // NB calls clear_all_start() & clear_all_end()
  void set_all_emit();
  void clear_all_emit();

  const char* element_descriptor() const { return "update flags"; }
  int  element_width() const { return 10; }
  void show_element (const bool& element, ostream& o) const { o << element; }
};

// Single_HMM_counts structure contains expected transition- and emission-usage counts
//
// NB:    Expected usage count for transition T   =   T * dP/dT    where P is the likelihood
//
struct Single_HMM_counts : Single_HMM<Prob>
{
  Loge   log_likelihood;    //  log(P)

  Single_HMM_counts (const Single_HMM_scores& hmm);   // constructor copies an existing HMM topology
  
  void   add_counts (const Single_HMM_counts& counts);
  void   subtract_counts (const Single_HMM_counts& counts);
  void   update_HMM_scores (Single_HMM_scores& hmm, const Single_HMM_mask& mask, bool sample = 0, double kT = 1) const;    // if sample==TRUE, probabilities are Dirichlet with pseudocounts scaled by 1/kT
  
  double dloglike_dx (const Single_HMM_scores& hmm, const Single_HMM_derivatives& deriv) const;         // deriv structure contains derivatives dT/dX of transition & emission probabilities w.r.t. X

  const char* element_descriptor() const { return "counts"; }
  int  element_width() const { return 10; }
  void show_element (const Prob& element, ostream& o) const { o << element; }

  // there are two versions of add_counts_from_state_path(), with & without metacounts
  void   add_counts_from_state_path (const Single_HMM_scores& hmm, const Named_profile& np, const vector<int>& state_path);
  void   add_counts_from_state_path (const Single_HMM_scores& hmm, const Named_profile& np, const vector<int>& state_path, vector<Metaprob>& metacounts);

  // overload read/write methods to include log_likelihood
  void write (ostream& out) const;
  void read (istream& in);
};


// template method code

template <class T>
void Single_HMM<T>::show (ostream& o) const
{
  int old_prec = o.precision(3);
  this->save_flags (o);
  this->right_align (o);

  const int w = this->element_width();

  this->show_transitions (o);

  o << "Emission profile " << this->element_descriptor() << ":\n";
  if (_alphabet)
    {
      o << "         ";
      for (int s = 0; s < _alphabet->size(); s++)
	{
	  sstring text;
	  text << (char) toupper (_alphabet->int2char(s)) << '=' << s;
	  o.width(w+1);
	  o << text;
	}
      o << "\n";
    }
  for (int i = 0; i < this->tm_states(); i++)
    {
      o << "State ";
      this->left_align (o);
      o.width(2);
      o << i << " (";
      this->right_align (o);
      for (int x = 0; x < (int) emit[i].size(); x++)
	{ o.width(w); show_element(emit[i][x],o); if (x < (int) emit[i].size()-1) o << " "; }
      o << ")";
      if (state_type[i] == Null) o << "  [null]";
      o << "\n";
    }
  this->restore_flags (o);
  o.precision (old_prec);
}

template <class T>
void Single_HMM<T>::write (ostream& out) const
{
  out << states() << " " << alphabet().size() << "\n";
  for (int s = 0; s < states(); ++s)
    out << this->start[s] << " ";
  out << "\n";
  for (int s = 0; s < states(); ++s)
    out << this->end[s] << " ";
  out << "\n";
  for (int s = 0; s < states(); ++s)
    {
      for (int d = 0; d < states(); ++d)
	out << this->transition(s,d) << " ";
      out << "\n";
    }
  for (int s = 0; s < states(); ++s)
    {
      for (int sym = 0; sym < alphabet().size(); ++sym)
	out << emit[s][sym] << " ";
      out << "\n";
    }
}

template <class T>
void Single_HMM<T>::read (istream& in)
{
  int in_states;
  int in_alph_sz;
  in >> in_states >> in_alph_sz;
  if (in_states != states() || in_alph_sz != alphabet().size())
    THROWEXPR ("Size mismatch while reading HMM");
  for (int s = 0; s < states(); ++s)
    in >> this->start[s];
  for (int s = 0; s < states(); ++s)
    in >> this->end[s];
  for (int s = 0; s < states(); ++s)
    {
      for (int d = 0; d < states(); ++d)
	in >> this->transition(s,d);
    }
  for (int s = 0; s < states(); ++s)
    {
      for (int sym = 0; sym < alphabet().size(); ++sym)
	in >> emit[s][sym];
    }
}

#endif
