#ifndef PAIRCFG_INCLUDED
#define PAIRCFG_INCLUDED

#include "hmm/singletmpl.h"
#include "hmm/pairhmm.h"
#include "seq/gff.h"
#include "scfg/foldenv.h"
#include "util/strsaver.h"

// No messing around here: #define alphabet stuff for speed
#define CFG_alphabet_size 4
#define CFG_alphabet      RNA_alphabet
// related params
#define CFG_alphabet_mask (CFG_alphabet_size - 1)  /* assumes alphabet size is power of 2 */
#define CFG_lex_size      (CFG_alphabet_size * CFG_alphabet_size * CFG_alphabet_size * CFG_alphabet_size)

// default names
#define CFG_default_name         "Pair_SCFG"
#define CFG_default_state_prefix "Pair_SCFG_state_"

// state type enumeration for pair CFG's
struct Pair_CFG_state_type_enum
{
  // state types
  enum State_type { EmitXL = 1, EmitXR = 2, EmitXLR = 3,
		    EmitYL = 4, EmitYR = 8, EmitYLR = 12,
		    EmitXLYL = 5, EmitXLYR = 9, EmitXLYLR = 13,
		    EmitXRYL = 6, EmitXRYR = 10, EmitXRYLR = 14,
		    EmitXLRYL = 7, EmitXLRYR = 11, EmitXLRYLR = 15,
		    Null = 0, Bifurc = 16, BifurcRevY = 32,
		    EmitStateTypes = 16,
		    Undefined = -1 };

  enum { ShiftXL = 0, ShiftXR = 1, ShiftYL = 2, ShiftYR = 3 };

  // query methods
  static bool is_emit_type (int t) { return t >= 0 && (t & 15) != 0; }
  static bool is_bifurc_type (int t) { return t >= 0 && (t & 48) != 0; }

  // display methods
  static const char* state_type_string (State_type t);

  // test methods
  static bool test_state_type_defined (int t) { return (t >= 0 && t <= 16) || t == 32; }

  // bifurcation connectivity descriptors
  struct Bifurcation
  {
    int l;
    int r;
    bool null() const { return l < 0 || r < 0; }
    bool operator== (const Bifurcation& b) const { return l == b.l && r == b.r; }
    Bifurcation () : l(-1), r(-1) { }
    Bifurcation (int l, int r) : l(l), r(r) { if (null()) THROWEXPR ("Tried to make illegal bifurcation"); }
  };
  struct Bifurcation_left_parent { int l; int p; Bifurcation_left_parent (int l, int p) : l(l), p(p) { } };
  struct Bifurcation_right_parent { int r; int p; Bifurcation_right_parent (int r, int p) : r(r), p(p) { } };
};
typedef Pair_CFG_state_type_enum::Bifurcation Bifurcation;

// RNA_pair_path holds the co-ordinates, alignment path & structures
// (but not the sequences) for a local alignment of two RNA sequences.
// The alignment is encoded as two gapped fold strings,
// e.g. xfold = "...<--<..>->.." is this alignment ****--****-***
//      yfold = "--.<<.<..>>>.-"                   --***********-
// plus structures for x(...<<..>>..) and y(.<<.<..>>>.)
struct RNA_pair_path : Fold_char_enum, Stream_saver
{
  // data
  int xseqlen, yseqlen;  // sequence lengths
  int xstart, ystart;    // coords of local alignment
  sstring xfold, yfold;  // gapped fold strings
  Score score;           // alignment score

  // ungapped fold string accessors
  sstring xfold_ungapped() const;
  sstring yfold_ungapped() const;
  int xfold_ungapped_len() const;
  int yfold_ungapped_len() const;

  // methods to return x & y Local_fold_string objects
  Local_fold_string local_xfold() const;
  Local_fold_string local_yfold() const;

  // method to get pairwise alignment
  Pairwise_path pairwise_path (bool global_pad = TRUE) const;

  // method to add or block out an alignment path from a fold envelope; calls Pair_envelope::add_pairwise_path()
  void add_pairwise_path (Pair_envelope& env, bool allow = TRUE) const;
};

// Pair_CFG_alignment with Named_profile's
struct Pair_CFG_alignment : RNA_pair_path, Pair_CFG_state_type_enum
{
  // data
  // aligned sequence references
  const Named_profile& npx;
  const Named_profile& npy;
  // State_type representation of parse tree (filled by Pair_CFG_branch)
  // Type Null is used to represent end nodes.
  vector<State_type> parse;

  // constructor
  Pair_CFG_alignment (const Named_profile& npx, const Named_profile& npy);

  // method to create Stockade local-alignment-with-subseqs for this alignment
  Stockade stockade (bool use_NSE = TRUE) const;  // if use_NSE is set, "Name/Start-End" name encoding is used

  // display method
  void show (ostream& o) const;
};

// local state path (single branch of parse tree, containing no bifurcations)
struct Pair_CFG_branch : Pair_CFG_state_type_enum, Fold_char_enum
{
  // coords of this subseq
  int xl;
  int xr;  // end + 1
  int yl;
  int yr;  // end + 1
  // tree connections
  int parent;
  int lchild;
  int rchild;
  bool revy;  // true if final state in path is of type BifurcRevY
  // state path
  vector<int> path;
  // constructor
  Pair_CFG_branch (int xl, int xr, int yl, int yr, int p);
  // accessors
  bool has_children() const { return lchild >= 0; }
  // recursive method to return alignment
  Pair_CFG_alignment alignment (const vector<Pair_CFG_branch>& parse_tree, const vector<State_type>& state_type, const Named_profile& npx, const Named_profile& npy) const;
  Pair_CFG_alignment alignment (const vector<Pair_CFG_branch>& parse_tree, const vector<State_type>& state_type, const set<int>& paired_states, const Named_profile& npx, const Named_profile& npy) const;
  // helper method
  static set<int> get_paired_states (const vector<State_type>& state_type);
};

// a parse tree is just a vector of branches
struct Pair_CFG_parse_tree : vector<Pair_CFG_branch>, Pair_CFG_state_type_enum
{
  // constructor
  Pair_CFG_parse_tree (int xlen, int ylen) { new_branch (0, xlen, 0, ylen, -1); }
  // build methods
  int new_branch (int xl, int xr, int yl, int yr, int p);
  int new_lchild (int xl, int xr, int yl, int yr, int p);
  int new_rchild (int xl, int xr, int yl, int yr, int p);
  // test methods
  bool test_connections() const;
  bool test_global (const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy) const;
  // method to return alignment
  Pair_CFG_alignment alignment (const vector<State_type>& state_type, const Named_profile& npx, const Named_profile& npy) const;
  Pair_CFG_alignment alignment (const vector<State_type>& state_type, const set<int>& paired_states, const Named_profile& npx, const Named_profile& npy) const;
  // debugging output method
  void show (ostream& out, const vector<State_type>* state_type = 0) const;
};

// local path structure... this is a bit too similar to a Pair_CFG_branch, but oh well
// used to return local tracebacks from CYK algorithm
struct Pair_CFG_local_path
{
  int xstart;
  int xlen;
  int ystart;
  int ylen;
  vector<int> path;
  // debugging output method
  void show (ostream& out) const;
  // constructors
  Pair_CFG_local_path();
  Pair_CFG_local_path (int xstart, int xlen, int ystart, int ylen, const vector<int>& path);
};

// state type container for pair CFG's
class Pair_CFG_state_typing : public Pair_CFG_state_type_enum
{
public:
  vector<State_type> state_type;
  vector<Bifurcation> bifurc;

  // constructors
  Pair_CFG_state_typing() : state_type(), bifurc() { }
  Pair_CFG_state_typing (int states) : state_type (states, (State_type) Undefined), bifurc (states) { }
  Pair_CFG_state_typing (const Pair_CFG_state_typing& base) : state_type (base.state_type), bifurc (base.bifurc) { }

  // assignment operator
  Pair_CFG_state_typing& operator= (const Pair_CFG_state_typing& base);
  void assign_state_typing (const Pair_CFG_state_typing& base);

  // accessors
  const inline Alphabet& alphabet() const { return CFG_alphabet; }
  inline int alphabet_size() const { return CFG_alphabet_size; }
  inline int states() const { return state_type.size(); }   // excludes start & end states

  // builders
  void init_bifurc (int state, int l, int r);

  // reset method
  void reset_state_types();

  // comparison method
  bool same_state_types (const Pair_CFG_state_typing& t) const;
  
  // sort methods
  vector<vector<Bifurcation_left_parent> > left_parent() const;
  vector<vector<Bifurcation_right_parent> > right_parent() const;

  // state path parsing methods
  Pair_CFG_parse_tree parse (const Pair_CFG_local_path& local_path) const;
  Pair_CFG_parse_tree parse (const vector<int>& state_path) const;
  Pair_CFG_parse_tree parse (const vector<int>& state_path, int xstart, int ystart) const;
  vector<int> unparse (const Pair_CFG_parse_tree& parse_tree) const;

  // test methods
  bool test_state_types_valid() const;
  bool test_branch_coords_consistent (const Pair_CFG_parse_tree& parse_tree) const;

  // method to check if this is a single or pair SCFG
  bool is_single_CFG() const;

  // indexing methods
  // for EmitXLRYLR, multipliers are XL=1, XR=4, YL=16, YR=64
  // see beginning of file "paircfg.cc" for details
  inline int emit_size (int type) const { return (type < 0 || type >= EmitStateTypes) ? 0 : _emit_sz[type]; }

  inline int emit_xl_mul (int type) const { return (type < 0 || type >= EmitStateTypes) ? 0 : _emit_xl_mul[type]; }
  inline int emit_xr_mul (int type) const { return (type < 0 || type >= EmitStateTypes) ? 0 : _emit_xr_mul[type]; }
  inline int emit_yl_mul (int type) const { return (type < 0 || type >= EmitStateTypes) ? 0 : _emit_yl_mul[type]; }
  inline int emit_yr_mul (int type) const { return (type < 0 || type >= EmitStateTypes) ? 0 : _emit_yr_mul[type]; }

  inline int emit_xl_idx (int type, int emit_idx) const { return emit_xl_mul(type)==0 ? -1 : ((emit_idx / emit_xl_mul(type)) & 3); }  // assumes CFG_alphabet_size == 4
  inline int emit_xr_idx (int type, int emit_idx) const { return emit_xr_mul(type)==0 ? -1 : ((emit_idx / emit_xr_mul(type)) & 3); }  // assumes CFG_alphabet_size == 4
  inline int emit_yl_idx (int type, int emit_idx) const { return emit_yl_mul(type)==0 ? -1 : ((emit_idx / emit_yl_mul(type)) & 3); }  // assumes CFG_alphabet_size == 4
  inline int emit_yr_idx (int type, int emit_idx) const { return emit_yr_mul(type)==0 ? -1 : ((emit_idx / emit_yr_mul(type)) & 3); }  // assumes CFG_alphabet_size == 4

  inline int emit_dxl (int type) const { return (type & EmitXL) >> ShiftXL; }
  inline int emit_dxr (int type) const { return (type & EmitXR) >> ShiftXR; }
  inline int emit_dyl (int type) const { return (type & EmitYL) >> ShiftYL; }
  inline int emit_dyr (int type) const { return (type & EmitYR) >> ShiftYR; }

  inline int emit_xylr_offset (int xl, int xr, int yl, int yr) const
  { return ((((yr&3)<<2) | ((yl&3))<<2) | ((xr&3))<<2) | (xl&3); }  // assumes CFG_alphabet_size == 4

  inline void decode_xylr_offset (int offset, int& xl, int& xr, int& yl, int& yr) const
  { xl=offset&3; offset>>=2; xr=offset&3; offset>>=2; yl=offset&3; offset>>=2; yr=offset&3; }
  
private:
  static int _emit_sz[EmitStateTypes];
  static int _emit_xl_mul[EmitStateTypes];
  static int _emit_xr_mul[EmitStateTypes];
  static int _emit_yl_mul[EmitStateTypes];
  static int _emit_yr_mul[EmitStateTypes];

  // recursive methods for building parse trees
  static void indent_xr (Pair_CFG_parse_tree& parse_tree, int node, int offset);
  static void indent_yr (Pair_CFG_parse_tree& parse_tree, int node, int offset);
  static void indent_yl (Pair_CFG_parse_tree& parse_tree, int node, int offset);
};

// template for pair stochastic context free grammars: transition matrix with emit profiles & bifurcations
template<class T>
class Pair_CFG : public Pair_CFG_state_typing, public Transition_matrix<T>
{
public:
  sstring name;
  vector<sstring> state_name;

  vector<vector<T> > emit;

  // metascore stuff
  vector<vector<int> > xlmeta_idx;   // xlmeta_idx[state] = index of xl-metascore vector used by state
  vector<vector<int> > xrmeta_idx;   // xrmeta_idx[state] = index of xr-metascore vector used by state
  vector<vector<int> > ylmeta_idx;   // ylmeta_idx[state] = index of yl-metascore vector used by state
  vector<vector<int> > yrmeta_idx;   // yrmeta_idx[state] = index of yr-metascore vector used by state

  void reset_meta();
  bool uses_metascores() const;
  void assert_no_metascores() const;
  int  max_xmeta_idx() const;
  int  max_ymeta_idx() const;
  int  max_xlmeta_idx() const;
  int  max_xrmeta_idx() const;
  int  max_ylmeta_idx() const;
  int  max_yrmeta_idx() const;
private:
  static int max_meta_idx (const vector<vector<int> >& meta_idx);
public:

  // constructors

  Pair_CFG();
  Pair_CFG (int states);
  Pair_CFG (int states, T t);
  Pair_CFG (const Pair_CFG<T>& cfg);

  template<class S>
  Pair_CFG (const Pair_CFG<S>& cfg, T t);

  // Single_HMM & Pair_HMM adaptor constructors
  // The resultant Pair SCFG is L-emitting.
  Pair_CFG (const Single_meta_HMM<T>& hmm);
  Pair_CFG (const Single_HMM<T>& hmm);
  Pair_CFG (const Pair_HMM<T>& hmm);
private:
  void init_from_HMM (const Single_HMM<T>& hmm);
public:

  // assignment operator
  Pair_CFG<T>& operator= (const Pair_CFG<T>& cfg);

  // builders
  void init_emit (int state, State_type type, const T& t);

  // emit accessor
  T& emit_by_chars (int state, const char* emit_string);

  // reset method
  void reset (const T& t);

  // output method
  void show (ostream& o) const;

  // test methods
  bool test_valid() const;

  // method to make a Pair_HMM from a Pair_CFG
  Pair_HMM<T> as_HMM() const;

protected:

  // dimensional comparison
  template<class S>
  bool same_dimensions (const Pair_CFG<S>& s) const;

  template<class S>
  void assert_same_dimensions (const Pair_CFG<S>& s) const;
};

// Pair_CFG_scores contains transition and emission scores
// It is used as a template by other Pair_CFG classes for dimensions, topology etc
//
class Pair_CFG_scores : public Pair_CFG<Score>
{
public:
  // constructors
  Pair_CFG_scores() : Pair_CFG<Score>() { }
  Pair_CFG_scores (int states) : Pair_CFG<Score> (states, -InfinityScore) { }
  Pair_CFG_scores (const Single_HMM_scores& hmm) : Pair_CFG<Score> (hmm) { }  // creates L-emitting SCFG
  Pair_CFG_scores (const Pair_HMM_scores& hmm) : Pair_CFG<Score> (hmm) { }  // creates L-emitting SCFG

  // templated copy constructor
  template<class S>
  Pair_CFG_scores (const Pair_CFG<S>& cfg);

  // set state type method
  void init_emit (int state, State_type type, Score sc = -InfinityScore);

  // path score calculation
  Score path_transition_score (const Pair_CFG_parse_tree& parse_tree) const;
  Score path_emit_score (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy) const;
  Score path_score (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy) const;

  // report/mask methods
  struct XYMask { vector<bool> xmask; vector<bool> ymask; };
  XYMask get_mask (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy, const set<int>& mask_states) const;
  GFF ymask2gff (const XYMask& mask, const Named_profile& xprof, const Named_profile& yprof) const;

  // sort methods
  vector<int> emit_states() const;
  vector<int> nonemit_states_unsorted() const;
  vector<int> nonemit_states() const;  // sorted topologically; throws exception if null cycle detected
  vector<vector<int> > incoming_states() const;
  vector<vector<int> > selected_outgoing_states (const vector<int>& selection) const;
  vector<vector<int> > selected_incoming_states (const vector<int>& selection) const;

  // Elimination of null bifurcations (transformation A->B) and null state paths (transformation B->C).
  // Notes:
  //  -- Caller must call Fold_envelope::remove_empty_bifurcations() to enforce non-empty constraint on DP algorithms when Pair SCFGs with eliminated null bifurcations are used.

  // Two independent transformations:
  // A->B: adds a new null states for each bifurcation, factoring in null-subtree bifurcations as direct transitions
  // method call on a_cfg returns b_cfg
  Pair_CFG_scores add_null_bifurc_transitions (vector<Prob>& empty_probs,  // a co-ordinates
					       vector<int>& b2a_state,  // b->a mapping
					       vector<int>& bifurc_prec_states) const;  // b co-ordinates

  // B->C: eliminates all null states, factoring paths through null states into direct transitions between emit/bifurcation states
  // method call on b_cfg returns c_cfg
  Pair_CFG_scores eliminate_null_states (Concrete_transition_probs& tp_orig,  // b co-ordinates
					 Concrete_transition_probs& tp_elim,  // b co-ordinates
					 const vector<int>& bifurc_prec_states,  // b co-ordinates
					 vector<int>& c2b_state) const;  // c->b mapping (state indices change due to removal of null states)

  // sample paths through null states (i.e. from a C-path, restore the B-path).
  // returned state path is in b co-ordinates
  static vector<int> sample_eliminated_states (const Transition_probs& tp_orig,  // b co-ords
					       const Transition_probs& tp_elim,  // b co-ords
					       const vector<int>& c2b_state,  // c->b
					       const Pair_CFG_scores& b_cfg,
					       const vector<int>& c_path);  // c co-ords

  // sample null (no-sequence emitting) subtrees (i.e. from a B-path, restore the A-path).
  // returned state path is in a co-ordinates
  static vector<int> sample_eliminated_bifurcations (const vector<Prob>& empty_probs,  // a co-ords
						     const vector<int>& b2a_state,  // b->a
						     const Pair_CFG_scores& a_cfg,
						     const Pair_CFG_scores& b_cfg,
						     const vector<int>& b_path);  // b co-ords

  // helper to return null (non-emitting) parse tree probabilities
  vector<Prob> empty_probs() const;  // null Inside probabilities, P(empty_sequence|state), indexed by state (excluding Start & End)

  // helper to sample from the distribution of null parse trees rooted at a particular state
  static vector<int> sample_null_subtree (const vector<Prob>& empty_probs,  // a co-ords
					  const Pair_CFG_scores& a_cfg,
					  int root_state);  // a co-ords

  // wrapper to combine transformations A->B and B->C into a single transformation: A->C
  Pair_CFG_scores eliminate_null_states_and_bifurcations (vector<Prob>& empty_probs,  // a co-ordinates
							  Concrete_transition_probs& tp_orig,  // b co-ordinates
							  Concrete_transition_probs& tp_elim,  // b co-ordinates
							  vector<int>& b2a_state,  // b->a mapping
							  vector<int>& c2b_state) const;  // c->b mapping

  // wrapper to restore null paths AND subtrees (i.e. from a C-path, restore the A-path).
  // returned state path is in a co-ordinates
  static vector<int> sample_eliminated_states_and_bifurcations (const vector<Prob>& empty_probs,  // a co-ords
								const Transition_probs& tp_orig,  // b co-ords
								const Transition_probs& tp_elim,  // b co-ords
								const vector<int>& b2a_state,  // b->a
								const vector<int>& c2b_state,  // c->b
								const Pair_CFG_scores& a_cfg,
								const Pair_CFG_scores& b_cfg,
								const vector<int>& c_path);  // c co-ords

  // End of null state & bifurcation elimination method prototypes


  // display methods
  const char* element_descriptor() const { return "scores"; }
  int  element_width() const { return 10; }
  void show_element (const Score& element, ostream& o) const { ShowScore (element, o); }

  // parse tree score breakdown display method
  void show_score_breakdown (const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy, ostream& out) const;

  // internal methods
private:
  void zero_outgoing_bifurcation_transitions();  // sets transitions from bifurc states to -InfinityScore
  void add_fake_bifurcation_transitions();  // adds two "fake" zero-score transitions from bifurc states, to fool topological sorter
};

// Pair_CFG_counts structure contains expected transition-usage counts
//
// NB:    Expected usage count for transition T   =   T * dP/dT    where P is the likelihood
//
struct Pair_CFG_counts : Pair_CFG<Prob>
{
  Loge log_likelihood;    // log(P)

  // constructors
  Pair_CFG_counts (const Pair_CFG_scores& cfg) : Pair_CFG<Prob> (cfg, 0.0), log_likelihood (0.0) { }
  Pair_CFG_counts (const Single_HMM_counts& hmm) : Pair_CFG<Prob> (hmm), log_likelihood (hmm.log_likelihood) { }
  Pair_CFG_counts (const Pair_HMM_counts& hmm) : Pair_CFG<Prob> (hmm), log_likelihood (hmm.log_likelihood) { }
  
  // display methods
  const char* element_descriptor() const { return "counts"; }
  int  element_width() const { return 10; }
  void show_element (const Prob& element, ostream& o) const { o << element; }

  // update methods
  void add_counts_from_parse_tree (const Pair_CFG_scores& cfg, const Pair_CFG_parse_tree& parse_tree, const Digitized_biosequence& dsqx, const Digitized_biosequence& dsqy);
};

// template method code

template<class T>
Pair_CFG<T>::Pair_CFG() :
  Pair_CFG_state_typing(),
  Transition_matrix<T>(),
  name (CFG_default_name),
  state_name(),
  emit(),
  xlmeta_idx(),
  xrmeta_idx(),
  ylmeta_idx(),
  yrmeta_idx()
{ }

template<class T>
Pair_CFG<T>::Pair_CFG (int states) :
  Pair_CFG_state_typing (states),
  Transition_matrix<T> (states),
  name (CFG_default_name),
  state_name (states, sstring(CFG_default_state_prefix)),
  emit (states),
  xlmeta_idx (states),
  xrmeta_idx (states),
  ylmeta_idx (states),
  yrmeta_idx (states)
{
  for (int s = 0; s < states; ++s)
    state_name[s] << (s + 1);
}
  
template<class T>
Pair_CFG<T>::Pair_CFG (int states, T t) :
  Pair_CFG_state_typing (states),
  Transition_matrix<T> (states, t),
  name (CFG_default_name),
  state_name (states, sstring(CFG_default_state_prefix)),
  emit (states),
  xlmeta_idx (states),
  xrmeta_idx (states),
  ylmeta_idx (states),
  yrmeta_idx (states)
{
  for (int s = 0; s < states; ++s)
    state_name[s] << (s + 1);
}

template<class T>
Pair_CFG<T>::Pair_CFG (const Pair_CFG<T>& cfg) :
  Pair_CFG_state_typing (cfg),
  Transition_matrix<T> (cfg),
  name (cfg.name),
  state_name (cfg.state_name),
  emit (cfg.emit),
  xlmeta_idx (cfg.xlmeta_idx),
  xrmeta_idx (cfg.xrmeta_idx),
  ylmeta_idx (cfg.ylmeta_idx),
  yrmeta_idx (cfg.yrmeta_idx)
{ }

template<class T>
template<class S>
Pair_CFG<T>::Pair_CFG (const Pair_CFG<S>& cfg, T t) :
  Pair_CFG_state_typing (cfg),
  Transition_matrix<T> (cfg.states()),
  name (cfg.name),
  state_name (cfg.state_name),
  xlmeta_idx (cfg.xlmeta_idx),
  xrmeta_idx (cfg.xrmeta_idx),
  ylmeta_idx (cfg.ylmeta_idx),
  yrmeta_idx (cfg.yrmeta_idx)
{
  reset (t);
}

template<class T>
Pair_CFG<T>& Pair_CFG<T>::operator= (const Pair_CFG<T>& cfg)
{
  assign_state_typing (cfg);
  Transition_matrix<T>::assign_transition_matrix (cfg);
  name = cfg.name;
  state_name = cfg.state_name;
  emit = cfg.emit;
  xlmeta_idx = cfg.xlmeta_idx;
  xrmeta_idx = cfg.xrmeta_idx;
  ylmeta_idx = cfg.ylmeta_idx;
  yrmeta_idx = cfg.yrmeta_idx;
  return *this;
}

template<class T>
  void Pair_CFG<T>::init_emit (int state, State_type type, const T& t)
  {
    state_type[state] = type;
    emit[state] = vector<T> (emit_size(type), t);
  }

template<class T>
  T& Pair_CFG<T>::emit_by_chars (int state, const char* emit_string)
{
  int emit_idx = 0;
  int pos = 0;
  int mul = 1;
  for (int b = 0; b < 4; ++b)
    if (state_type[state] & (1 << b))
      {
	emit_idx += mul * alphabet().char2int_strict (emit_string[pos++]);
	mul *= alphabet_size();
      }
  if (pos != (int) strlen (emit_string)) THROWEXPR ("Bad emit string: '" << emit_string << "'");
  return emit[state][emit_idx];
}

template<class T>
void Pair_CFG<T>::reset (const T& t)
{
  Transition_matrix<T>::reset_transitions (t);
  emit = vector<vector<T> > (states());
  for (int s = 0; s < states(); ++s) emit[s] = vector<T> (emit_size (state_type[s]), t);
  xlmeta_idx = vector<vector<int> > (states());
  xrmeta_idx = vector<vector<int> > (states());
  ylmeta_idx = vector<vector<int> > (states());
  yrmeta_idx = vector<vector<int> > (states());
}

template<class T>
template<class S>
  bool Pair_CFG<T>::same_dimensions (const Pair_CFG<S>& s) const
    {
      if (state_type != s.state_type) return 0;
      for (int i = 0; i < states(); ++i)
	if (emit[i].size() != s.emit[i].size()) return 0;
      return 1;
    }

template<class T>
template<class S>
  void Pair_CFG<T>::assert_same_dimensions (const Pair_CFG<S>& s) const
  {
    if (!same_dimensions(s)) THROW Standard_exception ("Pair CFG size mismatch");
  }

template <class T>
Pair_CFG<T>::Pair_CFG (const Single_HMM<T>& hmm) :
  Pair_CFG_state_typing (hmm.states()),
  Transition_matrix<T> (hmm),
  name (hmm.name),
  state_name (hmm.state_name),
  emit (hmm.states()),
  xlmeta_idx (hmm.states()),
  xrmeta_idx (hmm.states()),
  ylmeta_idx (hmm.states()),
  yrmeta_idx (hmm.states())
{
  init_from_HMM (hmm);
}

template <class T>
Pair_CFG<T>::Pair_CFG (const Single_meta_HMM<T>& hmm) :
  Pair_CFG_state_typing (hmm.states()),
  Transition_matrix<T> (hmm),
  name (hmm.name),
  state_name (hmm.state_name),
  emit (hmm.states()),
  xlmeta_idx (hmm.metascore_idx),
  xrmeta_idx (hmm.states()),
  ylmeta_idx (hmm.states()),
  yrmeta_idx (hmm.states())
{
  init_from_HMM (hmm);
}

template <class T>
void Pair_CFG<T>::init_from_HMM (const Single_HMM<T>& hmm)
{
  if (&alphabet() != &hmm.alphabet()) THROWEXPR ("Alphabet mismatch");
  for (int s = 0; s < states(); ++s)
    switch (hmm.state_type[s])
      {
      case Single_state_typing::Null:
	state_type[s] = Null;
	break;
      case Single_state_typing::Emit:
	state_type[s] = EmitXL;
	emit[s] = hmm.emit[s];
	break;
      default:
	THROWEXPR ("State type unknown");
	break;
      }
}

template <class T>
Pair_HMM<T> Pair_CFG<T>::as_HMM() const
{
  Pair_HMM<T> hmm (states(), &CFG_alphabet);
  hmm.name = name;
  hmm.assign_transition_matrix (*this);
  for (int s = 0; s < states(); ++s)
    {
      hmm.state_name[s] = state_name[s];
      const State_type t = state_type[s];
      switch (t)
	{
	case EmitXL:
	  hmm.init_emit (s, Pair_HMM_state_type_enum::EmitX, emit[s][0]);
	  for (int i = 0; i < CFG_alphabet_size; ++i)
	    hmm.single_emit[s][i] = emit[s][i];
	  break;
	case EmitYL:
	  hmm.init_emit (s, Pair_HMM_state_type_enum::EmitY, emit[s][0]);
	  for (int i = 0; i < CFG_alphabet_size; ++i)
	    hmm.single_emit[s][i] = emit[s][i];
	  break;
	case EmitXLYL:
	  hmm.init_emit (s, Pair_HMM_state_type_enum::EmitXY, emit[s][0]);
	  for (int i = 0; i < CFG_alphabet_size; ++i)
	    for (int j = 0; j < CFG_alphabet_size; ++j)
	      hmm.pair_emit[s](i,j) = emit[s][i*emit_xl_mul(t) + j*emit_yl_mul(t)];
	  break;
	case Null:
	  hmm.state_type[s] = Pair_HMM_state_type_enum::Null;
	  break;
	default:
	  THROWEXPR ("Can't convert states of type " << state_type_string (t) << " to HMM states");
      }
    }
  return hmm;
}
  
template <class T>
Pair_CFG<T>::Pair_CFG (const Pair_HMM<T>& hmm) :
  Pair_CFG_state_typing (hmm.states()),
  Transition_matrix<T> (hmm.states()),
  name (hmm.name),
  state_name (hmm.state_name),
  emit (hmm.states()),
  xlmeta_idx (hmm.states()),  // Pair_HMM's currently don't use metascores
  xrmeta_idx (hmm.states()),
  ylmeta_idx (hmm.states()),
  yrmeta_idx (hmm.states())
{
  if (&alphabet() != hmm.alphabet) THROWEXPR ("Alphabet mismatch");
  this->start_to_end() = hmm.start_to_end();
  for (int s = 0; s < states(); ++s)
    {
      this->start[s] = hmm.start[s];
      this->end[s] = hmm.end[s];
      for (int d = 0; d < states(); ++d)
	this->transition (s, d) = hmm.transition (s, d);
      switch (hmm.state_type[s])
	{
	case Pair_HMM<T>::EmitX:
	  state_type[s] = EmitXL;
	  emit[s] = hmm.single_emit[s];
	  break;
	case Pair_HMM<T>::EmitY:
	  state_type[s] = EmitYL;
	  emit[s] = hmm.single_emit[s];
	  break;
	case Pair_HMM<T>::EmitXY:
	  {
	    const State_type t = state_type[s] = EmitXLYL;
	    emit[s] = vector<T> (emit_size (t));
	    for (int i = 0; i < alphabet_size(); ++i)
	      for (int j = 0; j < alphabet_size(); ++j)
		emit[s] [i*emit_xl_mul(t) + j*emit_yl_mul(t)] = hmm.pair_emit[s] (i, j);
	    break;
	  }
	default:
	  THROWEXPR ("State type unknown");
	  break;
	}
    }
}

template <class T>
bool Pair_CFG<T>::uses_metascores() const
{
  for_const_contents (vector<vector<int> >, xlmeta_idx, m)
    if (m->size())
      return 1;
  for_const_contents (vector<vector<int> >, xrmeta_idx, m)
    if (m->size())
      return 1;
  for_const_contents (vector<vector<int> >, ylmeta_idx, m)
    if (m->size())
      return 1;
  for_const_contents (vector<vector<int> >, yrmeta_idx, m)
    if (m->size())
      return 1;
  return 0;
}

template <class T>
void Pair_CFG<T>::assert_no_metascores() const
{
  if (uses_metascores()) THROW Standard_exception ("Can't handle metascores");
}

template <class T>
int Pair_CFG<T>::max_meta_idx (const vector<vector<int> >& meta_idx)
{
  int max_i = -1;
  for_const_contents (vector<vector<int> >, meta_idx, i)
    if (i->size())
      {
	const int state_max = *(max_element (i->begin(), i->end()));
	if (state_max > max_i) max_i = state_max;
      }
  return max_i;
}

template <class T>
int Pair_CFG<T>::max_xmeta_idx() const
{ return max (max_xlmeta_idx(), max_xrmeta_idx()); }

template <class T>
int Pair_CFG<T>::max_ymeta_idx() const
{ return max (max_ylmeta_idx(), max_yrmeta_idx()); }

template <class T>
int Pair_CFG<T>::max_xlmeta_idx() const
{ return max_meta_idx (xlmeta_idx); }

template <class T>
int Pair_CFG<T>::max_xrmeta_idx() const
{ return max_meta_idx (xrmeta_idx); }

template <class T>
int Pair_CFG<T>::max_ylmeta_idx() const
{ return max_meta_idx (ylmeta_idx); }

template <class T>
int Pair_CFG<T>::max_yrmeta_idx() const
{ return max_meta_idx (yrmeta_idx); }

template <class T>
void Pair_CFG<T>::show (ostream& o) const
{
  int old_prec = o.precision(3);
  this->save_flags (o);
  this->right_align (o);

  this->show_transitions (o);

  o << "Emission profile " << this->element_descriptor() << ":\n";

  const int w = this->element_width();

  int max_sz = 0;
  for (int s = 0; s < this->tm_states(); ++s) max_sz = max (max_sz, emit_size(state_type[s]));
  for (; max_sz > 1; max_sz = max_sz / alphabet_size())
    {
      o << "                     ";
      for (int s = 0; s < max_sz; s++)
	{
	  sstring text;
	  for (int i = 1; i < max_sz; i *= alphabet_size())
	    text << (char) toupper (alphabet().int2char((s / i) % alphabet_size()));
	  text << '=' << s;
	  o.width(w+1);
	  o << text;
	}
      o << "\n";
    }
  for (int i = 0; i < this->tm_states(); i++)
    {
      sstring text;
      text << "State " << i << "[" << state_type_string (state_type[i]) << "]";
      this->left_align (o);
      o.width(20);
      o << text << " (";
      this->right_align (o);
      for (int x = 0; x < (int) emit[i].size(); x++) { o.width(w); Pair_CFG<T>::show_element(emit[i][x],o); if (x < (int) emit[i].size()-1) o << " "; }
      o << ")";
      for_const_contents (vector<int>, xlmeta_idx[i], m) o << " +xl[" << *m << "]";
      for_const_contents (vector<int>, xrmeta_idx[i], m) o << " +xr[" << *m << "]";
      for_const_contents (vector<int>, ylmeta_idx[i], m) o << " +yl[" << *m << "]";
      for_const_contents (vector<int>, yrmeta_idx[i], m) o << " +yr[" << *m << "]";
      if (state_type[i] == Bifurc || state_type[i] == BifurcRevY)
	o << " (bifurcates to " << bifurc[i].l << ", " << bifurc[i].r << ")";
      o << "\n";
    }

  this->restore_flags (o);
  o.precision (old_prec);
}

template <class T>
bool Pair_CFG<T>::test_valid() const
{
  if ((int) emit.size() != states()) return 0;
  if ((int) xlmeta_idx.size() != states()) return 0;
  if ((int) xrmeta_idx.size() != states()) return 0;
  if ((int) ylmeta_idx.size() != states()) return 0;
  if ((int) yrmeta_idx.size() != states()) return 0;
  if (!test_state_types_valid()) return 0;
  for (int s = 0; s < states(); ++s)
    if ((int) emit[s].size() != emit_size(state_type[s])) return 0;
  return 1;
}

template<class T>
void Pair_CFG<T>::reset_meta()
{
  xlmeta_idx = vector<vector<int> > (states(), vector<int>());
  xrmeta_idx = vector<vector<int> > (states(), vector<int>());
  ylmeta_idx = vector<vector<int> > (states(), vector<int>());
  yrmeta_idx = vector<vector<int> > (states(), vector<int>());
}

template<class S>
Pair_CFG_scores::Pair_CFG_scores (const Pair_CFG<S>& cfg) :
  Pair_CFG<Score> (cfg, -InfinityScore)
{ }

#endif /* PAIRCFG_INCLUDED */
