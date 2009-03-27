#ifndef TRIPLETSCFG_INCLUDED
#define TRIPLETSCFG_INCLUDED

#include "seq/biosequence.h"
#include "hmm/transmat.h"
#include "scfg/foldenv.h"
#include "scfg/paircfg.h"
#include "util/strsaver.h"

/// No messing around here: #define alphabet stuff for speed
#define SCFG_alphabet_size 4
#define SCFG_alphabet      RNA_alphabet

/// Enumeration of allowed SCFG state types.
struct SCFG_state_type_enum
{
  /// State types
  /*
   * State types are hashed by taking the decimal of (xl yl zl zr yr xr), where xl is a binary variable (1 bit)
   * such that 0 indicates no X left-emission and 1 indicates a X left-emission.
   * State types Start, End are handled separately as special cases by Grammar_state_enum (HMM_state_enum for paircfg case)
   * (see transmat.h).
   */
  enum State_type { Null = 0,
		    EmitXR = 1, EmitYR = 2, EmitXRYR = 3, 
		    EmitZR = 4, EmitXRZR = 5, EmitYRZR = 6, EmitXRYRZR = 7, 
		    EmitZL = 8, EmitXRZL = 9, EmitYRZL = 10, EmitXRYRZL = 11, 
		    EmitZLR = 12, EmitXRZLR = 13, EmitYRZLR = 14, EmitXRYRZLR = 15, 
		    EmitYL = 16, EmitXRYL = 17, EmitYLR = 18, EmitXRYLR = 19, 
		    EmitYLZR = 20, EmitXRYLZR = 21, EmitYLRZR = 22, EmitXRYLRZR = 23, 
		    EmitYLZL = 24, EmitXRYLZL = 25, EmitYLRZL = 26, EmitXRYLRZL = 27, 
		    EmitYLZLR = 28, EmitXRYLZLR = 29, EmitYLRZLR = 30, EmitXRYLRZLR = 31, 
		    EmitXL = 32, EmitXLR = 33, EmitXLYR = 34, EmitXLRYR = 35, 
		    EmitXLZR = 36, EmitXLRZR = 37, EmitXLYRZR = 38, EmitXLRYRZR = 39, 
		    EmitXLZL = 40, EmitXLRZL = 41, EmitXLYRZL = 42, EmitXLRYRZL = 43, 
		    EmitXLZLR = 44, EmitXLRZLR = 45, EmitXLYRZLR = 46, EmitXLRYRZLR = 47, 
		    EmitXLYL = 48, EmitXLRYL = 49, EmitXLYLR = 50, EmitXLRYLR = 51, 
		    EmitXLYLZR = 52, EmitXLRYLZR = 53, EmitXLYLRZR = 54, EmitXLRYLRZR = 55, 
		    EmitXLYLZL = 56, EmitXLRYLZL = 57, EmitXLYLRZL = 58, EmitXLRYLRZL = 59, 
		    EmitXLYLZLR = 60, EmitXLRYLZLR = 61, EmitXLYLRZLR = 62, EmitXLRYLRZLR = 63,
		    Bifurc = 64,
		    Undefined = -1 };

  enum State_type_ancestral { NullW = 0,
			      EmitWR = 1, EmitWL = 2, EmitWLR = 3,
			      BifurcW = 4 };

  // to do: why doesn't the below work?
  //const State_type& operator++ (State_type &t) { return t = static_cast<State_type> (t+1); }

  /// Is type emit?
  static bool is_emit_type (int t) { return (t > 0 && t < 64); }
  /// Is type bifurcation?
  static bool is_bifurc_type (int t) { return (t == Bifurc); }

  /// Get state type string.
  static const char* state_type_string (State_type t);

  /// Store the child states of a given bifurcation state.
  /*
   *  Used in SCFG_state_typing.
   */
  struct Bifurcation
  {
    int l; /// left child state
    int r; /// right child state
    bool null() const { return l < 0 || r < 0; }
    bool operator== (const Bifurcation& b) const { return l == b.l && r == b.r; }
    
    // constructors
    Bifurcation() : l (-1), r (-1) {}
    Bifurcation (int l, int r) : l (l), r (r) { if (null()) THROWEXPR ("Tried to make illegal bifurcation"); }
  };
  /// Store the parent and left sibling of a "right" bifurcation state.
  /*
   * To be used for precomputing handy things before doing DP loops.
   */
  struct Bifurcation_left_parent {
    int l; /// "left sibling" state
    int p; /// parent state
    Bifurcation_left_parent (int l, int p) : l(l), p(p) { }
  };
  /// Store the parent and right sibling of a "left" bifurcation state.
  /*
   * To be used for precomputing handy things before doing DP loops.
   */
  struct Bifurcation_right_parent {
    int r; /// "right sibling" state
    int p; /// parent state
    Bifurcation_right_parent (int r, int p) : r(r), p(p) { }
  };

};
typedef SCFG_state_type_enum::State_type State_type;
//typedef SCFG_state_type_enum::Bifurcation Bifurcation;


/// Store type and bifurcation information about all states.
/*
 * Primarily useful as a parent class of Triplet_SCFG.
 */
struct SCFG_state_typing : SCFG_state_type_enum
{
 public:
  /// State typing for all states.
  /* 
   * Used to store the state types of all states in a SCFG (Triplet_SCFG inherits from SCFG_state_typing).
   * The state type of state s is accessed as state_type[s].
   * We can exploit our hashing of emit state types to determine if a state s e.g. left-emits to X
   * with (state_type[s] & EmitXL != 0).
   */
  vector<State_type> state_type;

  // Ancestral emit typing
  typedef map<int, State_type_ancestral> State_type_ancestral_map;
  State_type_ancestral_map state_type_ancestral_map;

  /// Number of possible emissions from state s.
  // _emit_size is defined in the .cc file.  This is ugly but fast.
  inline int emit_size (int s) const { return (state_type[s] < 0 || state_type[s] >= 64) ? 0 : _emit_size[state_type[s]]; }

  /// Multipliers for mapping emission strings to integers for array access.
  /*
   * We hash emissions as (xl yl zl zr yr xr).  These functions allow us to 
   * map the set of emission strings of a particular state type into a continuous 
   * set of integers (starting at 0).
   * Also see the member 'emit' of Triplet_SCFG.
   */
  static inline int emit_xl_mul (int type) { return (type < 0 || type >= 64) ? 0 : _emit_xl_mul[type]; }
  static inline int emit_xr_mul (int type) { return (type < 0 || type >= 64) ? 0 : _emit_xr_mul[type]; }
  static inline int emit_yl_mul (int type) { return (type < 0 || type >= 64) ? 0 : _emit_yl_mul[type]; }
  static inline int emit_yr_mul (int type) { return (type < 0 || type >= 64) ? 0 : _emit_yr_mul[type]; }
  static inline int emit_zl_mul (int type) { return (type < 0 || type >= 64) ? 0 : _emit_zl_mul[type]; }
  static inline int emit_zr_mul (int type) { return (type < 0 || type >= 64) ? 0 : _emit_zr_mul[type]; }

  /// Helper to get the hash for a particular emission.
  static inline int emit_idx (int t, int xl, int xr, int yl, int yr, int zl, int zr);

  /// Maps an emit hash to the corresponding emit string (xl yl zl zr yr xr).
  /*
   * Emit hashes are a la those defined in Triplet_SCFG::path_emit_score() and testripletscfg.cc.
   */
  sstring emit_hash_to_string (int s, int emit_hash) const;


  /// List of the left and right children of bifurcation states.
  /*
   * bifurc[a] gets the Bifurcation struct corresponding to bifurcation state 'a'.
   */
  vector<Bifurcation> bifurc;

  // constructors
  SCFG_state_typing() : state_type(), state_type_ancestral_map(), bifurc() { } /// minimal constructor
  SCFG_state_typing (int num_states) : state_type (num_states, (State_type) Undefined), bifurc (num_states) { }
  SCFG_state_typing (const SCFG_state_typing& base) : state_type (base.state_type), state_type_ancestral_map (base.state_type_ancestral_map), bifurc (base.bifurc) { } /// copy constructor

  /// typing assignments
  SCFG_state_typing& operator= (const SCFG_state_typing& base);
  void assign_state_typing (const SCFG_state_typing& base);

  /// accessors
  const inline Alphabet& alphabet() const { return SCFG_alphabet; }
  inline int alphabet_size() const { return SCFG_alphabet_size; }
  inline int num_states() const { return state_type.size(); }   // return number of states, excluding start & end states

  // sorting
  vector<vector<Bifurcation_left_parent> > left_parent() const;
  vector<vector<Bifurcation_right_parent> > right_parent() const;

 private:
  static int _emit_size[64];
  static int _emit_xl_mul[64];
  static int _emit_xr_mul[64];
  static int _emit_yl_mul[64];
  static int _emit_yr_mul[64];
  static int _emit_zl_mul[64];
  static int _emit_zr_mul[64];

};


/// Holds a local alignment (coordinates and structures) of 3 RNA sequences.
/*
 * The alignment is encoded as two gapped fold strings,
 * e.g. xfold = "...<--<..>->.." is this alignment ****--****-***
 *      yfold = "--.<<.<..>>>.-"                   --***********-
 * plus structures for x(...<<..>>..) and y(.<<.<..>>>.)
 * (description taken from scfg/paircfg.h)
 */
struct RNA_triplet_path : Fold_char_enum, Stream_saver
{
  // data
  int seqlen_x, seqlen_y, seqlen_z;  /// sequence lengths
  int start_x, start_y, start_z; /// coordinates of local alignment
  sstring foldstring_w, foldstring_x, foldstring_y, foldstring_z; /// gapped fold strings
  Score score; /// alignment score

  /// Get ungapped X fold string.
  sstring foldstring_ungapped_x() const;
  sstring foldstring_ungapped_y() const;
  sstring foldstring_ungapped_z() const;
  int foldstring_ungapped_x_len() const;
  int foldstring_ungapped_y_len() const;
  int foldstring_ungapped_z_len() const;

  /// Convert the fold strings into a 3-way alignment.
  /*
   * This is the analog of the alignment accessor pairwise_path() in scfg/paircfg.h.
   * Here "alignment" means a matrix with 3 rows, one for each sequence, of 0s and 1s
   * A 1 indicates that the sequence cooresponding to that row has a non-gap character there in the 3-way alignment.
   * (This is how it's implemented in scfg/paircfg.cc.  Basically the routine pairwise_path()
   *  looks at the fold strings and reads off the alignment by determining which positions have gap characters. )
   * It's not implemented yet.  I'll do it when I need it.
   */
  Alignment_path alignment_path() const;

};

/// Holds a local alignment with Named_profile's.
/*
 * The alignment is stored by the X, Y, Z foldstrings, which 
 * are members of the parent class RNA_triplet_path.
 */

// IH, 8/7/2008:
// Having this inherit from SCFG_state_typing is slightly dangerous.
// Why? Because this means that there are member variables state_type[] and state_type_ancestral_map[] that are meaningless.
// Thus, this object can (for example) access state_type without compile-time errors, but the results will be undefined.
// It should inherit from SCFG_state_type_enum instead.
struct Triplet_SCFG_alignment : RNA_triplet_path, SCFG_state_typing
{

  /// aligned sequences
  const Named_profile np_x, np_y, np_z;

  /// State_type representation of the parse tree.
  vector<State_type> parse;

  /// constructor
  Triplet_SCFG_alignment (const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z);

  /// display
  void show (ostream& o) const;

};

/// Holds a local state path, meaning a single branch of the parse tree with no bifurcations.
/*
 * Each branch stores the coordinates of the X,Y,Z subsequences it corresponds to
 * as well as the state path along the branch.
 * It also contains the indices for connecting branches at bifurcations in a parse tree.
 */
// IH, 8/7/2008:
// Having this inherit from SCFG_state_typing is slightly dangerous.
// Why? Because this means that there are member variables state_type[] and state_type_ancestral_map[] that are meaningless.
// Thus, this object can (for example) access state_type without compile-time errors, but the results will be undefined.
// It should inherit from SCFG_state_type_enum instead.
struct Triplet_SCFG_branch : SCFG_state_typing, Fold_char_enum
{

  /// coordinates of the subsequence of X represented by the subtree rooted at this branch
  int xl, xr; // xr = end + 1
  int yl, yr;
  int zl, zr;

  /// parse tree connections
  /*
   * Used to store indices in the parse tree Triplet_SCFG_parse_tree for the corresponding parent and child branches.
   */
  int parent;
  int lchild, rchild;

  /// state path along the branch
  vector<int> path;

  /// Constructor
  Triplet_SCFG_branch (int xl, int xr, int yl, int yr, int zl, int zr, int p);

  /// Does this branch have child branches?
  bool has_children() const { return lchild >= 0; }

  /// Return the alignment for the parse subtree passed as an argument.
  /*
   * The alignment is stored as the X, Y, Z foldstrings in an instance of Triplet_SCFG_alignment.
   * This is a recursive method which recurses down the subtree rooted at this branch and is used by Triplet_SCFG_parse_tree::alignment().
   * state_type is a vector giving the state type of each state (so it will typically be passed as 'my_Triplet_SCFG.state_type').
   */
  Triplet_SCFG_alignment alignment (const vector<Triplet_SCFG_branch>& parse_tree, const vector<State_type>& state_type, const State_type_ancestral_map& state_type_ancestral_map, const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z) const;

};

/// We represent a parse tree as a vector of Triplet_SCFG_branch'es.
/*
 * Connections are held within the branches themselves.
 */
struct Triplet_SCFG_parse_tree : vector<Triplet_SCFG_branch>, SCFG_state_type_enum
{

  /// Constructor.
  /*
   * Create a root branch with subsequence indices [0,len_x], etc., and no parent.
   */
  Triplet_SCFG_parse_tree (int len_x, int len_y, int len_z) { new_branch (0, len_x, 0, len_y, 0, len_z, -1); } // -1 b/c no parent

  /// Create a new branch and return the index of the new branch in the Triplet_SCFG_parse_tree.
  int new_branch (int xl, int xr, int yl, int yr, int zl, int zr, int p);

  /// Create a new left child branch of a given parent branch and return its index.
  int new_lchild (int xl, int xr, int yl, int yr, int zl, int zr, int p);

  /// Creates a new right child branch of a given parent branch and return its index.
  int new_rchild (int xl, int xr, int yl, int yr, int zl, int zr, int p);

  /// Test that all of the branches of the parse tree are properly linked together.
  bool test_connections() const;

  /// Test that the subsequences corresponding to the branches of the parse tree span the entire sequences.
  bool test_global (const Digitized_biosequence& dsq_x, const Digitized_biosequence& dsq_y, const Digitized_biosequence& dsq_z) const;

  /// Return alignment.
  Triplet_SCFG_alignment alignment (const vector<State_type>& state_type, const SCFG_state_typing::State_type_ancestral_map& state_type_ancestral_map, const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z) const;

  /// Display parse tree.
  /*
   * state_type, if defined, is used to display the state type strings
   */
  void show (ostream& o, const vector<State_type>* state_type = 0) const;

};

/// Holds a single alignment and associated parse_tree's.
struct Triplet_alignment_with_tracebacks : vector<Triplet_SCFG_parse_tree>, SCFG_state_type_enum
{

  const Triplet_SCFG_alignment alignment;
  
  /// Constructor
  Triplet_alignment_with_tracebacks (const Triplet_SCFG_alignment& alignment);

};


/// Representation of a triplet SCFG.
/*
 */
struct Triplet_SCFG : SCFG_state_typing, Grammar_state_enum
{

 public:
 

  /// transition matrix
  /*
   * An instance of the class 'Concrete_transition_scores' (see hmm/transmat.h) has 
   * all of the members/methods of the class 'Transition_matrix', plus a couple of show methods.
   *
   * transition_scores.transition (i,j) returns a reference to the score of transition from i to j
   * The Start state is (i = -1) and the End state is (j = -2).
   */
   //
   // If 'trans' is an instance of 'Concrete_sparse_transition_scores', then
   //    transition (i,j)
   // returns a reference to the (i,j) entry of the member '_matrix' (an instance of 'array2d') of 'trans'.
   // This is implemented as follows: The class 'Transition_matrix' has a method 'transition()', defined in hmm/transmat.h as:
   //    inline T& transition (int source, int dest)
   //      { return _matrix (tm_index(source), tm_index(dest)); }
   // where '_matrix' is a member instance of 'array2d'.  'array2d' has a method which defines a new operator, 
   //   reference operator()(int x, int y)
   //      { return _data.get_xy (x, y, _xsize, _ysize); }
  Concrete_transition_scores transition_scores;

  /// emissions and corresponding scores for each state
  /*
   * Each emit state of a particular type can represent all possible emissions
   * corresponding to its emission profile.
   * emit[s][K] == score implies that state 's' emits the string corresponding to (integer) K with an associated score 'score'.
   * The integer K is obtained with the helper functions emit_xl_mul(), etc.
   * See SCFG_state_typing::emit_idx() for an example.
   */
  vector<vector<Score> > emit;

  // constructors
  Triplet_SCFG();  /// default constructor
  Triplet_SCFG (int num_states);  /// initializes internal data but doesn't populate
  Triplet_SCFG (const Pair_CFG_scores& paircfg); /// initialize from a Pair_CFG_scores
  Triplet_SCFG (const Triplet_SCFG& scfg); /// copy constructor

  /// assignment operator
  Triplet_SCFG& operator= (const Triplet_SCFG& base);

  /// Initialize an emit state.
  /*
   * Takes the state type as an argument and initializes the (empty) emission vector emit, etc.
   */
  void init_emit (int state, State_type type);

  /// Initialize a bifurcation state.
  void init_bifurc (int state, int l, int r);

  /// display methods
  void show_element (const Score& element, ostream& o) const { ShowScore (element, o); } // ShowScore() is defined in util/score.h
 
  // to do: still need to fix this!
  /// Show the SCFG.
  void show (ostream& o) const;

  // to do: fix reliance on zero_outgoing_bifurcation_transitions() and make this a method of SCFG_state_sorter
  /// Get all destination (outgoing) non-bifurcation states for each source (incoming) state in selection.
  vector<vector<int> > selected_outgoing_states (const vector<int>& selection) const;

  /// Get all source (incoming) non-bifurcation states for each dest (outgoing) state in selection.
  vector<vector<int> > selected_incoming_states (const vector<int>& selection) const;

  /// Build the parse_tree corresponding to a given state path.
  /*
   * state_path corresponds to a depth-first, left->right traversal of the parse tree.
   * xstart is the subsequence index [xstart,] for sequence X
   * parse() is a method of Triplet_SCFG rather than SCFG_state_typing (as it is in paircfg.h) because, well, I feel like it
   */
  Triplet_SCFG_parse_tree parse (const vector<int>& state_path, int xstart, int ystart, int zstart) const;

  Triplet_SCFG_parse_tree parse (const vector<int>& state_path) const;

  /// Get the transition score associated with the parse tree.
  Score path_transition_score (const Triplet_SCFG_parse_tree& parse_tree) const;

  /// Get the emit score associated with the parse tree.
  Score path_emit_score (const Triplet_SCFG_parse_tree& parse_tree, const Digitized_biosequence& dsq_x, const Digitized_biosequence& dsq_y, const Digitized_biosequence& dsq_z) const;

  /// Get the total (transition and emit) score associated with the parse tree.
  Score path_score (const Triplet_SCFG_parse_tree& parse_tree, const Digitized_biosequence& dsq_x, const Digitized_biosequence& dsq_y, const Digitized_biosequence& dsq_z) const;

  // to do: get rid of this!
  /// Sets transitions from bifurcation states to -InfinityScore.
  void zero_outgoing_bifurcation_transitions();

  /// Adds two "fake" zero-score transitions from bifurc states, to fool topological sorter.
  void add_fake_bifurcation_transitions();

 private:
  /// Indent all subsequence indices xr of the branches of the subtree rooted at 'node' by 'offset'.
  /*
   * Helper method for Triplet_SCFG::parse().
   */
  static void indent_xr (Triplet_SCFG_parse_tree& parse_tree, int node, int offset);
  static void indent_yr (Triplet_SCFG_parse_tree& parse_tree, int node, int offset);
  static void indent_zr (Triplet_SCFG_parse_tree& parse_tree, int node, int offset);

};


// NB: Don't forget to omit the static keyword when defining the function
// elsewhere in the .h file!
inline int SCFG_state_typing::emit_idx (int t, int xl, int xr, int yl, int yr, int zl, int zr)
{
  int emit_idx = 0; // map emit string to ints for array access
  if (t & EmitXL) emit_idx += emit_xl_mul (t) * xl;
  if (t & EmitXR) emit_idx += emit_xr_mul (t) * xr;
  if (t & EmitYL) emit_idx += emit_yl_mul (t) * yl;
  if (t & EmitYR) emit_idx += emit_yr_mul (t) * yr;
  if (t & EmitZL) emit_idx += emit_zl_mul (t) * zl;
  if (t & EmitZR) emit_idx += emit_zr_mul (t) * zr;
  return emit_idx;
}

#endif /* TRIPLETSCFG_INCLUDED */
