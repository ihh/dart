#ifndef PAIRCFGDP_INCLUDED
#define PAIRCFGDP_INCLUDED

#include "scfg/paircfg.h"
#include "scfg/foldenv.h"
#include "scfg/subseqiter.h"
#include "util/strsaver.h"

// Uncomment the following line to use dense(fast,fat) rather than sparse(slow,slim) DP matrices
// (this is now controlled by the Makefile - ihh, 1/27/2005)
// #define DART_USE_DENSE_CFGDP_ALLOCATORS

// (should probably template this, make a factory class & let program auto-detect optimal memory strategy)
// Switched to sparse matrices by default, since Rfam training was aborting due to lack of memory (October 1, 2003 -- ihh)
#ifdef DART_USE_DENSE_CFGDP_ALLOCATORS
#define CFG_DP_MATRIX_BASE_CLASS Pair_CFG_DP_dense_allocator
#else /* DART_USE_DENSE_CFGDP_ALLOCATORS */
#define CFG_DP_MATRIX_BASE_CLASS Pair_CFG_DP_sparse_allocator
#endif /* DART_USE_DENSE_CFGDP_ALLOCATORS */

/// A class for figuring out the accessible DP matrix cells for a given subseq pair, given the fold envelope.
/*
 * The "flag_pair" is a combination of two subsequence connection flags,
 * one for the X-fold-envelope (call this FX) and one for the Y-fold-envelope (FY).
 * It is defined as (FX | (FY << 3)).
 * (Note CFLAG_L=1, CFLAG_R=2, CFLAG_LR=4.)
 * 
 * Thus the flag_pair summarizes the directions of the available 
 * connections from a given subsequence-pair (i.e. cell) in the 
 * foldenv-constrained DP matrix.
 * 
 * The Pair_CFG_filter thus answers questions like this:
 * 
 * "Given that I am at subsequence coords (SX,SY) and in nonterminal N, and 
 * I am allowed to emit any combination of the following...
 * 
 *   { a basepair to X,
 *     a single leftward base to X,
 *     a single rightward base to X, and/or
 *     a basepair to Y },
 * 
 * ...which nonterminals in adjacent cells can I reach?"
 */
struct Pair_CFG_filter : Pair_CFG_state_type_enum, Connection_enum
{
  // the following enum defines the range of the flag_pair index
  enum { FlagShift = 3, FlagRange = 8,        // FlagShift = N_CONN, FlagRange = 1 << FlagShift
	 FlagPairRange = 64 };  // FlagPairRange = FlagRange^2
  // lists of accessible states for inside & outside fills
  vector<vector<int> > allowed_in_states;  // indexed by [flag_pair], sorted by inside fill order (i.e. emit states first, then nonemit states sorted reverse-topologically)
  vector<vector<int> > allowed_out_states;  // indexed by [flag_pair], sorted by outside fill order (i.e. emit states first, then nonemit states sorted topologically)
  // lists of accessible dest/source states from/to a given source/dest state, on inside/outside fills
  vector<vector<vector<int> > > allowed_dest_states;  // indexed by [dest flag_pair] [source state]
  vector<vector<vector<int> > > allowed_src_states;   // indexed by [src flag_pair] [dest state]
  // bifurcation outside-lookup
  vector<vector<Bifurcation_left_parent> > left_parent;
  vector<vector<Bifurcation_right_parent> > right_parent;
  // constructor
  Pair_CFG_filter (const Pair_CFG_scores& cfg);
  // flag_pair index calculation
  static inline int flag_pair (int xflag, int yflag);
};

// allocator base class -- holds dimensions, envelopes etc
struct Pair_CFG_DP_allocator_base
{
  // fold envelopes
  const Fold_envelope& xenv;
  const Fold_envelope& yenv;
  // alignment envelope
  Pair_envelope pair_env;
  // complete sequence data
  const Named_profile& npx;
  const Named_profile& npy;
  // digitized sequences
  const Digitized_biosequence& xdsq;
  const Digitized_biosequence& ydsq;
  // metascores
  const vector<Metascore>& xmeta;
  const vector<Metascore>& ymeta;
  // CFG stuff
  const Pair_CFG_scores& cfg;
  const Pair_CFG_filter filter;
  // subseqs
  const int xsubseqs;
  const int ysubseqs;
  // states
  const int states;
  // sequence lengths
  const int xlen;
  const int ylen;
  // total number of cells and bifurcations
  unsigned long total_cells;
  unsigned long total_bifs;
  // approx no. of bytes allocated (lower bound)
  unsigned int bytes_allocated;
  // constructor
  Pair_CFG_DP_allocator_base (const Named_profile& npx,
			      const Named_profile& npy,
			      const Fold_envelope& xenv,
			      const Fold_envelope& yenv,
			      const Pair_envelope& pair_env,
			      const Pair_CFG_scores& cfg);
  // method to count cells in pair envelope
  void count_cells();  // sets total_cells
  // helpers
  inline bool startpos_in_pair_envelope (int xstart, int ystart) const;
  inline bool endpos_in_pair_envelope (int xend, int yend) const;
};

// sparse DP matrix container
// slower but slimmer
class Pair_CFG_DP_sparse_allocator : public Pair_CFG_DP_allocator_base
{
private:
  // private structs & typedefs
  struct By_start_indices { unsigned long offset; int xsize; By_start_indices() : offset(0), xsize(0) { } };
  typedef vector<Subseq>::const_iterator Subseq_iterator;
  // private data
  vector<Score> _cell;
  array2d<By_start_indices> by_start_indices;
  const Score _dummy_inf_score;  // initialised to -InfinityScore
  const Subseq_iterator xenv_begin;
  const Subseq_iterator yenv_begin;
  // public methods
public:
  // accessors
  inline Score& cell (int state, int xsubseq_idx, int ysubseq_idx);  // does no bounds-checking, returns a reference
  inline const Score& cell (int state, int xsubseq_idx, int ysubseq_idx) const;  // does bounds-checking, returns a reference
  inline Score read_cell (int state, int xsubseq_idx, int ysubseq_idx) const;  // does bounds-checking, doesn't return a reference
  // constructor
  Pair_CFG_DP_sparse_allocator (const Named_profile& npx,
				const Named_profile& npy,
				const Fold_envelope& xenv,
				const Fold_envelope& yenv,
				const Pair_envelope& pair_env,
				const Pair_CFG_scores& cfg);
  // builder
  void alloc();
};

// dense DP matrix container
// faster and fatter
class Pair_CFG_DP_dense_allocator : public Pair_CFG_DP_allocator_base
{
  // private data
private:
  vector<Score> _cell;
  // public methods
public:
  // accessors
  inline Score& cell (int state, int xsubseq_idx, int ysubseq_idx);
  inline const Score& cell (int state, int xsubseq_idx, int ysubseq_idx) const;
  inline Score read_cell (int state, int xsubseq_idx, int ysubseq_idx) const;
  // constructor
  Pair_CFG_DP_dense_allocator (const Named_profile& npx,
			       const Named_profile& npy,
			       const Fold_envelope& xenv,
			       const Fold_envelope& yenv,
			       const Pair_envelope& pair_env,
			       const Pair_CFG_scores& cfg);
  // builder
  void alloc();
};

// DP matrix class
class Pair_CFG_DP_matrix : public Pair_CFG_state_type_enum, public CFG_DP_MATRIX_BASE_CLASS, protected Stream_saver
{
  // typedefs
public:
  struct Cell_coords { int state, xsubseq, ysubseq; inline Cell_coords (int s, int x, int y); };

  // data
public:

  // local alignment flag (this is now deprecated; use Fold_envelope::make_local() method instead)
  const bool local;

  // results of DP
  Score final_score;  // "final" score (actually, score for start state, since DP goes in reverse)

protected:
  bool _show_out;  // for display only, passed to fold envelope

  // methods
public:
  // constructors
  Pair_CFG_DP_matrix (const Named_profile& npx,
		      const Named_profile& npy,
		      const Fold_envelope& xenv,
		      const Fold_envelope& yenv,
		      const Pair_CFG_scores& cfg,
		      const Pair_envelope& pair_env,
		      bool local = FALSE);

  // iterators
  inline Inside_subseq_iterator in_subseq_begin() const;
  inline Inside_subseq_iterator in_subseq_end() const;
  inline Outside_subseq_iterator out_subseq_begin() const;
  inline Outside_subseq_iterator out_subseq_end() const;

  // helpers
  inline bool local_start_allowed (const Subseq& xsubseq, const Subseq& ysubseq) const;  // checks local & out_flags
  inline bool in_pair_envelope (const Subseq& xsubseq, const Subseq& ysubseq) const;  // wraps pair_env method
  inline const vector<int>& allowed_in_states (const Subseq& xsubseq, const Subseq& ysubseq) const;
  inline const vector<int>& allowed_dest_in_states (int src_state, const Subseq& dest_xsubseq, const Subseq& dest_ysubseq) const;
  inline const vector<int>& allowed_out_states (const Subseq& xsubseq, const Subseq& ysubseq) const;
  inline const vector<int>& allowed_src_out_states (int dest_state, const Subseq& src_xsubseq, const Subseq& src_ysubseq) const;
  inline int dest_xsubseq_idx (int xsubseq_idx, State_type t) const;
  inline int dest_ysubseq_idx (int ysubseq_idx, State_type t) const;
  inline int src_xsubseq_idx (int xsubseq_idx, State_type t) const;
  inline int src_ysubseq_idx (int ysubseq_idx, State_type t) const;
  inline int in_emit_idx (State_type t, const Subseq& xsubseq, const Subseq& ysubseq) const;
  inline int out_emit_idx (State_type t, const Subseq& xsubseq, const Subseq& ysubseq) const;
  inline Score in_meta_sc (int state, const Subseq& xsubseq, const Subseq& ysubseq) const;
  inline Score out_meta_sc (int state, const Subseq& xsubseq, const Subseq& ysubseq) const;
  static inline Score calc_meta_sc (const vector<int>& meta_idx, const vector<Metascore>& meta_sc, int pos);

  // display
  void show (ostream& o) const;
  void show_expanded_trace (const Pair_CFG_parse_tree& parse_tree, ostream& out, bool outside = FALSE) const;
  sstring cfg_dump() const;  // for debugging; dumps CFG + matrix
};

// Summary statistics for an ensemble of parse trees: specifically, subseq-pair hit counts
struct Pair_CFG_parse_stats : Connection_enum, Pair_CFG_state_type_enum, Grammar_state_enum
{
  // typedefs
  typedef pair<int,int> Subseq_index_pair;
  typedef map<Subseq_index_pair,int> Subseq_index_pair_count;
  // data
  const Pair_CFG_DP_matrix& matrix;
  Subseq_index_pair_count pair_count;
  // constructor
  Pair_CFG_parse_stats (const Pair_CFG_DP_matrix& matrix);
  // method to add counts for a particular parse tree
  void count_subseqs (const Pair_CFG_parse_tree& parse_tree);
};

// In the Inside matrix, cell(STATE,XSUBSEQ,YSUBSEQ) is the score for all parse subtrees rooted in STATE,
// where the subtree from STATE downwards has emitted XSUBSEQ and YSUBSEQ (including the emission from STATE).
struct Pair_inside_matrix : Pair_CFG_DP_matrix
{
  // constructors
  Pair_inside_matrix (const Named_profile& npx,
		      const Named_profile& npy,
		      const Fold_envelope& xenv,
		      const Fold_envelope& yenv,
		      const Pair_CFG_scores& cfg,
		      bool local = FALSE,
		      bool fill_now = TRUE);

  Pair_inside_matrix (const Named_profile& npx,
		      const Named_profile& npy,
		      const Fold_envelope& xenv,
		      const Fold_envelope& yenv,
		      const Pair_CFG_scores& cfg,
		      const Pair_envelope& pair_env,
		      bool local = FALSE,
		      bool fill_now = TRUE);

  // fill method
  void fill();

  // sampled traceback methods
  void traceback_from (int xsubseq_idx, int ysubseq_idx, int traceback_state, vector<int>& path, bool prepend_Start_state = TRUE) const;  // samples a state path from specified start coords
};

// Pair_inside_matrix cell sorter, for fast traceback
struct Pair_inside_cell_sorter
{
  // typedefs
  typedef Pair_CFG_DP_matrix::Cell_coords Cell_coords;
  // data
  const Pair_inside_matrix& inside;
  vector<Cell_coords> sorted_cells;  // cells sorted by descending probability
  vector<Prob> sorted_cells_prob;  // probability for each cell in sorted_cells array
  Prob total_prob;  // sum of sorted_cells_prob; may be different to Score2Prob(inside.final_score) due to rounding error
  // constructor
  Pair_inside_cell_sorter (const Pair_inside_matrix& inside);
  // sampled traceback methods
  Cell_coords sample_coords() const;

  // various amplified traceback methods
  Pair_CFG_local_path traceback_with_coords() const;  // calls Pair_inside_matrix::traceback_from()
  Pair_CFG_parse_tree parse_tree() const;  // calls traceback_with_coords()
  Pair_CFG_alignment alignment() const;  // calls parse_tree()
};

// In the Outside matrix, cell(STATE,XSUBSEQ,YSUBSEQ) is the score for all parse sub-trees ending in STATE,
// where the subtree from STATE upwards has emitted everything outside XSUBSEQ and YSUBSEQ
// (so the emission from STATE is also outside XSUBSEQ and YSUBSEQ).
struct Pair_outside_matrix : Pair_CFG_DP_matrix
{
  // data
  const Pair_inside_matrix& inside;
  Pair_CFG_counts  count;
  vector<Metaprob> xmetacount;
  vector<Metaprob> ymetacount;
  // constructor
  Pair_outside_matrix (const Pair_inside_matrix& inside, bool fill_now = TRUE);
  void fill();
  // metacount accumulation methods
  static inline void count_meta (Prob post_prob, const vector<int>& meta_idx, vector<Metaprob>& meta_count, int pos);
  inline void out_count_meta (Prob post_prob, int state, const Subseq& xsubseq, const Subseq& ysubseq);
};

// Inside-Outside data structure
struct Pair_inside_outside_matrix
{
  // data
  Pair_inside_matrix      inside;
  Pair_outside_matrix     outside;
  const Pair_CFG_counts&  count;
  const vector<Metaprob>& xmetacount;
  const vector<Metaprob>& ymetacount;
  // constructors
  Pair_inside_outside_matrix (const Named_profile& npx,
			      const Named_profile& npy,
			      const Fold_envelope& xenv,
			      const Fold_envelope& yenv,
			      const Pair_CFG_scores& cfg,
			      bool local = FALSE,
			      bool fill_now = TRUE);

  Pair_inside_outside_matrix (const Named_profile& npx,
			      const Named_profile& npy,
			      const Fold_envelope& xenv,
			      const Fold_envelope& yenv,
			      const Pair_CFG_scores& cfg,
			      const Pair_envelope& pair_env,
			      bool local = FALSE,
			      bool fill_now = TRUE);
  // display
  void show (ostream& o) const;

  // Posterior score calculation methods.
  // These give the posterior probability that a parse tree contains a subtree rooted in state DEST,
  // where the subtree from DEST downwards has emitted XSUBSEQ and YSUBSEQ (including the emission from DEST).
  // (In other words, the subseq-indexing convention here follows Inside/CYK, not Outside/KYC.)
  // In the case of post_transition_sc(), it's also required that the parent node of DEST must be SRC.
  Score post_transition_sc (int src_state, int dest_state, int xsubseq_idx, int ysubseq_idx) const;
  Score post_state_sc (int dest_state, int xsubseq_idx, int ysubseq_idx) const;  // = psum_{src} post_transition_sc(src,...)
};

// In the CYK matrix, cell(STATE,XSUBSEQ,YSUBSEQ) is the max score for any parse subtree rooted in STATE,
// where the subtree from STATE downwards has emitted XSUBSEQ and YSUBSEQ (including the emission from STATE).
struct Pair_CYK_matrix : Pair_CFG_DP_matrix
{
  // traceback start coords
  int final_xsubseq_idx;  // x coord of final cell, equal to (xsubseqs-1) for global matrices
  int final_ysubseq_idx;  // y coord of final cell, equal to (ysubseqs-1) for global matrices
  int final_state;  // state of final cell

  // constructors
  Pair_CYK_matrix (const Named_profile& npx,
		   const Named_profile& npy,
		   const Fold_envelope& xenv,
		   const Fold_envelope& yenv,
		   const Pair_CFG_scores& cfg,
		   bool local = FALSE,
		   bool fill_now = TRUE);

  Pair_CYK_matrix (const Named_profile& npx,
		   const Named_profile& npy,
		   const Fold_envelope& xenv,
		   const Fold_envelope& yenv,
		   const Pair_CFG_scores& cfg,
		   const Pair_envelope& pair_env,
		   bool local = FALSE,
		   bool fill_now = TRUE);

  void fill();
  // traceback method
  // if prepend_Start_state is set (which it is by most of the helper methods),
  // then a 'Start' state is prepended to the beginning of the path.
  void traceback_from (int xsubseq_idx, int ysubseq_idx, int traceback_state, vector<int>& path, bool prepend_Start_state = TRUE) const;

  // various amplified traceback methods
  vector<int> traceback() const;  // calls traceback_from(); useless without coords, unless matrix is global
  Pair_CFG_local_path traceback_with_coords() const;  // calls traceback()
  Pair_CFG_parse_tree parse_tree() const;  // calls traceback_with_coords()
  Pair_CFG_alignment alignment() const;  // calls parse_tree()
};

// The KYC matrix is to the CYK matrix as Outside is to Inside.
// In the KYC matrix, cell(STATE,XSUBSEQ,YSUBSEQ) is the max score for any parse subtree ending in STATE,
// where the subtree from STATE upwards has emitted everything outside XSUBSEQ and YSUBSEQ
// (so the emission from STATE is also outside XSUBSEQ and YSUBSEQ).
struct Pair_KYC_matrix : Pair_CFG_DP_matrix
{
  // data
  const Pair_CYK_matrix& cyk;
  // constructors
  Pair_KYC_matrix (const Pair_CYK_matrix& cyk, bool fill_now = TRUE);
  void fill();

  // Traceback method.
  // This yields the highest-scoring parse tree containing a subtree rooted in state DEST,
  // where the subtree from DEST downwards has emitted XSUBSEQ and YSUBSEQ (including the emission from DEST).
  // (In other words, the subseq-indexing convention here follows Inside/CYK, not Outside/KYC.)
  // In the case of traceback_from_transition(), the parent of DEST in the full parse tree must be SRC.
  void traceback_from (int xsubseq_idx, int ysubseq_idx, int dest_state, vector<int>& path) const;
  void traceback_from_transition (int xsubseq_idx, int ysubseq_idx, int dest_state, int src_state, vector<int>& path) const;
};

// CYK-KYC data structure
struct Pair_CYK_KYC_matrix
{
  // data
  Pair_CYK_matrix cyk;
  Pair_KYC_matrix kyc;
  // constructors
  Pair_CYK_KYC_matrix (const Named_profile& npx,
		       const Named_profile& npy,
		       const Fold_envelope& xenv,
		       const Fold_envelope& yenv,
		       const Pair_CFG_scores& cfg,
		       bool local = FALSE,
		       bool fill_now = TRUE);

  Pair_CYK_KYC_matrix (const Named_profile& npx,
		       const Named_profile& npy,
		       const Fold_envelope& xenv,
		       const Fold_envelope& yenv,
		       const Pair_CFG_scores& cfg,
		       const Pair_envelope& pair_env,
		       bool local = FALSE,
		       bool fill_now = TRUE);
  // display
  void show (ostream& o) const;

  // Max-score-through-cell methods.
  // These give the max score for any parse tree containing a subtree rooted in state DEST,
  // where the subtree from DEST downwards has emitted XSUBSEQ and YSUBSEQ (including the emission from DEST).
  // (In other words, the subseq-indexing convention here follows Inside/CYK, not Outside/KYC.)
  // In the case of max_transition_sc(), it's also required that the parent node of DEST must be SRC.
  Score max_transition_sc (int src_state, int dest_state, int xsubseq_idx, int ysubseq_idx) const;
  Score max_state_sc (int dest_state, int xsubseq_idx, int ysubseq_idx) const;  // = max_{src} max_transition_sc(src,...)
};

// CYK cell sorter.
struct Pair_CYK_KYC_cell_sorter
{
  // typedefs
  typedef Pair_CFG_DP_matrix::Cell_coords Cell_coords;
  // data
  const Pair_CYK_KYC_matrix& cyk_kyc;   // CYK-KYC matrix
  vector<Cell_coords>        sorted_cells;     // cells sorted by descending max traceback score
  // constructor
  Pair_CYK_KYC_cell_sorter (const Pair_CYK_KYC_matrix& cyk_kyc);
};

// inline method defs

// method to convert a pair of Connection_enum::Connection_flag's into a Pair_CFG_filter flag_pair index
inline int Pair_CFG_filter::flag_pair (int xflag, int yflag)
{ return xflag | (yflag << FlagShift); }

// Pair_CFG_DP_sparse_allocator methods: cell accessors
// note that mutable accessor doesn't check bounds, since mutable cells *must* be in the array
// to force const behaviour, use read_cell
inline Score& Pair_CFG_DP_sparse_allocator::cell (int state, int xsubseq_idx, int ysubseq_idx)
{
  const Subseq& xsubseq = xenv_begin[xsubseq_idx];
  const Subseq& ysubseq = yenv_begin[ysubseq_idx];
  const By_start_indices& ind = by_start_indices (xsubseq.start, ysubseq.start);
  return _cell [ind.offset + xsubseq.by_start_index + ind.xsize * (state + states * ysubseq.by_start_index)];
}

inline const Score& Pair_CFG_DP_sparse_allocator::cell (int state, int xsubseq_idx, int ysubseq_idx) const
{
  const Subseq& xsubseq = xenv_begin[xsubseq_idx];
  const Subseq& ysubseq = yenv_begin[ysubseq_idx];
  const By_start_indices& ind = by_start_indices (xsubseq.start, ysubseq.start);
  const int xs = ind.xsize;
  if (xs == 0) return _dummy_inf_score;
  return _cell [ind.offset + xsubseq.by_start_index + xs * (state + states * ysubseq.by_start_index)];
}

inline Score Pair_CFG_DP_sparse_allocator::read_cell (int state, int xsubseq_idx, int ysubseq_idx) const
{
  const Subseq& xsubseq = xenv_begin[xsubseq_idx];
  const Subseq& ysubseq = yenv_begin[ysubseq_idx];
  const By_start_indices& ind = by_start_indices (xsubseq.start, ysubseq.start);
  const int xs = ind.xsize;
  if (xs == 0) return -InfinityScore;
  return _cell [ind.offset + xsubseq.by_start_index + xs * (state + states * ysubseq.by_start_index)];
}

// Pair_CFG_DP_dense_allocator methods: cell accessors
inline Score& Pair_CFG_DP_dense_allocator::cell (int state, int xsubseq_idx, int ysubseq_idx)
{ return _cell [state + states * (xsubseq_idx + ysubseq_idx * xsubseqs)]; }

inline const Score& Pair_CFG_DP_dense_allocator::cell (int state, int xsubseq_idx, int ysubseq_idx) const
{ return _cell [state + states * (xsubseq_idx + ysubseq_idx * xsubseqs)]; }

inline Score Pair_CFG_DP_dense_allocator::read_cell (int state, int xsubseq_idx, int ysubseq_idx) const
{ return _cell [state + states * (xsubseq_idx + ysubseq_idx * xsubseqs)]; }

// Pair_CFG_DP_matrix::Cell_coords constructor
Pair_CFG_DP_matrix::Cell_coords::Cell_coords (int s, int x, int y)
  : state (s), xsubseq (x), ysubseq (y)
{ }

// Pair_CFG_DP_matrix methods
// iterators
inline Inside_subseq_iterator Pair_CFG_DP_matrix::in_subseq_begin() const
{ Inside_subseq_iterator iter (xenv, yenv, pair_env, 0, -1); return ++iter; }

inline Inside_subseq_iterator Pair_CFG_DP_matrix::in_subseq_end() const
{ return Inside_subseq_iterator (xenv, yenv, pair_env, xsubseqs, 0); }

inline Outside_subseq_iterator Pair_CFG_DP_matrix::out_subseq_begin() const
{ Outside_subseq_iterator iter (xenv, yenv, pair_env, xsubseqs, 0); return ++iter; }

inline Outside_subseq_iterator Pair_CFG_DP_matrix::out_subseq_end() const
{ return Outside_subseq_iterator (xenv, yenv, pair_env, 0, -1); }

// method to check if co-ords are in the pair envelope
inline bool Pair_CFG_DP_matrix::local_start_allowed (const Subseq& xsubseq, const Subseq& ysubseq) const
{ return local || (xsubseq.start_flag && ysubseq.start_flag); }

inline bool Pair_CFG_DP_matrix::in_pair_envelope (const Subseq& xsubseq, const Subseq& ysubseq) const
{ return Subseq_coords::allowed_by_pair_env (pair_env, xsubseq, ysubseq); }

inline bool Pair_CFG_DP_allocator_base::startpos_in_pair_envelope (int xstart, int ystart) const
{ return Subseq_coords::startpos_allowed_by_pair_env (pair_env, xstart, ystart); }

inline bool Pair_CFG_DP_allocator_base::endpos_in_pair_envelope (int xend, int yend) const
{ return Subseq_coords::endpos_allowed_by_pair_env (pair_env, xend, yend); }

// methods to check if given states are allowed at a particular position
inline const vector<int>& Pair_CFG_DP_matrix::allowed_in_states (const Subseq& xsubseq, const Subseq& ysubseq) const
{ return filter.allowed_in_states [Pair_CFG_filter::flag_pair (xsubseq.in_flags, ysubseq.in_flags)]; }

inline const vector<int>& Pair_CFG_DP_matrix::allowed_dest_in_states (int src_state, const Subseq& dest_xsubseq, const Subseq& dest_ysubseq) const
{ return filter.allowed_dest_states [Pair_CFG_filter::flag_pair (dest_xsubseq.in_flags, dest_ysubseq.in_flags)] [src_state]; }

inline const vector<int>& Pair_CFG_DP_matrix::allowed_out_states (const Subseq& xsubseq, const Subseq& ysubseq) const
{ return filter.allowed_out_states [Pair_CFG_filter::flag_pair (xsubseq.out_flags, ysubseq.out_flags)]; }

inline const vector<int>& Pair_CFG_DP_matrix::allowed_src_out_states (int dest_state, const Subseq& src_xsubseq, const Subseq& src_ysubseq) const
{ return filter.allowed_src_states [Pair_CFG_filter::flag_pair (src_xsubseq.out_flags, src_ysubseq.out_flags)] [dest_state]; }

// methods to find neighbour cells
inline int Pair_CFG_DP_matrix::dest_xsubseq_idx (int xsubseq_idx, State_type t) const
{
  switch (t & EmitXLR)
    {
    case EmitXL:  return xenv.subseq[xsubseq_idx].in_l(); break;
    case EmitXR:  return xenv.subseq[xsubseq_idx].in_r(); break;
    case EmitXLR: return xenv.subseq[xsubseq_idx].in_lr(); break;
    default: return xsubseq_idx; break;
    }
  return -1;  // unreachable
}

inline int Pair_CFG_DP_matrix::dest_ysubseq_idx (int ysubseq_idx, State_type t) const
{
  switch (t & EmitYLR)
    {
    case EmitYL:  return yenv.subseq[ysubseq_idx].in_l(); break;
    case EmitYR:  return yenv.subseq[ysubseq_idx].in_r(); break;
    case EmitYLR: return yenv.subseq[ysubseq_idx].in_lr(); break;
    default: return ysubseq_idx; break;
    }
  return -1;  // unreachable
}

inline int Pair_CFG_DP_matrix::src_xsubseq_idx (int xsubseq_idx, State_type t) const
{
  switch (t & EmitXLR)
    {
    case EmitXL:  return xenv.subseq[xsubseq_idx].out_l(); break;
    case EmitXR:  return xenv.subseq[xsubseq_idx].out_r(); break;
    case EmitXLR: return xenv.subseq[xsubseq_idx].out_lr(); break;
    default: return xsubseq_idx; break;
    }
  return -1;  // unreachable
}

inline int Pair_CFG_DP_matrix::src_ysubseq_idx (int ysubseq_idx, State_type t) const
{
  switch (t & EmitYLR)
    {
    case EmitYL:  return yenv.subseq[ysubseq_idx].out_l(); break;
    case EmitYR:  return yenv.subseq[ysubseq_idx].out_r(); break;
    case EmitYLR: return yenv.subseq[ysubseq_idx].out_lr(); break;
    default: return ysubseq_idx; break;
    }
  return -1;  // unreachable
}

// methods to find emit indices for a cell
inline int Pair_CFG_DP_matrix::in_emit_idx (State_type t, const Subseq& xsubseq, const Subseq& ysubseq) const
{
  int emit_idx = 0;
  if (t & EmitXL) emit_idx += cfg.emit_xl_mul(t) * xdsq [xsubseq.start];
  if (t & EmitXR) emit_idx += cfg.emit_xr_mul(t) * xdsq [xsubseq.start + xsubseq.len - 1];
  if (t & EmitYL) emit_idx += cfg.emit_yl_mul(t) * ydsq [ysubseq.start];
  if (t & EmitYR) emit_idx += cfg.emit_yr_mul(t) * ydsq [ysubseq.start + ysubseq.len - 1];
  return emit_idx;
}

inline int Pair_CFG_DP_matrix::out_emit_idx (State_type t, const Subseq& xsubseq, const Subseq& ysubseq) const
{
  int emit_idx = 0;
  if (t & EmitXL) emit_idx += cfg.emit_xl_mul(t) * xdsq [xsubseq.start - 1];
  if (t & EmitXR) emit_idx += cfg.emit_xr_mul(t) * xdsq [xsubseq.start + xsubseq.len];
  if (t & EmitYL) emit_idx += cfg.emit_yl_mul(t) * ydsq [ysubseq.start - 1];
  if (t & EmitYR) emit_idx += cfg.emit_yr_mul(t) * ydsq [ysubseq.start + ysubseq.len];
  return emit_idx;
}

// metascore methods
inline Score Pair_CFG_DP_matrix::in_meta_sc (int state, const Subseq& xsubseq, const Subseq& ysubseq) const
{
  const State_type t = cfg.state_type[state];
  const vector<int>& xlmeta_idx = cfg.xlmeta_idx[state];
  const vector<int>& xrmeta_idx = cfg.xrmeta_idx[state];
  const vector<int>& ylmeta_idx = cfg.ylmeta_idx[state];
  const vector<int>& yrmeta_idx = cfg.yrmeta_idx[state];
  Score sc = 0;
  if (t & EmitXL && xlmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (xlmeta_idx, xmeta, xsubseq.start));
  if (t & EmitXR && xrmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (xrmeta_idx, xmeta, xsubseq.start + xsubseq.len - 1));
  if (t & EmitYL && ylmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (ylmeta_idx, ymeta, ysubseq.start));
  if (t & EmitYR && yrmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (yrmeta_idx, ymeta, ysubseq.start + ysubseq.len - 1));
  return sc;
}

inline Score Pair_CFG_DP_matrix::out_meta_sc (int state, const Subseq& xsubseq, const Subseq& ysubseq) const
{
  const State_type t = cfg.state_type[state];
  const vector<int>& xlmeta_idx = cfg.xlmeta_idx[state];
  const vector<int>& xrmeta_idx = cfg.xrmeta_idx[state];
  const vector<int>& ylmeta_idx = cfg.ylmeta_idx[state];
  const vector<int>& yrmeta_idx = cfg.yrmeta_idx[state];
  Score sc = 0;
  if (t & EmitXL && xlmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (xlmeta_idx, xmeta, xsubseq.start - 1));
  if (t & EmitXR && xrmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (xrmeta_idx, xmeta, xsubseq.start + xsubseq.len));
  if (t & EmitYL && ylmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (ylmeta_idx, ymeta, ysubseq.start - 1));
  if (t & EmitYR && yrmeta_idx.size()) ScorePMulAcc (sc, calc_meta_sc (yrmeta_idx, ymeta, ysubseq.start + ysubseq.len));
  return sc;
}

inline Score Pair_CFG_DP_matrix::calc_meta_sc (const vector<int>& meta_idx, const vector<Metascore>& meta_sc, int pos)
{
  Score sc = 0;
  for_const_contents (vector<int>, meta_idx, mi)
    ScorePMulAcc (sc, meta_sc[*mi][pos]);
  return sc;
}

// Pair_outside_matrix
inline void Pair_outside_matrix::out_count_meta (Prob post_prob, int state, const Subseq& xsubseq, const Subseq& ysubseq)
{
  const State_type t = cfg.state_type[state];
  const vector<int>& xlmeta_idx = cfg.xlmeta_idx[state];
  const vector<int>& xrmeta_idx = cfg.xrmeta_idx[state];
  const vector<int>& ylmeta_idx = cfg.ylmeta_idx[state];
  const vector<int>& yrmeta_idx = cfg.yrmeta_idx[state];
  if (t & EmitXL && xlmeta_idx.size()) count_meta (post_prob, xlmeta_idx, xmetacount, xsubseq.start - 1);
  if (t & EmitXR && xrmeta_idx.size()) count_meta (post_prob, xrmeta_idx, xmetacount, xsubseq.start + xsubseq.len);
  if (t & EmitYL && ylmeta_idx.size()) count_meta (post_prob, ylmeta_idx, ymetacount, ysubseq.start - 1);
  if (t & EmitYR && yrmeta_idx.size()) count_meta (post_prob, yrmeta_idx, ymetacount, ysubseq.start + ysubseq.len);
}

inline void Pair_outside_matrix::count_meta (Prob post_prob, const vector<int>& meta_idx, vector<Metaprob>& meta_count, int pos)
{
  for_const_contents (vector<int>, meta_idx, mi)
    meta_count[*mi][pos] += post_prob;
}

#endif /* PAIRCFGDP_INCLUDED */
