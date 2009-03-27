#ifndef TRIPLETSCFGDP_INCLUDED
#define TRIPLETSCFGDP_INCLUDED

#include "scfg/foldenv.h"
#include "indiegram/tripletscfg.h"
#include "indiegram/statesorter.h"

// Choose whether to use dense (faster, fatter) or sparse (slower, slimmer) DP matrices.
// Uncomment the following line to use dense matrices.
// Note: Using sparse matrices fails for long (well, medium-length) sequences with an error
//   terminate called after throwing an instance of 'std::length_error'
//    what():  vector::_M_fill_insert
// presumably because we tried to create a vector whose length exceeded its maximum permitted size.

#define INDIEGRAM_USE_DENSE_DP

#ifdef INDIEGRAM_USE_DENSE_DP
#define TRIPLET_DP_MATRIX_BASE Triplet_DP_dense_allocator
#else /* INDIEGRAM_USE_DENSE_DP */
#define TRIPLET_DP_MATRIX_BASE Triplet_DP_sparse_allocator
#endif /* INDIEGRAM_USE_DENSE_DP */


// NB: Throughout the DP routines, we use a formulation where the "from" or "source" state
// (the state on the left-hand side of the production rule) emits symbols
// rather than the "to" or "destination" state.
// This matches up with Ian's formulation in e.g. Table 1 of Holmes & Rubin.
// I use the opposite formulation (where the destination state emits symbols) in my notes.

/// Class for figuring out which DP matrix cells are accessible for a given subseq triplet, given the fold envelopes.
/*
 * The "flag_triplet" is a combination of three subsequence connection flags,
 * one for the X fold envelope, one for the Y fold envelope and one for the Z fold envelope.
 * (Note CFLAG_L=1, CFLAG_R=2, CFLAG_LR=4.)
 * 
 * Thus the flag_pair summarizes the directions of the available 
 * connections from a given subsequence-triplet (i.e. cell) in the 
 * foldenv-constrained DP matrix.
 * 
 * The Triplet_SCFG_filter thus answers questions like this:
 * 
 * "Given that I am at subsequence coords (SX,SY,SZ) and in nonterminal N, and 
 * I am allowed to emit any combination of the following...
 * 
 *   { a basepair to X,
 *     a single leftward base to X,
 *     a single rightward base to X,
 *     a basepair to Y,
 *     and/or a basepair to Z}
 * 
 * ...which nonterminals in adjacent cells can I reach?"
 */
struct Triplet_SCFG_filter : SCFG_state_type_enum
{

  enum { FlagShift = 3, FlagRange = 8, // FlagShift = N_CONN, FlagRange = 1 << FlagShift
	 FlagTripletRange = 512 }; // FlagPairRange = FlagRange^3

  /// Vector of states (in Inside order) whose emissions correspond to a given connection flag triplet.
  /*
   * Sorted in reverse topological order (emit states last) for use in the Inside/CYK iteration.
   * For example: if the subsequence triplet has an X connection flag CFLAG_LR, then 
   * only emit states which LR-emit to X will appear in the corresponding vector
   * of allowed states.
   * Indexed as allowed_in_states[flag_triplet].
   */
  vector<vector<int> > allowed_in_states;

  /// Vector of states (in Outside order) whose emissions correspond to a given connection flag triplet.
  /*
   * Sorted in topological order for use in the Outside/KYC iteration.
   * For example: if the subsequence triplet has an X connection flag CFLAG_LR, then 
   * only emit states which LR-emit to X will appear in the corresponding vector
   * of allowed states.
   * Indexed as allowed_out_states[flag_triplet].
   */
  vector<vector<int> > allowed_out_states;

  /// Given a "source" state and "destination" connection flag triplet, get the vector of all possible destination states.
  /* 
   * Get vector of "destination" states which "source" state src_state has transitions into.
   * Allowed destination states are constrained by the specified connection flags.
   * Indexed as allowed_dest_states[destination flag_triplet][source state].
   * Used during the Inside/CYK fill.
   */
  vector<vector<vector<int> > > allowed_dest_states;

  /// Given a "dest" state and "source" connection flag triplet, get the vector of all possible source states.
  /* 
   * Get vector of "source" states which "dest" state dest_state has transitions from.
   * Allowed source states are constrained by the specified connection flags.
   * Indexed as allowed_src_states[source flag_triplet][dest state].
   * Used during the Outside fill.
   */
  vector<vector<vector<int> > > allowed_src_states;

  // bifurcation outside-lookup
  vector<vector<Bifurcation_left_parent> > left_parent;
  vector<vector<Bifurcation_right_parent> > right_parent;

  /// Constructor
  Triplet_SCFG_filter (const Triplet_SCFG& scfg, const SCFG_state_sorter& sorter);

  /// Convert a triplet of Connection_enum::Connection_flag's into a Triplet_SCFG_filter flag_triplet index.
  static inline int flag_triplet (int xflag, int yflag, int zflag);

};


// We currently only allow fold constraints, so we don't have an alignment envelope member.
/// Base class for Triplet_DP_allocator variants.
struct Triplet_DP_allocator_base
{

  /// Triplet SCFG
  const Triplet_SCFG& scfg;

  /// SCFG state sorter
  const SCFG_state_sorter sorter; // NB: Can't be a reference because will always be local to this object!

  /// Named_profiles for sequence data
  /**
   * A Named_profile holds name, "cruft" (i.e. the text following the name in FASTA files),
   * sequence data (in text, digitized, score and/or probability form) and Metascores.
   * Named_profile defined in seq/biosequence.h:
   */
  const Named_profile& np_x;
  const Named_profile& np_y;
  const Named_profile& np_z;

  /// Fold envelopes
  /*
   * These must be global; local alignment isn't implemented.
   */
  const Fold_envelope& foldenv_x;
  const Fold_envelope& foldenv_y;
  const Fold_envelope& foldenv_z;

  /// Digitized_biosequence's
  const Digitized_biosequence& dsq_x;
  const Digitized_biosequence& dsq_y;
  const Digitized_biosequence& dsq_z;
  
  /// Sequence lengths
  const int seqlen_x, seqlen_y, seqlen_z;  
  
  /// Filter
  const Triplet_SCFG_filter filter;

  /// Approx memory usage (lower bound)
  unsigned int megabytes_allocated;

  /// Number of cells allocated.
  unsigned long num_cells;

  /// Number of bifurcations.
  unsigned long num_bifs;


  /// Constructor.
  /*
   * Just copies in references to input vars.
   */
  Triplet_DP_allocator_base (const Triplet_SCFG& scfg,
			     const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
			     const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z);

  /// Count cells allocated and total number of bifurcations allowed.
  void count_cells();


};

/// Allocates (# states) * (# X subseqs) * (# X subseqs) * (# Z subseqs).
struct Triplet_DP_dense_allocator : public Triplet_DP_allocator_base
{

  Triplet_DP_dense_allocator (const Triplet_SCFG& scfg,
			      const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
			      const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z);


  /// Alloc method: sets up _cell[][][][]
  void alloc(); 
  /// cell accessor
  inline Score& cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z);
  /// Score accessor
  inline Score read_cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const;

 private:

  /// DP matrix
  typedef vector<Score> StateVec;
  typedef vector<StateVec> ZVec;
  typedef vector<ZVec> YVec;
  typedef vector<YVec> XVec;
  //  vector <vector <vector <vector <Score> > > > _cell;  // _cell[x][y][z][state] has type Score
  XVec _cell;  /// cell[x][y][z][state] has type Score

};

/// Sparser memory allocation.
/*
 * If a constraining Triplet_alignment is specified, then only allocates memory for valid positions.
 */
struct Triplet_DP_sparse_allocator : public Triplet_DP_allocator_base
{

 public:

  void alloc();
  inline Score& cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z);
  inline Score read_cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const;

  Triplet_DP_sparse_allocator (const Triplet_SCFG& scfg,
			       const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
			       const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z);

 private:

  // For each triplet of start indices (x,y,z) of subseqs of X,Y,Z, we need
  // a way to find the appropriate place in the vector<Score> _cell.  'offset' specifies
  // this.
  struct By_start_indices {
    unsigned long offset; // specifies the number of previous entries in _cell before the ones for
    int xsize; // the number of x subseqs for this 
    int zsize;
  By_start_indices() : offset (0), xsize (0), zsize (0) { }
  };


  // We need a nested structure because we have 3 seqs, not 2.  It goes like vector<array2d<By_start_indices> >,
  // structured for access as:
  // by_start_indices[zsubseq.start] (xsubseq.start, ysubseq.start)
  // The slightly weird indexing (z, then x and y) is done for ease of comparison with how 
  // it's done in scfg/paircfgdp.h.
  vector<array2d<By_start_indices> > by_start_indices;

  typedef vector<Subseq>::const_iterator Subseq_iterator;
  const Subseq_iterator foldenv_x_begin;
  const Subseq_iterator foldenv_y_begin;
  const Subseq_iterator foldenv_z_begin;

  vector<Score> _cell;

};


/// Base class for DP matrices.
struct Triplet_DP_matrix_base : public SCFG_state_type_enum, public TRIPLET_DP_MATRIX_BASE, protected Stream_saver
{

  /// Results of DP.
  /*
   * Final score for Start state.
   */
  Score final_sc;

  /// Constructor.
  /*
   * Just copies in references to input vars.
   * Note that Ian allows a (non-full) alignment envelope to be used as constraints via
   * Pair_CFG_DP_allocator_base::pair_env.
   * We currently only allow fold constraints, so we don't have an alignment envelope member.
   */
  Triplet_DP_matrix_base (const Triplet_SCFG& scfg,
			  const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
			  const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z);

  /// Show DP matrix.
  void show (ostream& o) const;

  /// Show DP matrix.
  void show_compact (ostream& o) const;

  /// For debugging: show scfg, fold envelopes and DP matrix.
  sstring scfg_dump() const;

  /// Given a subsequence triplet, get all states whose production rules match the triplet's (inside) connection flags.
  /*
   * Relies on Triplet_SCFG_filter::allowed_in_states.
   */
  inline const vector<int>& allowed_in_states (const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const;

  inline const vector<int>& allowed_out_states (const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const;


  /// Given a "source" state and a "destination" subsequence triplet, get all possible "destination" states.
  /*
   * For use in the Inside/CYK recursion.
   * Relies on Triplet_SCFG_filter::allowed_dest_states.
   */
  inline const vector<int>& allowed_dest_in_states (int src_state, const Subseq& dest_subseq_x, const Subseq& dest_subseq_y, const Subseq& dest_subseq_z) const;

  inline const vector<int>& allowed_src_out_states (int dest_state, const Subseq& src_subseq_x, const Subseq& src_subseq_y, const Subseq& src_subseq_z) const;



  /// Get the index for the (inside) subsequence connected to this subsequence (*this) by an emission.
  /*
   * t is the state type of the "source" state (in our formulation the source, not destination, state emits).
   * When using for e.g. the CYK fill, subseq_idx_x is the index for the X subseq of a cell and
   * t the state type of the state label for that cell.
   */
  inline int dest_subseq_idx_x (int subseq_idx_x, State_type t) const;
  inline int dest_subseq_idx_y (int subseq_idx_y, State_type t) const;
  inline int dest_subseq_idx_z (int subseq_idx_z, State_type t) const;

  inline int src_subseq_idx_x (int subseq_idx_x, State_type t) const;
  inline int src_subseq_idx_y (int subseq_idx_y, State_type t) const;
  inline int src_subseq_idx_z (int subseq_idx_z, State_type t) const;



  /// Helper to get the hash for a particular emission.
  /*
   * Maps an emit string to an int for Triplet_SCFG::emit[s][] array access.
   * Returns the emit_idx of the emission (xl yl zl zr yr xr), where the 
   * characters xl, etc. are read from the start and end of the Subseq's
   * (this is where the 'in' bit comes in).
   */
  inline int in_emit_idx (State_type t, const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const;

  inline int out_emit_idx (State_type t, const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const;

  /// Is a start transition allowed for this triplet of subsequences?
  /*
   * We currently only allow global alignment, so only the outermost subsequences have a start transition.
   */
  inline bool local_start_allowed (const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const;

 protected:

  /// Bool for display: show outside v. inside connection flags.
  /*
   * Set to false for the Inside and CYK matrices and true for the Outside matrix.
   */
  bool _show_out;

};


/// Inside DP matrix
struct Triplet_inside_matrix : Triplet_DP_matrix_base
{
  
  /// Constructor.
  /*
   * Just copies in references to input vars.
   */
  Triplet_inside_matrix (const Triplet_SCFG& scfg,
			 const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
			 const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z,
			 bool fill_now = true);

  /// fill method
  void fill();

};

// In the Outside matrix, cell(STATE,XSUBSEQ,YSUBSEQ) is the score for all parse sub-trees ending in STATE,
// where the subtree from STATE upwards has emitted everything outside XSUBSEQ and YSUBSEQ
// (so the emission from STATE is also outside XSUBSEQ and YSUBSEQ).
struct Triplet_outside_matrix : Triplet_DP_matrix_base
{
  // data
  const Triplet_inside_matrix& inside;

  // constructor
  Triplet_outside_matrix (const Triplet_inside_matrix& inside, bool fill_now = TRUE);

  /// fill method
  void fill();

};


/// CYK DP matrix
struct Triplet_CYK_matrix : Triplet_DP_matrix_base
{
  /// traceback start coordinates
  int final_subseq_idx_x; // x coordinate of final cell (should be foldenv_x.subseq.size()-1 for global matrices, which mine are)
  int final_subseq_idx_y;
  int final_subseq_idx_z;
  int final_state; /// the child state of Start in the CYK parse tree

  /// Constructor.
  /*
   * Just copies in references to input vars.
   */
  Triplet_CYK_matrix (const Triplet_SCFG& scfg,
		      const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
		      const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z,
		      bool fill_now = true);

  /// fill method
  void fill();

  /// traceback methods
  /// Perform CYK-traceback.
  /*
   * Returns the state path corresponding to a preorder traversal of the parse tree.
   */
  vector<int> traceback (bool prepend_Start_state = true) const;

  /// Perform traceback and return corresponding parse tree.
  /*
   * Relies on traceback ().
   */
  Triplet_SCFG_parse_tree parse_tree() const;
  
  /// Perform traceback and return corresponding alignment.
  /*
   * Relies on parse_tree().
   */
  Triplet_SCFG_alignment alignment() const;

};


// inline functions

inline int Triplet_SCFG_filter::flag_triplet (int xflag, int yflag, int zflag)
{ return xflag | (yflag << FlagShift) | ((zflag << FlagShift) << FlagShift); }

inline const vector<int>& Triplet_DP_matrix_base::allowed_in_states (const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const
{
  return filter.allowed_in_states[Triplet_SCFG_filter::flag_triplet (subseq_x.in_flags, subseq_y.in_flags, subseq_z.in_flags)];
}

inline const vector<int>& Triplet_DP_matrix_base::allowed_out_states (const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const
{
  return filter.allowed_out_states[Triplet_SCFG_filter::flag_triplet (subseq_x.out_flags, subseq_y.out_flags, subseq_z.out_flags)];
}

inline const vector<int>& Triplet_DP_matrix_base::allowed_dest_in_states (int src_state, const Subseq& dest_subseq_x, const Subseq& dest_subseq_y, const Subseq& dest_subseq_z) const
{
  // For use in the Inside/CYK recursion:
  //  use .in_flags to only return destination states whose emission profiles
  //  correspond to the flag_triplet of the destination subsequence in the Inside recursion
  return filter.allowed_dest_states[Triplet_SCFG_filter::flag_triplet (dest_subseq_x.in_flags, dest_subseq_y.in_flags, dest_subseq_z.in_flags)][src_state];
}


inline const vector<int>& Triplet_DP_matrix_base::allowed_src_out_states (int dest_state, const Subseq& src_subseq_x, const Subseq& src_subseq_y, const Subseq& src_subseq_z) const
{
  return filter.allowed_src_states[Triplet_SCFG_filter::flag_triplet (src_subseq_x.out_flags, src_subseq_y.out_flags, src_subseq_z.out_flags)][dest_state];
}

// NB: We're using a formulation where the "source" state
// (the state on the left-hand side of the production rule)
// emits symbols.
inline int Triplet_DP_matrix_base::dest_subseq_idx_x (int subseq_idx_x, State_type t) const
{
  switch (t & EmitXLR) // t is the state type of the state corresponding to the X subseq subseq_idx_x
    {
      // in_l() returns the index of the (inside) subsequence (a subsequence of foldenv_x.subseq[subseq_idx_x])
      // reachable from foldenv_x.subseq[subseq_idx_x] by a left emission
    case EmitXL:  return foldenv_x.subseq[subseq_idx_x].in_l(); break;
    case EmitXR:  return foldenv_x.subseq[subseq_idx_x].in_r(); break;
    case EmitXLR: return foldenv_x.subseq[subseq_idx_x].in_lr(); break;
    default: return subseq_idx_x; break;
    }
  return -1; // unreachable
}


inline int Triplet_DP_matrix_base::dest_subseq_idx_y (int subseq_idx_y, State_type t) const
{
  switch (t & EmitYLR)
    {
    case EmitYL:  return foldenv_y.subseq[subseq_idx_y].in_l(); break;
    case EmitYR:  return foldenv_y.subseq[subseq_idx_y].in_r(); break;
    case EmitYLR: return foldenv_y.subseq[subseq_idx_y].in_lr(); break;
    default: return subseq_idx_y; break;
    }
  return -1; // unreachable
}


inline int Triplet_DP_matrix_base::dest_subseq_idx_z (int subseq_idx_z, State_type t) const
{
  switch (t & EmitZLR)
    {
    case EmitZL:  return foldenv_z.subseq[subseq_idx_z].in_l(); break;
    case EmitZR:  return foldenv_z.subseq[subseq_idx_z].in_r(); break;
    case EmitZLR: return foldenv_z.subseq[subseq_idx_z].in_lr(); break;
    default: return subseq_idx_z; break;
    }
  return -1; // unreachable
}

inline int Triplet_DP_matrix_base::src_subseq_idx_x (int subseq_idx_x, State_type t) const
{
  switch (t & EmitXLR) // t is the state type of the state corresponding to the X subseq subseq_idx_x
    {
    case EmitXL:  return foldenv_x.subseq[subseq_idx_x].out_l(); break;
    case EmitXR:  return foldenv_x.subseq[subseq_idx_x].out_r(); break;
    case EmitXLR: return foldenv_x.subseq[subseq_idx_x].out_lr(); break;
    default: return subseq_idx_x; break;
    }
  return -1; // unreachable
}

inline int Triplet_DP_matrix_base::src_subseq_idx_y (int subseq_idx_y, State_type t) const
{
  switch (t & EmitYLR) // t is the state type of the state corresponding to the Y subseq subseq_idx_y
    {
    case EmitYL:  return foldenv_y.subseq[subseq_idx_y].out_l(); break;
    case EmitYR:  return foldenv_y.subseq[subseq_idx_y].out_r(); break;
    case EmitYLR: return foldenv_y.subseq[subseq_idx_y].out_lr(); break;
    default: return subseq_idx_y; break;
    }
  return -1; // unreachable
}

inline int Triplet_DP_matrix_base::src_subseq_idx_z (int subseq_idx_z, State_type t) const
{
  switch (t & EmitZLR) // t is the state type of the state corresponding to the Z subseq subseq_idx_z
    {
    case EmitZL:  return foldenv_z.subseq[subseq_idx_z].out_l(); break;
    case EmitZR:  return foldenv_z.subseq[subseq_idx_z].out_r(); break;
    case EmitZLR: return foldenv_z.subseq[subseq_idx_z].out_lr(); break;
    default: return subseq_idx_z; break;
    }
  return -1; // unreachable
}

// see e.g. Triplet_SCFG::path_emit_score () and Pair_CFG_DP_matrix::in_emit_idx ().
inline int Triplet_DP_matrix_base::in_emit_idx (State_type t,  const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const
{

  int emit_idx = 0; // map emit string to ints for array access
  if (t & EmitXL) emit_idx += scfg.emit_xl_mul (t) * dsq_x[subseq_x.start];
  if (t & EmitXR) emit_idx += scfg.emit_xr_mul (t) * dsq_x[subseq_x.start + subseq_x.len - 1];
  if (t & EmitYL) emit_idx += scfg.emit_yl_mul (t) * dsq_y[subseq_y.start];
  if (t & EmitYR) emit_idx += scfg.emit_yr_mul (t) * dsq_y[subseq_y.start + subseq_y.len - 1];
  if (t & EmitZL) emit_idx += scfg.emit_zl_mul (t) * dsq_z[subseq_z.start];
  if (t & EmitZR) emit_idx += scfg.emit_zr_mul (t) * dsq_z[subseq_z.start + subseq_z.len - 1];
  return emit_idx;

  // Note that this is equivalent to (this version is slightly more transparent but slower):
  //  int xl, xr, yl, yr, zl, zr;
  //  if (t & EmitXL) xl = dsq_x[subseq_x.start];
  //  if (t & EmitXR) xr = dsq_x[subseq_x.start + subseq_x.len - 1];
  //  if (t & EmitYL) yl = dsq_y[subseq_y.start];
  //  if (t & EmitYR) yr = dsq_y[subseq_y.start + subseq_y.len - 1];
  //  if (t & EmitZL) zl = dsq_z[subseq_z.start];
  //  if (t & EmitZR) zr = dsq_z[subseq_z.start + subseq_z.len - 1];
  //  return (*this).scfg.emit_idx (t, xl, xr, yl, yr, zl, zr);

  //  return (*this).scfg.emit_idx (t, dsq_x[subseq_x.start], dsq_x[subseq_x.start + subseq_x.len - 1], dsq_y[subseq_y.start], dsq_y[subseq_y.start + subseq_y.len - 1], dsq_z[subseq_z.start], dsq_z[subseq_z.start + subseq_z.len - 1]);
}

inline int Triplet_DP_matrix_base::out_emit_idx (State_type t,  const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const
{

  int emit_idx = 0; // map emit string to ints for array access
  if (t & EmitXL) emit_idx += scfg.emit_xl_mul (t) * dsq_x[subseq_x.start - 1];
  if (t & EmitXR) emit_idx += scfg.emit_xr_mul (t) * dsq_x[subseq_x.start + subseq_x.len];
  if (t & EmitYL) emit_idx += scfg.emit_yl_mul (t) * dsq_y[subseq_y.start - 1];
  if (t & EmitYR) emit_idx += scfg.emit_yr_mul (t) * dsq_y[subseq_y.start + subseq_y.len];
  if (t & EmitZL) emit_idx += scfg.emit_zl_mul (t) * dsq_z[subseq_z.start - 1];
  if (t & EmitZR) emit_idx += scfg.emit_zr_mul (t) * dsq_z[subseq_z.start + subseq_z.len];
  return emit_idx;

  // Equivalent to:
  //  int xl, xr, yl, yr, zl, zr;
  //  if (t & EmitXL) xl = dsq_x[subseq_x.start - 1];
  //  if (t & EmitXR) xr = dsq_x[subseq_x.start + subseq_x.len];
  //  if (t & EmitYL) yl = dsq_y[subseq_y.start - 1];
  //  if (t & EmitYR) yr = dsq_y[subseq_y.start + subseq_y.len];
  //  if (t & EmitZL) zl = dsq_z[subseq_z.start - 1];
  //  if (t & EmitZR) zr = dsq_z[subseq_z.start + subseq_z.len];
  //  return (*this).scfg.emit_idx (t, xl, xr, yl, yr, zl, zr);

  //  return (*this).scfg.emit_idx (t, dsq_x[subseq_x.start-1], dsq_x[subseq_x.start + subseq_x.len], dsq_y[subseq_y.start-1], dsq_y[subseq_y.start + subseq_y.len], dsq_z[subseq_z.start-1], dsq_z[subseq_z.start + subseq_z.len]);
}


inline Score& Triplet_DP_dense_allocator::cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z)
{
  return _cell[subseq_idx_x][subseq_idx_y][subseq_idx_z][state];
}


inline Score Triplet_DP_dense_allocator::read_cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const
{
  return _cell[subseq_idx_x][subseq_idx_y][subseq_idx_z][state];
}

// Here's how this works.  The basic idea is that we can store the contents
// of a 2D container (a matrix) as a 1D container (a vector)
// by interleaving the entries such that element (i,j) is accessed as
// i + (dimension of i)*j
// We do this repeatedly for n-dimensional containers.
// For example, this is done for 2 sequences X and Y as follows.
// For each pair (x,y) of subseq start coordinates,
// we access cells as
// xsubseq.by_start_index + ind.xsize * (state + scfg.num_states() * ysubseq.by_start_index)
// which in English is
// (index of X subseq) + (# of X subseqs) * (state + (# of states * index of Y subseq))
// For our case of 3 seqs X,Y,Z the container is structured as (Z (X (Y (states) ) ) ).
// The weird structure is for consistency with the X,Y case (see scfg/paircfgdp.h).
// NB: The 'offset' keeps track of where we are in the vector for each pair (x,y) of subseq start coordinates.
// See Triplet_DP_sparse_allocator::alloc() for details.
inline Score& Triplet_DP_sparse_allocator::cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z)
{
  const Subseq& xsubseq = foldenv_x_begin[subseq_idx_x];
  const Subseq& ysubseq = foldenv_y_begin[subseq_idx_y];
  const Subseq& zsubseq = foldenv_z_begin[subseq_idx_z];
  const By_start_indices& ind = by_start_indices[zsubseq.start] (xsubseq.start, ysubseq.start);
  return _cell[ind.offset + zsubseq.by_start_index + ind.zsize * (xsubseq.by_start_index + ind.xsize * (state + scfg.num_states() * ysubseq.by_start_index))];
}

inline Score Triplet_DP_sparse_allocator::read_cell (int state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const
{
  const Subseq& xsubseq = foldenv_x_begin[subseq_idx_x];
  const Subseq& ysubseq = foldenv_y_begin[subseq_idx_y];
  const Subseq& zsubseq = foldenv_z_begin[subseq_idx_z];
  const By_start_indices& ind = by_start_indices[zsubseq.start] (xsubseq.start, ysubseq.start);
  if (ind.xsize == 0) return -InfinityScore;
  return _cell[ind.offset + zsubseq.by_start_index + ind.zsize * (xsubseq.by_start_index + ind.xsize * (state + scfg.num_states() * ysubseq.by_start_index))];
}


inline bool Triplet_DP_matrix_base::local_start_allowed (const Subseq& subseq_x, const Subseq& subseq_y, const Subseq& subseq_z) const
{
  return (subseq_x.start_flag && subseq_y.start_flag && subseq_z.start_flag);
}


#endif /* TRIPLETSCFGDP_INCLUDED */
