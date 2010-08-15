#ifndef FOLDENV_INCLUDED
#define FOLDENV_INCLUDED

#include <iterator>

#include "util/sstring.h"
#include "seq/biosequence.h"
#include "seq/stockholm.h"
#include "seq/gff.h"
#include "util/strsaver.h"
#include "hmm/pairenv.h"

// Uncomment the following line to use bifurcation "pseudovector" iterators (slower)
//  rather than storing bifurcation vectors explicitly (takes O(L^3) memory)
// #define DART_USE_BIFURCATION_PSEUDOVECTORS

// (this is now controlled by the Makefile - ihh, 1/27/2005)
// (at last test, they didn't seem to be working anyway -- time to delete the code? hmmmm, could still be worth fixing one day -- ih, 12/19/2008)

// First, an enumeration of various constants indexing subsequence connections & flags.
// to do: figure out if I'm right about in_flag and out_flag before checking this in
/// A Connection_enum::Connection_flag is a flag that indicates the connectivity of a subsequence in a fold envelope.
/*
 * Each subsequence in the fold envelope for the sequence (accessed as Fold_envelope::subseq[subseq_idx])
 * contains "flag" information specifying which other cells (subsequences) it can be reached from,
 * ie, which are valid emissions to reach it.
 *
 * In more detail:
 * Suppose F is such a flag for a subsequence, S, with start & end co-ords 
 * (START,END). 
 * (More precisely, let F be the "in-flag" S.in_flags).
 * Then...
 * 
 * If (F & CFLAG_L) is nonzero,
 * S is connected to the subseq with coords (START+1,END)
 *  NB: This also implies that the out-flag F' of the subseq with coords (START+1,END)
 *  is such that (F' & CFLAG_L) is nonzero.
 * 
 * If (F & CFLAG_R) is nonzero,
 * S is connected to the subseq with coords (START,END-1)
 * 
 * If (F & CFLAG_LR) is nonzero,
 * S is connected to the subseq with coords (START+1,END-1)
 * 
 * Specifying these connections allows you to extend the idea of a "fold 
 * envelope", from constraining the set of subsequences that are allowed, 
 * to constraining the set of emissions that are allowed to connect those 
 * subsequences.
 * 
 * This is a more fine-grained level of constraint that is helpful if you 
 * want to precisely constrain the positions at which basepairs or single 
 * nucleotides are allowed to be emitted.
 *
 * NB: Flag information is automatically set when the fold envelope is constructed from e.g. a fold string
 * (see Fold_envelope::initialise_from_fold_string ()).  For example, a subsequence with paired nucleotides
 * at the ends can only be reached by a paired emission (set via Fold_envelope::connect_lr ()).
 */
struct Connection_enum
{
  // Connection flag indices (used as array indices & bit positions in flag masks)
  enum { CONN_L = 0, CONN_R = 1, CONN_LR = 2, N_CONN = 3 };
  // Flag masks
  enum Connection_flag { CFLAG_NONE = 0, CFLAG_L = 1, CFLAG_R = 2, CFLAG_LR = 4,
			 CFLAG_L_R = 3, CFLAG_L_R_LR = 7, CFLAG_TOTAL = 8 };
  // Flag encoding used by *_envelope_string() methods
  enum { CONN_BASE_CHAR = '@' };
};

// Subseq_coords describes a subsequence.
// For a sequence of length L, "start" can be [0,1,2...L] and "len" can be [0,1,2...L-start]
// e.g. if start==0 and len==1, then the subsequence is the first residue of the sequence
//      if start==1 and len==2, then the subsequence is the second and third residues... etc
struct Subseq_coords
{
  // data
  int start;
  int len;       // can be zero

  // constructor
  Subseq_coords (int start, int len);

  // accessors
  inline int end() const { return start + len; }

  // Sorting.
  // Ordering: "FOR subseq = 0 to N { ... }" implies "FOR end = 0 TO L { FOR start = end TO 0 { ... } }"
  // i.e. ascending order is ascending endpos then decreasing startpos.
  // This is the same order as used by Fold_envelope, and is an acceptable fill order for inside DP.
  bool operator< (const Subseq_coords& coords) const
    { return end() < coords.end() || (end() == coords.end() && start > coords.start); }

  // Pair_envelope helper methods
  static inline bool allowed_by_pair_env (const Pair_envelope& pair_env, const Subseq_coords& xsubseq, const Subseq_coords& ysubseq);
  static inline bool startpos_allowed_by_pair_env (const Pair_envelope& pair_env, int xstart, int ystart);
  static inline bool endpos_allowed_by_pair_env (const Pair_envelope& pair_env, int xend, int yend);
};

// Sets and maps involving Subseq_coords
typedef set<Subseq_coords>      Subseq_coords_set;
typedef map<Subseq_coords,int>  Subseq_coords_count;
typedef map<Subseq_coords,char> Subseq_coords_annot;

// A Subseq is an extension of Subseq_coords with info about connections to inside&outside emit&bifurcation neighbours.
// Because every bifurcation is stored, a fully connected Subseq takes O(L) memory.
struct Subseq : Subseq_coords, Fold_char_enum, Connection_enum, Stream_saver
{
  // indices of inside & outside connections
  int in[N_CONN];     // L, R, LR
  int out[N_CONN];    // L, R, LR

  // flags(CFLAG_*) for inside & outside (L,R,LR) connections
  int in_flags;
  int out_flags;

  // flag for outside start connection
  bool start_flag;

  // Bifurcation connection structures
  struct Bifurc_in { int l; int r; Bifurc_in (int l = -1, int r = -1) : l(l), r(r) { } };
  struct Bifurc_out_l { int out; int l; Bifurc_out_l (int out = -1, int l = -1) : out(out), l(l) { } };
  struct Bifurc_out_r { int out; int r; Bifurc_out_r (int out = -1, int r = -1) : out(out), r(r) { } };

  // Subseq_pointer and Subseq_index_pointer (vector<> iterators)

  // Subseq_pointer is an iterator through the fold envelope's "subseq" vector
  typedef vector<Subseq>::const_iterator Subseq_pointer;

  // Subseq_index_pointer is an iterator through the fold envelope's by_len[N] and by_start[N] vectors
  typedef vector<int>::const_iterator Subseq_index_pointer;

  // Subseq_index_rev_pointer is a reverse iterator through the fold envelope's by_len[N] and by_start[N] vectors
  typedef vector<int>::const_reverse_iterator Subseq_index_rev_pointer;

  // Fold_envelope base class.
  // A "fold envelope" explicitly lists an iteration over subsequences in a context-free grammar.
  // As such, it takes O(L^3) memory for the full CFG iteration, but less if constrained,
  // e.g. O(LM^2) if the maximum subsequence length is M.
  // (This complexity will be reduced when bifurcation iterators are finalized.)
  struct Fold_envelope_base
  {
    // The main Subseq vector is sorted by increasing endpos, then decreasing startpos.
    // e.g. in (start+length) notation:
    //  (0+0) (1+0) (0+1) (2+0) (1+1) (0+2) (3+0) ...
    //      in (start-end) notation:
    //  (0-0) (1-1) (0-1) (2-2) (1-2) (0-2) (3-3) ...
    // This is an inside->outside ordering, the same as defined by Subseq::operator<().
    // The builder methods use new_subseq() to create sequentially numbered Subseq's, carefully maintaining this order.
    vector<Subseq> subseq;

    // first_with_end[END] is the index of the first Subseq with endpos equal to END
    // first_after_end[END] is the index of the first Subseq with endpos greater than to END
    // Both are equal to zero if there are no Subseq's with endpos equal to END.
    vector<int> first_with_end, first_after_end;

    // We also store the indices of the Subseq's sorted by length & start pos.
    // by_len[LEN][] contains indices of subseqs with length LEN, sorted by increasing end (& start) pos.
    // NB by_len[0][POS] is the empty subseq at position POS.
    // This is a valid inside->outside ordering.
    vector<vector<int> > by_len;

    // by_start[START][] contains indices of subseqs with start pos START, sorted by increasing end pos.
    // Note this is not a valid inside->outside ordering.
    // by_start.size() = L+1, where L is the sequence length.
    vector<vector<int> > by_start;

    // accessors
    inline Subseq_pointer by_end_begin (int end) const;
    inline Subseq_pointer by_end_end (int end) const;

    // Subseq pseudovec initialiser
    void init_pseudovecs();
  };

  // Bifurcation typedefs
  typedef vector<Bifurc_in> Bifurc_in_pseudovec;
  typedef vector<Bifurc_out_l> Bifurc_outl_pseudovec;
  typedef vector<Bifurc_out_r> Bifurc_outr_pseudovec;

  // Bifurcation connections
  Bifurc_in_pseudovec bif_in;
  Bifurc_outl_pseudovec bif_out_l;
  Bifurc_outr_pseudovec bif_out_r;

  // the following variables hold the position of this Subseq within the Fold_envelope's by_start and by_len vectors
  int by_start_index;
  int by_len_index;

  // constructors
  Subseq (int start, int len);

  // comparison
  friend bool operator== (const Bifurc_in& b1, const Bifurc_in& b2);
  friend bool operator== (const Bifurc_out_l& b1, const Bifurc_out_l& b2);
  friend bool operator== (const Bifurc_out_r& b1, const Bifurc_out_r& b2);
  bool operator== (const Subseq& subseq) const;

  // accessors
  inline int in_l() const { return in[CONN_L]; }   /// Return the index for the (inside) subsequence connected to (*this) by a left emission.
  inline int in_r() const { return in[CONN_R]; }
  inline int in_lr() const { return in[CONN_LR]; }

  inline int out_l() const { return out[CONN_L]; }
  inline int out_r() const { return out[CONN_R]; }
  inline int out_lr() const { return out[CONN_LR]; }

  // test methods
  bool test_coords() const;
  bool test_flags() const;
  bool test_in_l (const Subseq& l_subseq) const;
  bool test_in_r (const Subseq& r_subseq) const;
  bool test_in_lr (const Subseq& lr_subseq) const;
  bool test_in_b (const Subseq& l_subseq, const Subseq& r_subseq) const;
  bool test_fold_string (const sstring& fold_string) const;  // tests in[] validity & ">" vs "<" balancing

  // display methods
  sstring terse_desc (int index, const Biosequence& seq, bool show_out = 0) const;
  sstring terser_desc() const;
};

// pull out the Subseq_pointer typedef
typedef Subseq::Subseq_pointer Subseq_pointer;

// A Local_fold_string is a fold_string with sequence length & startpos
struct Local_fold_string
{
  // data
  sstring fold;
  int start;
  int seqlen;
  // constructors
  Local_fold_string();
  Local_fold_string (const sstring& fold, int start, int seqlen);
  // accessors
  inline int end() const;
  bool is_global() const;
  // method to create global fold string
  sstring global_fold() const;
  // merge method
  void add (const Local_fold_string& local_fold_string);
};

// Local_fold_string_database indexed by sequence name
typedef map<sstring,Local_fold_string> Local_fold_string_database;

// A Fold_envelope explicitly lists an iteration over subsequences in a context-free grammar.
// As such, it takes O(L^3) memory for the full CFG iteration, but less if constrained,
// e.g. O(LM^2) if the maximum subsequence length is M.
// The class provides initialise_*() builder methods for constructing various kinds of fold envelope.
class Fold_envelope : public Subseq::Fold_envelope_base, public Fold_char_enum, public Connection_enum
{
public:
  // The default constructor calls initialise_full() and so creates a full envelope.
  // This can later be changed by calling an initialise_*() build method.
  // The default is for a sequence of length zero.
  Fold_envelope (int seqlen = 0);

  // copy constructor
  Fold_envelope (const Fold_envelope& env);

  // Assignment operator.
  Fold_envelope& operator= (const Fold_envelope& env);

  // swap method
  void swap (Fold_envelope& env);

  // Equality comparison operators.
  bool operator== (const Fold_envelope& env) const;
  bool operator!= (const Fold_envelope& env) const;

  // Accessors.
  inline int seqlen() const;
  inline int subseqs() const;
  int bifurcations() const;
  int find_subseq_idx (int start, int len) const;  // returns -1 if not found
  int first_start_subseq() const;  // returns index of first Subseq with start_flag==TRUE

  // Public build methods.
  // In general, building an envelope involves calling an initialise_*() method,
  // then possibly calling connect_all() and/or make_local() to add extra connections.
  // Most of the initialise_*() methods call make_global();
  // Exceptions: initialise_from_envelope_string() sets all flags itself
  //             initialise_local() calls make_local()
  //             initialise_local_fold_string_neighbourhood() sets flags itself

  // First, envelope_string serialisation, together with the to_envelope_string() method.
  sstring to_envelope_string() const;
  void initialise_from_envelope_string (const sstring& str);

  // initialise_local() makes a local version of an existing global envelope,
  // limiting the maximum subseq length & the minimum length for local start transitions into the envelope.
  // Set max_subseq_len to -1 to unlimit Subseq length (i.e. copy the whole envelope).
  void initialise_local (const Fold_envelope& global, int max_subseq_len, int min_start_len = 0);

  // initialise_local_fold_string_neighbourhood() is used to extend an existing structure.
  // Part of the sequence is constrained to fit a specified fold string.
  // Nearby Subseqs are also created, allowing structures to extend, given a minimum overlap with the existing fold.
  // This method can be used in progressive multiple alignment, to extend a hit.
  // Set max_subseq_len to -1 to unlimit Subseq length.
  // If band_3prime is TRUE, then 3' banding sequences from i..L will be created, where
  void initialise_local_fold_string_neighbourhood (const Local_fold_string& local_fold_string,
						   int max_subseq_len, int min_loop_len, int min_overlap_len,
						   bool band_3prime);

  // initialise_local_fold_string_constrained() calls initialise_local_fold_string_neighbourhood() followed by make_global().
  void initialise_local_fold_string_constrained (const Local_fold_string& local_fold_string, int max_subseq_len, int min_loop_len);

  // The following methods all call make_global(), so only the outermost subsequence has a start transition.
  // Envelopes can be made local (adding more start transitions) by calling make_local().

  // The initialise_3prime() method allows emulation of HMM DP iteration.
  // Only L-emit transitions, and no bifurcations, are added.
  void initialise_3prime (int seqlen);  // 3' envelope (L-emitting regular grammar == HMM)
  
  // initialise_5prime() emulates a "backwards" (R-emitting) HMM.
  // Only R-emit transitions, and no bifurcations, are added.
  void initialise_5prime (int seqlen);  // 5' envelope (R-emitting regular grammar)

  // initialise_full() gives the most permissive envelope.
  // All Subseqs, emissions & bifurcations are added.
  void initialise_full (int seqlen, int min_loop_len = 0);  // full envelope

  // initialise_3prime_banded() creates all Subseqs i..j and all Subseqs i..L,
  // where L is sequence length, 0<=i<=L and |j-i|<=M, where M is max_subseq_len.
  // LR emit transitions are added for all subseqs of length |j-i|>=K, where K is min_loop_len.
  // All bifurcations are added.
  void initialise_3prime_banded (int seqlen, int max_subseq_len, int min_loop_len = 0);  // 3' envelope plus all subseqs of length <= max_subseq_len (-1 to unlimit)

  // The following method initialises from a global fold string, resulting in a single-fold envelope.
  // Call the connect_all() method to add extra basepairs & bifurcations.
  void initialise_from_fold_string (const sstring& fold_string, int max_subseq_len = -1);

  // The following methods initialise from a set of subsequences, adding all emissions & bifurcations.
  // First, a Subseq_coords_set.
  void initialise_from_subseq_coords (int seqlen, const Subseq_coords_set& subseq_coords, int min_loop_len = 0);

  // Next, a Subseq_coords_set *with* another envelope that is used to constrain the connection flags
  // (this is used by stemloc for multiple sequence alignment, to avoid adding new connections when extending alignments)
  void initialise_from_subseq_coords_with_cflag_mask (const Subseq_coords_set& subseq_coords, const Fold_envelope& cflag_mask_env, int min_loop_len = 0);

  // Now, a list of GFF features.
  void initialise_from_gff_list (const GFF_list& gff_list, const sstring& seqname, int min_loop_len = 0);
  void to_gff (GFF_list& gff_list, const sstring& seqname) const;

  // The following methods can be called post_build, to allow local start transitions into the envelope.
  // By default, make_global() is called by many of the initialise_*() methods, so the outermost subseq is connected.
  // To allow local start transitions above a minimum subseq length, call make_local().
  // To allow local start transitions with a minimum overlap of a predefined fold, call make_local_overlap().
  // These methods set start_flag for various Subseqs.
  void make_global();  // automatically called by most build methods
  void make_local (int min_start_len = 0);  // called by initialise_local()
  void make_local_overlap (int fold_start, int fold_len, int min_overlap_len);

  // The following method fully connects an envelope.
  void connect_all (int min_loop_len = 0);

  // The following method removes empty bifurcations from an envelope.
  void remove_empty_bifurcations();

  // The following methods can be used to test envelope correctness in various ways.
  bool test_consistent() const;  // painfully slow (uses bifurcation pseudovecs inefficiently)
  bool test_fold_string (const sstring& fold_string) const;
  bool test_grounded() const;  // tests that every subseq has a path to an empty subseq
  bool test_bifs() const;  // tests consistency of bifurcation co-ordinates

  // accessor to check if an envelope is global
  bool is_global() const;

  // Output methods.
  // First a simple string dump, & stream wrapper object.
  void dump (ostream& o) const;
  struct Dump
  {
    const Fold_envelope& env;
    Dump (const Fold_envelope& env) : env (env) { }
    friend ostream& operator<< (ostream& o, const Dump& d) { d.env.dump (o); return o; }
  };

  // Now, beautiful ANSI color dotplots.
  static void render_dotplot_from_counts (ostream& out, const Subseq_coords_count& count, const Biosequence& seq,
					  int max_level = 1, bool use_ansi_color = TRUE,
					  const Subseq_coords_annot* annot = 0);
  // creates a Subseq_coords_count object, then calls render_dotplot_from_counts()
  void render_dotplot (ostream& out, const Biosequence& seq, bool use_ansi_color = TRUE) const;
  // annotates a dotplot with symbols representing connection flags, & colors by number of bifurcations
  void render_annotated_dotplot (ostream& out, const Biosequence& seq, bool use_ansi_color = TRUE) const;

protected:
  // Protected builder methods, called by initialise_*().
  // Methods to reset the envelope
  void clear();
  void resize (int size);

  // Method to create a new subsequence
  int new_subseq (int start, int len);  // throws an exception of a subsequence is created out of order

  // Methods to find indices of a subsequence's neighbours
  int find_l (int subseq_idx) const;
  int find_r (int subseq_idx) const;
  int find_lr (int subseq_idx) const;

  // Methods to connect a subsequence.
  // These methods set in_flags & out_flags bits and in[], out[], bif_in[], bif_out_l[] and bif_our_r[] index arrays.
  void connect_l (int outside, int inside_l);
  void connect_r (int outside, int inside_r);
  void connect_lr (int outside, int inside_lr);
  void connect_b (int outside, int inside_l, int inside_r);

  // connect_subseq() calls all flagged connect_*() methods, plus connect_b()
  void connect_subseq (int subseq_idx, Connection_flag conn = CFLAG_L_R_LR);

  // Methods to *disconnect* a subsequence.
  void disconnect_b (int outside, int inside_l, int inside_r);

  // Method to quickly calculate the index of a sequence in the full envelope
  static int full_subseq_idx (int start, int len) { return len + (start + len) * (start + len + 1) / 2; }

  // Method to calculate connection flags for local envelopes
  static Connection_flag cflag_mask (int len, int min_loop_len)
  { return (Connection_flag) (len >= min_loop_len + 2 ? CFLAG_L_R_LR : CFLAG_L_R); }
};

// Fold_envelope_database is indexed by sequence identifier
struct Fold_envelope_database : map<sstring,Fold_envelope>
{
  // shorthand methods
  void load (const char* filename);
  void load (istream& ins);
  void save (ostream& outs) const;
  // GFF methods
  void from_gff (const GFF_list& gff_list);
  void to_gff (GFF_list& gff_list) const;
  void from_gff (const char* gff_filename);
  void to_gff (const char* gff_filename) const;
};


// inline method defs

// Subseq_coords/Pair_envelope helpers
inline bool Subseq_coords::allowed_by_pair_env (const Pair_envelope& pair_env,
						const Subseq_coords& xsubseq,
						const Subseq_coords& ysubseq)
{
  // cut the following line to stop allowing all zero-length subseqs, regardless of pair envelope; October 1, 2003 -- ihh
  //  xsubseq.len == 0 || ysubseq.len == 0 ||
  return pair_env.allow_cut (xsubseq.start, ysubseq.start) && pair_env.allow_cut (xsubseq.end(), ysubseq.end());
}

inline bool Subseq_coords::startpos_allowed_by_pair_env (const Pair_envelope& pair_env, int xstart, int ystart)
{ return pair_env.allow_cut (xstart, ystart); }

inline bool Subseq_coords::endpos_allowed_by_pair_env (const Pair_envelope& pair_env, int xend, int yend)
{ return pair_env.allow_cut (xend, yend); }

// Local_fold_string
int Local_fold_string::end() const
{ return start + fold.size(); }

// Fold_envelope
int Fold_envelope::seqlen() const
{ return ((int) by_start.size()) - 1; }

int Fold_envelope::subseqs() const
{ return subseq.size(); }

Subseq::Subseq_pointer Subseq::Fold_envelope_base::by_end_begin (int end) const
{ return subseq.begin() + first_with_end[end]; }

Subseq::Subseq_pointer Subseq::Fold_envelope_base::by_end_end (int end) const
{ return subseq.begin() + first_after_end[end]; }

#endif /* FOLDENV_INCLUDED */
