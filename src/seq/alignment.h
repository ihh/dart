// multiple alignment structure

#ifndef ALIGNMENT_INCLUDED
#define ALIGNMENT_INCLUDED

#include <list>
#include <set>
#include <algorithm>
#include <numeric>
#include "seq/biosequence.h"
#include "util/nstring.h"
#include "util/strsaver.h"

#define DEFAULT_GAP_CHARS "-._"

class Alignment_path
{
 public:

  typedef vector<int>   Sequence_coords;
  typedef vector<int>   Row_index_set;
  typedef vector<int>   Column_index_set;
  typedef vector<bool>  Row;
  typedef pair<int,int> Row_pair;
  typedef map<Row_pair,Alignment_path> Decomposition;

 protected:
  // data
  vector<Row> path_data;     // path_data[R][C] = row R, column C

 public:
  // constructor
  Alignment_path (int rows = 0, int cols = 0);

  // equality comparison operator
  bool operator== (const Alignment_path& a) const;
  bool operator!= (const Alignment_path& a) const;

  // dimensions
  inline int rows() const;        // number of rows
  inline int columns() const;     // number of columns (assumes that all rows are same length)

  // path element accessors
  inline Row& row (int row);  // return reference to row #row
  inline const Row& row (int row) const;   // return const reference to row #row
  inline vector<bool> get_column (int col) const;    // create a new vector containing column #col
  inline void set_column (int col, const vector<bool>& col_data);        // set contents of column #col

  // path element accessors, operator shorthand
  inline bool operator() (int row, int col) const;  // shorthand; return value at (row,col)
  inline Row& operator[] (int row);                 // shorthand; return reference to row #row
  inline const Row& operator[] (int row) const;     // shorthand; return const reference to row #row

  // build methods
  // methods to reset the alignment path
  inline void clear();         // set number of rows to zero
  inline void reset_rows();    // set each row to zero length without changing the number of rows

  // methods to make the alignment path flush
  void make_flush();    // pad out rows with 0's until all rows same length
  void make_flush (int n_rows);   // ensures there are n_rows rows

  // methods to insert and delete rows and columns
  // methods to insert before row #row
  void insert_rows (int row, int n = 1);  // inserts n empty rows
  void insert_rows (int row, const Row& row_data, int n = 1);  // inserts n copies of row_data

  // methods to erase columns
  void erase_rows (int row, int n = 1);  // erase rows #row ... #(row+n-1)
  void erase_rows (const Row_index_set& rows);  // erases rows whose indices are specified by rows vector

  // methods to insert before column #col
  void insert_columns (int col, int n = 1);  // inserts n empty columns
  void insert_columns (int col, const vector<bool>& col_data, int n = 1);  // inserts n copies of col_data

  // method to erase columns
  void erase_columns (int col, int n = 1);  // erase columns #col ... #(col+n-1)
  void erase_empty_columns();  // erase every column that's all 0's
  bool column_is_empty (int col) const;  // TRUE if column contains all 0's

  // methods to append rows and columns
  inline void append_row (const Row& row_data);
  inline void append_column  (const vector<bool>& col_data);

  // swap method
  inline void swap_path (Alignment_path& p);

  // method to count the non-gap characters in a row
  inline int count_steps_in_row (int r) const;  // count non-zero entries in row #r

  // method to get the number of non-gap characters up to & including column #c in row #r
  // WARNING: this method takes O(r) steps to finish.
  // If using it on every column in the alignment, consider using create_seq_coords() and inc_seq_coords() instead
  inline int get_seq_pos (int r, int c) const;

  // method to test if an alignment is ungapped
  bool is_ungapped() const;

  // Sequence_coords methods
  inline Sequence_coords create_seq_coords() const;  // create a sequence co-ordinates vector, initialise to all 0's
  Sequence_coords seq_coords_begin() const;  // same as create_seq_coords()
  Sequence_coords seq_coords_end() const;    // calls create_seq_coords() then advances through whole alignment path
  void inc_seq_coords (Sequence_coords& seq_coords, int col) const;   // seq_coords += get_column(col)
  void dec_seq_coords (Sequence_coords& seq_coords, int col) const;   // seq_coords -= get_column(col)
  void zero_seq_coords (Sequence_coords& seq_coords) const;

  // methods to create empty & full Rows
  inline Row create_empty_row() const;  // return a row of 0's
  inline Row create_full_row() const;  // return a row of 1's
  
  // various helper methods to create dummies, masks & such
  inline Row create_dummy_sequence (int row) const;  // return this row, with all the 0's removed
  Row flag_nonempty_columns() const;  // return a row with 1's for each column in the alignment containing a 1

  // methods to explode & implode rows in an Alignment_path.
  // explode_row() replaces every 1 in row #row with a successive column from exploded_path
  void explode_row (int row, const Alignment_path& exploded_path);
  // explode() is like explode_row(), but does multiple rows at once, and re-orders them
  void explode (const vector<Row_index_set>& row_sets, const vector<Alignment_path>& exploded_paths);
  // implode() is the exact inverse of explode()
  vector<Alignment_path> implode (const vector<Row_index_set>& row_sets);

  // methods to interconvert an Alignment_path and a Decomposition.
  // compose() method builds a multiple alignment from a set of pairwise paths
  // if transducer_gap_ordering is true, then insertions from higher-numbered rows will be placed first
  // (c.f. Holmes 2003 transducer paper); otherwise, insertions from lower-numbered rows are placed first.
  void compose (const Decomposition& decomp, bool squash_columns, bool transducer_gap_ordering = false);
  // compose_and_log() method is same as compose(), but chatters to logfile if logging is turned on
  void compose_and_log (const Decomposition& decomp, bool squash_columns, bool transducer_gap_ordering = false);
  // decompose() method breaks up a multiple alignment into a set of pairwise paths
  Decomposition decompose (const vector<Row_pair>& pairs) const;
  // verify_decomposition() throws an error if Decomposition is invalid
  static void verify_decomposition (const Decomposition& decomp);

  // The following three realign_* methods can be used to disassemble an Alignment_path
  // according to a tree, adjust alignments on a branch or the neighbourhood of a node,
  // and reassemble the Alignment_path.

  // All three methods return a map from the old to the new Alignment_path as a pairwise path.
  // Every 1 in the first row of the map corresponds to a column in the old Alignment_path,
  // and every 1 in the second row corresponds to a column in the new Alignment_path.

  // realign_pair() splits the Alignment_path into pairwise paths, changes one path,
  // then puts everything back together again.
  Alignment_path realign_pair (const Row_pair& pair, const Alignment_path& new_align,
			       const vector<Row_pair>& dependencies);

  // realign_row() splits the Alignment_path into pairwise paths, changes one row,
  // modifies the adjacent pairwise paths, then puts everything back together again.
  // Note that the align_to_new_row_map specifies the relationship of the new row
  // to the whole existing Alignment_path, not just the old row.
  Alignment_path realign_row (int row, const Alignment_path& align_to_new_row_map, const vector<Row_pair>& dependencies);

  // realign_subpath() splits the Alignment_path in two row sets, changes one row set,
  // then puts the two row sets back together again.
  Alignment_path realign_subpath (const Row_index_set& subpath_rows, const Alignment_path& new_subpath,
				  const Alignment_path& old_to_new_subpath_map);

  // debugging output method
  void show (ostream& o, const vector<sstring>& row_names) const;

private:
  // Private member functions for constructing a multiple alignment from a tree of pairwise alignments.
  // As a self-contained algorithm, this should probably be shunted into a new class, but oh well.
  // These methods basically work by maintaining a cursor for each pairwise alignment.
  typedef map<Row_pair,int> Row_pair_cursor;

  // A nascent column is a map of row indices to bools.
  typedef map<int,bool> Nascent_column;

  // add_child_step_to_buffer():
  //  top-level function, initially called on the "root" row.
  //  advances the cursor for a particular row, calling advance_child_dependent_cursors() to propagate this down the tree,
  //  then adds a column to the alignment buffer, anchored on the new residue.
  //  returns FALSE if there are no more residues to be added for that row.
  bool add_child_step_to_buffer (const Row_pair& pair, const Decomposition& decomp,
				 const vector<Row_index_set>& children, bool squash_columns,
				 Row_pair_cursor& cursor);

  // advance_child_dependent_cursors():
  //  advances the cursors for all the pairwise alignments that involve a particular row
  //  (actually calls advance_cursor_to_next_parent_step() for each of these pairwise alignments),
  //  storing the growing column in col_data.
  //  returns FALSE if all of these cursors are at the end of their respective pairwise alignments.
  //
  bool advance_child_dependent_cursors (const Row_pair& pair, const Decomposition& decomp,
					const vector<Row_index_set>& children, bool squash_columns,
					Row_pair_cursor& cursor, Nascent_column& col_data);

  // advance_cursor_to_next_parent_step():
  //  moves the cursor for a pairwise alignment to the next position where the alignment has a '1' in the parent row,
  //  then adds the appropriate column information to col_data.
  //  calls add_child_step_to_buffer() and advance_child_dependent_cursors() as appropriate.
  //  returns FALSE if cursor is at end of pairwise alignment (i.e. no more parent steps were found).
  //
  bool advance_cursor_to_next_parent_step (const Row_pair& pair, const Decomposition& decomp,
					   const vector<Row_index_set>& children, bool squash_columns,
					   Row_pair_cursor& cursor, Nascent_column& col_data);

public:
  // Output_mask holds co-ords & printable columns for local (i.e. padding gap columns trimmed) or global alignment
  struct Output_mask
  {
    Sequence_coords first_match;  // residue index of first non-ignorable (i.e. aligned or annotated) residue in each sequence, or -1 if none exist
    Sequence_coords last_match;   // residue index of final non-ignorable (i.e. aligned or annotated) residue in each sequence, or -1 if none exist
    Column_index_set printable_columns;  // list of columns to print
    // build methods
    void initialise_local (const Alignment_path& path);  // initialise first_match and last_match based on aligned residues with other rows (so unaligned sequences at end will be trimmed)
    void initialise_global (const Alignment_path& path);  // initialise first_match and last_match globally
    void update_printable_columns (const Alignment_path& path);  // called by initialise_local and initialise_global, also by external classes (notably Stockholm) that mess with first_match and last_match (e.g. using annotations)
  };
};

// A Pairwise_path is a subclass of Alignment_path with some trivial helper methods.
struct Pairwise_path : Alignment_path
{
  // constructors
  Pairwise_path();
  Pairwise_path (const Alignment_path& path);
  Pairwise_path (const Alignment_path& path, int row0, int row1, bool collapse_empty_columns);

  // accessors
  inline Row& parent();
  inline Row& child();
  inline const Row& parent() const;
  inline const Row& child()  const;

  // method to swap parent & child
  inline void swap_parent_child();

  // an overloaded append_column() is provided, along with the original
  inline void append_column (bool row0, bool row1);
  inline void append_column (const vector<bool>& col_data);

  // methods to count match columns & overlap
  int match_columns() const;  // returns total number of match columns
  int match_overlap (const Pairwise_path& p) const;  // returns number of overlapping match columns for two paths through the same (lengths of) sequences

  // not-very-useful method for finding mean displacement of match columns from central diagonal
  // (used for ordered over-relaxation attempt)
  double mean_match_displacement() const;

  // method for sorting columns such that insertions are collected before deletions
  void sort_indels();
  
  // shorthand composition operators
  Pairwise_path& operator*= (const Pairwise_path& p2);
  Pairwise_path  operator*  (const Pairwise_path& p2) const;
};

// Subalignment_path is another subclass of Alignment_path for creating subalignments
struct Subalignment_path : Alignment_path
{
  // constructors
  Subalignment_path() : Alignment_path() {}
  Subalignment_path (const Alignment_path& path, const Row_index_set& rows, bool collapse_empty_columns);
  // the following constructor collapses empty columns by definition,
  // returning the alignment-to-subalignment map in align_subalign_map:
  Subalignment_path (const Alignment_path& path, const Row_index_set& rows, Pairwise_path& align_subalign_map);
};

// Named_rows has one member: row_name[], a vector of row names in an alignment.
struct Named_rows
{
  // typedefs
  typedef Alignment_path::Row_index_set    Row_index_set;
  typedef Alignment_path::Column_index_set Column_index_set;
  // data
  vector<sstring> row_name;
  // constructor
  Named_rows (int rows = 0);
  // method to return number of rows
  inline int rows() const;
  // methods to find rows by name
  // find_rows_by_name() returns a set of matching row indices
  Row_index_set find_rows_by_name (const char* name) const;
  // find_row_by_name() throws an exception unless the name matches one, and only one, row
  int find_row_by_name (const char* name) const;
  // methods to assert all names in this alignment are unique
  bool names_unique() const;
  void assert_names_unique() const;
};

// Alignment extends Named_rows, providing in addition to row_name[]:
//  path, an Alignment_path;
//  prof, a vector of pointers to Score_profile's (sequence data).
// Note: read() method ignores Stockholm lines ("#=GF" etc); use Stockholm subclass for that.
struct Alignment : Named_rows, Stream_saver
{
  // typedefs
  typedef Alignment_path::Output_mask Output_mask;

  // static data for gap characters
  // these should be regarded as private; use accessors to read/write
  static char primary_gap_char;
  static vector<int> char_is_gap;

  // static methods returning gap and space characters
  inline static char gap_char() { return primary_gap_char; }
  inline static bool is_gap_char (char c) { return char_is_gap[c]; }
  inline static bool is_gapspace_char (char c) { return is_gap_char(c) || c == ' '; }
  inline static bool is_gapspacenull_char (char c) { return is_gap_char(c) || c == ' ' || c == '\0'; }

  // static accessor to write gap characters
  static void set_gap_chars (const sstring& gap_chars);

  // data
  Alignment_path path;
  vector<const Score_profile*> prof;  // if prof[R] = 0, then row R has no associated sequence data

  // constructors
  Alignment (int rows = 0);
  Alignment (const Alignment_path& path);
  Alignment (const Pairwise_path& path, const Named_profile& xseq, const Named_profile& yseq);
  Alignment (const Alignment& align, const Row_index_set& row_set, bool collapse_empty_columns);
  Alignment (const Sequence_database_index& seqdb_index);  // sets up prof, but all rows have zero columns

  // accessors
  // number of columns
  inline int columns() const;

  // wrapper methods for asking what Symbol_score_map (or gap) is at a particular row & column
  inline bool not_gap (int row, int col) const;
  inline bool is_gap (int row, int col) const;
  const Symbol_score_map& get_ssm (int row, int col) const;

  // reset methods
  void reset (int rows = 0, int cols = 0);
  void clear() { reset(0); }

  // test method, throws an exception if the alignment row lengths don't match the sequence lengths
  void assert_rows_fit_sequences() const;

  // Primitive I/O; Stockholm methods are more up-to-date.
  // read_MUL() needs a Sequence_database to store the sequences that are read in
  void read_MUL (istream& in, Sequence_database& db);
  // write_MUL() needs an Alphabet
  void write_MUL (ostream& out, const Alphabet& alphabet, const sstring* column_labels = 0, bool print_empty_rows = false) const;
  
  // discard_wild_sequences() throws out Score_profiles that contain only wildcards
  // (should also erase them from Sequence_database; potential memory leak)
  void discard_wild_sequences (const Alphabet& alphabet);

  // method to return sorted rows
  const Alignment::Row_index_set sorted_rows() const;

  // method to swap rows
  void swap_rows (int row1, int row2);

  // Various analysis methods.
  // null_emit_score() method returns emit score under a null model.
  Score null_emit_score (const vector<Score>& null_model) const;

  // Methods to compare alignments by counting shared residue/indel pairs or columns.
  // Residue_pair_set methods use Pairform_* classes, which represent a pairwise alignment
  // as a set of pairs of aligned residue indices in the two sequences.
  typedef pair<int,int>       Pairform_match;
  typedef set<Pairform_match> Pairform_alignment;
  typedef set<int> Gap_coord_set;
  typedef pair<sstring,sstring> Row_name_pair;
  typedef map<Row_name_pair,Gap_coord_set> Pairwise_gap_set;
  typedef vector<int> Col_coord_vector;
  typedef set<Col_coord_vector> Col_set;
  // A Residue_pair_set is a structure for holding all residue pairs in an alignment
  typedef map<Row_name_pair,Pairform_alignment> Residue_pair_set;
  // method to find the set of residue pairs
  Residue_pair_set residue_pair_set() const;
  Pairwise_gap_set pairwise_deletions() const;  // map from row-pairs to h_D
  Pairwise_gap_set pairwise_insertions() const;  // map from row-pairs to h_I  
  //method returning set of columns in alignment
  Col_set col_set() const;
  // method returning size of residue_pair_set()
  int residue_pairs() const;
  //methods returning # of insertions and deletions
  int count_pairwise_insertions() const;
  int count_pairwise_deletions() const;
  // comparison methods
  // residue_pair_overlap() computes number of residue pairs shared between two alignments
  int residue_pair_overlap (const Alignment& align) const;
  static int residue_pair_overlap (const Residue_pair_set& set1, const Residue_pair_set& set2);
  int pairwise_deletion_overlap (const Alignment& align) const;
  int pairwise_insertion_overlap (const Alignment& align) const; 
  static int pairwise_gap_overlap (const Pairwise_gap_set& set1, const Pairwise_gap_set& set2);
  // method to compute the number of columns exactly shared between two alignments
  int column_overlap (const Alignment& align) const;
  static int column_overlap (const Col_set& set1, const Col_set& set2);
  // method to compute the normalized AMA similarity score given two alignments
  float ama_similarity (const Alignment& align) const;
};

// dummy singleton object to set Alignment static data
struct Alignment_initializer
{
  Alignment_initializer();
};
extern Alignment_initializer singleton_alignment_initializer;

// Aligned_score_profile is a quick lookup table of pointers to Symbol_score_profile's
// for each (row,col) position in an alignment.
// Optimized for large genomic alignments, Aligned_score_profile avoids the need to create new Symbol_score_map's,
// using pointers to Symbol_score_map's owned by the Alphabet instead.
// To do this, the Aligned_score_profile constructor must be supplied with a vector of pointers to Named_profile's,
// so it can access the original Biosequence objects.
//
// Suppose that asp is an Aligned_score_profile; then:
//  asp.xsize() == number of rows in alignment
//  asp.ysize() == number of columns in alignment
//     asp(X,Y) == pointer to Symbol_score_map for row X, column Y
// [NB this is a little confusing since "row number" suggests a Y-coordinate.]

struct Aligned_score_profile : array2d<const Symbol_score_map*> {
  // accessors; these override (& swap around) the superclass methods. yes this is grim.
  inline int rows() const { return xsize(); }
  inline int columns() const { return ysize(); }
  // initialiser
  void init (const Alignment& align, const vector<Named_profile*>& np, const Alphabet& alph);
  // sliding-window method
  Aligned_score_profile window (int start, int len) const;
  // output
  void show (ostream& out) const;
};

// inline method defs
// Alignment_path

vector<bool> Alignment_path::get_column (int col) const
{
  vector<bool> col_data (rows());
  for (int row = 0; row < rows(); ++row)
    col_data[row] = (*this) (row, col);
  return col_data;
}

void Alignment_path::set_column (int col, const vector<bool>& col_data)
{
  for (int row = 0; row < rows(); row++) path_data[row][col] = col_data[row];
}

bool Alignment_path::operator() (int row, int col) const
{ return path_data[row][col]; }

Alignment_path::Row& Alignment_path::operator[] (int row)
{ return path_data[row]; }

const Alignment_path::Row& Alignment_path::operator[] (int row) const
{ return path_data[row]; }

int Alignment_path::rows() const
{ return path_data.size(); }

int Alignment_path::columns() const
{ return rows() ? path_data[0].size() : 0; }

void Alignment_path::clear()
{ path_data.clear(); }

void Alignment_path::reset_rows()
{ for_contents (vector<Row>, path_data, i) (*i).clear(); }

Alignment_path::Row& Alignment_path::row (int row)
{ return path_data[row]; }

const Alignment_path::Row& Alignment_path::row (int row) const
{ return path_data[row]; }

void Alignment_path::append_row (const Row& row_data)
{ path_data.push_back (row_data); }

void Alignment_path::append_column  (const vector<bool>& col_data)
{ insert_columns (columns(), col_data, 1); }

void Alignment_path::swap_path (Alignment_path& p)
{ path_data.swap (p.path_data); }

int Alignment_path::count_steps_in_row (int r) const
{ return accumulate (row(r).begin(), row(r).end(), 0); }

int Alignment_path::get_seq_pos (int r, int c) const
{ return accumulate (row(r).begin(), row(r).begin() + c, 0); }

Alignment_path::Sequence_coords Alignment_path::create_seq_coords() const
{ return vector<int> (rows(), (int) 0); }

Alignment_path::Row Alignment_path::create_empty_row() const
{ return vector<bool> (columns(), (bool) 0); }

Alignment_path::Row Alignment_path::create_full_row() const
{ return vector<bool> (columns(), (bool) 1); }

Alignment_path::Row Alignment_path::create_dummy_sequence (int row) const
{ return vector<bool> (count_steps_in_row (row), (bool) 1); }

// Pairwise_path

Alignment_path::Row& Pairwise_path::parent() { return row(0); }
Alignment_path::Row& Pairwise_path::child()  { return row(1); }
const Alignment_path::Row& Pairwise_path::parent() const { return row(0); }
const Alignment_path::Row& Pairwise_path::child()  const { return row(1); }
  
void Pairwise_path::swap_parent_child() { parent().swap(child()); }
  
void Pairwise_path::append_column (const vector<bool>& col_data)
{ ((Alignment_path*)this)->append_column (col_data); }

void Pairwise_path::append_column (bool row0, bool row1)
{ path_data[0].push_back (row0); path_data[1].push_back (row1); }

// Named_rows
inline int Named_rows::rows() const
{ return row_name.size(); }

// Alignment
inline int Alignment::columns() const
{ return path.columns(); }

bool Alignment::not_gap (int row, int col) const
{ return path(row,col); }

bool Alignment::is_gap (int row, int col) const
{ return !path(row,col); }

#endif
