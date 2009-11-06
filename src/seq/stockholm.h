// pk 3/05 accept long input lines

#ifndef STOCKHOLM_INCLUDED
#define STOCKHOLM_INCLUDED

#include "seq/alignment.h"
#include "seq/gff.h"
#include "util/phonebook.h"

// Stockholm tags

#define Stockholm_pound_space                       "# "
#define Stockholm_format_identifier                 "STOCKHOLM"
#define Stockholm_version_identifier                "1.0"
#define Stockholm_alignment_separator               "//"

#define Stockholm_header Stockholm_pound_space Stockholm_format_identifier " " Stockholm_version_identifier "\n"
#define Stockholm_footer Stockholm_alignment_separator "\n"

#define Stockholm_file_annotation                   "#=GF"
#define Stockholm_column_annotation                 "#=GC"
#define Stockholm_sequence_annotation               "#=GS"
#define Stockholm_sequence_column_annotation        "#=GR"

#define Stockholm_consensus_suffix                  "_cons"
#define StockholmConsensus(STOCKHOLM_TAG)           STOCKHOLM_TAG Stockholm_consensus_suffix

// #=GR and #=GC tags
#define Stockholm_primary_sequence_tag              "PS"
#define Stockholm_secondary_structure_tag           "SS"
#define Stockholm_surface_accessibility_tag         "SA"
#define Stockholm_transmembrane_tag                 "TM"
#define Stockholm_posterior_probability_tag         "PP"
#define Stockholm_ligand_binding_tag                "LI"
#define Stockholm_active_site_tag                   "AS"
#define Stockholm_intron_tag                        "IN"

// #=GF tags
#define Stockholm_comment_tag                       "CC"
#define Stockholm_identifier_tag                    "ID"
#define Stockholm_accession_tag                     "AC"
#define Stockholm_New_Hampshire_tag                 "NH"
#define Stockholm_bit_score_tag                     "SC"
#define Stockholm_Dart_alignment_weight_tag         "WT"
#define Stockholm_GFF_tag                           "GFF"

// other stuff
#define Default_alignment_wildcard_char             '*'
#define Default_annotation_wildcard_char            '.'

// Fold string chars for specifying RNA secondary structures & alignments
struct Fold_char_enum
{
  enum { FOLD_LCHAR = '<', FOLD_RCHAR = '>', NOFOLD_CHAR = '.', GAP_CHAR = '-', NONWC_LCHAR = '(', NONWC_RCHAR = ')', CONTIGUOUS_UNPAIRED_CHAR = '~' };
  static bool is_lchar (char c) { return c == '<' || c == '(' || c == '{' || c == '['; }
  static bool is_rchar (char c) { return c == '>' || c == ')' || c == '}' || c == ']'; }
  // is_contiguous_digram returns true if a bifurcation is prohibited between two adjacent fold string chars
  static bool is_contiguous_digram (char lc, char rc) { return lc == '~' && rc == '~'; }
};

// Stockholm-format multiple alignment with annotations
struct Stockholm : Alignment, Fold_char_enum
{
  // typedefs
  typedef pair<sstring,sstring> Tag_value;   // #=GF line
  typedef map<sstring,set<int> > Tag_index_map;  // map from #=GF tags (e.g. "AC", "CC", "RT", etc) to indices into gf_annot vector
  typedef map<sstring,sstring> Annotation;  // annotation, keyed by tag (e.g. "SS")
  typedef map<sstring,Annotation> Row_annotation;  // keyed by sequence identifier
  // data
  Phonebook row_index;  // mapping from row names to row indices; update after adding rows, by calling update_row_index()
  vector<Tag_value> gf_annot;  // free-form annotation for whole alignment
  Tag_index_map gf_index;
  Annotation gc_annot;  // column annotation for whole alignment
  Row_annotation gs_annot;  // free-form annotation for each sequence
  Row_annotation gr_annot;  // column annotation for each sequence
  vector<Named_profile*> np;   // this is messy, as we have already have pointers to constituent Score_profile's and names, but certain objects need the sequences as Named_profile's

  // constructor
  Stockholm (int rows = 0, int cols = 0);

  // reset method
  void reset (int rows = 0, int cols = 0);
  void clear() { reset(0); }
  // update method
  void update_row_index();
  // I/O
  void read_Stockholm (istream& in, Sequence_database& db, int max_line = 0);		// pk accept long input lines
  // The write_Stockholm method displays the entire alignment.
  // The write_Stockholm_NSE method trims off columns that do not contain any matches or annotations (except for wildcard annotations, i.e. ".").
  // Thus write_Stockholm_NSE effectively displays the interesting "local" part of the alignment.
  // The co-ordinates of the local part are munged artlessly into the sequence name, using NSE (Name/Start-End) notation. (eugh)
  void write_Stockholm (ostream& out) const;
  void write_Stockholm_NSE (ostream& out, char annotation_wildcard_char = Default_annotation_wildcard_char) const;  // writes a local alignment with NSE encoding

  // lower level I/O
  static void write_Stockholm_header (ostream& out);   // writes format & version identifiers
  static void write_Stockholm_separator (ostream& out);   // writes separator line
  void write_Stockholm_body (ostream& out, const Output_mask& out_mask,
			     bool use_NSE_names) const;  // writes alignment

  // helper methods to add rows & set all sequence pointers/refs.
  // after calling add_row(), caller needs to call update_row_index() to update the row_index lookup.
  int add_row();  // appends a new (empty) row
  int add_row (const vector<bool>& new_row, const sstring& new_row_name, Named_profile* new_np = 0);  // appends a new row
  void set_np (int row, const Named_profile& row_np);

  // override Alignment::discard_wild_sequences to update np
  void discard_wild_sequences (const Alphabet& alphabet);

  // const annotation accessors
  sstring get_gf_annot (const sstring& tag) const;  // GF accessor (concatenates multiline annotations)
  sstring get_gc_annot (const sstring& tag) const;  // GC accessor
  sstring get_gs_annot (const sstring& seqname, const sstring& tag) const;  // GS accessor
  sstring get_gr_annot (const sstring& seqname, const sstring& tag) const;  // GR accessor

  // nonconst annotation accessors
  // GF is problematic, can be split over multiple lines (well, all of them can really, but we ignore this)
  void clear_gf_annot();
  void clear_gf_annot (const sstring& tag);
  void add_gf_annot (const sstring& tag, const sstring& val);  // GF accessor; can only add lines
  void set_gc_annot (const sstring& tag, const sstring& val);  // GC accessor
  void set_gs_annot (const sstring& seqname, const sstring& tag, const sstring& val);  // GS accessor
  void set_gr_annot (const sstring& seqname, const sstring& tag, const sstring& val);  // GR accessor

  // method to erase empty columns & update #=GC/#=GR annotations
  void erase_empty_columns();
  void assert_flush() const;  // check that all rows & annotations have same length

  // accessors for specific annotations
  double get_alignment_weight() const;  // returns "#=GF WT" tag
  sstring get_fold_string (const sstring& seqname) const;  // returns "#=GR SS" tag (with gaps removed) or empty string
  void set_fold_string (const sstring& seqname, const sstring& fold_string);  // sets "#=GR SS" tag (adding gaps)
  void add_gff (const GFF& gff);
  void add_gff (const GFF_list& gff_list);
  sstring get_ID() const;
  sstring get_AC() const;
  sstring get_name() const;  // returns ID (if it exists), AC (if it exists & ID doesn't) or empty string

  // helpers for specific annotations
  void propagate_consensus_fold (bool override_row_folds = FALSE);  // maps "#=GC SS_cons" to "#=GR SS" for each row
  void make_consensus_fold();  // maps "#=GR SS" to "#=GC SS_cons"
};

// Stockade class: owns a Stockholm & Named_profiles
struct Stockade
{
  // data
  Stockholm align;
  vector<Named_profile> np;
  // constructors
  Stockade (int rows = 0, int cols = 0);
  Stockade (const Stockholm& s);
  // copy constructor
  Stockade (const Stockade& s);
  // assignment operator
  Stockade& operator= (const Stockade& s);
  // method to add a row (NB alignment path is left empty)
  void add_np (const Named_profile& np);
};

// Stockholm multiple alignment database (e.g. Pfam, Rfam)
class Stockholm_database
{
public:
  // data
  list<Stockholm> align;
  vector<Stockholm*> align_index;
  // sequence name -> (nAlign,nRow) map
  struct Align_row { int nAlign; int nRow; Align_row (int n = -1, int r = -1) : nAlign(n), nRow(r) { } };
  map<sstring,Align_row> align_row;
  // accessors
  int size() const { return align.size(); }

  // I/O methods
  void read (istream& in, Sequence_database& seq_db, int max_line = 0);		// increase max_line to accept long input lines
  void write (ostream& out) const;

  // builder method
  void add (Stockholm& stock);

  // the following method autodetects FASTA/Stockholm (if FASTA, creates a separate alignment for each sequence)
  void read_Stockholm_or_FASTA (istream& in, Sequence_database& seq_db);

  // helpers for various annotations
  void propagate_consensus_folds (bool override_row_folds = FALSE);  // calls propagate_consensus_fold() on each alignment

private:
  // method to update align_row & align_index: called by read(), read_Stockholm_or_FASTA()
  void update_index();
};

#endif /* STOCKHOLM_INCLUDED */
