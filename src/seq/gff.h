#ifndef GFF_INCLUDED
#define GFF_INCLUDED

#include <list>
#include "util/sstring.h"
#include "util/Regexp.h"
#include "seq/biosequence.h"

// special keys for group field (GFF3)
#define GFF_ID_tag     "ID"
#define GFF_Parent_tag "Parent"

// enumerated values for strand and frame fields
struct GFF_enum
{
  // enumerations
  enum Strand { PlusStrand = 1, MinusStrand = -1, NoStrand = 0 };
  enum Frame { FrameZero = 0, FrameOne = 1, FrameTwo = 2, NoFrame = -1 };
};

// GFF feature
struct GFF : GFF_enum, NSE
{
  // GFF regexp
  static Regexp regexp;

  // fields
  sstring source;
  sstring feature;  // aka "type"
  double  score;
  Strand  strand;
  Frame   frame;
  sstring group;

  // constructor
  GFF() : NSE(), score(0.), strand(NoStrand), frame(NoFrame) { }
  GFF (const NSE& nse) : NSE(nse), score(0.), strand(NoStrand), frame(NoFrame) { }

  // helpers
  inline int length() const { return end + 1 - start; }
  static Strand string2strand (const sstring& strand_str);
  static Frame string2frame (const sstring& frame_str);

  // method to parse a line of a GFF file
  void initialise (const sstring& s);

  // methods to munge key=value pairs in and out of the group field
  static Regexp key_value_regexp;
  sstring get_value (const char* key) const;
  map<sstring,sstring> get_key_value_map() const;
  void set_value (const char* key, const char* val);
  void set_values (const map<sstring,sstring>& key_val);

  // specific methods for masking
  static const char* mask_tens_key;
  static const char* mask_units_key;
  static const char* seq_key;
  void get_mask (Metascore& meta, bool invert = 1) const;  // copies mask into meta
  void set_mask (const Metascore& meta, bool invert = 1);  // sucks mask out of meta

  // I/O
  friend ostream& operator<< (ostream& outs, const GFF& gff);
  friend istream& operator>> (istream& ins, GFF& gff);
};

struct GFF_list : list<GFF>, GFF_enum
{
  // basic IO
  friend ostream& operator<< (ostream& out, const GFF_list& gff_list);
  friend istream& operator>> (istream& in, GFF_list& gff_list);

  // file IO
  void load (const char* filename);
  void save (const char* filename) const;

  // hacky way of creating a unique ID
  sstring create_unique_id() const
  {
    sstring unique_id;
    unique_id << size() + 1;
    return unique_id;
  }

  // idiosyncratic mask methods
  void acquire_mask (const Metaprob& meta, Prob min_prob, const char* seqname, const char* seq = 0, const char* source = "", const char* feature = "", bool invert = 0);
  void apply_mask (Sequence_database_index& index, int mask_metascore_idx, bool invert = 1) const;
  // the following two methods are untested... doh
  void write_n (ostream& out) const;  // as operator<<, but preceded by number of lines
  void read_n (istream& in);  // reads number of lines first; does not clear
};

#endif
