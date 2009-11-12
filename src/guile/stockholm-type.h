#ifndef STOCKHOLM_TYPE_INCLUDED
#define STOCKHOLM_TYPE_INCLUDED

#include "guile/guile-defs.h"

// guile functions
extern scm_t_bits stockholm_tag;
void init_stockholm_type (void);
SCM make_stockholm_smob (const Stockholm& stock);

// guile smobs
struct Stockholm_smob {
  // data
  Stockholm stock;
  Sequence_database seqdb;

  // constructors
  Stockholm_smob() { }
  Stockholm_smob (Stockholm s) : stock(s) { }

  // methods (just wrappers for Stockholm methods)
  void read_from_file (const char* filename)
  {
    ifstream infile(filename);
    stock.read_Stockholm (infile, seqdb);
  }

  void write_to_file (const char* filename)
  {
    ofstream outfile(filename);
    stock.write_Stockholm (outfile);
  }

  void read_from_string (const char* s)
  {
    stringstream ss (s, stringstream::in);
    stock.read_Stockholm (ss, seqdb);
  }

  void write_to_string (sstring& s)
  {
    stock.write_Stockholm(s);
  }

  // cast method
  static Stockholm_smob* cast_from_scm (SCM stock_smob);
};

#endif /* STOCKHOLM_TYPE_INCLUDED */
