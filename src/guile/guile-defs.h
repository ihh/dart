#include <libguile.h>
#include "seq/stockholm.h"
#include "util/sexpr.h"
#include "ecfg/ecfg.h"

extern scm_t_bits stockholm_tag;

struct Stockholm_smob {
  // data
  Stockholm stock;
  Sequence_database seqdb;

  // constructors
  Stockholm_smob() { }
  Stockholm_smob (Stockholm s) : stock(s) { }

  // methods (just wrappers for Stockholm methods)
  void read (const char* filename)
  {
    ifstream infile(filename);
    stock.read_Stockholm (infile, seqdb);
  }

  void write (const char* filename)
  {
    ofstream outfile(filename);
    stock.write_Stockholm (outfile);
  }

  // cast method
  static Stockholm_smob* cast_from_scm (SCM stock_smob);
};

void init_stockholm_type (void);
void init_xrate_primitives (void);

SCM make_stockholm_smob (const Stockholm& stock);

SExpr* scm_to_new_sexpr (SCM scm);
SCM string_to_scm (const char* s);
SCM sexpr_to_scm (SExpr* sexpr);
SCM ecfg_to_scm (const ECFG_scores& ecfg, const ECFG_counts* counts = 0);
