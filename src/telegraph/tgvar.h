#ifndef TGVAR_INCLUDED
#define TGVAR_INCLUDED

#include "seq/pvar.h"
#include "seq/dirichlet.h"

// I/O adaptor for PScores object
struct Telegraph_PScores_adaptor
{
  // reference to PScores
  PScores& pscores;
  // set of PGroups to save/load
  set<int> io_groups;
  // constructors
  Telegraph_PScores_adaptor (PScores& pscores);  // puts all PGroups into io_groups
  Telegraph_PScores_adaptor (Dirichlet_prior& prior);  // only puts the PGroups from the prior into io_groups
  // method to return the Telegraph parameter name for group g, var v
  sstring telegraph_varname (int g, int v) const;
  // I/O methods
  void read (const char* filename, bool ignore_unknown_params = false, bool ignore_undefined_params = false);
  void read (istream& in, bool ignore_unknown_params = false, bool ignore_undefined_params = false);
  void read_from_string (const char* param_str, bool ignore_unknown_params = false, bool ignore_undefined_params = false);
  void write (const char* filename) const;
  void write (ostream& out) const;
};

#endif /* TGVAR_INCLUDED */
