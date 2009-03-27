#ifndef EXPR_INCLUDED
#define EXPR_INCLUDED

#include "seq/biosequence.h"
#include "util/opts_list.h"
#include "seq/distmat.h"

// a table of expression data
struct Expr_table
{
  // data, labels
  vector<vector<double> > profile;
  vector<sstring>         probe_name;
  vector<sstring>         expt_name;

  // indices
  Phonebook probe_idx;
  Phonebook expt_idx;

  // constructors
  Expr_table();
  Expr_table (const char* filename);

  // input methods
  void read (istream& in);
  void read (const char* filename);

  // accessors
  int probes() const { return probe_name.size(); }
  int experiments() const { return expt_name.size(); }
  
  const vector<double>& profile_by_name (const sstring& name) const;

  // statistics
  struct Stats { vector<double> mean; vector<double> variance; };
  Stats experiment_stats() const;

  // format help method
  static bool format_help (Opts_list* ol);
};

// distance function for expression profiles: sum-of-squares, weighted by expression variance
struct Expr_dist_func : Dist_func
{
  const Expr_table& expr_tab;
  const Expr_table::Stats& expr_stats;
  Expr_dist_func (const Expr_table& expr_tab, const Expr_table::Stats& expr_stats);
  double operator() (int i, int j);
};

// avoid using this class if possible--- it is legacy code used by "gp.h"--- stores an offset & scale for each profile
// anyway, centering your data will send you to a frequentists' Hell full of correlation coefficients and t-tests
struct Scaled_expr_table : Expr_table
{
  // centering stuff, set up by "normalise()"; avoid if possible
  vector<double> probe_offset;
  vector<double> probe_scaling;
  vector<double> experiment_offset;
  vector<double> experiment_scaling;
  // methods
  void reset_scale();  // sets offsets to 0 & scalings to 1
  void read (istream& in);  // calls superclass, then calls reset_scale().... ugh.... overriding is error-prone
  void normalise();  // the standard, hacky, normalising of expression data. avoid if possible
  // constructor
  Scaled_expr_table (const char* filename);
};

#endif
