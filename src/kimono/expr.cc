#include <stdlib.h>
#include "kimono/expr.h"

const vector<double>& Expr_table::profile_by_name (const sstring& name) const
{
  const int seq_idx = probe_idx.lookup (name);
  if (seq_idx < 0) THROWEXPR ("No expression profile for '" << name << "'");
  return profile[seq_idx];
}

void Expr_table::read (istream& in)
{
  sstring s;
  s.getline (in);
  s.chomp();
  expt_name = s.split ("\t", 0);
  expt_name.erase (expt_name.begin(), expt_name.begin() + 1);
  int experiments = expt_name.size();

  probe_name.clear();
  profile.clear();
  while (1)
    {
      s.getline (in);
      s.chomp();
      if (s.size() == 0) break;
      vector<sstring> field = s.split ("\t", 0);
      if ((int) field.size() != experiments + 1) THROW Format_exception (in, "Wrong number of experiments\n");
      probe_name.push_back (field[0]);
      profile.push_back (vector<double> (experiments));
      for (int i = 0; i < experiments; ++i) profile.back()[i] = atof (field[i+1].c_str());
    }

  probe_idx.clear();
  expt_idx.clear();
  probe_idx.make_index (probe_name.begin(), probe_name.end());
  expt_idx.make_index (expt_name.begin(), expt_name.end());
}

void Expr_table::read (const char* filename)
{
  CLOG(7) << "Reading expression profiles from file '" << filename << "'\n";
  ifstream expr_stream (filename);
  if (!expr_stream) { sstring e; e << "Expression data file '" << filename << "' not found"; THROW Standard_exception (e); }
  read (expr_stream);
}

bool Expr_table::format_help (Opts_list* ol)
{
  CLOGERR << "\n";
  CLOGERR << "Input format for expression profiles is tabular, with rows corresponding to probes and columns to experiments.\n";
  CLOGERR << "First line should contain a tab character (0x09) followed by tab-separated experiment labels.\n";
  CLOGERR << "Subsequent lines should start with the probe label, followed by a tab, followed by tab-separated expression ratios.\n";
  CLOGERR << "\n";
  return 0;
}

Expr_table::Expr_table() { }
Expr_table::Expr_table (const char* filename) { read (filename); }

Expr_table::Stats Expr_table::experiment_stats() const
{
  vector<double> sum1 (experiments(), 0.0);
  vector<double> sum2 (experiments(), 0.0);
  for (int i = 0; i < probes(); ++i)
    for (int j = 0; j < experiments(); ++j)
      {
	const double d = profile[i][j];
	sum1[j] += d;
	sum2[j] += d*d;
      }
  for (int j = 0; j < experiments(); ++j)
    {
      sum1[j] = sum1[j] / (double) experiments();
      sum2[j] = (sum2[j] / (double) experiments()) - sum1[j] * sum1[j];
    }
  Stats stats;
  stats.mean.swap (sum1);
  stats.variance.swap (sum2);
  return stats;
}

Expr_dist_func::Expr_dist_func (const Expr_table& expr_tab, const Expr_table::Stats& expr_stats)
  : expr_tab (expr_tab), expr_stats (expr_stats)
{ }

double Expr_dist_func::operator() (int i, int j)
{
  double ssq = 0;
  for (int k = 0; k < expr_tab.experiments(); ++k)
    {
      const double diff = expr_tab.profile[i][k] - expr_tab.profile[j][k];
      ssq += diff * diff / expr_stats.variance[k];
    }
  return sqrt (ssq);
}

void Scaled_expr_table::reset_scale()
{
  probe_offset = vector<double> (probes(),  (double) 0);
  probe_scaling = vector<double> (probes(),  (double) 1);
  experiment_offset = vector<double> (experiments(), (double) 0);
  experiment_scaling = vector<double> (experiments(), (double) 1);
}

void Scaled_expr_table::read (istream& in)
{
  Expr_table::read (in);
  reset_scale();
}

void Scaled_expr_table::normalise()
{
  probe_offset  = vector<double> (probes(),  (double) 0);
  probe_scaling = vector<double> (probes(),  (double) 1);
  experiment_offset    = vector<double> (experiments(), (double) 0);
  experiment_scaling   = vector<double> (experiments(), (double) 1);

  for (int n = 0; n < probes(); ++n)
    {
      double sum = 0, sq_sum = 0;
      for_const_contents (vector<double>, profile[n], x) { sum += *x; sq_sum += *x * *x; }
      double mean = sum / experiments();
      double sd = sqrt (sq_sum / experiments() - mean * mean);
      probe_offset[n]  = -mean;
      probe_scaling[n] = sd == 0 ? 1 : 1/sd;
    }
}

Scaled_expr_table::Scaled_expr_table (const char* filename) : Expr_table (filename)
{
  reset_scale();
}
