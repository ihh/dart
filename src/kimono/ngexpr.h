#ifndef NGEXPR_INCLUDED
#define NGEXPR_INCLUDED

#include "kimono/expr.h"
#include "util/strsaver.h"

struct Normal_Gamma_expr_model : Stream_saver
{
  // parameters
  vector<double> mu;  // mean
  vector<double> tau;  // precision
  // hyperparameters
  // expected mean = mu0
  // expected precision = alpha/beta
  double alpha;  // roughly, half size of prior pseudo-training-set for tau
  vector<double> beta;  // by experiment
  vector<double> mu0;  // by experiment
  double lambda;  // roughly, size of prior pseudo-training-set for mu
  // constructor
  Normal_Gamma_expr_model (const Expr_table::Stats& expr_stats, double alpha, double lambda);
  // accessors
  int experiments() const { return mu.size(); }
  // method to calculate likelihood
  Loge loglike (const vector<double>& data) const;
  // sufficient stats for EM update
  struct Update_stats
  {
    double m0;  // zeroth moment
    vector<double> m1;  // first moment
    vector<double> m2;  // second moment
    Update_stats (const Normal_Gamma_expr_model& ng_expr);
    void clear();
    void add_profile (const vector<double>& profile, double weight = 1.0);
  };
  // EM update
  void optimise (const Update_stats& stats);
  // seed methods
  void default_seed();
  void seed (const vector<double>& profile);
  // display method
  void display (ostream& out) const;
};

#endif /* NGEXPR_INCLUDED */
