#include "kimono/ngexpr.h"

Normal_Gamma_expr_model::Normal_Gamma_expr_model (const Expr_table::Stats& expr_stats, double alpha, double lambda)
  : mu (expr_stats.mean.size(), 0.0),
    tau (expr_stats.mean.size(), 1.0),
    alpha (alpha),
    beta (expr_stats.mean.size()),
    mu0 (expr_stats.mean),
    lambda (lambda)
{
  for (int i = 0; i < experiments(); ++i)
    beta[i] = expr_stats.variance[i] * alpha;
  default_seed();
}

Loge Normal_Gamma_expr_model::loglike (const vector<double>& data) const
{
  const Loge log2pi = Prob2Nats (2*Pi);
  Loge L = 0;
  for (int i = 0; i < experiments(); ++i)
    {
      const double d = data[i] - mu[i];
      L += Prob2Nats(tau[i]) - log2pi;
      L -= tau[i] * d * d;
    }
  return L / 2;
}

Normal_Gamma_expr_model::Update_stats::Update_stats (const Normal_Gamma_expr_model& ng_expr)
  : m0 (0.0),
    m1 (ng_expr.experiments(), 0.0),
    m2 (ng_expr.experiments(), 0.0)
{ }

void Normal_Gamma_expr_model::Update_stats::clear()
{
  m0 = 0;
  m1 = vector<double> (m1.size(), 0.0);
  m2 = vector<double> (m2.size(), 0.0);
}

void Normal_Gamma_expr_model::Update_stats::add_profile (const vector<double>& profile, double weight)
{
  m0 += weight;
  for (int i = 0; i < (int) profile.size(); ++i)
    {
      m1[i] += weight * profile[i];
      m2[i] += weight * profile[i] * profile[i];
    }
}

void Normal_Gamma_expr_model::optimise (const Update_stats& stats)
{
  const double lambda_new = lambda + stats.m0;
  const double alpha_new = alpha + stats.m0 / 2;
  for (int i = 0; i < experiments(); ++i)
    {
      const double num = lambda * mu0[i] + stats.m1[i];
      const double beta_new = beta[i] + (lambda*mu0[i]*mu0[i] + stats.m2[i] - num*num/lambda_new) / 2;
      const double mu0_new = num / lambda_new;
      mu[i] = mu0_new;
      tau[i] = alpha_new / beta_new;
    }
}

void Normal_Gamma_expr_model::seed (const vector<double>& profile)
{
  Update_stats stats (*this);
  stats.add_profile (profile, 1.0);
  optimise (stats);
}

void Normal_Gamma_expr_model::default_seed()
{
  Update_stats stats (*this);
  optimise (stats);
}

void Normal_Gamma_expr_model::display (ostream& out) const
{
  save_flags (out);
  left_align (out);

  out << "Expression profile:\n";

  out << "Mean ";
  for (int i = 0; i < experiments(); ++i)
    {
      out.width (10);
      out << mu[i];
    }
  out << "\n";

  out << "  SD ";
  for (int i = 0; i < experiments(); ++i)
    {
      out.width (10);
      out << sqrt (1 / tau[i]);
    }
  out << "\n";

  restore_flags (out);
}
