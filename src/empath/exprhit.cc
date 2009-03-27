#include "empath/exprhit.h"
#include "util/vector_output.h"

Expr_hitter::Expr_hitter (Hitter& base_hitter, Normal_Gamma_expr_model& ng_expr, const Expr_table& expr_tab)
  : base_hitter (base_hitter),
    expr_tab (expr_tab),
    ng_expr (ng_expr),
    ng_stats (ng_expr)
{ }

void Expr_hitter::seed (const Named_profile& np, int pos)
{
  base_hitter.seed (np, pos);
  if (expr_tab.probe_idx.contains (np.name))
    ng_expr.seed (expr_tab.profile_by_name (np.name));
  else
    ng_expr.default_seed();
}

void Expr_hitter::initialise()
{
  base_hitter.initialise();
  ng_stats.clear();
}

Loge Expr_hitter::calc_loglike (const Named_profile& np)
{
  Loge L = base_hitter.calc_loglike (np);
  if (expr_tab.probe_idx.contains (np.name))
    NatsPMulAcc (L, ng_expr.loglike (expr_tab.profile_by_name (np.name)));
  return L;
}

void Expr_hitter::optimise()
{
  base_hitter.optimise();
  ng_expr.optimise (ng_stats);
}

vector<double> Expr_hitter::get_params()
{
  vector<double> params = base_hitter.get_params();
  params.reserve (params.size() + 2 * ng_expr.experiments());
  params.insert (params.end(), ng_expr.mu.begin(), ng_expr.mu.end());
  params.insert (params.end(), ng_expr.tau.begin(), ng_expr.tau.end());
  return params;
}

void Expr_hitter::set_params (const vector<double>& params)
{
  const int offset = ((int) params.size()) - 2 * ng_expr.experiments();
  const vector<double> base_params (params.begin(), params.begin() + offset);
  base_hitter.set_params (base_params);
  copy (params.begin() + offset, params.begin() + offset + ng_expr.experiments(), ng_expr.mu.begin());
  copy (params.begin() + offset + ng_expr.experiments(), params.end(), ng_expr.tau.begin());
}

void Expr_hitter::send_update (ostream& out, const Named_profile& np, Prob weight)
{
  base_hitter.send_update (out, np, weight);
  if (expr_tab.probe_idx.contains (np.name))
    out << expr_tab.profile_by_name (np.name) << " ";
  else
    for (int i = 0; i < ng_expr.experiments(); ++i)
      out << "0 ";
  out << weight << "\n";
}

void Expr_hitter::receive_update (istream& in)
{
  base_hitter.receive_update (in);
  vector<double> profile (ng_expr.experiments());
  double weight;
  for (int i = 0; i < ng_expr.experiments(); ++i)
    in >> profile[i];
  in >> weight;
  ng_stats.add_profile (profile, weight);
}

void Expr_hitter::display (ostream& out)
{
  base_hitter.display (out);
  ng_expr.display (out);
}
