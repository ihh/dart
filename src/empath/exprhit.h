#ifndef EXPRHIT_INCLUDED
#define EXPRHIT_INCLUDED

#include "empath/multihit.h"
#include "kimono/ngexpr.h"

struct Expr_hitter : Hitter
{
  // data
  Hitter& base_hitter;  // Hitter to which we delegate
  const Expr_table& expr_tab;  // the expression data
  Normal_Gamma_expr_model& ng_expr;  // the model
  Normal_Gamma_expr_model::Update_stats ng_stats;  // the model's update statistics
  // constructor
  Expr_hitter (Hitter& base_hitter, Normal_Gamma_expr_model& ng_expr, const Expr_table& expr_tab);
  // virtual methods
  void seed (const Named_profile& np, int pos);
  void initialise();  // resets counts, updates params
  Loge calc_loglike (const Named_profile& np);
  void optimise();
  vector<double> get_params();
  void set_params (const vector<double>& params);
  // parallelisable update read/write methods
  void send_update (ostream& out, const Named_profile& np, Prob weight);
  void receive_update (istream& in);
  // display method
  void display (ostream& out);
};


#endif /* EXPRHIT_INCLUDED */
