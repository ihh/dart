#ifndef PAIRCFG_XFOLD_INCLUDED
#define PAIRCFG_XFOLD_INCLUDED

// pk 8/05 remove dependence on evoldoer directory -
// copy the appropriate structs, defines here
// from evoldoer/tkfstructuretree.h, .cc
// also use PFOLD_loop_matrix, PFOLD_stem_matrix

#include "ecfg/pfold.h"

// duplicated code from handel/tkfparams.h, to avoid handel dependency.... yuck....

struct TKF_params_redundant
{
  // data
  Substitution_matrix_factory& submat_factory;
  double                       birth_rate;
  double                       death_rate;
  double&                      lambda;  // reference to birth_rate
  double&                      mu;      // reference to death_rate

  // methods
  // constructors
  TKF_params_redundant (Substitution_matrix_factory& submat_factory, double birth_rate, double death_rate) :
    submat_factory(submat_factory), birth_rate(birth_rate), death_rate(death_rate), lambda(this->birth_rate), mu(this->death_rate) {}

  TKF_params_redundant (const TKF_params_redundant& p) :
    submat_factory(p.submat_factory), birth_rate(p.birth_rate), death_rate(p.death_rate), lambda(birth_rate), mu(death_rate) {}

  // operators
  TKF_params_redundant& operator= (const TKF_params_redundant& p)
    { submat_factory = p.submat_factory; birth_rate = p.birth_rate; death_rate = p.death_rate; return *this; }

  // helpers
  double equilibrium_length_prob () const
    { return lambda / mu; }

  double expected_sequence_length () const
    { return equilibrium_length_prob() / (1 - equilibrium_length_prob()); }

  double ancestral_survival_prob (double time) const
    { return exp (-mu * time); }

  double more_descendants_prob (double time) const
    {
      double u = lambda * (1 - exp ((lambda - mu) * time));
      double v = mu - lambda * exp ((lambda - mu) * time);
      return u / v;
    }

  double orphan_survival_prob (double time) const
    {
      double u = -mu * (1 - exp ((lambda - mu) * time));
      double v = (1 - exp (-mu * time)) * (mu - lambda * exp ((lambda - mu) * time));
      return 1 + u / v;
    }


  // Partial derivatives wrt time
  // Thanks to MuPAD for checking these ... :^)

  double danc_dt (double time) const
    { return -mu * exp (-mu * time); }

  double ddes_dt (double time) const
    {
      double u = lambda * (1 - exp ((lambda - mu) * time));
      double v = mu - lambda * exp ((lambda - mu) * time);
      double du_dt = lambda * (- (lambda - mu) * exp ((lambda - mu) * time));
      double dv_dt = - lambda * (lambda - mu) * exp ((lambda - mu) * time);
      return (v * du_dt - u * dv_dt) / (v*v);
    }

  double dorp_dt (double time) const
    {
      double u = -mu * (1 - exp ((lambda - mu) * time));
      double v = (1 - exp (-mu * time)) * (mu - lambda * exp ((lambda - mu) * time));
      double du_dt = -mu * (-(lambda - mu) * exp ((lambda - mu) * time));
      double dv_dt = mu * exp (-mu * time) * (mu - lambda * exp ((lambda - mu) * time)) + (1 - exp (-mu * time)) * (-lambda * (lambda - mu) * exp ((lambda - mu) * time));
      return (v * du_dt - u * dv_dt) / (v*v);
    }

  // Partial derivatives wrt lambda

  double deqm_dlambda () const
    { return 1 / mu; }
  
  double ddes_dlambda (double time) const
    {
      double u = lambda * (1 - exp ((lambda - mu) * time));
      double v = mu - lambda * exp ((lambda - mu) * time);
      double du_dlambda = 1 - exp ((lambda - mu) * time) * (1 + lambda * time);
      double dv_dlambda =  -exp ((lambda - mu) * time) * (1 + lambda * time);
      return (v * du_dlambda - u * dv_dlambda) / (v*v);
    }

  double dorp_dlambda (double time) const
    {
      double u = -mu * (1 - exp ((lambda - mu) * time));
      double v = (1 - exp (-mu * time)) * (mu - lambda * exp ((lambda - mu) * time));
      double du_dlambda = -mu * (-exp ((lambda - mu) * time) * time);
      double dv_dlambda = (1 - exp(-mu * time)) * (- exp((lambda - mu) * time) - lambda * time * exp((lambda - mu) * time));
      return (v * du_dlambda - u * dv_dlambda) / (v*v);
    }

  // Partial derivatives wrt mu

  double deqm_dmu () const
    { return -lambda / (mu * mu); }
  
  double danc_dmu (double time) const
    { return -time * exp(-mu * time); }

  double ddes_dmu (double time) const
    {
      double u = lambda * (1 - exp ((lambda - mu) * time));
      double v = mu - lambda * exp ((lambda - mu) * time);
      double du_dmu = lambda * time * exp((lambda - mu) * time);
      double dv_dmu = lambda * time * exp((lambda - mu) * time) + 1;
      return (v * du_dmu - u * dv_dmu) / (v*v);
    }

  double dorp_dmu (double time) const
    {
      double u = -mu * (1 - exp ((lambda - mu) * time));
      double v = (1 - exp (-mu * time)) * (mu - lambda * exp ((lambda - mu) * time));
      double du_dmu = exp((lambda - mu) * time) - mu * time * exp((lambda - mu) * time) - 1;
      double dv_dmu = time * exp(-mu * time) * (mu - lambda * exp((lambda - mu) * time)) + (1 - exp(-mu * time)) * (lambda * time * exp((lambda - mu) * time) + 1);
      return (v * du_dmu - u * dv_dmu) / (v*v);
    }


  // Shorthands

  double alpha (double time) const { return ancestral_survival_prob (time); }
  double beta (double time) const { return more_descendants_prob (time); }
  double gamma (double time) const { return orphan_survival_prob (time); }
  double kappa () const { return equilibrium_length_prob(); }

  double eqm () const { return equilibrium_length_prob(); }
  double anc (double time) const { return ancestral_survival_prob (time); }
  double des (double time) const { return more_descendants_prob (time); }
  double orp (double time) const { return orphan_survival_prob (time); }

  double neqm () const { return 1 - eqm(); }
  double nanc (double time) const { return 1 - anc(time); }
  double ndes (double time) const { return 1 - des(time); }
  double norp (double time) const { return 1 - orp(time); }

  double dnanc_dt (double time) const { return -danc_dt(time); }
  double dndes_dt (double time) const { return -ddes_dt(time); }
  double dnorp_dt (double time) const { return -dorp_dt(time); }
  
  double dnanc_dlambda (double time) const { return 0; }
  double dndes_dlambda (double time) const { return -ddes_dlambda(time); }
  double dnorp_dlambda (double time) const { return -dorp_dlambda(time); }
  
  double dnanc_dmu (double time) const { return -danc_dmu(time); }
  double dndes_dmu (double time) const { return -ddes_dmu(time); }
  double dnorp_dmu (double time) const { return -dorp_dmu(time); }
  
};

struct PFOLD_TKFST_params
{
  // data
  PFOLD_loop_matrix  loop_mat;
  PFOLD_stem_matrix  stem_mat;
  TKF_params_redundant loop;
  TKF_params_redundant stem;
  double     stem_prob;

  // constructor
  PFOLD_TKFST_params();
};

struct PFOLD_pair_CFG : Pair_CFG_scores, PFOLD_state_enum
{
  // data
  PFOLD_TKFST_params& params;
  double time;
  // states taken from PFOLD_state_enum

  // constructor
  PFOLD_pair_CFG (PFOLD_TKFST_params& params, double time);
};

#endif /* PAIRCFG_XFOLD_INCLUDED */
