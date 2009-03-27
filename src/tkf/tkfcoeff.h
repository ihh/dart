#ifndef TKFCOEFF_INCLUDED
#define TKFCOEFF_INCLUDED

#include "tkf/tkfparams.h"
#include "seq/alignment.h"

/* Structure definitions */

class TKF_transition_counts {
public:
  double n0;   /*       D->(not I) transitions (DM,DD,DE) */
  double n1;   /* (not D)->(not I) transitions (SM,SD,SE,MM,MD,ME,IM,ID,IE) */
  double n01;  /*       D->I       transition  (DI) */
  double n11;  /* (not D)->I       transitions (SI,MI,II) */
  // methods
  TKF_transition_counts();
  void add_path_counts (const Pairwise_path& path);
  void show (ostream& o) const;

  // private stuff
private:
  enum { MAT = 0, DEL = 1, INS = 2 };
  int TKF_pair_state (const Pairwise_path& path, int pos);
};

struct TKF_counts {
  double T;  /* total time */
  double Z;  /* total expected site-time */
  double I;  /* total expected insertion events */
  double D;  /* total expected deletion events */
  // methods
  TKF_counts();
  void reset_counts();
  void set_params (TKF_params& params);  /* the M step of EM */
  void show (ostream& o) const;
};

struct TKF_coeff {
  /* Time parameters */
  double dt; /* timestep */
  int N;     /* number of timesteps */
  double T;  /* max time; equal to (N-1)*T */
  /* Pair HMM parameters */
  double g;  /* lambda / mu */
  vector<double> alpha;
  vector<double> beta;
  vector<double> gamma;
  /* Effective insertion rates */
  vector<double> lambdaD;
  vector<double> lambdaI;
  vector<double> lambdaE;
  /* Expected insertions by zone type */
  vector<double> I0;
  vector<double> I1;
  vector<double> I01;
  vector<double> I11;
  /* Expected site time by zone type */
  vector<double> Z0;
  vector<double> Z1;
  vector<double> Z01;
  vector<double> Z11;
  // methods
  TKF_coeff (double dt, double T);
  void initialise (const TKF_params& params);
  void accumulate_counts (const TKF_transition_counts& trans, double T, TKF_counts& counts);  /* the E step of EM */
  void accumulate_naive_counts (const TKF_transition_counts& trans, double T, TKF_counts& counts);  /* naive (biased) E-step */
};

#endif /* TKFCOEFF_INCLUDED */
