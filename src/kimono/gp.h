#ifndef GP_INCLUDED
#define GP_INCLUDED

#include "kimono/expr.h"
#include "newmat/newmat.h"
#include "util/opts_list.h"
#include "util/strsaver.h"

struct Gaussian_cluster_params
{
  const Scaled_expr_table& table;

  ColumnVector    Emu;             // expected prior value of mu, the centre of the cluster
  SymmetricMatrix C;               // covariance matrix for cluster centres
  SymmetricMatrix N;               // covariance matrix for within-cluster noise
  
  double          log_invN_norm;   // log(|invN| / sqrt((2*Pi)^d))  = normalising coefficient for exp(-x*invN*x/2)
  double          log_invC_norm;   // log(|invC| / sqrt((2*Pi)^d))  = normalising coefficient for exp(-mu*invC*mu/2)

  SymmetricMatrix invN;            // inverse of N
  SymmetricMatrix invC;            // inverse of C

  static void read_ColumnVector (istream& in, ColumnVector& v, const Phonebook& index, int size);
  static void read_SymmetricMatrix (istream& in, SymmetricMatrix& m, const Phonebook& index, int size);

  void calculate_inverses();
  void read (const sstring& Emu_file, const sstring& C_file, const sstring& N_file, double N_scale = 1);   // calls calculate_inverses()

  void debug_display() const;  // displays C, Emu, N
  void debug_display_inv() const;  // displays invC, invN

  static bool format_help (Opts_list* ol)
    {
      CLOGERR <<
	"\n"
	"Input format for expression profile vectors is on two lines:\n"
	" (1) one tab character (0x09) followed by a tab-separated list of data labels;\n"
	" (2) a tab-separated list of corresponding values.\n"
	"Covariance matrices are symmetric and so only the lower diagonal need be specified.\n"
	"The first line of the input is, again, a tab character followed by a tab-separated\n"
	"list of column labels; successive lines correspond to rows of the matrix, with the\n"
	"first field being the row label, e.g.:\n"
	"\tt=0\tt=1\tt=2\n"
	"t=0\t1.0\n"
	"t=1\t0.5\t1.0\n"
	"t=2\t0\t0.5\t1.0\n"
	"\n";
      return 0;
    }
  
  Gaussian_cluster_params (const Scaled_expr_table& table) : table (table) { }
};

struct Gaussian_cluster_null_model
{
  const Scaled_expr_table& table;
  const Gaussian_cluster_params&  params;

  vector<ColumnVector*> profile;
  vector<double>        profile_null_exponent;      // null log-likelihoods (sans normalising constant) for each probe

  int profile_null_score (int i, double kappa = 1) const { return Nats2Score (params.log_invC_norm - 0.5*log(kappa) + (1/kappa) * profile_null_exponent[i]); }

  Gaussian_cluster_null_model (const Scaled_expr_table& table, const Gaussian_cluster_params& params);
  ~Gaussian_cluster_null_model() { for_contents (vector<ColumnVector*>, profile, p) delete *p; profile.clear(); }
};

struct Gaussian_cluster : Stream_saver
{
  // external parameters
  
  const Gaussian_cluster_params&      params;
  const Gaussian_cluster_null_model&  null_model;
  const vector<ColumnVector*>&        x;            // x[i] is the i'th expression profile vector

  // our variables

  ColumnVector     mu;           // mean expression profile (centre of cluster)
  vector<double>   weight;       // probabilities of cluster membership for x[i]

  SymmetricMatrix  K;            // K = invC + effective_vectors() * invN
  SymmetricMatrix  invK;
  SymmetricMatrix  sqrtK;

  bool             invK_valid;   // cleared every time membership probabilities are changed
  bool             sqrtK_valid;  // cleared every time membership probabilities are changed

  // methods

  int    experiments() const { return x.size() == 0 ? -1 : x[0]->Nrows(); }
  int    vectors() const { return x.size(); }
  double effective_vectors() const { return accumulate (weight.begin(), weight.end(), (double) 0); }

  void   dirty_K() { invK_valid = sqrtK_valid = 0; }    // call this when membership probabilities or probe/experiment scalings are changed
  void   calculate_invK();       // calculates K and invK - called by calculate_mu_mean() after dirty_K() has been called
  void   calculate_sqrtK();      // calculates sqrtK - called by sample_parameters()

  void   calculate_mu_mean (ColumnVector& mu_mean);   // mu_mean = invK * (invC * Emu + invN * weighted_x_sum)
  
  Gaussian_cluster (const Gaussian_cluster_null_model& null_model);
  
  void   set_membership_probability (int point, int score);           // NB this calls dirty_K() [slow] so try to set membership probabilities all at once
  int    log_likelihood_sc (int point, double kappa = 1) const;       // kappa is the variance scaling parameter
  int    log_parameter_prior_sc (double kappa = 1) const;             // returns the log of the prior probability of the model parameters
  void   sample_parameters (double kT = 1);
  void   optimise_parameters();
  
  void   display (ostream& o);
  void   display_scores (ostream& o, double kappa = 1, double threshold = .01);
  static void display_all_scores (ostream& o, const Gaussian_cluster_null_model& null, const vector<Gaussian_cluster*>& model, double kappa = 1);
};

struct Gaussian_cluster_factory
{
  const Gaussian_cluster_params&     params;
  const Gaussian_cluster_null_model& null;

  Gaussian_cluster* new_model() const;

  Gaussian_cluster_factory (const Gaussian_cluster_params& params, const Gaussian_cluster_null_model& null);
};

#endif

