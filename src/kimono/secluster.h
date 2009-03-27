#ifndef SECLUSTER_INCLUDED
#define SECLUSTER_INCLUDED

#include "kimono/gibbs.h"
#include "kimono/gp.h"
#include "util/strsaver.h"

struct SE_cluster_model : Stream_saver
{
  const Gibbs_alignment_factory&  gibbs_factory;
  const Gaussian_cluster_factory& gauss_factory;

  int k;                            // number of models  (excluding null model)
  vector<Gibbs_alignment*> gibbs;   // sequence models   (excluding null model)
  vector<Gaussian_cluster*> gauss;  // expression models (excluding null model)

  vector<int>          model_prior_sc;  // model_prior_sc[m] = prior score for model m                             (m==0 => null model)
  vector<vector<int> > multi_sc;        // multi_sc[g][m] = posterior score that gene g came from joint SE model m (m==0 => null model)
  vector<vector<int> > gibbs_sc;        // gibbs_sc[g][m] = likelihood score for gene g by sequence model m        (m==0 => null model)
  vector<vector<int> > gauss_sc;        // gauss_sc[g][m] = likelihood score for gene g by expression model m      (m==0 => null model)

  vector<int>  gibbs_param_sc;  // gibbs_param_sc[m] = parameter prior for sequence model m (excluding null model)
  vector<int>  gauss_param_sc;  // gauss_param_sc[m] = parameter prior for expression model m (excluding null model)

  SE_cluster_model (const Gibbs_alignment_factory& gibbs_factory, const Gaussian_cluster_factory& gauss_factory, int k);
  SE_cluster_model (const SE_cluster_model& model);
  ~SE_cluster_model();

  SE_cluster_model& operator= (const SE_cluster_model& model);

  int genes() const { return gauss_factory.params.table.probes(); }

  void read_model_prior_from_file (const char* model_prior_file);

  // for the following four methods, m=0 => null model
  //
  double log_gene_model_likelihood (int g, int m) const;  // log P(g|m,theta_m)
  double log_gene_likelihood (int g) const;               // log P(g|theta) = log \sum_{m} (P(m) . P(g|m,theta_m))
  double log_likelihood() const;                          // log \prod_{g} P(g|theta)
  double log_parameter_prior() const;                     // log P(theta) = log \prod_{m} P(theta_m)
  double log_probability() const;                         // log P(g,theta) = log P(g|theta) . P(theta)
  double log_gene_model_posterior (int g, int m) const;   // log P(m|g,theta) = log P(m) . P(g|m,theta_m) / P(g|theta)

  void assign_gene_to_model (int g, int m);
  void assign_genes_to_random_models();     // need to call update_model_multiplicities() after this

  void read_gff_motifs (const char* gff_filename);

  void calculate_expected_multiplicities (double& rms_delta_ret); // need to call update_model_multiplicities() after this
  void update_model_multiplicities();       // calls Gibbs & Gauss set_membership_probability() methods

  void sample_alignments (double kT = 1, double gibbs_resample_period = 0.0);
  void sample_gibbs_indent (double kT = 1, const char* log_prefix = "");
  void optimise_gibbs_constrained();
  void optimise_gauss (const char* log_prefix = "");

  void update_gibbs_scores();
  void update_gauss_scores (double kappa = 1);
  void update_scores (double kappa = 1) { update_gibbs_scores(); update_gauss_scores (kappa); }

  void greedy_gauss (const char* log_prefix = "", double kappa = 1, double min_rms_delta = .01, int max_iter = 20);

  void display_verbose (ostream& o);
  void display_terse (ostream& o);
  void display_scores (ostream& o, double kappa = 1);           // displays gibbs_sc & gauss_sc arrays
  void display_mixture_consensi (ostream& o);
  void display_alphabet_consensi (ostream& o);
  void display_multiplicities (ostream& o);   // displays multi_sc array

  void save (ostream& o);
  void load (istream& i);
};


#endif


