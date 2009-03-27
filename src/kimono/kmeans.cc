#include <iomanip>
#include <algorithm>
#include "kimono/trues.h"
#include "util/logfile.h"
#include "util/vector_output.h"

// main

int main(int argc, char* argv[])
{
  Opts_list opts (argc, argv);
  opts.short_description = "k-means clustering for microarray data";
  opts.syntax = "[options] <expression file>";
  Rnd::add_opts (opts);
  opts.newline();
  Log_stream::add_opts (opts);

  sstring mu_file;
  sstring sigma_file;
  sstring cov_file;
  double  noise_scale;
  bool    normalise;
  int     k;
  double  min_improvement;
  bool    sanity;
  
  opts.newline();
  opts.add ("signalmean", mu_file = "",               "file containing expectation vector for each cluster's mean expression profile (default is zero)", 0);
  opts.add ("signalcov",  sigma_file = "",            "file containing covariance matrix for each cluster's mean expression profile (default is identity matrix)", 0);
  opts.add ("noisecov",   cov_file = "",              "file containing covariance matrix for deviations from the mean expression profile within each cluster (default is identity matrix)", 0);
  opts.add ("nscale",     noise_scale = 1,            "\tscaling factor for noise covariance matrix");
  opts.add ("normalise",  normalise = 1,              "\tnormalise expression profile vectors");
  opts.add ("k",          k = 10,                     "\tnumber of models (the 'k' in k-means clustering)");
  opts.add ("mininc",     min_improvement = .001,     "\tminimum fractional score increment for EM");
  opts.add ("sanity",     sanity = 0,                 "\tdo sanity test for means calculation");
  
  opts.newline();
  opts.add ("exprformathelp", &Expr_table::format_help,              "\tdescribe expression file format");
  opts.add ("gpformathelp",   &Gaussian_cluster_params::format_help, "\t\tdescribe Gaussian prior file format");

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 1) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // get parameters

      sstring expr_file = opts.args[0];

      // check we're not being asked for sanity tests in inappropriate circumstances

      if (sanity && (mu_file != "" || sigma_file != "" || cov_file != "")) {
	CLOGERR << "Oops; can't do sanity tests for fancy covariance matrices yet\n";
	sanity = 0;
      }

      // read in the expression vectors
      
      Scaled_expr_table expr_tab (expr_file.c_str());
      if (normalise) expr_tab.normalise();

      // set up the Gaussian_cluster_factory

      Gaussian_cluster_params gauss_params (expr_tab);
      gauss_params.read (mu_file, sigma_file, cov_file, noise_scale);

      Gaussian_cluster_null_model gauss_null (expr_tab, gauss_params);
      Gaussian_cluster_factory    gauss_factory (gauss_params, gauss_null);

      // set up array of Gaussian_cluster's
      
      vector<Gaussian_cluster*> gauss (k, (Gaussian_cluster*) 0);
      for (int c = 0; c < k; ++c) gauss[c] = gauss_factory.new_model();

      vector<int> which_cluster (expr_tab.probes(), (int) -1);

      // assign each gene to a random cluster

      CLOG(7) << "Assigning genes to random clusters\n";

      for (int i = 0; i < expr_tab.probes(); ++i) {
	for (int c = 0; c < k; ++c) gauss[c]->set_membership_probability (i, -InfinityScore);
	int my_c = Rnd::rnd_int (k);
	which_cluster[i] = my_c;
	gauss[my_c]->set_membership_probability (i, 0);
      }

      // main loop

      int last_sc = -InfinityScore;
      for (int iter = 0; ; ++iter) {

	vector<int> last_cluster = which_cluster;

	// find the means

	for (int c = 0; c < k; ++c) {
	  gauss[c]->optimise_parameters();

	  // sanity check
	  if (sanity) {
	    double sq_err = 0, sq_mag = 0, sq_disp = 0;
	    for (int j = 0; j < expr_tab.experiments(); ++j) {
	      const double mu_mean = gauss[c]->mu (j+1);
	      double sum = 0, sq_dev = 0, n = 0;
	      for (int i = 0; i < expr_tab.probes(); ++i) {
		if (which_cluster[i] == c) {
		  double d = (*gauss[c]->x[i]) (j+1);
		  d += expr_tab.probe_offset[i] + expr_tab.experiment_offset[j];
		  d *= expr_tab.probe_scaling[i] * expr_tab.experiment_scaling[j];
		  sum += d;
		  sq_dev += (d - mu_mean) * (d - mu_mean);
		  n++;
		}
	      }
	      if (n == 0) n = 1;
	      const double mean = sum / n;
	      sq_err += (mean - mu_mean) * (mean - mu_mean);
	      sq_mag += mu_mean * mu_mean;
	      sq_disp += sq_dev / n;
	    }
	    CLOG(6) << "Cluster_" << c << ": sanity check square-error= " << sq_err << ", sq-mag= " << sq_mag << ", sq-dispersal= " << sq_disp << "\n";
	  }
	}

	int tmp_sc = 0;
	for (int i = 0; i < expr_tab.probes(); ++i) {
	  tmp_sc += gauss[which_cluster[i]]->log_likelihood_sc (i);
	}
	CLOG(7) << "Iteration " << iter << ": score after M-step = " << Score2Bits (tmp_sc) << " bits\n";

	// assign each gene to closest mean

	int total_sc = 0;
	for (int i = 0; i < expr_tab.probes(); ++i) {
	  int best_sc = -InfinityScore;
	  int my_c = -1;
	  for (int c = 0; c < k; ++c) {
	    int sc = gauss[c]->log_likelihood_sc (i);
	    if (sc > best_sc) {
	      best_sc = sc;
	      my_c = c;
	    }
	    gauss[c]->set_membership_probability (i, -InfinityScore);
	  }
	  gauss[my_c]->set_membership_probability (i, 0);
	  which_cluster[i] = my_c;
	  total_sc += best_sc;
	}

	CLOG(7) << "Iteration " << iter << ": score after E-step = " << Score2Bits (total_sc) << " bits\n";

	// check for termination
	if (abs (((double) (total_sc - last_sc)) / (double) total_sc) < min_improvement) break;
	last_sc = total_sc;

	bool done = 1;
	for (int i = 0; i < expr_tab.probes(); ++i) {
	  if (which_cluster[i] != last_cluster[i]) { done = 0; break; }
	}
	if (done) break;
      }

      // sort by magnitude of cluster vectors

      vector<double> mu_mag (k, (double) 0);
      for (int c = 0; c < k; ++c) {
	double& mag = mu_mag[c] = 0;
	for (int i = 0; i < expr_tab.experiments(); ++i)
	  mag += gauss[c]->mu(i+1) * gauss[c]->mu(i+1);
      }
      
      vector<int> cluster_order (k);
      for (int c = 0; c < k; ++c) cluster_order[c] = c;

      Schwartzian<double> by_mu_mag (mu_mag);
      sort (cluster_order.begin(),
	    cluster_order.end(),
	    by_mu_mag);
      
      // print out the clusters

      for (int c = 0; c < k; ++c) {
	cout << "Cluster_" << c+1 << ":";
	for (int i = 0; i < expr_tab.probes(); ++i) {
	  if (which_cluster[i] == cluster_order[c]) {
	    cout << " " << expr_tab.probe_name[i];
	  }
	}
	cout << "\n";
      }
      

    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}

