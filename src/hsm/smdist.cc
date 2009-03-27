#include "util/logfile.h"
#include "util/rnd.h"
#include "tree/pam.h"
#include "hsm/em_matrix.h"
#include "seq/alignment.h"
#include "util/vector_output.h"

struct Submat_cache
{
  Substitution_matrix_factory* factory;
  vector<array2d<Score> > submat;
  vector<int> computed;
  double tmul;  // scaled by expected substitution rate
  Submat_cache (Substitution_matrix_factory* factory, int n_times, double tmul)
    : factory (factory), submat (n_times), computed (n_times, (int) 0), tmul (tmul) { }
  Score pairwise_score (const Alignment& align, int ti, int row1, int row2)
  {
    if (!computed[ti])
      {
	CLOG(3) << "Computing substitution matrix #" << ti << "\n";
	const array2d<Prob> sm_prob = factory->create_joint_substitution_matrix (ti * tmul);
	submat[ti] = Prob2ScoreArray2d (sm_prob);
	computed[ti] = 1;
      }
    const array2d<Score>& sm = submat[ti];
    Score align_sc = 0;
    vector<int> seq_coords = align.path.create_seq_coords();
    for (int c = 0; c < align.columns(); align.path.inc_seq_coords (seq_coords, c++))
      if (align.not_gap(row1,c) && align.not_gap(row2,c))
	{
	  Score col_sc = -InfinityScore;
	  for_const_contents (Symbol_score_map, (*align.prof[row1])[seq_coords[row1]], ss1)
	    for_const_contents (Symbol_score_map, (*align.prof[row2])[seq_coords[row2]], ss2)
	    ScorePSumAcc (col_sc, ScorePMul3 (ss1->second, ss2->second, sm (ss1->first, ss2->first)));
	  ScorePMulAcc (align_sc, col_sc);
	}
    return align_sc;
  }
};

int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  Opts_list opts (argc, argv);
  opts.short_description = "calculate a substitution distance matrix for aligned sequences\n";
  opts.syntax = "[options] <alignment>";
  Rnd::add_opts (opts);

  opts.newline();
  Log_stream::add_opts (opts);

  double tmax;
  double tres;
  int max_fork;
  sstring hsm_file;

  opts.newline();
  opts.add ("tmax", tmax = 10,  "maximum separation time between two taxa");
  opts.add ("tres", tres = .01, "fractional resolution of separation times");
  opts.add ("hsm", hsm_file = "", "\tload hidden substitution matrix from file", 0);
  opts.add ("fork", max_fork = 2, "\tnumber of processes to fork");

  // parse the command line
  //
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
      // get args
      const char* align_filename = opts.args[0].c_str();

      // display initial logging messages
      Score_fns::describe_scoring_scheme (CLOG(8));
      
      // get substitution matrix factory
      Substitution_matrix_factory_with_expected* subst_mat_factory;
      if (hsm_file.size())
	{
	  const Alphabet& alphabet = Protein_alphabet;
	  EM_matrix* em_matrix = new EM_matrix (1, alphabet.size(), max_fork);
	  ifstream hsm_stream (hsm_file.c_str());
	  if (!hsm_stream) THROWEXPR ("Couldn't find substitution matrix file '" << hsm_file << "'");
	  em_matrix->read (hsm_stream);
	  em_matrix->init_alphabet (alphabet);
	  subst_mat_factory = em_matrix;
	}
      else
	subst_mat_factory = new PAM_factory;  // always PAM, for now

      // calibrate the matrix
      const double expected_sub_rate = subst_mat_factory->expected_substitution_rate();
      CLOG(6) << "Calibrating: expected substitution rate for model is " << expected_sub_rate << "\n";

      // initialise submat cache
      const int n_times = (int) (tmax / tres + 1.5);
      Submat_cache cache (subst_mat_factory, n_times, tres / expected_sub_rate);

      // read in sequences
      CLOG(7) << "Reading alignment\n";
      Sequence_database seqs;
      Alignment align;
      ifstream align_file (align_filename);
      if (!align_file) THROWEXPR ("Alignment file '" << align_filename << "' not found");
      align.read_MUL (align_file, seqs);
      const Alphabet& alphabet = subst_mat_factory->alphabet();
      seqs.seqs2scores (alphabet);

      // calculate distance matrix entries
      int nseqs = seqs.size();
      array2d<double> best_time (nseqs, nseqs);

      for (int row = 0; row < nseqs; ++row)
	{
	  for (int col = 0; col < row; ++col)
	    {
	      CLOG(6) << "Estimating divergence time of " << align.row_name[row] << " and " << align.row_name[col] << "\n";
	      
	      vector<Score> log_likelihood (n_times, -InfinityScore);
	      // do initial sweep
	      int step = 32;
	      int n_to_keep = (n_times / step) / 4;
	      double zoom = 1.5;
	      set<int> t_done;
	      set<int> t_next;
	      multiset<Score> sorted_scores;
	      for (int i = step / 2; i <= n_times - step / 2; i += step)
		t_next.insert (i);
	      while (1)
		{
		  // calculate scores
		  for_const_contents (set<int>, t_next, t)
		    if (t_done.find (*t) == t_done.end())
		      {
			const Score sc = cache.pairwise_score (align, *t, row, col);
			log_likelihood[*t] = sc;
			sorted_scores.insert (sc);
			t_done.insert (*t);
		      }
		  // check for termination
		  step /= 2;
		  if (step < 1) break;
		  // zoom in on top (n_to_keep) timepoints
		  multiset<int>::iterator ss_iter = sorted_scores.end();
		  for (int k = 0; k < n_to_keep; ++k) --ss_iter;
		  const Score min_sc = *ss_iter;
		  n_to_keep = (int) (n_to_keep * zoom);  // NB grows much slower than (1/step)
		  // add neighbours of kept set to queue
		  t_next.clear();
		  for_const_contents (set<int>, t_done, t)
		    if (log_likelihood[*t] >= min_sc)
		      {
			if (*t - step >= 0) t_next.insert (*t - step);
			if (*t + step < n_times) t_next.insert (*t + step);
		      }
		}
	      // get best score
	      const int max_i = max_element (log_likelihood.begin(), log_likelihood.end()) - log_likelihood.begin();
	      best_time(row,col) = best_time(col,row) = max_i * tres;
	    }
	  best_time(row,row) = 0;
	}
      
      cout << nseqs << "\n";
      const int old_prec = cout.precision(6);
      for (int row = 0; row < nseqs; ++row)
	{
	  cout << align.row_name[row];
	  for (int col = 0; col < nseqs; ++col)
	    {
	      if (col % 9 == 0 && col != 0) cout << "\n";
	      cout << "  ";
	      cout.width(8);
	      cout << best_time(row,col);
	    }
	  cout << "\n";
	}
      cout.precision (old_prec);

      delete subst_mat_factory;
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
