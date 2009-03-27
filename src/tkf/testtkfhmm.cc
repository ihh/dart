#include "tkf/tkfhmm.h"
#include "tkf/tkfopts.h"
#include "util/rnd.h"
#include "util/logfile.h"
#include "util/maximise.h"
#include "util/checkderiv.h"

// function objects for extracting emissions & transitions of an HMM
template<class F>
struct HMM_entry
{
  typedef F                               func_type;
  typedef typename func_type::result_type hmm_type;
  typedef typename hmm_type::entry_type   entry_type;
  typedef double                          argument_type;
  typedef entry_type                      result_type;
  func_type& f;
  HMM_entry (func_type& f) : f(f) { }
};

template<class F>
struct HMM_single_emit : HMM_entry<F>
{
  int state;
  int res;
  HMM_single_emit (typename HMM_entry<F>::func_type& f, int state, int res) : HMM_entry<F> (f), state (state), res (res) { }
  typename HMM_entry<F>::result_type operator() (typename HMM_entry<F>::argument_type t) { return f(t).single_emit[state][res]; }
};
template<class F>
HMM_single_emit<F> hmm_single_emit (F& f, int state, int res)
{ return HMM_single_emit<F> (f, state, res); }

template<class F>
struct HMM_pair_emit : HMM_entry<F>
{
  int state;
  int xres;
  int yres;
  HMM_pair_emit (typename HMM_entry<F>::func_type& f, int state, int xres, int yres) : HMM_entry<F> (f), state (state), xres (xres), yres (yres) { }
  typename HMM_entry<F>::result_type operator() (typename HMM_entry<F>::argument_type t) { return f(t).pair_emit[state](xres,yres); }
};
template<class F>
HMM_pair_emit<F> hmm_pair_emit (F& f, int state, int xres, int yres)
{ return HMM_pair_emit<F> (f, state, xres, yres); }

template<class F>
struct HMM_transition : HMM_entry<F>
{
  int src;
  int dest;
  HMM_transition (typename HMM_entry<F>::func_type& f, int src, int dest) : HMM_entry<F> (f), src (src), dest (dest) { }
  typename HMM_entry<F>::result_type operator() (typename HMM_entry<F>::argument_type t) { return f(t).transition(src,dest); }
};
template<class F>
HMM_transition<F> hmm_transition (F& f, int src, int dest)
{ return HMM_transition<F> (f, src, dest); }

template<class F>
struct HMM_start : HMM_entry<F>
{
  int state;
  HMM_start (typename HMM_entry<F>::func_type& f, int state) : HMM_entry<F> (f), state (state) { }
  typename HMM_entry<F>::result_type operator() (typename HMM_entry<F>::argument_type t) { return f(t).start[state]; }
};
template<class F>
HMM_start<F> hmm_start (F& f, int state)
{ return HMM_start<F> (f, state); }

template<class F>
struct HMM_end : HMM_entry<F>
{
  int state;
  HMM_end (typename HMM_entry<F>::func_type& f, int state) : HMM_entry<F> (f), state (state) { }
  typename HMM_entry<F>::result_type operator() (typename HMM_entry<F>::argument_type t) { return f(t).end[state]; }
};
template<class F>
HMM_end<F> hmm_end (F& f, int state)
{ return HMM_end<F> (f, state); }

// Score2Prob function adaptor
template<class F>
struct Score2Prob_adaptor : unary_function<typename F::argument_type,Prob>
{
  F f;
  Score2Prob_adaptor (F f) : f(f) { }
  Prob operator() (typename F::argument_type x) { return Score2Prob(f(x)); }
};
template<class F>
Score2Prob_adaptor<F> s2p (F f) { return Score2Prob_adaptor<F> (f); }

// main
int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  TKF_opts opts (argc, argv, 1);
  opts.short_description = "test TKF HMM objects\n";
  opts.syntax = "[options] <sequence or alignment file>";
  Rnd::add_opts (opts);

  opts.newline();
  Log_stream::add_opts (opts);

  double tmin;
  double tmax;
  double tres;
  bool   fbtime;

  double pairwise_divergence_estimate = 1;

  opts.newline();
  opts.add ("tmin",   tmin = .03,  "\tmin time");
  opts.add ("tmax",   tmax = 3,    "\tmax time");
  opts.add ("tres",   tres = .01, "\ttime step");
  opts.add ("fbtime", fbtime = 0, "\t\tuse forward-backward algorithm");

  opts.add_align_detection();

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
      // display initial logging messages

      Score_fns::describe_scoring_scheme (CLOG(8));
      
      // read in sequences; autodetect aligned/unaligned and DNA/protein
      
      Sequence_database seqs;
      Alignment         align;

      ifstream seq_file (opts.args[0].c_str());
      if (!seq_file) THROW String_exception ("Sequence file not found: ", opts.args[0].c_str());
      
      bool prealigned = opts.detect_aligned(seq_file);
      if (prealigned)
	{
	  CLOG(7) << "Reading alignment\n";
	  align.read_MUL (seq_file, seqs);
	}
      else
	{
	  CLOG(7) << "Reading sequences\n";
	  seqs.read_FASTA (seq_file);
	}

      opts.detect_alphabet (seqs);

      // set up parameters
      TKF_params params = opts.params();
      const Alphabet& alphabet = params.submat_factory.alphabet();
      seqs.seqs2scores (alphabet);

      // iterate over sequence pairs
      int nseqs = seqs.size();
      int row;
      Sequence_database::const_iterator row_prof;
      bool first_row_col = 1;
      for (row = 0, row_prof = seqs.begin(); row < nseqs; ++row, ++row_prof)
	{
	  int col;
	  Sequence_database::const_iterator col_prof;
	  for (col = 0, col_prof = seqs.begin(); col < row; ++col, ++col_prof)
	    {
	      CLOG(6) << "Estimating divergence time of " << (*row_prof).name << " and " << (*col_prof).name << "\n";

	      TKF_counts_function* f;
	      if (prealigned)
		f = new TKF_aligned_counts_function (align, row, col);
	      else
		{
		  if (fbtime)
		    f = new TKF_unaligned_counts_function ((*row_prof).prof_sc, (*col_prof).prof_sc, &params.submat_factory.alphabet());
		  else
		    {
		      // estimate a pairwise alignment then use it to guess at divergence time
		      TKF_joint_pair_HMM_scores tkf_pair_hmm (params, pairwise_divergence_estimate);
		      Pair_Viterbi_DP_matrix tkf_pair_matrix (tkf_pair_hmm, (*row_prof).prof_sc, (*col_prof).prof_sc);
		      vector<int> vtraceback = tkf_pair_matrix.optimal_state_path();
		      Pairwise_path vpath = tkf_pair_hmm.convert_state_path_to_alignment (vtraceback);
		      f = new TKF_aligned_counts_function ((*row_prof).prof_sc, (*col_prof).prof_sc, vpath);
		    }
		}

	      TKF_functions funcs (params, *f, TRUE, TRUE, tres);

	      // first time round, check time derivatives of HMM transitions and emissions
	      if (first_row_col)
		{
		  first_row_col = 0;
		  cout << "Checking time derivatives of HMM transitions\n";
		  const int states = 3;
		  for (int src = 0; src < states; ++src)
		    {
		      // emissions
		      if (src < 2)   // single emit
			for (int res = 0; res < alphabet.size(); ++res)
			  {
			    cout << "State " << src << " single emit (" << res << ")\n";
			    check_deriv (s2p (hmm_single_emit (funcs.pair_scores_function, src, res)),
					 hmm_single_emit (funcs.pair_dt_function, src, res),
					 tmin, tmax, tres, cout);
			  }
		      else           // pair emit
			for (int xres = 0; xres < alphabet.size(); ++xres)
			  for (int yres = 0; yres < alphabet.size(); ++yres)
			    {
			      cout << "State " << src << " pair emit (" << xres << "," << yres << ")\n";
			      check_deriv (s2p (hmm_pair_emit (funcs.pair_scores_function, src, xres, yres)),
					   hmm_pair_emit (funcs.pair_dt_function, src, xres, yres),
					   tmin, tmax, tres, cout);
			    }
		      // transitions
		      for (int dest = 0; dest < states; ++dest)
			{
			  cout << "Transition from " << src << " to " << dest << "\n";
			  check_deriv (s2p (hmm_transition (funcs.pair_scores_function, src, dest)),
				       hmm_transition (funcs.pair_dt_function, src, dest),
				       tmin, tmax, tres, cout);
			}
		      cout << "State " << src << " start transition\n";
		      check_deriv (s2p (hmm_start (funcs.pair_scores_function, src)),
				   hmm_start (funcs.pair_dt_function, src),
				   tmin, tmax, tres, cout);
		      cout << "State " << src << " end transition\n";
		      check_deriv (s2p (hmm_end (funcs.pair_scores_function, src)),
				   hmm_end (funcs.pair_dt_function, src),
				   tmin, tmax, tres, cout);
		    }
		}

	      // check the derivatives of the function object
	      check_deriv (funcs.log_like_function, funcs.log_like_dt_function, tmin, tmax, tres, cout);
	    }
	}
   }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
