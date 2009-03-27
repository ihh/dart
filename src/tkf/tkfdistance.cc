#include "tkf/tkfhmm.h"
#include "tkf/tkfopts.h"
#include "util/rnd.h"
#include "util/logfile.h"
#include "util/maximise.h"
#include "util/checkderiv.h"

int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  TKF_opts opts (argc, argv, 1);
  opts.short_description = "calculate a time (distance) matrix using the TKF model\n";
  opts.syntax = "[options] <sequence or alignment file>";

  double tmax;
  double tres;
  bool   fbtime;
  bool   use_subst;
  bool   use_indels;

  double pairwise_divergence_estimate = 1;

  opts.newline();
  opts.print ("Pairwise distance estimation options\n");
  opts.print ("------------------------------------\n");

  opts.add ("tmax",    tmax = 10,  "\tmaximum separation time between two taxa");
  opts.add ("tres",    tres = .01, "\tfractional resolution of separation times");
  opts.add ("fbtime",  fbtime = FALSE, "\t\tuse forward-backward algorithm to get ML estimate of pairwise divergence times");
  opts.add ("cs -countsubst", use_subst = TRUE, "count substitutions when estimating distances");
  opts.add ("ci -countindel", use_indels = TRUE, "count indels when estimating distances");

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
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  try
    {
      // display initial logging messages

      Score_fns::describe_scoring_scheme (CLOG(8));
      
      // read in sequences; autodetect aligned/unaligned and DNA/protein
      
      Sequence_database seqs;
      Stockholm         align;

      ifstream seq_file (opts.args[0].c_str());
      if (!seq_file) THROW String_exception ("Sequence file not found: ", opts.args[0].c_str());
      
      bool prealigned = opts.detect_aligned(seq_file);
      if (prealigned)
	{
	  CLOG(7) << "Reading alignment\n";
	  align.read_Stockholm (seq_file, seqs);
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

      // calculate distance matrix

      int nseqs = seqs.size();
      array2d<double> best_time (nseqs, nseqs);

      int row;
      Sequence_database::const_iterator row_prof;
      for (row = 0, row_prof = seqs.begin(); row < nseqs; ++row, ++row_prof)
	{
	  int col;
	  Sequence_database::const_iterator col_prof;
	  for (col = 0, col_prof = seqs.begin(); col < row; ++col, ++col_prof)
	    {
	      if (CLOGGING(6))
		{
		  CL << "Estimating divergence time of ";
		  CL << (*row_prof).name << " (" << (*row_prof).size() << " residues) and ";
		  CL << (*col_prof).name << " (" << (*col_prof).size() << " residues)\n";
		}

	      TKF_counts_function* f;
	      if (prealigned)
		f = new TKF_aligned_counts_function (align, row, col, &alphabet);
	      else
		{
		  if (fbtime)
		    f = new TKF_unaligned_counts_function ((*row_prof).prof_sc, (*col_prof).prof_sc, &alphabet);
		  else
		    {
		      // estimate a pairwise alignment then use it to guess at divergence time
		      TKF_joint_pair_HMM_scores tkf_pair_hmm (params, pairwise_divergence_estimate);
		      Pair_Viterbi_DP_matrix tkf_pair_matrix (tkf_pair_hmm, (*row_prof).prof_sc, (*col_prof).prof_sc);
		      vector<int> vtraceback = tkf_pair_matrix.optimal_state_path();
		      Pairwise_path vpath = tkf_pair_hmm.convert_state_path_to_alignment (vtraceback);
		      f = new TKF_aligned_counts_function ((*row_prof).prof_sc, (*col_prof).prof_sc, vpath, &alphabet);
		    }
		}

	      TKF_functions funcs (params, *f, use_indels, use_subst, tres);

	      double t1, t2, t3, tbest, fbest;
	      bracket_maximum (funcs.log_like, t1, t2, t3, 0., tmax);
	      brent_deriv (funcs.log_like, funcs.log_like_dt, t1, t2, t3, tres, tbest, fbest);

	      best_time(row,col) = best_time(col,row) = tbest;
	    }
	  best_time(row,row) = 0;
	}

      cout << nseqs << "\n";
      int old_prec = cout.precision(6);
      for (row = 0, row_prof = seqs.begin(); row < nseqs; ++row, ++row_prof)
	{
	  cout << (*row_prof).name;
	  for (int col = 0; col < nseqs; ++col)
	    {
	      if (col % 9 == 0 && col != 0) cout << "\n";
	      cout << "  ";
	      // formatting lines commented out as they mysteriously introduced '\0's on my SuSE linux box [ihh]
	      //	      cout.width(8);
	      //	      cout.fill(' ');
	      cout << best_time(row,col);
	    }
	  cout << "\n";
	}
      cout.precision (old_prec);
   }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      exit(1);
    }

  return 0;
}
