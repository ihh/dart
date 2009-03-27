// this program strips the poly-A tails off sequences, tolerant to occasional misreads in the tail.
// an almost trivial piece of DP to write, but illustrates basic functionality of the single-sequence HMM classes.

#include "hmm/singlehmm.h"

int main (int argc, char** argv)
{
  // create an options handler
  Opts_list opts (argc, argv);
  opts.short_description = "poly-A tail stripper";
  opts.syntax = "[options] <sequence file>";

  // add in error logging parameters for the options handler
  opts.newline();
  Log_stream::add_opts (opts);

  // add hooks for command-line parameters
  double p_misread;

  opts.newline();
  opts.add ("misread", p_misread = 1e-2, "probability of misreads in the tail (# of A's to balance this is approx log(misread)/log(p[A]))");

  // parse the command line
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

  // do the work
  try
    {
      sstring seq_file = opts.args[0];
      // we're *all about* deoxyribose (uhm, even though these are mRNAs)
      const Alphabet& alphabet = DNA_alphabet;
      const int dna_size = alphabet.size();               // this is 4, for anyone who's asleep

      // read in sequences
      // the following line creates a sequence "database", digitised and name-indexed:
      FASTA_sequence_database seq_db (seq_file.c_str(), &alphabet);

      // create an HMM with 2 states
      Single_HMM_scores hmm (2, alphabet);

      // set up emit profile for pre-poly-A state
      const int PREPOLY = 0;   // index of state
      vector<Prob> null_model = seq_db.get_null_model (dna_size);
      // Prob2ScoreVec converts a vector of probabilities into a vector of scores:
      hmm.emit[PREPOLY] = Prob2ScoreVec (null_model);
      
      // set up emit profile for poly-A state
      const int POLY = 1;   // index of state
      hmm.emit[POLY] = vector<Score> (dna_size);
      for (int sym = 0; sym < dna_size; ++sym)
	hmm.emit[POLY][sym] = Prob2Score (null_model[sym] * p_misread);
      // ScorePSumAcc(A,B) evaluates to: A = log (exp(A) + exp(B))
      ScorePSumAcc (hmm.emit[POLY][alphabet.char2int('a')], Prob2Score (1.0 - p_misread));

      // add transitions with score 0
      // not *strictly* probabilistic, but this won't make it any worse at poly-A-tail stripping.
      hmm.start[PREPOLY] = hmm.end[POLY]
	= hmm.transition(PREPOLY,PREPOLY) = hmm.transition(PREPOLY,POLY) = hmm.transition(POLY,POLY) = 0;

      // output the HMM scores for debugging
      if (CLOGGING(5)) hmm.show(CL);

      // do the search
      // for_contents is a cute macro, defined in "dart/macros.h",
      // that declares an iterator and loops it over a container.
      // in this case, the container is 'seq_db' of type 'Sequence_database'
      // and the iterator is 'np' of type 'Sequence_database::iterator',
      // i.e. effectively of type 'Named_profile*'
      for_contents (Sequence_database, seq_db, np)
	{
	  Single_fast_Viterbi_matrix vit (hmm, *np);	// construct the DP matrix
	  vector<int> path = vit.optimal_state_path();  // get the traceback path
	  // inspect the path to see how long it spends in the poly-A state
	  for (int i = path.size() - 1; i >= 0; --i)
	    // strip off one base for every POLY state in the path
	    if (path[i] == POLY) (*np).seq.chop();
	    else break;

	  // output a friendly log message
	  CLOG(6) << "Chopped " << path.size() - (*np).seq.size() << " bases off sequence '" << (*np).name << "'\n";

	  // output the stripped sequence.
	  // since we messed with it, we have to call an update method
	  // or write_FASTA will complain. (yes i know this is uncool)
	  (*np).seq_update (alphabet);
	  (*np).write_FASTA (cout);   // output the stripped sequence
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
