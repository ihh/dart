#include "psw/pswmain.h"
#include "seq/stockholm.h"

PSW_aligner::PSW_aligner (Opts_list& opts, Simple_PSW_params& params)
  : opts(opts), params(params), alphabet(params.alphabet)
{
  // register command-line options with Opts_list object
  opts.newline();
  opts.add ("go -gapopen", params.gap_open, "\tgap open probability");
  opts.add ("gx -gapextend", params.gap_extend, "gap extend probability");
  opts.add ("fo -flankopen", params.flank_open, "flank open probability (set to 0 for global alignment)");
  opts.add ("ax -alignextend", params.align_extend, "alignment extension probability");
  opts.add ("nx -nullextend", params.null_extend, "null model extension probability");

  opts.newline();
  opts.add ("or -oddsratio", odds_ratio = TRUE, "\t\tuse odds-ratios instead of likelihoods");
  opts.add ("oa -optacc", opt_acc = FALSE, "\t\tuse \"optimal accuracy\" posterior decoding algorithm instead of Viterbi");

  opts.newline();
  opts.add ("st -subtab", subtab_prefix = "", "\tsave substitution-probability table(s) to file with this prefix", FALSE);
  opts.add ("pt -posttab", posttab_prefix = "", "\tsave posterior-probability table(s) to file with this prefix", FALSE);

  opts.newline();
  opts.add ("sm -submat", submat_filename = "", "\tread substitution matrix from file with this name [format: null emit probs followed by joint probs]", FALSE);
  opts.add ("p -params", params_filename = "", "\tread all parameters, including substitution matrix, from file with this name", FALSE);
  opts.add ("fh -formathelp", &Simple_PSW_params::format_help, "\tprint description of substitution matrix & parameter file formats");
  
  add_simple_PSW_references (opts);  // add some help text
}

void PSW_aligner::run()
{
  opts.parse_or_die();  // parse the command-line options
  try
    {
      // get sequence filename
      const sstring seq_filename = opts.args[0];

      // read in sequences
      FASTA_sequence_database seq_db (seq_filename.c_str(), &alphabet);

      // read in parameters, if asked
      if (params_filename.size())
	{
	  ifstream infile (params_filename.c_str());
	  if (!infile) THROWEXPR ("Could not open params file: " << params_filename);
	  params.read_params (infile);
	}
      else if (submat_filename.size())
	{
	  ifstream infile (submat_filename.c_str());
	  if (!infile) THROWEXPR ("Could not open submat file: " << params_filename);
	  params.read_submat (infile);
	}
      
      // create necessary objects for PSW Pair HMM
      // this is a bit convoluted, because we are using the DART PGroup machinery;
      // the advantage of doing things this way is that we can automate the Baum-Welch training.
      const Simple_PSW_pscores pscores (params);  // create symbolic representations of params; convert probabilities to scores
      const Simple_PSW_PHMM symbolic_hmm (pscores, odds_ratio);  // create symbolic representation of Pair HMM
      const Pair_HMM_scores hmm = symbolic_hmm.eval_hmm_scores (pscores);  // create the HMM itself

      // loop over all pairs of sequences
      for (int x = 0; x < seq_db.size(); ++x)
	for (int y = x + 1; y < seq_db.size(); ++y)
	  {
	    // get the two sequences
	    const Named_profile& xseq = seq_db.get_seq(x);
	    const Named_profile& yseq = seq_db.get_seq(y);
	    // extract the Score_profiles
	    const Score_profile& xprofile = xseq.prof_sc;
	    const Score_profile& yprofile = yseq.prof_sc;

	    // find the best state path through the Pair HMM (Viterbi version)
	    if (!opt_acc)
	      {
		CTAG(7,PSW) << "Finding Viterbi alignment\n";
		const Pair_Viterbi_DP_matrix dp_matrix (hmm, xprofile, yprofile);  // fill the DP matrix
		const vector<int> state_path = dp_matrix.optimal_state_path();  // do traceback
		const Pairwise_path pairwise_path = hmm.convert_state_path_to_alignment (state_path); // convert the Pair HMM state path into an alignment path
		const Alignment align (pairwise_path, xseq, yseq);  // label the alignment
		write_alignment (align, alphabet, Score2Bits(dp_matrix.final_score), "bits");  // and output
	      }

	    // create the substitution table, if asked
	    if (subtab_prefix.size())
	      {
		CTAG(7,PSW) << "Making substitution table\n";
		const array2d<Prob> joint_submat = params.joint_submat();
		const Submat_ProbArray2d submat_matrix (xseq.dsq, yseq.dsq, joint_submat);
		sstring subtab_name;
		subtab_name << subtab_prefix << '-' << xseq.name << '-' << yseq.name;
		ofstream subtab_file (subtab_name.c_str());
		subtab_file << submat_matrix;
	      }
	    
	    // calculate Forward-Backward matrix if necessary
	    if (posttab_prefix.size() || opt_acc)
	      {
		CTAG(7,PSW) << "Finding posterior probabilities\n";
		const Pair_forward_backward_DP_matrix fb (hmm, xprofile, yprofile);  // fill Forward-Backward DP matrix
		const Post_pair_HMM post_matrix (fb);  // get the posterior match probabilities

		// render & save posterior table, if asked
		if (posttab_prefix.size())
		  {
		    CTAG(7,PSW) << "Making posterior table\n";
		    sstring posttab_name;
		    posttab_name << posttab_prefix << '-' << xseq.name << '-' << yseq.name;
		    ofstream posttab_file (posttab_name.c_str());
		    posttab_file << post_matrix;
		  }

		// do posterior decoding alignment (optimal accuracy algorithm)
		if (opt_acc)
		  {
		    CTAG(7,PSW) << "Finding optimal accuracy alignment by posterior decoding\n";
		    const Optimal_accuracy_DP_matrix opt_acc_matrix (post_matrix);  // fill the posterior decoding matrix
		    const Pairwise_path opt_acc_path = opt_acc_matrix.traceback();  // trace back the best path
		    const Alignment align (opt_acc_path, xseq, yseq);  // label the alignment
		    write_alignment (align, alphabet, opt_acc_matrix.final_score(), "expected matches");  // and output
		  }
	      }
	  }
    }
  catch (const Dart_exception& e)  // exception; bail out gracefully
    {
      CLOGERR << e.what();
      exit(1);
    }
}

void PSW_aligner::write_alignment (const Alignment& align, const Alphabet& alphabet, double score, const char* score_units)
{
  // TO DO: trim off the gaps at the end, to reflect the fact that alignment is local
  cout << Stockholm_header;
  align.write_MUL (cout, alphabet);
  cout << Stockholm_file_annotation << " " Stockholm_bit_score_tag << " " << score << " " << score_units << "\n";
  cout << Stockholm_footer;
}


PSW_trainer::PSW_trainer (Opts_list& opts, Simple_PSW_params& params)
  : opts(opts), params(params), alphabet(params.alphabet)
{
  opts.newline();
  opts.add ("mi -max_iter", max_iter = 10, "\tmaximum no. of Baum-Welch iterations");
  opts.add ("sym -symmetrise", symmetrise = TRUE, "\tsymmetrise model by flipping alignments");
  opts.add ("ia -ignorealign", ignore_align = FALSE, "\tignore alignments provided; sum over all alignments instead");

  add_simple_PSW_references (opts);  // add some help text
}

void PSW_trainer::run()
{
  opts.parse_or_die();  // parse the command-line options
  try
    {
      // get alignment-list filename
      const sstring align_list_filename = opts.args[0];
      
      // open alignment-list file
      ifstream align_list_file (align_list_filename.c_str());
      if (!align_list_file) THROWEXPR ("Couldn't open alignment list file " << align_list_filename);

      // read in alignments
      Sequence_database seq_db;
      vector<Alignment> align_list;
      while (align_list_file)
	{
	  sstring align_name;
	  align_name.getline (align_list_file);
	  align_name.chomp();
	  if (align_name.size())
	    {
	      // read in the alignment
	      ifstream align_file (align_name.c_str());
	      if (!align_file) THROWEXPR ("Couldn't open alignment file: " << align_name);
	      Alignment align;
	      align.read_MUL (align_file, seq_db);
	      // check it's a pairwise alignment
	      if (align.rows() != 2)
		THROWEXPR ("Alignment '" << align_name << "' is not a pairwise alignment\n");
	      // add the alignment to the list
	      align_list.push_back (align);
	      // flip parent-child rows & add to list (ensures model is reversible)
	      if (symmetrise)
		{
		  align.swap_rows (0, 1);
		  align_list.push_back (align);
		}
	    }
	}
      seq_db.seqs_update (alphabet);  // convert text sequences into internal representations

      // optimise null model from sequence database, with Laplace +1 pseudocounts
      params.null_emit = seq_db.get_null_model (alphabet.size(), 1.);
      params.null_extend = seq_db.get_null_extend (1.);
      
      // create necessary objects for PSW HMM
      // this is a bit convoluted, because we are using the DART PGroup machinery;
      // the advantage of doing things this way is that we can do automated Baum-Welch training
      Simple_PSW_pscores pscores (params);  // create symbolic representations of params; convert probabilities to scores
      const Simple_PSW_PHMM symbolic_hmm (pscores, FALSE);  // create symbolic representation of Pair HMM (NB not odds-ratio)

      // do Baum-Welch
      Loge best_loglike = -InfinityLoge;
      int iteration = 0;
      while (1)
	{
	  // update iteration counter
	  ++iteration;

	  // print iteration counter
	  CTAG(7,PSWTRAIN) << "Iteration #" << iteration;
	  if (iteration > 1)
	    CL << "; best log-likelihood " << Nats2Bits(best_loglike) << " bits";
	  CL << "\n";

	  // create an HMM from the current PSW param values
	  const Pair_HMM_scores hmm = symbolic_hmm.eval_hmm_scores (pscores);

	  // zero Baum-Welch counts for HMM transitions and emissions
	  Pair_HMM_counts hmm_counts (hmm);

	  // loop over alignments, collecting Baum-Welch counts for HMM
	  Loge loglike = 0.;
	  for_const_contents (vector<Alignment>, align_list, align)
	    {
	      Score align_sc;
	      if (ignore_align)
		align_sc = hmm_counts.add_counts_from_unaligned_sequences (hmm, *align->prof[0], *align->prof[1]);
	      else
		align_sc = hmm_counts.add_counts_from_aligned_sequences (hmm, *align->prof[0], *align->prof[1], align->path);
	      NatsPMulAcc (loglike, Score2Nats (align_sc));
	    }

	  // convert HMM transition/emission counts into PSW parameter counts
	  PCounts pcounts (pscores);
	  symbolic_hmm.inc_var_counts (pcounts, pscores, hmm_counts);

	  // optimise PSW parameters
	  pscores.optimise_from_counts (pcounts);
	  pscores.update_params (params);

	  // check if log-likelihood did not increase; if so, exit Baum-Welch loop
	  if (loglike <= best_loglike)
	    {
	      CTAG(7,PSWTRAIN) << "Log-likelihood went from " << Nats2Bits(best_loglike) << " bits to " << Nats2Bits(loglike) << "; stopping\n";
	      break;
	    }
	  best_loglike = loglike;

	  // check if exceeded max no. of iterations; if so, exit gracefully
	  if (iteration >= max_iter)
	    {
	      CTAG(7,PSWTRAIN) << "Passed " << max_iter << " iterations; stopping\n";
	      break;
	    }
	}
      
      // output optimal PSW parameters
      params.write_params (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
