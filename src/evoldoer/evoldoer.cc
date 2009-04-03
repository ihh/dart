#include "evoldoer/tkfstructuretree.h"
#include "tkf/tkfhmm.h"
#include "scfg/postenv.h"
#include "scfg/paircfgdp.h"
#include "scfg/cfgdotplot.h"
#include "scfg/cfgsexpr.h"
#include "ecfg/ecfgsexpr.h"

int main (int argc, char** argv)
{
  try
    {
      // create default parameters
      TKFST_default_params params;

      // create opts
      INIT_OPTS_LIST (opts, argc, argv, -1, "[options] <sequence file(s)>",
		      "pairwise RNA structural alignment using TKF structure tree");

      sstring loop_subst_model_filename;
      sstring stem_subst_model_filename;
      sstring grammar_filename;

      opts.newline();
      opts.print_title ("TKF structure tree options");

      double t;

      opts.add ("t -time", t = 1, "evolutionary time");
      opts.add ("ll -loop-lambda", params.loop.lambda, "insertion rate in loops");
      opts.add ("lm -loop-mu",     params.loop.mu,     "deletion rate in loops");
      opts.add ("ls -loop-subst",  loop_subst_model_filename, "single-nucleotide substitution model for loops (xrate chain format)", false);
      opts.add ("sl -stem-lambda", params.stem.lambda, "insertion rate in stems");
      opts.add ("sm -stem-mu",     params.stem.mu,     "deletion rate in stems");
      opts.add ("ss -stem-subst",  stem_subst_model_filename, "covariant dinucleotide substitution model for stems (xrate chain format)", false);
      opts.add ("sp -stem-prob",   params.stem_prob,   "stem probability");
      opts.add ("g -grammar",      grammar_filename,   "save TKFST grammar to file (stemloc pair-SCFG format)", false);

      opts.newline();
      opts.print_title ("Envelope generation options");

      int nfold;
      int nalign;
      int max_subseq_len;
      int min_loop_len;
      int align_band;

      opts.add ("nf -nfold", nfold = -1, "number of folds to sample; -1 to unlimit");
      opts.add ("na -nalign", nalign = 100, "number of alignments to sample; -1 to unlimit");

      opts.add ("len -maxstemlen", max_subseq_len = -1, "maximum length of envelope subsequences; -1 to unlimit");
      opts.add ("loop -minlooplen", min_loop_len = 0, "minimum length of loop regions");
      opts.add ("band -bandalign", align_band = -1, "band size for pairwise alignment; -1 to unlimit");

      opts.newline();
      opts.print_title ("Posterior probability matrix options");

      sstring align_prefix, fold_prefix, palign_prefix, pfold_prefix;

      opts.add ("af -alignfile", align_prefix = "", "filename prefix for pre-alignment post.prob. matrices", false);
      opts.add ("ff -foldfile", fold_prefix = "", "filename prefix for pre-fold post.prob. matrices", false);
      opts.add ("paf -postalignfile", palign_prefix = "", "filename prefix for alignment post.prob. matrices", false);
      opts.add ("pff -postfoldfile", pfold_prefix = "", "filename prefix for fold post.prob. matrices", false);

      opts.newline();
      opts.print_title ("Common shorthand options");

      opts.newline();
      opts.add ("verbose", "-log 5 -log ALLOC -log CFGDP", "show progress of dynamic programming recursions");

      opts.newline();
      opts.newline();
      opts.print ("Input sequence file can be in either Stockholm or FASTA format.\n");
      opts.print ("If in Stockholm format, the corresponding multiple alignment is taken as fixed.");
      opts.newline();

      opts.newline();
      opts.print_title ("References");
      opts.print ("Holmes I, 2004. A probabilistic model for the evolution of RNA structure.\n");
      opts.print ("BMC Bioinformatics 5:166.   http://www.biomedcentral.com/1471-2105/5/166\n");

      // parse opts
      opts.parse_or_die();

      // get sequence filenames
      const vector<sstring>& seq_filename = opts.args;

      // create databases & alphabet
      const Alphabet& alphabet (CFG_alphabet);
      Stockholm_database align_db;
      FASTA_sequence_database seq_db;

      // check for at least one sequence file
      if (seq_filename.size() == 0 && grammar_filename.size() == 0)
	THROWEXPR (opts.short_help() << "You need to specify at least one sequence file");

      // read everything in as a Stockholm alignment
      for_const_contents (vector<sstring>, seq_filename, sf)
	{
	  ifstream in (sf->c_str());
	  if (!in) THROWEXPR ("Sequence file not found: '" << *sf << "'");
	  align_db.read_Stockholm_or_FASTA (in, seq_db);
	}
      // check names are unique
      seq_db.index.assert_names_unique (true);
      // conservatively propagate consensus folds
      align_db.propagate_consensus_folds (false);
      // update the sequence database index
      seq_db.update_index();
      seq_db.update_alphabet (&alphabet);
      seq_db.seqs_update (alphabet, (Profile_flags_enum::Profile_flags) (Profile_flags_enum::DSQ | Profile_flags_enum::SCORE));

      // read loop subst model
      if (loop_subst_model_filename.size())
	{
	  SExpr_file param_sexpr (loop_subst_model_filename.c_str());
	  ECFG_builder::init_chain_given_alphabet (params.loop_mat, alphabet, param_sexpr.sexpr, 1);
	}

      // read stem subst model
      if (stem_subst_model_filename.size())
	{
	  SExpr_file param_sexpr (stem_subst_model_filename.c_str());
	  ECFG_builder::init_chain_given_alphabet (params.stem_mat, alphabet, param_sexpr.sexpr, 2);
	}

      // create a TKF Pair HMM
      TKF_joint_pair_HMM_scores tkf_pair_hmm (params.loop, t, true, true);
      tkf_pair_hmm.alphabet = &CFG_alphabet;

      // get non-pad states
      set<int> tkf_pair_hmm_states;
      for (int s = 0; s < tkf_pair_hmm.states(); ++s)
	tkf_pair_hmm_states.insert (s);

      // create a TKF Pair HMM-emulating Pair_CFG_scores
      const Pair_CFG_scores tkf_pair_hmm_cfg (tkf_pair_hmm);

      // create a TKF structure tree Single SCFG
      const TKFST_single_CFG tkfst_single_cfg (params, true);

      // create a TKF structure tree Pair SCFG
      const TKFST_pair_CFG tkfst_pair_cfg (params, t, false, true);

      // save TKFST pair SCFG to file, if asked
      if (grammar_filename.size())
	{
	  const Pair_PCFG tkfst_pair_pcfg (tkfst_pair_cfg);
	  PScores dummy_pscores;
	  set<int> dummy_mutable_pgroups;
	  ofstream gfstream (grammar_filename.c_str());
	  PCFG_builder::grammar2stream (gfstream, tkfst_pair_pcfg, dummy_pscores, dummy_mutable_pgroups);
	}

      // prepare fold envelope database
      Fold_envelope_database foldenv_db;
      for (int i = 0; i < seq_db.size(); ++i)
	{
	  // get Named_profile
	  const Named_profile& np = *seq_db.index.profile[i];

	  // create Fold_envelope
	  Fold_envelope env;

	  // keep a flag saying whether we've got an envelope yet
	  bool got_foldenv = false;

	  // check if a structure was specified for this sequence; if so, don't try to pre-fold it
	  if (!got_foldenv)
	    {
	      const Stockholm& stock = *align_db.align_index[align_db.align_row[np.name].nAlign];
	      const sstring fold_string = stock.get_fold_string (np.name);
	      const bool got_ss = fold_string.size() > 0;
	      if (got_ss)
		{
		  // initialise envelope from fold string
		  env.initialise_from_fold_string (fold_string);
		  // show dotplot
		  if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PREFOLD))
		    {
		      CL << "Pre-specified fold envelope for sequence '" << np.name << "':\n";
		      env.render_dotplot (CL, np.seq, true);
		    }
		  got_foldenv = true;
		}
	    }

	  // if no structure was specified, do the pre-folding thing
	  if (!got_foldenv)
	    {
	      // initialise full, banded fold envelope
	      env.initialise_3prime_banded (np.size(), max_subseq_len, min_loop_len);
	      // show dotplot
	      if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_FOLDENV STEMLOC_BANDED_FOLDENV))
		{
		  CL << "Banded folding pre-envelope for '" << np.name << "':\n";
		  env.render_dotplot (CL, np.seq, true);
		}

	      // create dummy y-sequence & envelope
	      const Fold_envelope dummy_yenv;
	      const Named_profile dummy_npy;

	      // if unlimited structures, stick with banded envelope; otherwise, sample structures
	      if (nfold < 0)
		CTAG(6,STEMLOC_PREFOLD) << "Using banded fold envelope for sequence '" << np.name << "'\n";
	      else
		{
		  // if doing dotplots, create Inside-Outside matrix & print to file
		  if (fold_prefix.size())
		    {
		      const Pair_inside_outside_matrix in_out (np, dummy_npy, env, dummy_yenv, tkfst_single_cfg, false);
		      const PairCFG_fold_dotplot dotplot (in_out, 0);
		      dotplot.write_dotplot (fold_prefix, np.name);
		    }

		  // find the sampled fold envelope
		  CTAG(6,STEMLOC_PREFOLD) << "Pre-folding sequence '" << np.name << "'\n";
		  Sampled_fold_envelope sample_env;
		  Subseq_coords_count subseq_count;
		  // create CYK-KYC DP matrix
		  const Pair_CYK_KYC_matrix cyk_kyc (np, dummy_npy, env, dummy_yenv, tkfst_single_cfg, false);
		  if (CTAGGING(1,STEMLOC_DP))
		    cyk_kyc.show (CL);
		  // sample the envelope
		  subseq_count = sample_env.initialise_best (cyk_kyc, nfold, min_loop_len);
		  // assign the envelope
		  env.swap (sample_env);
		  // show dotplot
		  if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_FOLDENV STEMLOC_SAMPLED_FOLDENV))
		    {
		      CL << "Sampled fold envelope for '" << np.name << "':\n";
		      Fold_envelope::render_dotplot_from_counts (CL, subseq_count, np.seq, nfold, true);
		    }
		}
	    }

	  // print no. of subseqs to log
	  if (CTAGGING(6,STEMLOC_FOLDENV STEMLOC_ENVSIZE))
	    CL << "Fold envelope for " << np.name << " has " << env.subseqs() << " subseqs\n";

	  // store in database
	  foldenv_db[np.name] = env;  // another costly assignment; there are better ways to do this
	}

      // loop through sequence pairs
      for (int i = 0; i < seq_db.size(); ++i)
	for (int j = i + 1; j < seq_db.size(); ++j)
	  {
	    // get sequences
	    const Named_profile& npx = *seq_db.index.profile[i];
	    const Named_profile& npy = *seq_db.index.profile[j];

	    // retrieve fold envelopes
	    const Fold_envelope& envx = foldenv_db[npx.name];
	    const Fold_envelope& envy = foldenv_db[npy.name];

	    // create the sampled alignment envelope
	    Pair_envelope pair_env;

	    // print log message
	    CTAG(6,STEMLOC) << "Finding alignment envelope for sequences '" << npx.name << "' and '" << npy.name << "'\n";
	    
	    // keep a flag saying whether we've got the alignment yet
	    bool got_align = false;

	    // if an alignment between the two sequences was prespecified, then don't bother pre-aligning
	    if (!got_align)
	      {
		const Stockholm_database::Align_row& arx = align_db.align_row[npx.name];
		const Stockholm_database::Align_row& ary = align_db.align_row[npy.name];
		if (arx.nAlign == ary.nAlign)
		  {
		    const Pairwise_path path (align_db.align_index[arx.nAlign]->path, arx.nRow, ary.nRow, true);
		    pair_env = Pair_envelope (npx.size(), npy.size(), 0);
		    pair_env.add_pairwise_path (path);
		    got_align = true;
		    if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PREALIGN))
		      {
			CL << "Using pre-specified alignment for '" << npx.name << "' vs '" << npy.name << "':\n";
			pair_env.render_dotplot (CL, npx.seq, npy.seq, 1, true);
		      }
		  }
	      }
	    
	    // if no alignment envelope has been found yet, do the pre-alignment thing
	    if (!got_align)
	      {
		// initialise full, banded alignment envelope
		Pair_envelope banded_pair_env;
		banded_pair_env.initialise_banded (npx.size(), npy.size(), align_band);
		if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PAIRENV STEMLOC_BANDED_PAIRENV))
		  {
		    CL << "Banded alignment pre-envelope for '" << npx.name << "' vs '" << npy.name << "':\n";
		    banded_pair_env.render_dotplot (CL, npx.seq, npy.seq, 1, true);
		  }
		
		// if unlimited structures, copy banded envelope to sampled envelope; otherwise, sample structures
		if (nalign < 0)
		  {
		    CTAG(6,STEMLOC_PREALIGN) << "Using banded alignment envelope for '" << npx.name << "' and '" << npy.name << "'\n";
		    pair_env = banded_pair_env;
		  }
		else
		  {
		    // create HMM-emulating fold envelopes
		    Fold_envelope hmm_xenv, hmm_yenv;
		    hmm_xenv.initialise_3prime (npx.size());
		    hmm_yenv.initialise_3prime (npy.size());

		    // if doing dotplots, create Forward-Backward matrix & print to file
		    if (align_prefix.size())
		      {
			const Pair_inside_outside_matrix in_out (npx, npy, hmm_xenv, hmm_yenv, tkf_pair_hmm_cfg, banded_pair_env, false);
			const PairCFG_alignment_dotplot dotplot (in_out);
			dotplot.write_dotplot (align_prefix, npx.name, npy.name);
		      }

		    // find the sampled alignment envelope
		    CTAG(6,STEMLOC_PREALIGN) << "Pre-aligning '" << npx.name << "' and '" << npy.name << "'\n";

		    // create CYK-KYC DP matrix
		    Pair_CYK_KYC_matrix cyk_kyc (npx, npy, hmm_xenv, hmm_yenv, tkf_pair_hmm_cfg, banded_pair_env, false, true);

		    // sample the envelope
		    Sampled_pair_envelope sampled_pair_env (npx, npy, hmm_xenv, hmm_yenv);
		    sampled_pair_env.initialise_best (cyk_kyc, tkf_pair_hmm, nalign, tkf_pair_hmm_states);
		    pair_env = sampled_pair_env;

		    // show dotplot
		    if (CTAGGING(4,DOTPLOT STEMLOC_DOTPLOT STEMLOC_PAIRENV STEMLOC_SAMPLED_PAIRENV))
		      {
			CL << "Sampled alignment envelope for '" << npx.name << "' vs '" << npy.name << "':\n";
			pair_env.render_dotplot (CL, npx.seq, npy.seq, nalign, true);
		      }
		  }
	      }
	    if (CTAGGING(6,STEMLOC_PAIRENV STEMLOC_ENVSIZE))
	      CL << "Alignment envelope for '" << npx.name << "' vs '" << npy.name << "' has " << pair_env.cells() << " cells\n";

	    // if doing dotplots, create Inside-Outside matrix & print to file
	    if (palign_prefix.size() || pfold_prefix.size())
	      {
		const Pair_inside_outside_matrix in_out (npx, npy, envx, envy, tkfst_pair_cfg, pair_env, false);
		if (palign_prefix.size())
		  {
		    const PairCFG_alignment_dotplot align_dotplot (in_out);
		    align_dotplot.write_dotplot (palign_prefix, npx.name, npy.name);
		  }
		if (pfold_prefix.size())
		  {
		    const PairCFG_fold_dotplot xfold_dotplot (in_out, 0);
		    xfold_dotplot.write_dotplot (pfold_prefix, npx.name);

		    const PairCFG_fold_dotplot yfold_dotplot (in_out, 1);
		    yfold_dotplot.write_dotplot (pfold_prefix, npy.name);
		  }
	      }

	    // print log message
	    CTAG(6,STEMLOC) << "Aligning & folding '" << npx.name << "' and '" << npy.name << "'\n";

	    // do CYK
	    // build the DP matrix
	    const Pair_CYK_matrix cyk (npx, npy, envx, envy, tkfst_pair_cfg, pair_env, false);
	    if (CTAGGING(1,STEMLOC_DP))
	      CL << cyk.cfg_dump();
	    
	    // test for -inf score
	    if (cyk.final_score == -InfinityScore)
	      {
		CLOGERR << "Warning: alignment score for '" << npx.name << "' and '" << npy.name << "' is -infinity; skipping.\n";
		CLOGERR << " (this usually means there is no valid alignment of the sequences, which often means that the precomputed structures or alignments are too restrictive;\n";
		CLOGERR << "  try setting -nf or -na higher?)\n";
		break;
	      }

	    // get alignment
	    const Pair_CFG_alignment cyk_alignment = cyk.alignment();

	    // output
	    cyk_alignment.show (cout);
	  }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      return 1;
    }
  return 0;
}
