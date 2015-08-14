#include "scfg/cfgsexpr.h"
#include "stemloc/quickalign.h"
#include "seq/pkeywords.h"

#define REFLOC_HMM_TAG "hmm"
#define REFLOC_NULL_EMIT_PREFIX "null"

int main (int argc, char** argv) {
  try {
    INIT_OPTS_LIST (opts, argc, argv, 2,
		    "[options] <refs.fasta> <reads.fasta>",
		    "configurable pairwise alignment");

    opts.newline();
    opts.print_title ("General options");
    
    sstring model_filename, trained_model_filename, untrained_model_filename;
    bool do_alignment, use_revcomp;

    opts.newline();
    opts.add ("m -model", model_filename = "", "read model & parameters from file", false);
    opts.add ("t -train", trained_model_filename = "", "do training, output trained model & parameters to file", false);
    opts.add ("d -dump", untrained_model_filename = "", "before doing any training, save model & parameters to file", false);
    opts.add ("a -align", do_alignment = true, "align sequences & output alignment", true);
    opts.add ("r -revcomp", use_revcomp = false, "also use reverse complement of refs", true);

    opts.parse_or_die();
    const sstring& refs_filename = opts.args[0];
    const sstring& reads_filename = opts.args[1];

    Alphabet& alphabet = CFG_alphabet;  // hardwired for now
    Pair_PHMM hmm (0, alphabet);
    PScores pscores;
    Alphabet_group null_emit;
    set<int> pad_states;
    set<int> mutable_pgroups;

    Quick_align* quick_align = NULL;

    if (model_filename.size()) {
      SExpr_file model_file (model_filename.c_str());
      hmm = PCFG_builder::init_hmm (model_file.sexpr[0], pscores, null_emit, pad_states, mutable_pgroups, false);

    } else {
      null_emit = pscores.new_alphabet_group (alphabet, REFLOC_NULL_EMIT_PREFIX);
      for (int i = 0; i < alphabet.size(); ++i) {
	pscores.group_suffix[null_emit.group_idx][i].clear();
	pscores.group_suffix[null_emit.group_idx][i] << REFLOC_NULL_EMIT_PREFIX << alphabet.int2char_uc(i);
      }
      quick_align = new Quick_align (pscores, null_emit);
      mutable_pgroups = quick_align->mutable_pgroups;
      set<int> non_pad = quick_align->non_pad_states();
      for (int s = 0; s < quick_align->states(); ++s)
	if (non_pad.find(s) == non_pad.end())
	  pad_states.insert (s);
      hmm = *quick_align;

      Dirichlet_prior prior (pscores);
      prior.assign_Laplace (pscores, mutable_pgroups);
      prior.assign_Laplace (null_emit);
      prior.initialise();
    }

    FASTA_sequence_database refs_db (refs_filename.c_str(), &alphabet);
    FASTA_sequence_database reads_db (reads_filename.c_str(), &alphabet);

    if (use_revcomp)
      refs_db.append_sequence_database (refs_db.revcomp (alphabet));
    
    if (untrained_model_filename.size()) {
      ofstream out (untrained_model_filename.c_str());
      PCFG_builder::hmm2stream (out, hmm, pscores, mutable_pgroups, REFLOC_HMM_TAG, "", null_emit.group_idx, &pad_states);
    }
    
    if (trained_model_filename.size()) {
      // do unsupervised training
      Dirichlet_prior prior (pscores);
      prior.assign_Laplace (pscores, mutable_pgroups);
      
      Score prev_sc = -InfinityScore;
      for (int round = 0; ; ++round) {
	CTAG(5,REFLOC) << "Beginning EM round " << round+1 << "\n";
	// create Pair_HMM_scores for current parameters
	const Pair_HMM_scores hmm_scores = hmm.eval_hmm_scores (pscores);
	if (CTAGGING(4,REFLOC))
	  hmm_scores.show (CL);
	// accumulate counts for whole database
	PCounts pcounts (pscores);
	Score total_sc = 0;
	for (int n_read = 0; n_read < reads_db.size(); ++n_read) {
	  const Score_profile& read_prof = reads_db.get_seq(n_read).prof_sc;
	  vector<Score> ref_sc;
	  vector<Pair_HMM_counts> ref_counts;
	  Score read_sc = -InfinityScore;
	  for (int n_ref = 0; n_ref < refs_db.size(); ++n_ref) {
	    const Score_profile& ref_prof = refs_db.get_seq(n_ref).prof_sc;
	    Pair_forward_backward_DP_matrix mx (hmm_scores, ref_prof, read_prof);
	    ref_sc.push_back (mx.forward_score);
	    ref_counts.push_back (mx.counts);
	    ScorePSumAcc (read_sc, mx.forward_score);
	    CTAG(3,REFLOC) << "Read " << reads_db.index.name[n_read] << ", reference " << refs_db.index.name[n_ref] << ": log-likelihood " << Score2Bits(mx.forward_score) << " bits\n";
	  }
	  for (int n_ref = 0; n_ref < refs_db.size(); ++n_ref)
	    hmm.inc_var_counts (pcounts, pscores, ref_counts[n_ref], Score2Prob (ScorePMul (ref_sc[n_ref], -read_sc)));
	  ScorePMulAcc (total_sc, read_sc);
	  CTAG(4,REFLOC) << "Read " << reads_db.index.name[n_read] << ": log-likelihood " << Score2Bits(read_sc) << " bits\n";
	}
	// strictly speaking we should include the Dirichlet prior in the EM score
	// but this is a bug woven deep into DART at this point. hmmmmm
	CTAG(6,REFLOC) << "EM iteration " << round+1 << ": log-likelihood " << Score2Bits(total_sc) << " bits\n";
	if (round > 0 && total_sc <= prev_sc)
	  break;
	prev_sc = total_sc;
	prior.optimise (pcounts);
      }

      ofstream out (trained_model_filename.c_str());
      PCFG_builder::hmm2stream (out, hmm, pscores, mutable_pgroups, REFLOC_HMM_TAG, " ", null_emit.group_idx, &pad_states);
    }
  
    if (do_alignment) {
	const Pair_HMM_scores hmm_scores = hmm.eval_hmm_scores (pscores);
	for (int n_read = 0; n_read < reads_db.size(); ++n_read) {
	  const Score_profile& read_prof = reads_db.get_seq(n_read).prof_sc;
	  Score best_sc = -InfinityScore;
	  int n_best_ref = -1;
	  vector<int> best_state_path;
	  for (int n_ref = 0; n_ref < refs_db.size(); ++n_ref) {
	    const Score_profile& ref_prof = refs_db.get_seq(n_ref).prof_sc;
	    Pair_Viterbi_DP_matrix mx (hmm_scores, ref_prof, read_prof);
	    if (n_ref == 0 || mx.final_score > best_sc) {
	      best_sc = mx.final_score;
	      n_best_ref = n_ref;
	      best_state_path = mx.optimal_state_path();
	    }
	  }
	  // output best alignment for this read
	  if (n_best_ref >= 0) {
	    const Pairwise_path path = hmm_scores.convert_state_path_to_alignment (best_state_path);
	    const Alignment align (path, refs_db.get_seq(n_best_ref), reads_db.get_seq(n_read));  // label the alignment
	    align.write_MUL (cout, alphabet);
	    cout << "#=SC " << Score2Bits(best_sc) << " bits\n";
	    cout << Stockholm_alignment_separator << "\n";
	  }
	}
    }

    if (quick_align)
      delete quick_align;
  
  } catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
