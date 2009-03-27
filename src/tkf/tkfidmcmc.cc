#include "tkf/tkfsequence.h"
#include "tkf/tkfopts.h"
#include "tkf/tkfdata.h"

#include "util/logfile.h"
#include "util/rnd.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  TKF_opts opts (argc, argv, 1);
  opts.syntax = "[options] <index file>";
  opts.short_description = "estimate TKF model indel rates using MCMC";
  opts.expect_args = 1;

  double mu_ceiling;
  int rounds;
  double seqlen_sample_rate;
  double mu_hop;
  double seqlen_hop;

  opts.newline();
  opts.print ("TKF EM parameters (Thorne et al, 1991)\n");
  opts.print ("--------------------------------------\n");
  opts.add ("md -maxdel", mu_ceiling = 1, "maximum deletion rate");
  opts.add ("r -rounds", rounds = 1000, "number of rounds of MCMC");
  opts.add ("sl -samplelength", seqlen_sample_rate = .5, "rate with which sequence length parameter (=insrate/delrate) is sampled");
  opts.add ("dh -delhop", mu_hop = .1, "maximum 'hop' size for deletion rate");
  opts.add ("slh -seqlenhop", seqlen_hop = .01, "maximum 'hop' size for sequence length parameter");

  opts.parse_or_die();
  try
    {
      // get args
      const char* index_filename = opts.args[0].c_str();

      // load database
      Sequence_database seq_db;
      Tree_alignment_database align_db (seq_db, index_filename);

      // digitise sequences
      const Alphabet& alph = seq_db.detect_alphabet();
      opts.use_pam = &alph == &Protein_alphabet;  // ensures that TKF_params.submat_factory has the correct alphabet
      seq_db.seqs2scores (alph);

      // get initial TKF_params
      TKF_params params = opts.params();

      // create TKF_align objects
      vector<TKF_align*> tkf_vec;
      Score old_score = 0;
      for (int i = 0; i < align_db.size(); ++i)
	{
	  // print name of alignment
	  CLOG(7) << "Building alignment: " << align_db.name[i] << "\n";
	  const Tree_alignment& tree_align = *align_db.tree_align[i];
	  if (CLOGGING(3)) tree_align.align.write_MUL (CL, alph);

	  // create TKF_align object
	  TKF_align* tkf = new TKF_align (params);
	  tkf->set_tree (tree_align.tree);
	  tkf->set_alignment (tree_align.align);
	  tkf->build_maps_from_names();
	  tkf->optimise_missing_nodes (FALSE);   // estimate sequences at internal nodes

	  // add to list
	  tkf_vec.push_back (tkf);

	  // get score
	  ScorePMulAcc (old_score, tkf->alignment_path_score());
	}

      // main loop
      for (int iter = 1; iter <= rounds; ++iter)
	{
	  // print log message
	  CLOG(7) << "Beginning MCMC round #" << iter << "\n";

	  // generate new candidate parameters
	  TKF_params old_params (params);
	  TKF_params new_params (params);
	  double mul;
	  if (Rnd::decide (seqlen_sample_rate))
	    {
	      CLOG(4) << "Resampling sequence length parameter\n";
	      const Prob old_seqlen = new_params.lambda / new_params.mu;

	      const Prob seqlen_min = max (old_seqlen - seqlen_hop, 0.);
	      const Prob seqlen_max = min (old_seqlen + seqlen_hop, 1.);
	      const Prob new_seqlen = Rnd::prob() * (seqlen_max - seqlen_min) + seqlen_min;

	      new_params.lambda = new_params.mu * new_seqlen;

	      const Prob new_seqlen_min = max (new_seqlen - seqlen_hop, 0.);
	      const Prob new_seqlen_max = min (new_seqlen + seqlen_hop, 1.);
	      mul = (seqlen_max - seqlen_min) / (new_seqlen_max - new_seqlen_min);
	    }

	  else
	    {
	      CLOG(4) << "Resampling deletion rate parameter\n";
	      const Prob old_seqlen = new_params.lambda / new_params.mu;

	      const Prob mu_min = max (old_params.mu - mu_hop, 0.);
	      const Prob mu_max = min (old_params.mu + mu_hop, mu_ceiling);
	      const Prob new_mu = Rnd::prob() * (mu_max - mu_min) + mu_min;

	      new_params.mu = new_mu;
	      new_params.lambda = old_seqlen * new_mu;

	      const Prob new_mu_min = max (new_mu - mu_hop, 0.);
	      const Prob new_mu_max = min (new_mu + mu_hop, mu_ceiling);
	      mul = (mu_max - mu_min) / (new_mu_max - new_mu_min);
	    }
	  CLOG(6) << "Proposed new parameters: insrate = " << new_params.lambda << ", delrate = " << new_params.mu << " (ins/del = " << new_params.lambda / new_params.mu << ")\n";

	  // get new score
	  params = new_params;
	  Score new_score = 0;
	  for (int i = 0; i < align_db.size(); ++i)
	    {
	      // reset the score cache in the TKF_align object
	      tkf_vec[i]->update_seq_scores();
	      tkf_vec[i]->tree_changed();
	      // calculate new score
	      const Score align_score = tkf_vec[i]->alignment_path_score();
	      // accumulate and print
	      ScorePMulAcc (new_score, align_score);
	      CLOG(5) << "Alignment " << align_db.name[i] << " has score " << Score2Bits(align_score) << " bits\n";
	    }
	  CLOG(6) << "New score for alignment database is " << Score2Bits(new_score) << " bits (old score " << Score2Bits(old_score) << " bits)\n";

	  // reject?
	  const Prob accept_prob = mul * Score2Prob (ScorePMul (new_score, -old_score));
	  if (accept_prob < 1.)
	    {
	      if (Rnd::decide (1. - accept_prob))
		{
		  params = old_params;
		  CLOG(6) << "Move rejected; P(accept) = " << accept_prob << "\n";
		}
	      else
		{
		  old_score = new_score;
		  CLOG(6) << "Move accepted with probability " << accept_prob << "\n";
		}
	    }
	  else
	    {
	      old_score = new_score;
	      CLOG(6) << "Move accepted with probability 1\n";
	    }

	  // print parameters & score to standard output
	  cout << params.lambda << " " << params.mu << " " << Score2Bits(old_score) << "\n";
	}

      // free TKF_align objects
      for_contents (vector<TKF_align*>, tkf_vec, tkf)
	delete *tkf;
      }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
