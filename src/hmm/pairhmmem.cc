#include "hmm/pairhmmem.h"
#include "hmm/pairenv.h"

Pair_HMM_trainer::Pair_HMM_trainer (const Stockholm_database& stock, Pair_PHMM& hmm, Dirichlet_prior& prior, PScores& pscore)
  : stock (stock), hmm (hmm), prior (prior), pscore (pscore),
    symmetrise (FALSE)
{ }

void Pair_HMM_trainer::train (double min_inc, int max_rounds, int forgive)
{
  // check alignment database isn't empty
  if (stock.size() == 0)
    THROWEXPR ("Alignment database is empty");
  // do EM
  Loge prev_loglike = -InfinityLoge;
  int dec = 0;
  for (int round = 0; ; ++round)
    {
      // check if max no. of rounds reached
      if (max_rounds > 0 && round >= max_rounds)
	{
	  CTAG(5,HMMEM) << "Stopping after " << round << " rounds of EM\n";
	  break;
	}
      CTAG(5,HMMEM) << "Beginning EM round " << round+1 << "\n";
      // create Pair_HMM_scores for current parameters
      const Pair_HMM_scores hmm_scores = hmm.eval_hmm_scores (pscore);
      // accumulate counts for whole database
      PCounts pcounts (pscore);
      Loge total_loglike = 0;
      int n_align = 0;
      for_const_contents (list<Stockholm>, stock.align, align)
	{
	  // update counter
	  ++n_align;
	  // count rows
	  const int rows = align->rows();
	  if (rows < 2)
	    {
	      CTAG(5,HMMEM) << "Warning: alignment " << n_align << " has " << rows << " row(s); skipping\n";
	      continue;
	    }
	  // get alignment weight
	  const double align_weight = align->get_alignment_weight() / (double) (rows - 1);
	  // loop through row pairs
	  for (int i = 0; i < rows; ++i)
	    for (int j = symmetrise ? 0 : i+1; j < rows; ++j)
	      if (i != j)
		{
		  // get row names
		  const sstring& iname = align->row_name[i];
		  const sstring& jname = align->row_name[j];
		  // print log message
		  CTAG(5,HMMEM) << "Preparing pair envelope for alignment " << n_align << ": row " << i+1 << " '" << iname << "' vs row " << j+1 << " '" << jname << "'\n";
		  // create pairwise envelope
		  const Pairwise_path pair_path (align->path, i, j, TRUE);
		  Pair_envelope pair_env (pair_path);
		  // fill forward-backward matrix; add final score; accumulate counts for this alignment
		  Pair_forward_backward_DP_matrix fwdback (hmm_scores,
							   *align->prof[i], *align->prof[j],
							   pair_env.allow_cut);
		  const Score sc = fwdback.forward_score;
		  if (sc > -InfinityScore)
		    {
		      CTAG(5,HMMEM) << "Adding counts for alignment " << n_align << ": row " << i+1 << " vs row " << j+1 << " (log-likelihood " << Score2Bits(sc) << " bits; alignment weight " << align_weight << ")\n";
		      NatsPMulAcc (total_loglike, Score2Nats(sc) * align_weight);
		      hmm.inc_var_counts (pcounts, pscore, fwdback.counts, align_weight);
		      if (CTAGGING(-1,HMM_CUMULATIVE_PCOUNTS))
			{
			  CL << "Cumulative PVar counts:\n";
			  pcounts.show (CL);
			}
		    }
		  else
		    CTAG(5,HMMEM) << "Skipping alignment " << n_align << " (log-likelihood was -infinity)\n";
		    
		}
	}
      // report counts & score; check for improvement
      CTAG(5,HMMEM) << "Total (weighted) log-likelihood is " << Nats2Bits(total_loglike) << " bits\n";
      if (CTAGGING(4,HMMEM)) { CL << "Parameter counts:\n"; pcounts.show (CL); }
      const double inc = abs ((total_loglike - prev_loglike) / (abs(prev_loglike) < TINY ? 1. : prev_loglike));
      if (total_loglike <= prev_loglike || inc < min_inc)
	{
	  if (total_loglike < prev_loglike)
	    CTAG(7,CFGEM) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_loglike)
			  << " to " << Nats2Bits(total_loglike) << " bits during EM\n";

	  if (++dec > forgive)
	    {
	      CTAG(7,HMMEM) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
	      break;
	    }
	}
      else
	dec = 0;
      prev_loglike = total_loglike;
      // use Dirichlet prior to update PScores
      prior.pscores = &pscore;
      prior.optimise (pcounts);
    }
}
