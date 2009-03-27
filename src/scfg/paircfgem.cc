#include "scfg/paircfgem.h"

Pair_CFG_trainer::Pair_CFG_trainer (const Stockholm_database& stock, Pair_PCFG& cfg, Dirichlet_prior& prior, PScores& pscore, bool local)
  : stock(stock), cfg(cfg), pscore(pscore), prior(prior), local(local)
{ }

Pair_CFG_trainer::Pair_CFG_trainer (const Stockholm_database& stock, Trainable_PCFG& cfg, bool local)
  : stock(stock), cfg(cfg), pscore(cfg.pscore), prior(cfg.default_prior()), local(local)
{ }

void Pair_CFG_trainer::train_pairwise (double min_inc, int max_rounds, int forgive)
{
  // print initial log message
  if (CTAGGING(3,CFGEM)) {
    CL << "Training the following pair CFG by EM:\n";
    cfg.show (CL);
    CL << "PScores before EM:\n";
    pscore.show (CL);
  }
  // constants
  const sstring ss_annot_string (Stockholm_secondary_structure_tag);  // "SS"
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
	  CTAG(5,CFGEM) << "Stopping after " << round << " rounds of EM\n";
	  break;
	}
      CTAG(5,CFGEM) << "Beginning EM round " << round+1 << "\n";
      // create Pair_CFG_scores for current parameters
      const Pair_CFG_scores cfg_scores = cfg.eval_cfg_scores (pscore);
      if (CTAGGING(4,CFGEM)) {
	CL << "Pair CFG for EM round " << round+1 << ":\n";
	cfg_scores.show (CL);
      }
      // accumulate counts for whole database
      PCounts pcounts (pscore);
      Loge total_loglike = 0;
      int n_align = 0;
      for_const_contents (list<Stockholm>, stock.align, align)
	{
	  ++n_align;  // update counter
	  // check that alignment is pairwise
	  if (align->rows() < 2)
	    {
	      CTAG(5,CFGEM) << "Warning: alignment " << n_align << " has " << align->rows() << " rows; skipping\n";
	      continue;
	    }
	  else if (align->rows() > 2)
	    CTAG(5,CFGEM) << "Warning: alignment " << n_align << " has " << align->rows() << " rows; I will only use the first two\n";
	  // print log message
	  CTAG(5,CFGEM) << "Preparing envelopes for alignment " << n_align << ": '" << align->row_name[0] << "' vs '" << align->row_name[1] << "'\n";
	  // create pairwise envelope
	  const Pairwise_path pair_path (align->path);
	  Pair_envelope pair_env (pair_path);
	  // create fold envelopes
	  vector<Fold_envelope> env (2);
	  for (int r = 0; r < 2; ++r)
	    {
	      const sstring fold_string = ((Stockholm&) *align).get_fold_string (align->row_name[r]);
	      if (fold_string.size())
		env[r].initialise_from_fold_string (fold_string);
	      else  // no SS, so make full fold envelope
		{
		  CTAG(5,CFGEM) << "Warning: found no secondary structure for alignment " << n_align << ", sequence " << align->row_name[r] << "\n";
		  env[r].initialise_full (align->path.count_steps_in_row (r));
		}
	    }
	  // fill inside-outside matrix; add final score; accumulate counts for this alignment
	  Pair_inside_outside_matrix inout (*align->np[0], *align->np[1],
					    env[0], env[1],
					    cfg_scores, pair_env,
					    local);
	  const Score sc = inout.inside.final_score;
	  const double weight = align->get_alignment_weight();
	  if (sc > -InfinityScore)
	    {
	      CTAG(5,CFGEM) << "Adding counts for alignment " << n_align << ": '" << align->np[0]->name << "' vs '" << align->np[1]->name << "' (log-likelihood " << Score2Bits(sc) << " bits; alignment weight is " << weight << ")\n";
	      NatsPMulAcc (total_loglike, Score2Nats(sc) * weight);
	      cfg.inc_var_counts (pcounts, pscore, inout.count, weight);
	      if (CTAGGING(-1,CFG_CUMULATIVE_PCOUNTS))
		{
		  CL << "Cumulative PVar counts:\n";
		  pcounts.show (CL);
		}
	    }
	}
      // report counts & score; check for improvement
      CTAG(5,CFGEM) << "Total (weighted) log-likelihood is " << Nats2Bits(total_loglike) << " bits\n";
      if (CTAGGING(4,CFGEM)) { CL << "Parameter counts:\n"; pcounts.show (CL); }
      const double inc = abs ((total_loglike - prev_loglike) / (abs(prev_loglike) < TINY ? 1. : prev_loglike));
      if (total_loglike <= prev_loglike || inc < min_inc)
	{
	  if (total_loglike < prev_loglike)
	    CTAG(7,CFGEM) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_loglike)
			  << " to " << Nats2Bits(total_loglike) << " bits during EM\n";

	  if (++dec > forgive)
	    {
	      CTAG(7,CFGEM) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
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

void Pair_CFG_trainer::train_single (double min_inc, int max_rounds, int forgive)
{
  // print initial log message
  if (CTAGGING(3,CFGEM)) {
    CL << "Training the following single-sequence CFG by EM:\n";
    cfg.show (CL);
    CL << "PScores before EM:\n";
    pscore.show (CL);
  }
  // check it's a single SCFG
  if (!cfg.is_single_CFG())
    THROWEXPR ("Attempt to use a single-sequence training method on a pair CFG");
  // check alignment database isn't empty
  if (stock.size() == 0)
    THROWEXPR ("Alignment database is empty");
  // constants
  const sstring ss_annot_string (Stockholm_secondary_structure_tag);  // "SS"
  // dummy Y sequence & fold envelope
  const Named_profile dummy_np;
  const Fold_envelope dummy_env;
  // do EM
  Loge prev_loglike = -InfinityLoge;
  int dec = 0;
  for (int round = 0; ; ++round)
    {
      // check if max no. of rounds reached
      if (max_rounds > 0 && round >= max_rounds)
	{
	  CTAG(5,CFGEM) << "Stopping after " << round << " rounds of EM\n";
	  break;
	}
      CTAG(5,CFGEM) << "Beginning EM round " << round+1 << "\n";
      // create Pair_CFG_scores for current parameters
      const Pair_CFG_scores cfg_scores = cfg.eval_cfg_scores (pscore);
      if (CTAGGING(4,CFGEM)) {
	CL << "Pair CFG for EM round " << round+1 << ":\n";
	cfg_scores.show (CL);
      }
      // accumulate counts for whole database
      PCounts pcounts (pscore);
      Loge total_loglike = 0;
      int n_align = 0;
      for_const_contents (list<Stockholm>, stock.align, align)
	{
	  ++n_align;  // update counter
	  // loop through each sequence in alignment
	  for (int r = 0; r < align->rows(); ++r)
	    {
	      // print log message
	      const sstring& seq_name = align->row_name[r];
	      const int seq_len = align->path.count_steps_in_row (r);
	      CTAG(5,CFGEM) << "Preparing envelopes for alignment " << n_align << ", sequence " << seq_name << "\n";
	      // create pairwise envelope
	      const Pair_envelope pair_env (seq_len, 0, 1);
	      // create fold envelope
	      Fold_envelope env;
	      const sstring fold_string = ((Stockholm&) *align).get_fold_string (align->row_name[r]);
	      if (fold_string.size())
		env.initialise_from_fold_string (fold_string);
	      else  // no SS, so make full fold envelope
		{
		  CTAG(5,CFGEM) << "Warning: found no secondary structure for alignment " << n_align << ", sequence " << seq_name << "\n";
		  env.initialise_full (seq_len);
		}
	      // fill inside-outside matrix; add final score; accumulate counts for this alignment
	      Pair_inside_outside_matrix inout (*align->np[r], dummy_np,
						env, dummy_env,
						cfg_scores, pair_env,
						local);
	      const Score sc = inout.inside.final_score;
	      const double weight = align->get_alignment_weight();
	      if (sc > -InfinityScore)
		{
		  CTAG(5,CFGEM) << "Adding counts for alignment " << n_align << ", row " << seq_name << " (log-likelihood " << Score2Bits(sc) << " bits; alignment weight is " << weight << ")\n";
		  NatsPMulAcc (total_loglike, Score2Nats(sc) * weight);
		  cfg.inc_var_counts (pcounts, pscore, inout.count, weight);
		  if (CTAGGING(-1,CFG_CUMULATIVE_PCOUNTS))
		    {
		      CL << "Cumulative PVar counts:\n";
		      pcounts.show (CL);
		    }
		}
	    }
	}
      // report counts & score; check for improvement
      CTAG(5,CFGEM) << "Total (weighted) log-likelihood is " << Nats2Bits(total_loglike) << " bits\n";
      if (CTAGGING(4,CFGEM)) { CL << "Parameter counts:\n"; pcounts.show (CL); }
      const double inc = abs ((total_loglike - prev_loglike) / (abs(prev_loglike) < TINY ? 1. : prev_loglike));
      if (total_loglike <= prev_loglike || inc < min_inc)
	{
	  if (total_loglike < prev_loglike)
	    CTAG(7,CFGEM) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_loglike)
			  << " to " << Nats2Bits(total_loglike) << " bits during EM\n";

	  if (++dec > forgive)
	    {
	      CTAG(7,CFGEM) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
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
