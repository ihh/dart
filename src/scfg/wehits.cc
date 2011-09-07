#include "scfg/wehits.h"

void WE_hits::get_hits (const Pair_CFG_scores& cfg, const Named_profile& np, bool local, int max_subseq_len, Score min_score, int max_hits)
{
  const set<int> paired_states = Pair_CFG_branch::get_paired_states (cfg.state_type);
  get_hits (cfg, paired_states, np, local, max_subseq_len, min_score, max_hits);
}

void WE_hits::get_hits (const Pair_CFG_scores& cfg, const set<int>& paired_states, const Named_profile& np, bool local, int max_subseq_len, Score min_score, int max_hits)
{
  CTAG(6,WEHITS) << "Folding sequence '" << np.name << "'\n";

  // create the dummy y-sequence
  const Fold_envelope dummy_yenv;
  const Named_profile dummy_npy = Named_profile();

  // create the subseq set
  Subseq_coords_set subseqs;
  for (int end = 0; end <= np.size(); ++end)
    for (int start = end; start >= 0 && ((end - start) < max_subseq_len || max_subseq_len < 0); --start)
      subseqs.insert (Subseq_coords (start, end - start));

  // loop through all hits
  for (int hit = 0; hit < max_hits || max_hits < 0; ++hit)
    {
      CTAG(5,WEHITS) << "Doing DP to find hit #" << hit+1 << "\n";

      // create fold envelope
      Fold_envelope env;
      env.initialise_from_subseq_coords (np.size(), subseqs);

      // create the matrix
      const Pair_CYK_matrix cyk (np, dummy_npy, env, dummy_yenv, cfg, local);
      if (CTAGGING(2,WEHITS_MATRIX))
	{
	  CL << "DP matrix for hit #" << hit+1 << ":\n";
	  CL << cyk.cfg_dump();
	}
      CTAG(5,WEHITS) << "DP matrix for hit #" << hit+1 << ": score " << cyk.final_score << "\n";

      // bail out if score too low
      if (cyk.final_score < min_score)
	{
	  CTAG(5,WEHITS) << "Score too low; bailing\n";
	  break;
	}

      // get the traceback
      const Pair_CFG_local_path lp = cyk.traceback_with_coords();
      const Pair_CFG_parse_tree parse_tree = cfg.parse (lp);

      // get subseqs
      Pair_CFG_parse_stats stats (cyk);
      stats.count_subseqs (parse_tree);

      // remove the traced-back subseqs from the fold envelope
      for_const_contents (Pair_CFG_parse_stats::Subseq_index_pair_count, stats.pair_count, sip_count)
	{
	  const Pair_CFG_parse_stats::Subseq_index_pair& subseq_index_pair = sip_count->first;
	  const Subseq_coords sc = cyk.xenv.subseq[subseq_index_pair.first];
	  const Subseq_coords_set::iterator p = subseqs.find (sc);
	  if (p == subseqs.end()) THROWEXPR ("Subseq (" << sc.start << "," << sc.end() << ") not found in envelope");
	  subseqs.erase (p);
	}

      // get start & end coords for this hit
      const Subseq& xsubseq = cyk.xenv.subseq[cyk.final_xsubseq_idx];

      // get the parse tree & alignment for this hit
      const Pair_CFG_alignment align = parse_tree.alignment (cfg.state_type, paired_states, np, dummy_npy);

      // create a GFF for this hit
      GFF gff;
      gff.seqname = np.name;
      gff.source = "DART";
      gff.feature = align.xfold;
      for (int i = 0; i < (int) gff.feature.size(); ++i)
	if (gff.feature[i] == NOFOLD_CHAR)
	  gff.feature[i] = np.seq[i + xsubseq.start];
      gff.start = xsubseq.start + 1;
      gff.end = xsubseq.end();
      gff.score = cyk.final_score;
      gff.strand = PlusStrand;
      gff.frame = NoFrame;

      const sstring xtmp (np.seq.begin() + xsubseq.start, np.seq.begin() + xsubseq.end());
      gff.group << hit + 1 << " " << xtmp;

      // add the GFF
      push_back (gff);
    }
}
