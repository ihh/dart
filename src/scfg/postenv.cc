#include "scfg/postenv.h"

// Pair_CFG_state_set
Pair_CFG_state_set::Pair_CFG_state_set (const Pair_CFG_state_typing& typing) : typing(typing)
{ }

void Pair_CFG_state_set::initialise_full()
{
  states.clear();
  for (int s = 0; s < typing.states(); ++s) states.insert (s);
}

// Sampled_fold_envelope
Sampled_fold_envelope::Sampled_fold_envelope() : Fold_envelope(0)
{ }

void Sampled_fold_envelope::initialise_CYK (const Pair_CYK_matrix& cyk, int min_loop_len)
{
  CTAG(5,POSTENV POSTENV_FOLD) << "Finding best RNA secondary structure for fold envelope\n";

  // do the traceback
  const Pair_CFG_parse_tree parse_tree = cyk.parse_tree();

  // count the subseqs
  Pair_CFG_parse_stats stats (cyk);
  stats.count_subseqs (parse_tree);

  // convert into single-sequence counts
  Subseq_coords_set scs;
  for_const_contents (Pair_CFG_parse_stats::Subseq_index_pair_count, stats.pair_count, sipc)
    scs.insert (cyk.xenv.subseq[sipc->first.first]);

  // make the envelope
  initialise_from_subseq_coords_with_cflag_mask (scs, cyk.xenv, min_loop_len);
}

Subseq_coords_count Sampled_fold_envelope::initialise_sampled (const Pair_inside_matrix& imx, int n_samples, int min_loop_len)
{
  CTAG(6,POSTENV POSTENV_FOLD) << "Sampling " << n_samples << " RNA secondary structures for fold envelope\n";

  // make the cell sorter
  Pair_inside_cell_sorter cell_sorter (imx);

  // get ready to collect subseq stats
  Pair_CFG_parse_stats stats (imx);

  // print log message
  CTAG(4,POSTENV POSTENV_FOLD) << "Sequence   " << imx.npx.seq << "\n";

  // do the tracebacks
  for (int n = 0; n < n_samples; ++n)
    {
      // get parse tree and count it
      const Pair_CFG_parse_tree parse_tree = cell_sorter.parse_tree();
      stats.count_subseqs (parse_tree);

      // log it
      if (CTAGGING(4,POSTENV POSTENV_FOLD))
	{
	  const Pair_CFG_alignment align = parse_tree.alignment (imx.cfg.state_type, imx.npx, imx.npy);
	  CL << "Fold ";
	  CL.width (5);
	  CL << n+1 << ' ';
	  for (int i = 0; i < align.xstart; ++i) CL << ' ';
	  CL << align.xfold << '\n';
	}
    }

  // count subseqs
  Subseq_coords_set scs;
  Subseq_coords_count subseq_count;
  for_const_contents (Pair_CFG_parse_stats::Subseq_index_pair_count, stats.pair_count, sipc)
    {
      const Subseq_coords& coords = imx.xenv.subseq[sipc->first.first];
      scs.insert (coords);
      subseq_count[coords] = sipc->second;
    }

  // make the envelope
  initialise_from_subseq_coords_with_cflag_mask (scs, imx.xenv, min_loop_len);

  // return the counts
  return subseq_count;
}

Subseq_coords_count Sampled_fold_envelope::initialise_best (const Pair_CYK_KYC_matrix& cyk_kyc, int n_best, int min_loop_len)
{
  CTAG(6,POSTENV POSTENV_FOLD) << "Finding " << n_best << " best RNA secondary structures for fold envelope\n";

  // get references to CYK & KYC matrices
  const Pair_CYK_matrix& cyk = cyk_kyc.cyk;
  const Pair_KYC_matrix& kyc = cyk_kyc.kyc;

  // make the cell sorter
  Pair_CYK_KYC_cell_sorter cell_sorter (cyk_kyc);

  // get ready to collect subseq stats
  Pair_CFG_parse_stats stats (cyk);

  // print log message
  CTAG(4,POSTENV POSTENV_FOLD) << "Sequence   " << cyk.npx.seq << "\n";

  // do the tracebacks
  int n_tracebacks = 0;
  Score last_score = InfinityScore;
  for_const_contents (vector<Pair_CFG_DP_matrix::Cell_coords>, cell_sorter.sorted_cells, cell)
    {
      // describe cell
      if (CTAGGING(-1,POSTENV_CELLS))
	CL << "(considering state " << cell->state << ", xsubseq " << cell->xsubseq << ")\n";

      // skip this cell if we've already seen it
      const Pair_CFG_parse_stats::Subseq_index_pair sip (cell->xsubseq, cell->ysubseq);
      if (stats.pair_count.find (sip) != stats.pair_count.end())
	continue;

      // bail out if score drops too low & we've done enough tracebacks
      const Score sc = cyk_kyc.max_state_sc (cell->state, cell->xsubseq, cell->ysubseq);
      if (sc <= -InfinityScore || (n_tracebacks >= n_best && sc < last_score))
	break;

      // get the traceback
      vector<int> path;
      kyc.traceback_from (cell->xsubseq, cell->ysubseq, cell->state, path);

      // build the parse tree
      const Pair_CFG_local_path local_path (0, cyk.xlen, 0, cyk.ylen, path);
      const Pair_CFG_parse_tree parse_tree = cyk.cfg.parse (local_path);

      // count it
      stats.count_subseqs (parse_tree);

      // log it
      if (CTAGGING(4,POSTENV POSTENV_FOLD))
	{
	  const Pair_CFG_alignment align = parse_tree.alignment (cyk.cfg.state_type, cyk.npx, cyk.npy);
	  CL << "Fold ";
	  CL.width (5);
	  CL << n_tracebacks+1 << ' ';
	  for (int i = 0; i < align.xstart; ++i) CL << ' ';
	  CL << align.xfold;

	  // it is poor practise to embed a test in working code like this, but what the hell
	  const Score test_sc = cyk.cfg.path_score (parse_tree, cyk.npx.dsq, cyk.npy.dsq);
	  CL << " (score " << Score2Bits(sc) << " bits)\n";
	  if (sc != test_sc)
	    {
	      CL << "WARNING: Pair_CFG says score is " << Score2Bits(test_sc) << " bits. Traceback was from state " << cell->state << ", ";
	      CL << "xsubseq " << cell->xsubseq << " (" << cyk.xenv.subseq[cell->xsubseq].start << '+' << cyk.xenv.subseq[cell->xsubseq].len << "), ";
	      CL << "ysubseq " << cell->ysubseq << " (" << cyk.yenv.subseq[cell->ysubseq].start << '+' << cyk.yenv.subseq[cell->ysubseq].len << ").\n";
	      CL << "CYK expanded trace:\n";
	      cyk.show_expanded_trace (parse_tree, CL, FALSE);
	      CL << "KYC expanded trace:\n";
	      kyc.show_expanded_trace (parse_tree, CL, TRUE);
	      CL << "CFG score breakdown:\n";
	      cyk.cfg.show_score_breakdown (parse_tree, cyk.npx.dsq, cyk.npy.dsq, CL);
	    }
	}

      // update the number of tracebacks done & the last score
      ++n_tracebacks;
      last_score = sc;
    }

  // count subseqs
  Subseq_coords_set scs;
  Subseq_coords_count subseq_count;
  for_const_contents (Pair_CFG_parse_stats::Subseq_index_pair_count, stats.pair_count, sipc)
    {
      const Subseq_coords& coords = cyk.xenv.subseq[sipc->first.first];
      scs.insert (coords);
      subseq_count[coords] = sipc->second;
    }

  // make the envelope
  initialise_from_subseq_coords_with_cflag_mask (scs, cyk.xenv, min_loop_len);

  // return the counts
  return subseq_count;
}

// Sampled_pair_envelope
Sampled_pair_envelope::Sampled_pair_envelope()
  : Pair_envelope(0,0,0), npx(0), npy(0), xenv(0), yenv(0),
    paths_added (0),
    cells_needed (0),
    cells_available (0)
{ }

Sampled_pair_envelope::Sampled_pair_envelope (const Named_profile& npx, const Named_profile& npy, const Fold_envelope& xenv, const Fold_envelope& yenv)
  : Pair_envelope(npx.size(),npy.size(),0),
    npx(&npx), npy(&npy), xenv(&xenv), yenv(&yenv),
    paths_added (0),
    cells_needed (0),
    cells_available (0)
{ }

void Sampled_pair_envelope::initialise_sampled (const Pair_forward_DP_matrix& fwd, int n_samples, const set<int>& states)
{
  CTAG(6,POSTENV POSTENV_ALIGN) << "Sampling " << n_samples << " pairwise alignments for alignment envelope\n";

  // resize
  initialise_full (fwd.xsize-1, fwd.ysize-1, 0);

  // do the tracebacks, make the envelope
  for (int n = 0; n < n_samples; ++n)
    {
      // get the alignment, incorporate it
      const vector<int> state_path = fwd.sample_state_path();
      if (!add_state_path (state_path, fwd.hmm, states))
	{
	  CTAG(4,POSTENV POSTENV_ALIGN) << "Ran out of memory while adding alignment " << n+1 << " to envelope\n";
	  break;
	}
      // log it
      if (CTAGGING(4,POSTENV POSTENV_ALIGN))
	{
	  CL << "Alignment " << n+1 << ":\n";
	  const Pairwise_path path = fwd.hmm.convert_state_path_to_alignment (state_path);
	  if (npx != 0 && npy != 0)
	    {
	      const Alignment align (path, *npx, *npy);
	      align.write_MUL (CL, *fwd.alphabet);
	    }
	  else
	    {
	      vector<sstring> row_name (2);
	      row_name[0] = "x";
	      row_name[0] = "y";
	      path.show (CL, row_name);
	    }
	}
    }
}

void Sampled_pair_envelope::initialise_best (const Pair_CYK_KYC_matrix& cyk_kyc, const Pair_HMM_scores& hmm, int n_best, const set<int>& states)
{
  CTAG(6,POSTENV POSTENV_ALIGN) << "Finding " << n_best << " best pairwise alignments for alignment envelope\n";

  // get references to CYK & KYC matrices
  const Pair_CYK_matrix& cyk = cyk_kyc.cyk;
  const Pair_KYC_matrix& kyc = cyk_kyc.kyc;

  // resize
  initialise_full (cyk.xlen, cyk.ylen, 0);

  // make the cell sorter
  Pair_CYK_KYC_cell_sorter cell_sorter (cyk_kyc);

  // get ready to collect subseq stats
  Pair_CFG_parse_stats stats (cyk);

  // do the tracebacks
  int n_tracebacks = 0;
  Score last_score = InfinityScore;
  for_const_contents (vector<Pair_CFG_DP_matrix::Cell_coords>, cell_sorter.sorted_cells, cell)
    {
      // describe cell
      if (CTAGGING(-1,POSTENV_CELLS))
	CL << "(considering state " << cell->state << ", xsubseq " << cell->xsubseq << ", ysubseq " << cell->ysubseq << ")\n";

      // skip this cell if we've already seen it
      const Pair_CFG_parse_stats::Subseq_index_pair sip (cell->xsubseq, cell->ysubseq);
      if (stats.pair_count.find (sip) != stats.pair_count.end())
	continue;

      // bail out if score drops too low & we've done enough tracebacks
      const Score sc = cyk_kyc.max_state_sc (cell->state, cell->xsubseq, cell->ysubseq);
      if (n_tracebacks >= n_best && sc < last_score)
	break;

      // get the traceback
      vector<int> path;
      kyc.traceback_from (cell->xsubseq, cell->ysubseq, cell->state, path);

      // incorporate this path
      if (!add_state_path (path, hmm, states))
	{
	  CTAG(4,POSTENV POSTENV_ALIGN) << "Ran out of memory while adding alignment " << n_tracebacks+1 << " to envelope\n";
	  break;
	}

      // record the subseq coords
      const Pair_CFG_local_path local_path (0, cyk.xlen, 0, cyk.ylen, path);
      const Pair_CFG_parse_tree parse_tree = cyk.cfg.parse (local_path);
      stats.count_subseqs (parse_tree);

      // log it
      if (CTAGGING(4,POSTENV POSTENV_ALIGN))
	{
	  const Score test_sc = cyk.cfg.path_score (parse_tree, cyk.npx.dsq, cyk.npy.dsq);
	  CL << "Alignment " << n_tracebacks+1 << " (score " << Score2Bits(sc) << " bits):\n";
	  const Pair_CFG_alignment cfg_align = parse_tree.alignment (cyk.cfg.state_type, cyk.npx, cyk.npy);
	  const Pairwise_path& pairwise_path = cfg_align.pairwise_path (TRUE);
	  const Alignment align (pairwise_path, cyk.npx, cyk.npy);
	  align.write_MUL (CL, CFG_alphabet);
	  if (sc != test_sc)
	    {
	      CL << "WARNING: Pair_CFG says score is " << Score2Bits(test_sc) << " bits. Traceback was from state " << cell->state << ", ";
	      CL << "xsubseq " << cell->xsubseq << " (" << cyk.xenv.subseq[cell->xsubseq].start << '+' << cyk.xenv.subseq[cell->xsubseq].len << "), ";
	      CL << "ysubseq " << cell->ysubseq << " (" << cyk.yenv.subseq[cell->ysubseq].start << '+' << cyk.yenv.subseq[cell->ysubseq].len << ").\n";
	      CL << "CYK expanded trace:\n";
	      cyk.show_expanded_trace (parse_tree, CL, FALSE);
	      CL << "KYC expanded trace:\n";
	      kyc.show_expanded_trace (parse_tree, CL, TRUE);
	      CL << "CFG score breakdown:\n";
	      cyk.cfg.show_score_breakdown (parse_tree, cyk.npx.dsq, cyk.npy.dsq, CL);
	    }
	}

      // update the number of tracebacks done & the last score
      ++n_tracebacks;
      last_score = sc;
    }
}

bool Sampled_pair_envelope::add_state_path (const vector<int>& path, const Pair_HMM_scores& hmm, const set<int>& states)
{
  CTAG(4,POSTENV) << "Adding state path to alignment envelope\n";

  unsigned long new_cells = 0;
  int x = 0;
  int y = 0;
  vector<int> new_x, new_y;
  for (int i = 0; i <= (int) path.size(); ++i)
    {
      bool add_cut = false;
      if (!allow_cut(x,y))
	{
	  if (i < (int) path.size())
	    {
	      const int s = path[i];
	      if (states.find (s) != states.end())
		add_cut = true;  // allow this cut if current state is in path
	      else if (i > 0)
		if (states.find (path[i-1]) != states.end())
		  add_cut = true;  // allow cut if previous state was in path
	    }
	  else if (path.size())
	    if (states.find (path.back()) != states.end())
	      add_cut = true;

	  if (add_cut)
	    {
	      allow_cut(x,y) = paths_added + 1;
	      new_x.push_back (x);
	      new_y.push_back (y);
	      if (xenv && yenv)
		new_cells += xenv->by_start[x].size() * yenv->by_start[y].size();
	    }
	}

      if (i == (int) path.size())
	break;

      const int s = path[i];
      if (s >= 0)
	{
	  x += hmm.dx (s);
	  y += hmm.dy (s);
	}
    }

  const unsigned long new_cells_needed = cells_needed + new_cells;
  const bool out_of_memory = cells_available > 0 && new_cells_needed > cells_available;

  CTAG(3,POSTENV) << "Added path #" << paths_added + 1 << " to alignment envelope; " << new_cells_needed << " cells needed, " << cells_available << " available\n";

  if (out_of_memory)
    for (int n = 0; n < (int) new_x.size(); ++n)
      allow_cut (new_x[n], new_y[n]) = 0;
 else
   {
     cells_needed = new_cells_needed;
     ++paths_added;
   }

 return !out_of_memory;
}
