#include "scfg/cfgdotplot.h"

#define CFGDOTPLOT_REPORT_INTERVAL 10000  /* number of cells between logfile messages during DP */

PairCFG_fold_dotplot::PairCFG_fold_dotplot (const Pair_inside_outside_matrix& in_out, int seq_index)
  : Dotplot (seq_index == 0 ? in_out.inside.npx.seq : in_out.inside.npy.seq,
		     seq_index == 0 ? in_out.inside.npx.seq : in_out.inside.npy.seq)
{
  const Fold_envelope& envx (in_out.inside.xenv);
  const Fold_envelope& envy (in_out.inside.yenv);
  const Pair_CFG_scores& cfg (in_out.inside.cfg);

  const Fold_envelope& env (seq_index == 0 ? envx : envy);
  const Fold_envelope& other_env (seq_index == 0 ? envy : envx);
  const State_type lr_emit_flags = (seq_index == 0 ? EmitXLR : EmitYLR);

  const int seqlen = xseq.size();
  const bool logging = CTAGGING(1,DOTPLOT_PROBS);

  const int total_cells = seqlen * (seqlen + 1) / 2;
  int cells_done = 0, last_report = 0;

  for (int start = 0; start < seqlen; ++start)
    for (int end = start + 1; end <= seqlen; ++end) {
      const int subseq = env.find_subseq_idx (start, end - start);
      if (subseq >= 0)
	{
	  Prob& p = (*this) (end - 1, start);
	  for (int s = 0; s < cfg.states(); ++s)
	    if ((cfg.state_type[s] & lr_emit_flags) == lr_emit_flags)
	      for (int other_subseq = 0; other_subseq < other_env.subseqs(); ++other_subseq) {
		int xsubseq, ysubseq;
		if (seq_index == 0) {
		  xsubseq = subseq;
		  ysubseq = other_subseq;
		} else {
		  xsubseq = other_subseq;
		  ysubseq = subseq;
		}
		const Prob q = Score2Prob (in_out.post_state_sc (s, xsubseq, ysubseq));
		p += q;
		if (logging)
		  CTAG(1,DOTPLOT_PROBS) << "P(Emit" << (seq_index==0 ? "X" : "Y") << "LR,state=" << s << ",xsubseq=[" << envx.subseq[xsubseq].start << ".." << envx.subseq[xsubseq].end() << "],ysubseq=[" << envy.subseq[ysubseq].start << ".." << envy.subseq[ysubseq].end() << "]) = " << q << "\n";
	      }

	  // print this entry to log
	  if (logging)
	    CTAG(2,DOTPLOT_PROBS) << "P(EmitLR,start=" << start << ",end=" << end << ")=" << p << "\n";

	  ++cells_done;
	  if (cells_done - last_report >= CFGDOTPLOT_REPORT_INTERVAL
	      || (cells_done == total_cells && cells_done > last_report))
	    {
	      last_report = cells_done;
	      CTAG(1,CFGDP) << "PairCFG_fold_dotplot: finished " << cells_done << " basepairs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%)\n";
	    }
	}
    }
}

PairCFG_alignment_dotplot::PairCFG_alignment_dotplot (const Pair_inside_outside_matrix& in_out)
  : Dotplot (in_out.inside.npx.seq, in_out.inside.npy.seq)
{
  const Pair_CFG_scores& cfg (in_out.inside.cfg);
  const Fold_envelope& envx (in_out.inside.xenv);
  const Fold_envelope& envy (in_out.inside.yenv);
  const bool logging = CTAGGING(1,DOTPLOT_PROBS);

  const int total_cells = xseq.size() * yseq.size();
  int cells_done = 0, last_report = 0;

  for (int y = 1; y <= (int) yseq.size(); ++y)
    for (int x = 1; x <= (int) xseq.size(); ++x) {
      Prob& p = (*this) (x-1, y-1);

      // subseqs with length>2 ending at (x,y)
      for (int y0 = 0; y0 < y-1; ++y0) {
	const int ysubseq = envy.find_subseq_idx (y0, y - y0);
	if (ysubseq >= 0)
	  for (int x0 = 0; x0 < x-1; ++x0) {
	    const int xsubseq = envx.find_subseq_idx (x0, x - x0);
	    if (xsubseq >= 0)
	      for (int s = 0; s < cfg.states(); ++s) {
		if ((cfg.state_type[s] & EmitXRYR) == EmitXRYR) {
		  const Prob q = Score2Prob (in_out.post_state_sc (s, xsubseq, ysubseq));
		  p += q;
		  if (logging)
		    CTAG(1,DOTPLOT_PROBS) << "P(EmitXRYR,state=" << s << ",xsubseq=[" << envx.subseq[xsubseq].start << ".." << envx.subseq[xsubseq].end() << "],ysubseq=[" << envy.subseq[ysubseq].start << ".." << envy.subseq[ysubseq].end() << "]) = " << q << "\n";
		}
	      }
	  }
      }

      // subseqs starting at (x-1,y-1)
      for_const_contents (vector<int>, envy.by_start[y-1], ysubseq)
	for_const_contents (vector<int>, envx.by_start[x-1], xsubseq)
	for (int s = 0; s < cfg.states(); ++s) {
	  if ((cfg.state_type[s] & EmitXLYL) == EmitXLYL) {
	    const Prob q = Score2Prob (in_out.post_state_sc (s, *xsubseq, *ysubseq));
	    p += q;
	    if (logging)
	      CTAG(1,DOTPLOT_PROBS) << "P(EmitXLYL,state=" << s << ",xsubseq=[" << envx.subseq[*xsubseq].start << ".." << envx.subseq[*xsubseq].end() << "],ysubseq=[" << envy.subseq[*ysubseq].start << ".." << envy.subseq[*ysubseq].end() << "]) = " << q << "\n";
	  }
	}

      // print this entry to log
      if (logging)
	CTAG(2,DOTPLOT_PROBS) << "P(EmitXY,x=" << x << ",y=" << y << ")=" << p << "\n";

      ++cells_done;
      if (cells_done - last_report >= CFGDOTPLOT_REPORT_INTERVAL
	  || (cells_done == total_cells && cells_done > last_report))
	{
	  last_report = cells_done;
	  CTAG(1,CFGDP) << "PairCFG_alignment_dotplot: finished " << cells_done << " residue pairs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%)\n";
	}
    }

}
