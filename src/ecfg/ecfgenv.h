#ifndef ECFG_FOLD_ENVELOPE_INCLUDED
#define ECFG_FOLD_ENVELOPE_INCLUDED

#include "scfg/foldenv.h"

// Fold envelopes for ECFG's are less flexible than for Pair_CFG's
// For parameters (seqlen=S, max_subseq_len=L) the following (start,len) Subseq_coords are allowed:
//   (s,l)   where 0 <= s <= S-l and 0 <= l <= L
//   (s,S-s) where 0 <= s < S-L
struct ECFG_envelope
{
  // data
  int seqlen, max_subseq_len;
  vector<Subseq_coords> subseq;

  // constructor
  ECFG_envelope (int seqlen = 0, int max_subseq_len = 0);

  // initialiser
  void init (int seqlen, int max_subseq_len);

  // accessors
  inline int subseqs() const { return subseq.size(); }

  // method to get index of a subseq; returns -1 if outside envelope
  inline int find_subseq_idx (int start, int len) const
  {
    int idx = -1;
    if (start >= 0 && len >= 0 && start + len <= seqlen)
      {
	if (len <= max_subseq_len)
	  idx = len * (2*seqlen + 3 - len) / 2 + start;   // = start + sum_{k=0}^{len-1} (seqlen + 1 - k)
	else if (start == 0)
	  idx = (max_subseq_len + 1) * (2*seqlen - max_subseq_len - 2) / 2 + len * 2;  // = 2*(len - max_subseq_len - 1) + sum_{k=0}^{max_subseq_len} (seqlen + 1 - k)
	else if (start + len == seqlen)
	  idx = (max_subseq_len + 1) * (2*seqlen - max_subseq_len - 2) / 2 + len * 2 + 1;  // = 1 + 2*(len - max_subseq_len - 1) + sum_{k=0}^{max_subseq_len} (seqlen + 1 - k)
      }
#ifdef DART_DEBUG
    if (idx >= 0)
      if (subseq[idx].start != start || subseq[idx].len != len)
	THROWEXPR ("Incorrectly calculated subsequence (" << start << '+' << len << ')');
#endif
    return idx;
  }

  // methods to get indices of bifurcation connections
  inline void get_bif_in (const Subseq_coords& out, vector<Subseq::Bifurc_in>& bif_in) const
  {
    bif_in.clear();
    const int start = out.start;
    const int end = out.end();
    const int min_r_start = end == seqlen ? start : (end - max_subseq_len);
    const int max_l_end = start == 0 ? end : (start + max_subseq_len);
    for (int mid = max_l_end; mid >= min_r_start; --mid)
      {
	const int l = find_subseq_idx (start, mid - start);
	if (l >= 0)
	  {
	    const int r = find_subseq_idx (mid, end - mid);
	    if (r >= 0)
	      bif_in.push_back (Subseq::Bifurc_in (l, r));
	  }
      }
  }

  inline void get_bif_outl (const Subseq_coords& r, vector<Subseq::Bifurc_out_l>& bif_outl) const
  {
    bif_outl.clear();
    const int mid = r.start;
    const int end = r.end();
    // TODO: calculate min_out_start and min_l_start here
    for (int start = 0; start <= mid; ++start)
      {
	const int out = find_subseq_idx (start, end - start);
	if (out >= 0)
	  {
	    const int l = find_subseq_idx (start, mid - start);
	    if (l >= 0)
	      bif_outl.push_back (Subseq::Bifurc_out_l (out, l));
	  }
      }
  }

  inline void get_bif_outr (const Subseq_coords& l, vector<Subseq::Bifurc_out_r>& bif_outr) const
  {
    bif_outr.clear();
    const int start = l.start;
    const int mid = l.end();
    // TODO: calculate max_out_end and max_r_end here
    for (int end = seqlen; end >= mid; --end)
      {
	const int out = find_subseq_idx (start, end - start);
	if (out >= 0)
	  {
	    const int r = find_subseq_idx (mid, end - mid);
	    if (r >= 0)
	      bif_outr.push_back (Subseq::Bifurc_out_r (out, r));
	  }
      }
  }

  // methods to count bifurcations
  inline double n_bif (int subseq_idx) const
  {
    const int len = subseq[subseq_idx].len;
    return len > max_subseq_len ? (max_subseq_len + 2) : len + 1;
  }
  inline double n_bif() const
  {
    double b = 0.;
    for (int i = 0; i < subseqs(); ++i)
      b += n_bif (i);
    return b;
  }
};

#endif /* ECFG_FOLD_ENVELOPE_INCLUDED */
