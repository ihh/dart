#include "ecfgenv.h"

ECFG_envelope::ECFG_envelope (int seqlen, int max_subseq_len)
  : use_foldenv(false)
{
  init (seqlen, max_subseq_len);
}

void ECFG_envelope::init_from_fold_string (const sstring& fold_string, int max_subseq_len)
{
  foldenv.initialise_from_fold_string (fold_string, max_subseq_len);
  use_foldenv = true;

  // copy Subseq_coords and seqlen
  subseq.clear();
  subseq.reserve (foldenv.subseq.size());
  for_const_contents (vector<Subseq>, foldenv.subseq, ss)
    subseq.push_back (*ss);

  seqlen = foldenv.seqlen();

  // initialize subseq_lookup
  subseq_lookup = vector<map<int,int> > (seqlen + 1);
  for (int start = 0; start <= seqlen; ++start)
    {
      const vector<int>& by_start = foldenv.by_start[start];
      subseq_lookup[start].resize (by_start.size());
      for_const_contents (vector<int>, by_start, subseq_idx)
	{
	  const int len = foldenv.subseq[*subseq_idx].len;
	  subseq_lookup[start][len] = *subseq_idx;
	}
    }

  // log
  if (CTAGGING(1,FOLDENV))
    foldenv.dump (CL);
}

void ECFG_envelope::init (int S, int L)
{
  CTAG(4,FOLDENV) << "Initialising fold envelope: " << S << " columns, max subseq length = " << L << '\n';
  if (S > 0 && L == 0)
    CL << "NB max subseq length is 0. This is only suitable for left- or right-regular grammars.\n";

  use_foldenv = false;
  foldenv = Fold_envelope();  // empty dummy

  seqlen = S;
  max_subseq_len = L >= 0 ? L : S;
  subseq.clear();
  for (int len = 0; len <= min (max_subseq_len, seqlen); ++len)
    for (int start = 0; start <= seqlen - len; ++start)
      subseq.push_back (Subseq_coords (start, len));
  for (int len = max_subseq_len + 1; len <= seqlen; ++len)
    {
      subseq.push_back (Subseq_coords (0, len));
      if (len < seqlen)
	subseq.push_back (Subseq_coords (seqlen - len, len));
    }

  // test... should be in a separate test program but whatthehell
  if (CTAGGING(-1,FOLDENV_TEST))
    {
      CL << "Running ECFG_envelope tests: seqlen = " << seqlen << ", max_subseq_len = " << max_subseq_len << "\n";
      CL << "Testing ECFG_envelope::find_subseq_idx\n";
      for (int idx = 0; idx < subseqs(); ++idx)
	if (find_subseq_idx (subseq[idx].start, subseq[idx].len) != idx)
	  THROWEXPR ("Subsequence index error: find_subseq_idx(start=" << subseq[idx].start << ",len=" << subseq[idx].len
		     << ") = " << find_subseq_idx (subseq[idx].start, subseq[idx].len) << ", expected " << idx
		     << " (seqlen=" << seqlen << ", max_subseq_len=" << max_subseq_len << ")");

      CL << "Testing ECFG_envelope::get_bif_in\n";
      vector<Subseq::Bifurc_in> bif_in;
      typedef map<int,Subseq::Bifurc_in> Bifurc_in_map;
      Bifurc_in_map test_bif_in;
      for (int idx = 0; idx < subseqs(); ++idx)
	{
	  const Subseq_coords& ss = subseq[idx];
	  get_bif_in (ss, bif_in);
	  int l, r;
	  test_bif_in.clear();
	  for (int mid = ss.start; mid <= ss.end(); ++mid)
	    if ((l = find_subseq_idx (ss.start, mid - ss.start)) >= 0
		&& (r = find_subseq_idx (mid, ss.end() - mid)) >= 0)
	      test_bif_in[mid] = Subseq::Bifurc_in (l, r);

	  bool fail = false;
	  if (test_bif_in.size() != bif_in.size())
	    fail = true;
	  for_const_contents (vector<Subseq::Bifurc_in>, bif_in, bi)
	    {
	      Bifurc_in_map::const_iterator tbi = test_bif_in.find (subseq[bi->r].start);
	      if (tbi == test_bif_in.end())
		fail = true;
	      else if (tbi->second.l != bi->l || tbi->second.r != bi->r)
		fail = true;
	    }

	  if (fail)
	    {
	      CL << "bif_in:";
	      for_const_contents (vector<Subseq::Bifurc_in>, bif_in, bi)
		CL << " (" << subseq[bi->l].start << ',' << subseq[bi->l].end()
		   << ")(" << subseq[bi->r].start << ',' << subseq[bi->r].end()
		   << ")";
	      CL << "\ntest_bif_in:";
	      for_const_contents (Bifurc_in_map, test_bif_in, tbi)
		CL << " (" << subseq[tbi->second.l].start << ',' << subseq[tbi->second.l].end()
		   << ")(" << subseq[tbi->second.r].start << ',' << subseq[tbi->second.r].end()
		   << ")";
	      CL << "\n";
	      THROWEXPR ("bif_in comparison failed");
	    }
	}

      // should test bif_outl and bif_outr too, but that's a lot of mind-numbing code... (yes I know that's a bad excuse)

      CL << "Passed ECFG_envelope tests\n";
    }
}
