#include <iomanip>
#include <stack>
#include "scfg/foldenv.h"
#include "util/logfile.h"

// periodicity (in columns) of log message when constructing a fold envelope
#define INIT_FOLD_ENV_LOG_COLUMNS 100

Subseq_coords::Subseq_coords (int start, int len) :
  start (start),
  len (len)
{ }

Subseq::Subseq (int start, int len) :
  Subseq_coords (start, len),
  in_flags (CFLAG_NONE),
  out_flags (CFLAG_NONE),
  start_flag (FALSE),
  bif_in(),
  bif_out_l(),
  bif_out_r(),
  by_start_index(-1),
  by_len_index(-1)
{
  for (int i = 0; i < N_CONN; ++i)
    in[i] = out[i] = -1;
}

bool operator== (const Subseq::Bifurc_in& b1, const Subseq::Bifurc_in& b2)
{
  return b1.l == b2.l && b1.r == b2.r;
}

bool operator== (const Subseq::Bifurc_out_l& b1, const Subseq::Bifurc_out_l& b2)
{
  return b1.out == b2.out && b1.l == b2.l;
}

bool operator== (const Subseq::Bifurc_out_r& b1, const Subseq::Bifurc_out_r& b2)
{
  return b1.out == b2.out && b1.r == b2.r;
}

bool Subseq::operator== (const Subseq& subseq) const
{
  return start == subseq.start
    && len == subseq.len
    && in_flags == subseq.in_flags
    && out_flags == subseq.out_flags
    && start_flag == subseq.start_flag
    && by_start_index == subseq.by_start_index
    && by_len_index == subseq.by_len_index;
}

bool Subseq::test_coords() const
{
  return start >= 0 && len >= 0;
}

bool Subseq::test_flags() const
{
  const int in_test = (in_l() < 0 ? 0 : CFLAG_L) | (in_r() < 0 ? 0 : CFLAG_R) | (in_lr() < 0 ? 0 : CFLAG_LR);
  const int out_test = (out_l() < 0 ? 0 : CFLAG_L) | (out_r() < 0 ? 0 : CFLAG_R) | (out_lr() < 0 ? 0 : CFLAG_LR);
  return in_test == (in_flags & CFLAG_L_R_LR) && out_test == (out_flags & CFLAG_L_R_LR);
}

bool Subseq::test_in_l (const Subseq& l_subseq) const
{
  return start == l_subseq.start - 1 && end() == l_subseq.end();
}

bool Subseq::test_in_r (const Subseq& r_subseq) const
{
  return start == r_subseq.start && end() == r_subseq.end() + 1;
}

bool Subseq::test_in_lr (const Subseq& lr_subseq) const
{
  return start == lr_subseq.start - 1 && end() == lr_subseq.end() + 1;
}

bool Subseq::test_in_b (const Subseq& l_subseq, const Subseq& r_subseq) const
{
  return start == l_subseq.start && end() == r_subseq.end() && l_subseq.end() == r_subseq.start;
}

bool Subseq::test_fold_string (const sstring& fold_string) const
{
  if (start > (int) fold_string.size() || end() > (int) fold_string.size()) return 0;
  if (len > 0)
    {
      const char lchar = fold_string[start];
      const char rchar = fold_string[end() - 1];
      if (in_l() >= 0 && (is_lchar(lchar) || is_rchar(lchar))) return 0;
      if (in_r() >= 0 && (is_lchar(rchar) || is_rchar(rchar))) return 0;
      if (in_lr() >= 0 && (!is_lchar(lchar) || !is_rchar(rchar))) return 0;
      int nest_lev = 0;
      for (int i = start; i < end(); ++i)
	{
	  const char c = fold_string[i];
	  if (is_lchar(c))
	    ++nest_lev;
	  else if (is_rchar(c))
	    if (--nest_lev < 0)
	      THROWEXPR ("Too many " << (char) FOLD_RCHAR << "'s in fold string '" << fold_string << "' (subseq start " << start << ", length " << len << ")");
	}
      if (nest_lev > 0) THROWEXPR ("Too many " << (char) FOLD_LCHAR << "'s in fold string '" << fold_string << "' (subseq start " << start << ", length " << len << ")");
    }
  return 1;
}

sstring Subseq::terse_desc (int index, const Biosequence& seq, bool show_out) const
{
  sstring td;

  td.width (4);
  left_align (td);
  td << index << ' ';

  td.width (4);
  right_align (td);
  td << start << '+';

  td.width (4);
  left_align (td);
  td << len;
  td << (start_flag ? '(' : ' ');

  if (len < 4)
    {
      td.width (4);
      right_align (td);
      td << seq.substr (start, len);
    }
  else
    td << seq[start] << ".." << seq[start + len - 1];
  td << (start_flag ? ')' : ' ');
  
  td.width (4);
  left_align (td);
  if (show_out)
    if (out_flags & CFLAG_L) td << out_l(); else td << "    ";
  else
    if (in_flags & CFLAG_L) td << in_l(); else td << "    ";
  td << ' ';

  td.width (4);
  if (show_out)
    if (out_flags & CFLAG_R) td << out_r(); else td << "    ";
  else
    if (in_flags & CFLAG_R) td << in_r(); else td << "    ";
  td << ' ';

  td.width (4);
  if (show_out)
    if (out_flags & CFLAG_LR) td << out_lr(); else td << "    ";
  else
    if (in_flags & CFLAG_LR) td << in_lr(); else td << "    ";
  td << ' ';

  return td;
}

sstring Subseq::terser_desc() const
{
  sstring td;
  td << '(' << start << '+' << len << ')';
  return td;
}

Local_fold_string::Local_fold_string()
{ }

Local_fold_string::Local_fold_string (const sstring& fold, int start, int seqlen)
  : fold (fold), start (start), seqlen (seqlen)
{ }

bool Local_fold_string::is_global() const
{
  return start == 0 && seqlen == (int) fold.size();
}

sstring Local_fold_string::global_fold() const
{
  sstring gf;
  gf << sstring (start, '.') << fold << sstring (seqlen - end(), '.');
  return gf;
}

void Local_fold_string::add (const Local_fold_string& local_fold_string)
{
  // check overlap
  if (local_fold_string.seqlen != seqlen)
    THROWEXPR ("Local_fold_strings have different lengths");
  if (local_fold_string.start > end() || local_fold_string.end() < start)
    THROWEXPR ("Local_fold_strings do not overlap");
  // get new coords & fold string
  const int new_start = min (local_fold_string.start, start);
  const int new_end = max (local_fold_string.end(), end());
  const int new_len = new_end - new_start;
  sstring new_fold (new_len);
  // copy this fold string
  for (int j = 0; j < (int) fold.size(); ++j)
    {
      // get character to write
      const char c = fold[j];
      const int pos = j + start;
      // write new fold string
      const int k = pos - new_start;
      new_fold[k] = c;
    }
  // copy other fold string & check overlaps
  for (int i = 0; i < (int) local_fold_string.fold.size(); ++i)
    {
      // get character to write
      const char c = local_fold_string.fold[i];
      const int pos = i + local_fold_string.start;
      // check if fold strings overlap; if so, check they agree, otherwise write character to fold string
      const int j = pos - start;
      if (j >= 0 && j < (int) fold.size())
	{
	  if (fold[j] != c)
	    THROWEXPR ("Local_fold_strings do not agree: character '" << c << "' at index " << i << " of '[" << local_fold_string.start << ']' << local_fold_string.fold << "' doesn't match '" << fold[j] << "' at index " << j << " of '[" << start << ']' << fold << "'");
	}
      else
	{
	  // write new fold string
	  const int k = pos - new_start;
	  new_fold[k] = c;
	}
    }
  // update
  ((char_string::basic_string&)fold).swap (new_fold);
  start = new_start;
}

void Subseq::Fold_envelope_base::init_pseudovecs()
{
#ifdef DART_USE_BIFURCATION_PSEUDOVECTORS
  for (int subseq_idx = 0; subseq_idx < (int) subseq.size(); ++subseq_idx)
    {
      subseq[subseq_idx].bif_in.init (*this, subseq_idx);
      subseq[subseq_idx].bif_out_l.init (*this, subseq_idx);
      subseq[subseq_idx].bif_out_r.init (*this, subseq_idx);
    }
#endif /* DART_USE_BIFURCATION_PSEUDOVECTORS */
}

Fold_envelope::Fold_envelope (int seqlen)
{
  initialise_full (seqlen);
}

Fold_envelope::Fold_envelope (const Fold_envelope& env)
{
  *this = env;
}

int Fold_envelope::bifurcations() const
{
  int bifurcs = 0;
  for (int i = 0; i < (int) subseq.size(); ++i) bifurcs += subseq[i].bif_in.size();
  return bifurcs;
}

void Fold_envelope::initialise_from_fold_string (const sstring& fold_string, int max_subseq_len)
{
  CTAG(6,FOLDENV) << "Initialising fold envelope from fold string '" << fold_string << "'\n";

  const int fold_size = fold_string.size();
  resize (fold_size);

  // initialize startpos_stack with an (empty) vector<int>
  stack<vector<int> > startpos_stack;
  startpos_stack.push (vector<int>());

  for (int end = 0; end <= fold_size; ++end)
    {
      // log
      if (end > 0 && end % INIT_FOLD_ENV_LOG_COLUMNS == 0)
	CTAG(6,FOLDENV) << "Still initialising fold envelope; finished " << end << " of " << fold_size << " columns (subseqs so far: " << subseq.size() << ")\n";

      // create empty Subseq
      const int empty_subseq_idx = new_subseq (end, 0);
      connect_b (empty_subseq_idx, empty_subseq_idx, empty_subseq_idx);
      // parse next char
      if (end > 0)
	{
	  const char fold_char = fold_string [end - 1];
	  if (is_lchar (fold_char))
	    {
	      startpos_stack.top().push_back (end - 1);
	      startpos_stack.push (vector<int>());
	    }
	  else
	    {
	      if (is_rchar (fold_char))
		{
		  startpos_stack.pop();
		  if (startpos_stack.empty())
		    THROWEXPR ("Too many <'s\nFold string: " << fold_string);
		  const int start = startpos_stack.top().back();
		  const int new_subseq_len = end - start;
		  if (max_subseq_len >= 0 && new_subseq_len > max_subseq_len && !(start == 0 || end == fold_size))
		    THROWEXPR ("Fold envelope contains a basepair between columns " << start+1 << " and " << end+1
			       << ", which implies a subseq of length " << new_subseq_len << ", but maximum allowed subseq length is " << max_subseq_len);
		  // create the new subseq
		  const int new_subseq_idx = new_subseq (start, new_subseq_len);
		  // note that connect_lr () automatically sets subsequence connectivity:
		  // here we implicitly set
		  //  (*this).subseq[new_subseq_idx].inflags = (*this).subseq[lr_subseq_idx].outflags = Subseq::CFLAG_LR
		  const int lr_subseq_idx = find_lr (new_subseq_idx);
		  if (lr_subseq_idx >= 0)
		    connect_lr (new_subseq_idx, lr_subseq_idx);
		  // dummy bifurcations involving empty subseqs at start & end
		  connect_b (new_subseq_idx, by_len[0][start], new_subseq_idx);
		  connect_b (new_subseq_idx, new_subseq_idx, by_len[0][end]);
		}
	      else
		startpos_stack.top().push_back (end - 1);

	      // create all subseqs (start,end) for each start in startpos
	      vector<int>& startpos = startpos_stack.top();
	      const int max_start_idx = startpos.size() - (is_rchar(fold_char) ? 2 : 1);  // if fold_char is '>', we have already created the smallest subseq, so don't create it again
	      for (int start_idx = max_start_idx; start_idx >= 0; --start_idx)
		{
		  const int start = startpos[start_idx];
		  const int len = end - start;
		  const vector<int>& start_subseqs = by_start[start];
		  if (max_subseq_len >= 0 && len > max_subseq_len && !(start == 0 || end == fold_size))
		    continue;
		  // create the new subseq
		  const int new_subseq_idx = new_subseq (start, len); // NB after this call, start_subseqs includes new subseq
		  // connect the new subseq, implicitly setting the Subseq::CFLAG_L and Subseq::CFLAG_R bits of in_flags/out_flags in the Subseq objects
		  if (!is_lchar (fold_string[start]))
		    {
		      const int l_subseq_idx = find_l (new_subseq_idx);
		      if (l_subseq_idx >= 0)
			connect_l (new_subseq_idx, l_subseq_idx);
		    }
		  if (!is_rchar (fold_char))
		    {
		      const int r_subseq_idx = find_r (new_subseq_idx);
		      if (r_subseq_idx >= 0)
			connect_r (new_subseq_idx, r_subseq_idx);
		    }
		  // quick sanity test...
		  if ((max_subseq_len < 0 || max_subseq_len >= fold_size) && (int) start_subseqs.size() != new_subseq_idx + 1 - empty_subseq_idx)
		    THROWEXPR ("Should be one end_subseq for every start_subseq");
		  // add every bifurcation point not flanked by CONTIGUOUS_UNPAIRED_CHAR characters
		  for (int bif_idx = 0; bif_idx < (int) start_subseqs.size(); ++bif_idx)
		    {
		      const int bif_l_idx = start_subseqs[bif_idx];
		      const int bif_r_start = subseq[bif_l_idx].end();
		      if ((start == 0 || end == fold_size) ? true : !is_contiguous_digram (fold_string[bif_r_start-1], fold_string[bif_r_start]))
			{
			  const int bif_r_len = end - bif_r_start;
			  const int bif_r_idx = find_subseq_idx (bif_r_start, bif_r_len);
			  if (bif_r_idx >= 0)
			    connect_b (new_subseq_idx, bif_l_idx, bif_r_idx);
			}
		    }
		}
	    }
	}
    }

  // expect to have exactly one vector<int> left on startpos_stack
  if (startpos_stack.empty())
    THROWEXPR ("Too many <'s\nFold string: " << fold_string);

  startpos_stack.pop();

  if (!startpos_stack.empty())
    THROWEXPR ("Too many >'s\nFold string: " << fold_string);

  make_global();
  init_pseudovecs();
}

void Fold_envelope::initialise_local_fold_string_neighbourhood (const Local_fold_string& local_fold_string, int max_subseq_len, int min_loop_len, int min_overlap_len, bool band_3prime)
{
  // FIXME: this function can potentially create a lot of inaccessible Subseqs.
  // While this does not change the result, it is inefficient, in terms of both time & memory.
  // It would be better to prune out inaccessible subseqs at the end of the function.

  // get stuff from Local_fold_string
  const sstring& fold_string = local_fold_string.fold;
  const int fold_string_start = local_fold_string.start;
  const int new_seqlen = local_fold_string.seqlen;

  // resize
  resize (new_seqlen);

  // get coords of fold_string
  const int fold_string_len = fold_string.size();
  const int fold_string_end = fold_string_start + fold_string_len;

  // -1 sets max_seqlen equal to sequence length
  if (max_subseq_len < 0 || max_subseq_len > new_seqlen)
    max_subseq_len = new_seqlen;
  else if (max_subseq_len < fold_string_len)
    max_subseq_len = fold_string_len;

  // check bounds
  THROWASSERT (fold_string_len > 0);
  THROWASSERT (fold_string_end <= seqlen());
  THROWASSERT (min_overlap_len <= fold_string_len);

  // print log message
  CTAG(6,FOLDENV) << "Seeding fold envelope on subseq " << fold_string_start << ".." << fold_string_end << ", fold string '" << fold_string << "', min overlap " << min_overlap_len << ", min loop len " << min_loop_len << ", max subseq len " << max_subseq_len << "\n";

  // create temporary global Fold_envelope for fold string
  Fold_envelope global_env;
  global_env.initialise_from_fold_string (fold_string);

  // Copy global envelope, adding the following sets of extension Subseqs:
  //   (i) all Subseq's overlapping fold_string by at least min_overlap_len, up to a maximum length of max_subseq_len;
  //  (ii) all Subseq's not overlapping fold_string, that are contained inside a member of (i);
  // (iii) all Subseq's i..L, where L = sequence length.
  // Set (i) has L, R and bifurcation connections.
  // Sets (ii) & (iii) are fully connected, within limits of min_loop_len.
  // start_flag is set for all Subseqs (including those copied from global_env) subject to min_overlap_len.
  //
  // Subseqs (S,E) sorted by increasing end order then decreasing start order, bounded by A <= S <= E <= B.
  // Note A <= X <= U <= V <= Y <= B.
  const int X = fold_string_start;
  const int Y = fold_string_end;
  const int U = X + min_overlap_len;  // endpoint of leftmost Subseq with start_flag==TRUE
  const int V = Y - min_overlap_len;  // startpoint of rightmost Subseq with start_flag==TRUE
  const int A = max (U - max_subseq_len, 0);         // leftmost startpoint in envelope
  const int B = min (V + max_subseq_len, seqlen());  // rightmost endpoint in envelope
  // Defined inside E loop:
  //        C = max (E - max_subseq_len, A)   = leftmost startpoint, for given E
  //        D = min (E, V) - min_overlap_len  = rightmost startpoint with start_flag==TRUE, for given E >= U

  // print these values to log
  CTAG(1,FOLDENV) << "In initialise_fold_string_neighbourhood(): A=" << A << ", B=" << B << ", X=" << X << ", Y=" << Y << ", U=" << U << ", V=" << V << "\n";
  
  // Creation order for Subseqs:
  //  (1)     For A <= E < X:                                               [E counting up]
  //  (1.1)     For C <= S <= E:                                            [S counting down]
  //  (1.1.1)     Create fully-connected (S,E), set (ii).                     (left of global_env)
  for (int E = A; E < X; ++E)
    {
      const int C = max (E - max_subseq_len, A);   // leftmost startpoint, for given E
      for (int S = E; S >= C; --S)
	connect_subseq (new_subseq (S, E-S), cflag_mask (E-S, min_loop_len));
    }

  //  (2)     For X <= E <= Y:                                              [E counting up]
  //  (2.1)     For X <= S <= E:                                            [S counting down]
  //  (2.1.1)     Copy (S,E) connections from global_env, if it exists.       (enclosed by global_env)
  //  (2.2)     For C <= S < X:                                             [S counting down]
  //  (2.2.1)     Create L- and bif-connected (S,E), set (i).                 (left extension of global_env)
  for (int E = X; E <= Y; ++E)
    {
      const int C = max (E - max_subseq_len, A);   // leftmost startpoint, for given E
      for (int S = E; S >= X; --S)
	{
	  // check if global_env contains this Subseq
	  const int idx = global_env.find_subseq_idx (S-X, E-S);
	  if (idx >= 0)
	    {
	      // get global Subseq
	      const Subseq& global_ss = global_env.subseq[idx];
	      // create local Subseq
	      const int new_idx = new_subseq (S, E-S);
	      // copy bifurcations exactly
	      for_const_contents (Subseq::Bifurc_in_pseudovec, global_ss.bif_in, ssbif)
		{
		  const Subseq& global_l = global_env.subseq[ssbif->l];
		  const Subseq& global_r = global_env.subseq[ssbif->r];
		  const int l_idx = find_subseq_idx (global_l.start + X, global_l.len);
		  const int r_idx = find_subseq_idx (global_r.start + X, global_r.len);
		  connect_b (new_idx, l_idx, r_idx);
		}
	      // copy emit connections manually
	      const int cflags = global_ss.in_flags;
	      if (cflags & CFLAG_L) connect_l (new_idx, find_l (new_idx));
	      if (cflags & CFLAG_R) connect_r (new_idx, find_r (new_idx));
	      if (cflags & CFLAG_LR) connect_lr (new_idx, find_lr (new_idx));
	    }
	}
      for (int S = X-1; S >= C; --S)
	connect_subseq (new_subseq (S, E-S), CFLAG_L);
    }
  
  //  (3)     For Y < E <= B:                                               [E counting up]
  //  (3.1)     For Y <= S <= E:                                            [S counting down]
  //  (3.1.1)     Create fully-connected (S,E), set (ii).                     (right of global_env)
  //  (3.2)     For X <= S < Y:                                             [S counting down]
  //  (3.2.1)     Create R- and bif-connected (S,E), set (i).                 (right extension of global_env)
  //  (3.3)     For C <= S < X:                                             [S counting down]
  //  (3.3.1)     Create fully-connected (S,E), set (i).                      (encloses global_env)
  for (int E = Y+1; E <= B; ++E)
    {
      const int C = max (E - max_subseq_len, A);   // leftmost startpoint, for given E
      for (int S = E; S >= Y; --S)
	connect_subseq (new_subseq (S, E-S), cflag_mask (E-S, min_loop_len));
      for (int S = Y-1; S >= X; --S)
	connect_subseq (new_subseq (S, E-S), CFLAG_R);
      for (int S = X-1; S >= C; --S)
	connect_subseq (new_subseq (S, E-S), cflag_mask (E-S, min_loop_len));
    }
  
  //  (4)     For U <= E <= B:                                              [E counting up]
  //  (4.1)     For C <= S <= D:                                            [S counting down]
  //  (4.1.1)     Set start_flag for (S,E).
  for (int E = U; E <= B; ++E)
    {
      const int C = max (E - max_subseq_len, A);   // leftmost startpoint, for given E
      const int D = min (E, Y) - min_overlap_len;  // rightmost startpoint with start_flag==TRUE, for given E >= U
      for (int S = D; S >= C; --S)
	{
	  const int idx = find_subseq_idx (S, E-S);
	  if (idx >= 0)
	    subseq[idx].start_flag = TRUE;
      }
    }

  // Create 3' banding sequences
  if (band_3prime)
    {
      const int E = seqlen();
      const int max_S = (B < E || subseq.size() == 0) ? E : subseq.back().start - 1;
      for (int S = max_S; S >= 0; --S)
	{
	  const int idx = new_subseq (S, E-S);
	  connect_subseq (idx, cflag_mask (E-S, min_loop_len));
	  // set start_flag if overlap is large enough.
	  // this is a bit futile, because if band_3prime is set we're almost certainly going to call make_global(),
	  // but at least it's consistent.
	  const int overlap = fold_string_end - max (fold_string_start, S);  // can be negative
	  if (overlap >= min_overlap_len)
	    subseq[idx].start_flag = TRUE;
	}
    }
  init_pseudovecs();

  if (CTAGGING(1,FOLDENV_DEBUG))
    {
      CL << "Exiting initialise_fold_string_neighbourhood() with the following fold envelope:\n";
      dump (CL);
    }
}

void Fold_envelope::initialise_local_fold_string_constrained (const Local_fold_string& local_fold_string, int max_subseq_len, int min_loop_len)
{
  initialise_local_fold_string_neighbourhood (local_fold_string, max_subseq_len, min_loop_len, 0, TRUE);
  make_global();
  init_pseudovecs();
}

void Fold_envelope::initialise_full (int seqlen, int min_loop_len)
{
  // print log message, estimating size & time
  const int L = seqlen;
  if (seqlen > 0)
    {
      const double subs = L*L/2;
      const double bifs = L*L*L/6;
      CTAG(6,FOLDENV) << "Initialising full fold envelope for " << seqlen << "-base sequence: approx " << setprecision(2) << subs << " subsequences & " << setprecision(2) << bifs << " bifurcations\n";
    }

  // do it
  resize (seqlen);
  for (int end = 0; end <= seqlen; ++end)
    {
      for (int start = end; start >= 0; --start)
	{
	  const int len = end - start;
	  const int new_subseq_idx = new_subseq (start, len);
	  if (new_subseq_idx != full_subseq_idx (start, len))
	    THROWEXPR ("Full subseq index miscalculation");

	  if (len > 0)
	    {
	      connect_l (new_subseq_idx, full_subseq_idx (start + 1, len - 1));
	      connect_r (new_subseq_idx, full_subseq_idx (start, len - 1));
	      if (len >= min_loop_len + 2)
		connect_lr (new_subseq_idx, full_subseq_idx (start + 1, len - 2));
	    }
	  
	  for (int mid = start; mid <= end; ++mid)
	    {
	      const int lsub_idx = full_subseq_idx (start, mid - start);
	      const int rsub_idx = full_subseq_idx (mid, end - mid);
	      connect_b (new_subseq_idx, lsub_idx, rsub_idx);
	    }
	}
    }
  if (seqlen >= 0)
    make_global();
  init_pseudovecs();
}

void Fold_envelope::initialise_5prime (int seqlen)
{
  CTAG(6,FOLDENV) << "Initialising 5' fold envelope (reverse HMM) for " << seqlen << "-base sequence\n";

  resize (seqlen);
  for (int end = 0; end <= seqlen; ++end)
    {
      const int start = 0;
      const int len = end - start;
      const int new_subseq_idx = new_subseq (start, len);
      if (new_subseq_idx != end)
	THROWEXPR ("5' subseq index miscalculation");
      
      if (len > 0)
	connect_r (new_subseq_idx, new_subseq_idx - 1);
    }
  make_global();
  init_pseudovecs();
}

void Fold_envelope::initialise_3prime (int seqlen)
{
  CTAG(6,FOLDENV) << "Initialising 3' fold envelope (HMM) for " << seqlen << "-base sequence\n";

  resize (seqlen);
  const int end = seqlen;
  for (int start = end; start >= 0; --start)
    {
      const int len = end - start;
      const int new_subseq_idx = new_subseq (start, len);
      if (new_subseq_idx != len)
	THROWEXPR ("3' subseq index miscalculation");
      
      if (len > 0)
	connect_l (new_subseq_idx, new_subseq_idx - 1);
    }
  make_global();
  init_pseudovecs();
}

void Fold_envelope::initialise_3prime_banded (int seqlen, int max_subseq_len, int min_loop_len)
{
  // print log message, estimating size & time
  const int L = seqlen;
  const int M = max_subseq_len;
  const int old_prec = CL.precision();
  if (max_subseq_len < 0)
    {
      const double subs = L*L/2;
      const double bifs = L*L*L/6;
      CTAG(6,FOLDENV) << "Initialising full fold envelope for " << seqlen << "-base sequence: approx " << setprecision(2) << subs << " subsequences & " << setprecision(2) << bifs << " bifurcations\n" << setprecision(old_prec);
    }
  else
    {
      const double subs = L*M;
      const double bifs = L*M*M/2;
      CTAG(6,FOLDENV) << "Initialising 3'/banded fold envelope for " << seqlen << "-base sequence, max subseq len " << max_subseq_len << ": approx " << setprecision(2) << subs << " subsequences & " << setprecision(2) << bifs << " bifurcations\n" << setprecision(old_prec);
    }

  // do it
  resize (seqlen);
  for (int end = 0; end < seqlen; ++end)
    for (int len = 0; len <= end && (len <= max_subseq_len || max_subseq_len < 0); ++len)
      {
	const int start = end - len;
	connect_subseq (new_subseq (start, len), cflag_mask (len, min_loop_len));
      }
  const int end = seqlen;
  for (int start = end; start >= 0; --start)
    {
      const int len = end - start;
      connect_subseq (new_subseq (start, len), cflag_mask (len, min_loop_len));
    }
  make_global();
  init_pseudovecs();
}

void Fold_envelope::initialise_local (const Fold_envelope& global, int max_subseq_len, int min_start_len)
{
  if (max_subseq_len < 0)
    CTAG(6,FOLDENV) << "Copying fold envelope for " << global.seqlen() << "-base sequence\n";
  else
    CTAG(6,FOLDENV) << "Making local version of global fold envelope for " << global.seqlen() << "-base sequence, max subseq len " << max_subseq_len << "\n";

  resize (global.seqlen());
  for (int i = 0; i < global.subseqs(); ++i)
    {
      const Subseq& ss = global.subseq[i];
      if (ss.len <= max_subseq_len || max_subseq_len < 0)
	{
	  const int n = new_subseq (ss.start, ss.len);
	  connect_subseq (n, (Connection_flag) ss.in_flags);
	}
    }
  make_local (min_start_len);
  init_pseudovecs();
}

void Fold_envelope::make_local (int min_start_len)
{
  CTAG(2,FOLDENV FOLDENV_START) << "Making fold envelope local, min start length " << min_start_len << ", seqlen " << seqlen() << ", subseqs " << subseq.size() << "\n";
  for_contents (vector<Subseq>, subseq, ss)
    if (ss->len >= min_start_len)
      ss->start_flag = TRUE;
}

void Fold_envelope::make_global()
{
  // don't print log message if envelope size is zero;
  // this avoids calls to uninitialised CLOGSTREAM when Fold_envelope's are declared globally in tests
  if (seqlen())
    CTAG(2,FOLDENV FOLDENV_START) << "Making fold envelope global, seqlen " << seqlen() << ", subseqs " << subseq.size() << "\n";
  if (subseq.back().len != seqlen())
    THROWEXPR ("Can't make envelope global if longest subseq is absent");
  for (int i = 0; i < (int) subseq.size() - 1; ++i)
    subseq[i].start_flag = FALSE;
  subseq.back().start_flag = TRUE;
}

void Fold_envelope::make_local_overlap (int fold_start, int fold_len, int min_overlap_len)
{
  CTAG(2,FOLDENV FOLDENV_START) << "Making fold envelope local for min overlap " << min_overlap_len << " with subseq " << fold_start << "+" << fold_len << ", seqlen " << seqlen() << ", subseqs " << subseq.size() << "\n";
  const int fold_end = fold_start + fold_len;
  for_contents (vector<Subseq>, subseq, ss)
    {
      const int overlap = min (ss->end(), fold_end) - max (ss->start, fold_start);  // can be negative
      ss->start_flag = overlap >= min_overlap_len;
    }
}

bool Fold_envelope::is_global() const
{
  for (int i = 0; i < (int) subseq.size() - 1; ++i)
    if (subseq[i].start_flag)
      return FALSE;
  if (subseq.back().len != seqlen())
    return FALSE;
  return subseq.back().start_flag;
}

void Fold_envelope::initialise_from_subseq_coords (int xlen, const Subseq_coords_set& subseq_coords, int min_loop_len)
{
  CTAG(6,FOLDENV) << "Initialising fold envelope for " << xlen << "-base sequence from subsequence co-ordinate set\n";

  resize (xlen);
  for_const_contents (Subseq_coords_set, subseq_coords, subseq)
    {
      const int n = new_subseq (subseq->start, subseq->len);
      connect_subseq (n, cflag_mask (subseq->len, min_loop_len));
    }
  if (subseq_coords.size())
    make_global();
  else
    CTAG(6,FOLDENV) << "Empty subsequence co-ordinate set; creating null fold envelope\n";
  init_pseudovecs();
}

void Fold_envelope::initialise_from_subseq_coords_with_cflag_mask (const Subseq_coords_set& subseq_coords, const Fold_envelope& cflag_mask_env, int min_loop_len)
{
  CTAG(6,FOLDENV) << "Initialising fold envelope for " << cflag_mask_env.seqlen() << "-base sequence from subsequence co-ordinate set, using only connections from previous fold envelope\n";

  resize (cflag_mask_env.seqlen());
  for_const_contents (Subseq_coords_set, subseq_coords, subseq)
    {
      const int n = new_subseq (subseq->start, subseq->len);
      const int mask_n = cflag_mask_env.find_subseq_idx (subseq->start, subseq->len);
      int env_mask =
	mask_n < 0
	? (int) CFLAG_L_R_LR
	: (int) cflag_mask_env.subseq[mask_n].in_flags;
      connect_subseq (n, (Connection_flag) (env_mask & (int) cflag_mask (subseq->len, min_loop_len)));
    }
  if (subseq_coords.size())
    make_global();
  else
    CTAG(6,FOLDENV) << "Empty subsequence co-ordinate set; creating null fold envelope\n";
  init_pseudovecs();
}

void Fold_envelope::clear()
{
  subseq.clear();
  by_start.clear();
  by_len.clear();
  first_with_end.clear();
  first_after_end.clear();
}

void Fold_envelope::resize (int seqlen)
{
  clear();
  first_with_end = vector<int> (seqlen + 1, 0);
  first_after_end = vector<int> (seqlen + 1, 0);
  by_start = vector<vector<int> > (seqlen + 1, vector<int>());
  by_len = vector<vector<int> > (seqlen + 1, vector<int>());
}

int Fold_envelope::new_subseq (int start, int len)
{
  const int new_subseq_idx = subseq.size();
  subseq.push_back (Subseq (start, len));
  Subseq& new_subseq = subseq[new_subseq_idx];
  if (new_subseq_idx > 0)
    {
      const int prev_subseq_idx = new_subseq_idx - 1;
      const Subseq& prev_subseq = subseq[prev_subseq_idx];
      if (!(prev_subseq < new_subseq))
	THROWEXPR("Subsequence (" << start  << "+" << len << ") added to fold envelope after subsequence (" << prev_subseq.start << "+" << prev_subseq.len << "), inconsistently with inside->outside ordering");
    }
  // update first_with_end & first_after_end
  const int end = new_subseq.end();
  if (first_with_end[end] == 0)
    first_with_end[end] = new_subseq_idx;
  first_after_end[end] = new_subseq_idx + 1;
  // update by_start & by_len
  by_start[start].push_back (new_subseq_idx);
  by_len[len].push_back (new_subseq_idx);
  new_subseq.by_start_index = by_start[start].size() - 1;
  new_subseq.by_len_index = by_len[len].size() - 1;
  // return
  return new_subseq_idx;
}

void Fold_envelope::connect_l (int outside, int inside_l)
{
  Subseq& outsub = subseq[outside];
  Subseq& insub = subseq[inside_l];
  outsub.in[Subseq::CONN_L] = inside_l;
  insub.out[Subseq::CONN_L] = outside;
  outsub.in_flags |= Subseq::CFLAG_L;
  insub.out_flags |= Subseq::CFLAG_L;
}

void Fold_envelope::connect_r (int outside, int inside_r)
{
  Subseq& outsub = subseq[outside];
  Subseq& insub = subseq[inside_r];
  outsub.in[Subseq::CONN_R] = inside_r;
  insub.out[Subseq::CONN_R] = outside;
  outsub.in_flags |= Subseq::CFLAG_R;
  insub.out_flags |= Subseq::CFLAG_R;
}

void Fold_envelope::connect_lr (int outside, int inside_lr)
{
  Subseq& outsub = subseq[outside];
  Subseq& insub = subseq[inside_lr];
  outsub.in[Subseq::CONN_LR] = inside_lr;
  insub.out[Subseq::CONN_LR] = outside;
  outsub.in_flags |= Subseq::CFLAG_LR;
  insub.out_flags |= Subseq::CFLAG_LR;
}

void Fold_envelope::connect_b (int outside, int inside_l, int inside_r)
{
#ifndef DART_USE_BIFURCATION_PSEUDOVECTORS
  Subseq& outsub = subseq[outside];
  Subseq& lsub = subseq[inside_l];
  Subseq& rsub = subseq[inside_r];
  outsub.bif_in.push_back (Subseq::Bifurc_in (inside_l, inside_r));
  lsub.bif_out_r.push_back (Subseq::Bifurc_out_r (outside, inside_r));
  rsub.bif_out_l.push_back (Subseq::Bifurc_out_l (outside, inside_l));
#endif /* DART_USE_BIFURCATION_PSEUDOVECTORS */
}

void Fold_envelope::disconnect_b (int outside, int inside_l, int inside_r)
{
#ifdef DART_USE_BIFURCATION_PSEUDOVECTORS
  THROWEXPR ("disconnect_b unimplemented for bifurcation pseudovectors");
#else /* DART_USE_BIFURCATION_PSEUDOVECTORS */
  Subseq& outsub = subseq[outside];
  Subseq& lsub = subseq[inside_l];
  Subseq& rsub = subseq[inside_r];

  for (int b = 0; b < (int) outsub.bif_in.size(); ++b)
    if (outsub.bif_in[b].l == inside_l && outsub.bif_in[b].r == inside_r)
      {
	outsub.bif_in.erase (outsub.bif_in.begin() + b);
	break;
      }

  for (int b = 0; b < (int) lsub.bif_out_r.size(); ++b)
    if (lsub.bif_out_r[b].out == outside && lsub.bif_out_r[b].r == inside_r)
      {
	lsub.bif_out_r.erase (lsub.bif_out_r.begin() + b);
	break;
      }

  for (int b = 0; b < (int) rsub.bif_out_l.size(); ++b)
    if (rsub.bif_out_l[b].out == outside && rsub.bif_out_l[b].l == inside_l)
      {
	rsub.bif_out_l.erase (rsub.bif_out_l.begin() + b);
	break;
      }

#endif /* DART_USE_BIFURCATION_PSEUDOVECTORS */
}

bool Fold_envelope::test_consistent() const
{
  // return affirmatively if envelope is null

  if (subseq.size() == 0)
    {
      if (by_start.size() || by_len.size()) return 0;
      return 1;
    }

  // test subseq vector

  for (int i = 0; i < (int) subseq.size(); ++i)
    {
      const Subseq& ss = subseq[i];
      if (!ss.test_coords()) return 0;
      if (!ss.test_flags()) return 0;

      // for each connection, test (1) reciprocity (2) co-ordinates

      if (ss.in_l() >= 0)
	{
	  const Subseq& ssl = subseq [ss.in_l()];
	  if (ssl.out_l() != i) return 0;
	  if (!ss.test_in_l (ssl)) return 0;
	}

      if (ss.in_r() >= 0)
	{
	  const Subseq& ssr = subseq [ss.in_r()];
	  if (ssr.out_r() != i) return 0;
	  if (!ss.test_in_r (ssr)) return 0;
	}

      if (ss.in_lr() >= 0)
	{
	  const Subseq& sslr = subseq [ss.in_lr()];
	  if (sslr.out_lr() != i) return 0;
	  if (!ss.test_in_lr (sslr)) return 0;
	}

      for (int j = 0; j < (int) ss.bif_in.size(); ++j)
	{
	  const int l_idx = ss.bif_in[j].l;
	  const int r_idx = ss.bif_in[j].r;
	  const Subseq& ssl = subseq [l_idx];
	  const Subseq& ssr = subseq [r_idx];
	  bool found_recip_l = 0;
	  for (int k = 0; k < (int) ssr.bif_out_l.size(); ++k)
	    if (ssr.bif_out_l[k].out == i)
	      {
		if (ssr.bif_out_l[k].l == l_idx)
		  {
		    if (found_recip_l) return 0;
		    found_recip_l = 1;
		  }
		else return 0;
	      }
	  if (!found_recip_l) return 0;
	  bool found_recip_r = 0;
	  for (int k = 0; k < (int) ssl.bif_out_r.size(); ++k)
	    if (ssl.bif_out_r[k].out == i)
	      {
		if (ssl.bif_out_r[k].r == r_idx)
		  {
		    if (found_recip_r) return 0;
		    found_recip_r = 1;
		  }
		else return 0;
	      }
	  if (!found_recip_r) return 0;
	  if (!ss.test_in_b (ssl, ssr)) return 0;
	}

      if (ss.out_l() >= 0)
	{
	  const Subseq& ssl = subseq [ss.out_l()];
	  if (ssl.in_l() != i) return 0;
	  if (!ssl.test_in_l (ss)) return 0;
	}

      if (ss.out_r() >= 0)
	{
	  const Subseq& ssr = subseq [ss.out_r()];
	  if (ssr.in_r() != i) return 0;
	  if (!ssr.test_in_r (ss)) return 0;
	}

      if (ss.out_lr() >= 0)
	{
	  const Subseq& sslr = subseq [ss.out_lr()];
	  if (sslr.in_lr() != i) return 0;
	  if (!sslr.test_in_lr (ss)) return 0;
	}

      for (int j = 0; j < (int) ss.bif_out_l.size(); ++j)
	{
	  const int l_idx = ss.bif_out_l[j].l;
	  const int out_idx = ss.bif_out_l[j].out;
	  const Subseq& ssl = subseq [l_idx];
	  const Subseq& ssout = subseq [out_idx];
	  bool found_reciprocal = 0;
	  for (int k = 0; k < (int) ssout.bif_in.size(); ++k)
	    if (ssout.bif_in[k].r == i)
	      {
		if (ssout.bif_in[k].l == l_idx) found_reciprocal = 1;
		else return 0;
	      }
	  if (!found_reciprocal) return 0;
	  if (!ssout.test_in_b (ssl, ss)) return 0;
	}

      for (int j = 0; j < (int) ss.bif_out_r.size(); ++j)
	{
	  const int r_idx = ss.bif_out_r[j].r;
	  const int out_idx = ss.bif_out_r[j].out;
	  const Subseq& ssr = subseq [r_idx];
	  const Subseq& ssout = subseq [out_idx];
	  bool found_reciprocal = 0;
	  for (int k = 0; k < (int) ssout.bif_in.size(); ++k)
	    if (ssout.bif_in[k].l == i)
	      {
		if (ssout.bif_in[k].r == r_idx) found_reciprocal = 1;
		else return 0;
	      }
	  if (!found_reciprocal) return 0;
	  if (!ssout.test_in_b (ss, ssr)) return 0;
	}

      // test sort order

      if (i > 0)
	{
	  const Subseq& ssprev = subseq[i - 1];
	  if (ss.end() < ssprev.end() || (ss.end() == ssprev.end() && ss.start > ssprev.start)) return 0;
	}
    }

  if (subseq.front().start != 0 || subseq.back().start != 0) return 0;
  
  // test by_start & by_len vectors

  const int seqlen = subseq.back().len;
  if ((int) by_start.size() != seqlen + 1 || (int) by_len.size() != seqlen + 1) return 0;
  int by_start_count = 0;
  int by_len_count = 0;
  for (int i = 0; i <= seqlen; ++i)
    {
      for (int j = 0; j < (int) by_start[i].size(); ++j)
	if (subseq[by_start[i][j]].start != i || (j > 0 && subseq[by_start[i][j]].len <= subseq[by_start[i][j-1]].len) || subseq[by_start[i][j]].by_start_index != j)
	  return 0;
      for (int j = 0; j < (int) by_len[i].size(); ++j)
	if (subseq[by_len[i][j]].len != i || (j > 0 && subseq[by_len[i][j]].start <= subseq[by_len[i][j-1]].start) || subseq[by_len[i][j]].by_len_index != j)
	  return 0;
      by_start_count += by_start[i].size();
      by_len_count += by_len[i].size();
    }
  if (by_start_count != (int) subseq.size() || by_len_count != (int) subseq.size()) return 0;

  // all done

  return 1;
}

bool Fold_envelope::test_fold_string (const sstring& fold_string) const
{
  if (subseq.back().len != (int) fold_string.size()) return 0;
  for (int i = 0; i < (int) subseq.size(); ++i)
    if (!subseq[i].test_fold_string (fold_string))
      return 0;
  return 1;
}

bool Fold_envelope::test_grounded() const
{
  THROWEXPR ("Unimplemented");
  return 0;
}

bool Fold_envelope::test_bifs() const
{
  for_const_contents (vector<Subseq>, subseq, ss)
    {
      for_const_contents (Subseq::Bifurc_in_pseudovec, ss->bif_in, b)
	if (subseq[(*b).l].start != ss->start
	    || subseq[(*b).r].end() != ss->end()
	    || subseq[(*b).l].end() != subseq[(*b).r].start)
	  {
	    CLOGERR << "Failed bifurcation check; ss=" << ss->terser_desc()
		    << " l=" << subseq[(*b).l].terser_desc()
		    << " r=" << subseq[(*b).r].terser_desc()
		    << "\n";
	    return false;
	  }

      for_const_contents (Subseq::Bifurc_outl_pseudovec, ss->bif_out_l, b)
	if (subseq[(*b).l].end() != ss->start
	    || subseq[(*b).out].end() != ss->end()
	    || subseq[(*b).out].start != subseq[(*b).l].start)
	  {
	    CLOGERR << "Failed bifurcation check; ss=" << ss->terser_desc()
		    << " out=" << subseq[(*b).out].terser_desc()
		    << " l=" << subseq[(*b).l].terser_desc()
		    << "\n";
	    return false;
	  }

      for_const_contents (Subseq::Bifurc_outr_pseudovec, ss->bif_out_r, b)
	if (subseq[(*b).out].start != ss->start
	    || subseq[(*b).r].start != ss->end()
	    || subseq[(*b).out].end() != subseq[(*b).r].end())
	  {
	    CLOGERR << "Failed bifurcation check; ss=" << ss->terser_desc()
		    << " out=" << subseq[(*b).out].terser_desc()
		    << " r=" << subseq[(*b).r].terser_desc()
		    << "\n";
	    return false;
	  }
    }
  return true;
}

void Fold_envelope::dump (ostream& o) const
{
  for (int i = 0; i < (int) subseq.size(); ++i)
    {
      const Subseq& ss = subseq[i];
      o << i << ": (" << ss.start << "+" << ss.len << ") ";
      o << "CIN(L:" << ss.in_l() << ",R:" << ss.in_r() << ",LR:" << ss.in_lr() << ") ";
      o << "COUT(L:" << ss.out_l() << ",R:" << ss.out_r() << ",LR:" << ss.out_lr() << ") ";
      o << "BIN(";
      for_const_contents (Subseq::Bifurc_in_pseudovec, ss.bif_in, b)
	o << (*b).l << '/' << (*b).r << ' ';
      o << ") BOUTL(";
      for_const_contents (Subseq::Bifurc_outl_pseudovec, ss.bif_out_l, b)
	o << (*b).out << '/' << (*b).l << ' ';
      o << ") BOUTR(";
      for_const_contents (Subseq::Bifurc_outr_pseudovec, ss.bif_out_r, b)
	o << (*b).out << '/' << (*b).r << ' ';
      o << ")\n";
    }
}

int Fold_envelope::find_subseq_idx (int start, int len) const
{
  const Subseq_coords coords (start, len);
  const int end = coords.end(), sl = seqlen();
  if (start < 0 || len < 0 || start > sl || end > sl)
    return -1;
  const vector<Subseq>::const_iterator
    b = subseq.begin() + first_with_end[end],
    e = subseq.begin() + first_after_end[end];
  const vector<Subseq>::const_iterator iter = lower_bound (b, e, coords);
  if (iter == e)
    return -1;
  if (iter->start != start || iter->len != len)
    return -1;
  return iter - subseq.begin();
}

int Fold_envelope::first_start_subseq() const
{
  for (int i = 0; i < (int) subseq.size(); ++i)
    if (subseq[i].start_flag)
      return i;
  THROWEXPR ("No start transition into fold envelope");
  return (int) subseq.size();
}

int Fold_envelope::find_l (int subseq_idx) const
{
  const Subseq& s = subseq[subseq_idx];
  return find_subseq_idx (s.start + 1, s.len - 1);
}

int Fold_envelope::find_r (int subseq_idx) const
{
  const Subseq& s = subseq[subseq_idx];
  return find_subseq_idx (s.start, s.len - 1);
}

int Fold_envelope::find_lr (int subseq_idx) const
{
  const Subseq& s = subseq[subseq_idx];
  return find_subseq_idx (s.start + 1, s.len - 2);
}

void Fold_envelope::connect_subseq (int subseq_idx, Connection_flag conn)
{
  const Subseq& s = subseq[subseq_idx];
  const int start = s.start;
  const int end = s.end();
  const int l_idx = find_l (subseq_idx);
  const int r_idx = find_r (subseq_idx);
  const int lr_idx = find_lr (subseq_idx);
  if ((conn & CFLAG_L) && l_idx >= 0) connect_l (subseq_idx, l_idx);
  if ((conn & CFLAG_R) && r_idx >= 0) connect_r (subseq_idx, r_idx);
  if ((conn & CFLAG_LR) && lr_idx >= 0) connect_lr (subseq_idx, lr_idx);
  for (int bpos = start; bpos <= end; ++bpos)
    {
      const int bl_idx = find_subseq_idx (start, bpos - start);
      const int br_idx = find_subseq_idx (bpos, end - bpos);
      if (bl_idx >= 0 && br_idx >= 0)
	connect_b (subseq_idx, bl_idx, br_idx);
    }
}

void Fold_envelope::connect_all (int min_loop_len)
{
  for (int subseq_idx = 0; subseq_idx < (int) subseq.size(); ++subseq_idx)
    connect_subseq (subseq_idx, cflag_mask (subseq[subseq_idx].len, min_loop_len));
}

void Fold_envelope::remove_empty_bifurcations()
{
  for (int subseq_idx = 0; subseq_idx < (int) subseq.size(); ++subseq_idx)
    {
      const Subseq& s = subseq[subseq_idx];
      disconnect_b (subseq_idx, by_len[0][s.start], subseq_idx);
      disconnect_b (subseq_idx, subseq_idx, by_len[0][s.end()]);
    }
}

Fold_envelope& Fold_envelope::operator= (const Fold_envelope& env)
{
  subseq = env.subseq;
  by_start = env.by_start;
  by_len = env.by_len;
  first_with_end = env.first_with_end;
  first_after_end = env.first_after_end;
  init_pseudovecs();
  return *this;
}

void Fold_envelope::swap (Fold_envelope& env)
{
  subseq.swap (env.subseq);
  by_len.swap (env.by_len);
  by_start.swap (env.by_start);
  first_with_end.swap (env.first_with_end);
  first_after_end.swap (env.first_after_end);
}

bool Fold_envelope::operator== (const Fold_envelope& env) const
{
  return subseq == env.subseq;
}

bool Fold_envelope::operator!= (const Fold_envelope& env) const
{
  return !(*this == env);
}

sstring Fold_envelope::to_envelope_string() const
{
  sstring s;
  s << by_start.size();
  for_const_contents (vector<Subseq>, subseq, ss)
    {
      const int munged_flags = ss->in_flags | (ss->start_flag ? CFLAG_TOTAL : 0);
      s << ' ' << ss->start << (char) (munged_flags + CONN_BASE_CHAR) << ss->len;
    }
  return s;
}

Regexp foldenv_re ("([0123456789]+)([@ABCDEFGHIJKLMNO])([0123456789]+)");
void Fold_envelope::initialise_from_envelope_string (const sstring& str)
{
  const vector<sstring> f = str.split();
  const int seqlen_plus1 = f[0].to_int();
  clear();
  if (seqlen_plus1 > 0)
    {
      resize (seqlen_plus1 - 1);
      subseq.reserve (f.size() - 1);
      for (int i = 1; i < (int) f.size(); ++i)
	{
	  if (!foldenv_re.Match (f[i].c_str()))
	    THROWEXPR ("Bad field in envelope string: " << f[i]);
	  const int start = foldenv_re[1].to_int();
	  const int munged_flags = foldenv_re[2][0] - CONN_BASE_CHAR;
	  const int len = foldenv_re[3].to_int();
	  const int n = new_subseq (start, len);
	  connect_subseq (n, (Connection_flag) (munged_flags & CFLAG_L_R_LR));
	  subseq[n].start_flag = munged_flags & CFLAG_TOTAL;
	}
    }
}

typedef map<int,set<int> > StartByLen;
void Fold_envelope::initialise_from_gff_list (const GFF_list& gff_list, const sstring& seqname, int min_loop_len)
{
  // StartByLen structure avoids duplicates
  StartByLen bylen;
  for_const_contents (GFF_list, gff_list, gff)
    if (gff->seqname == seqname)
      {
	const int start = gff->start - 1;
	const int len = gff->end + 1 - gff->start;
	bylen[len].insert (start);
	if (len == 1)
	  bylen[0].insert (start);
      }
  // build the fold envelope
  clear();
  for_const_contents (StartByLen, bylen, len_startset)
    {
      const int len = (*len_startset).first;
      for_const_contents (set<int>, (*len_startset).second, start)
	{
	  const int n = new_subseq (*start, len);
	  connect_subseq (n, cflag_mask (len, min_loop_len));
	}
    }
  make_global();
}

void Fold_envelope::to_gff (GFF_list& gff_list, const sstring& seqname) const
{
  GFF gff; 
  gff.seqname = seqname;
  for_const_contents (vector<Subseq>, subseq, s)
    if (s->len > 0)
      {
	gff.start = s->start + 1;
	gff.end = s->start + s->len;
	gff_list.push_back (gff);
      }
}

void Fold_envelope::render_dotplot_from_counts (ostream& out, const Subseq_coords_count& count, const Biosequence& seq, int max_level, bool use_ansi_color, const Subseq_coords_annot* annot)
{
  // prepare chars
  const int cols = 8;
  sstring esc (2);
  esc[0] = 27;
  esc[1] = '[';
  const double scale = (double) (cols-1) / (double) max_level;
  sstring inverse_text;
  sstring normal_text;
  vector<sstring> color_text (cols);
  vector<sstring> cell_text (cols);
  if (use_ansi_color)
    {
      const int colseq[] = { 31, 33, 32, 34, 36, 35, 37 };  // rainbow
      for (int k = 1; k < cols; ++k)
	color_text[k] << esc << colseq[k-1] << 'm';
      if (!annot)
	for (int k = 0; k < cols; ++k)
	  cell_text[k] = " ";
      inverse_text << esc << "7m";
      normal_text << esc << "27m";
    }
  else
    {
      cell_text[0] = " ";
      for (int k = 1; k < cols; ++k)
	cell_text[k] << k;
    }
  // render, with axes
  out << "Vertical axis = start, horizontal axis = end\n";
  const bool double_width = (annot != 0) && !use_ansi_color;
  if (double_width)
    {
      out << " *";
      for (int i = 0; i < (int) seq.size(); ++i)
	out << ' ' << seq[i];
      out << '\n';
    }
  else
    out << " *" << seq << '\n';
  const sstring pad_text (double_width ? 2 : 1, '\\');
  for (int j = 0; j <= (int) seq.size(); ++j)
    {
      out << (char) (j==0 ? '*' : seq[j-1]);
      for (int i = 0; i < j-1; ++i)
	out << pad_text;
      if (j > 0)
	out << pad_text;
      int last_col = -1;
      for (int i = j; i <= (int) seq.size(); ++i)
	{
	  // make Subseq_coords
	  const Subseq_coords my_coords (j, i-j);
	  // determine color
	  Subseq_coords_count::const_iterator val_iter = count.find (my_coords);
	  const int val = val_iter == count.end() ? 0 : val_iter->second;
	  const int col = val<=0 ? 0 : min ((int) ((scale * (double) val) + 1), (int) (cols-1));
	  // render cell
	  if (i == 0 || col != last_col)
	    {
	      if (last_col >= 0) out << color_text[last_col];
	      if (last_col <= 0 && col > 0) out << normal_text << inverse_text;
	      else if (last_col > 0 && col <= 0) out << inverse_text << normal_text;
	      out << color_text[col];
	      last_col = col;
	    }
	  out << cell_text[col];
	  // print annotation
	  if (annot)
	    {
	      Subseq_coords_annot::const_iterator annot_iter = annot->find (my_coords);
	      if (annot_iter == annot->end())
		out << ' ';
	      else
		out << annot_iter->second;
	    }
	}
      if (last_col >= 0) out << color_text[last_col] << inverse_text << normal_text;
      out << (char) (j==0 ? '*' : seq[j-1]);
      out << '\n';
    }
  out << " *" << seq << '\n';
  out << "Color code: [0] " << normal_text << inverse_text;
  for (int i = 1; i < cols; ++i)
    out << color_text[i] << ' ' << color_text[i];
  out << inverse_text << normal_text << "[" << max_level << "]\n";
}

void Fold_envelope::render_dotplot (ostream& out, const Biosequence& seq, bool use_ansi_color) const
{
  Subseq_coords_count count;
  for_const_contents (vector<Subseq>, subseq, ss)
    count[*ss] = 1;
  render_dotplot_from_counts (out, count, seq, 1, use_ansi_color, 0);
}

void Fold_envelope::render_annotated_dotplot (ostream& out, const Biosequence& seq, bool use_ansi_color) const
{
  Subseq_coords_count count;
  Subseq_coords_annot annot;
  int max_count = 0;
  for_const_contents (vector<Subseq>, subseq, ss)
    {
      max_count = max (max_count, count[*ss] = ss->bif_in.size() + 1);
      annot[*ss] = (char) (ss->in_flags + (ss->start_flag ? 'A' : 'a'));
    }
  render_dotplot_from_counts (out, count, seq, max_count, use_ansi_color, &annot);
  out << "Connections (colour=#bifurcs+1): ";
  vector<sstring> flags;
  for (int i = 0; i < CFLAG_TOTAL*2; ++i)
    {
      flags.clear();
      if (i & CFLAG_L) flags.push_back (sstring ("L"));
      if (i & CFLAG_R) flags.push_back (sstring ("R"));
      if (i & CFLAG_LR) flags.push_back (sstring ("LR"));
      if (i & CFLAG_TOTAL) flags.push_back (sstring ("S"));
      if (i > 0) out << ", ";
      out << (char) ((i&CFLAG_L_R_LR) + (i&CFLAG_TOTAL?'A':'a')) << '=' << sstring::join (flags, "+");
      if (i == 0) out << "none";
    }
  out << "\n";
}

typedef map<sstring,GFF_list> GFF_index;
void Fold_envelope_database::from_gff (const GFF_list& gff_list)
{
  GFF_index by_seqname;
  for_const_contents (GFF_list, gff_list, gff)
    by_seqname[gff->seqname].push_back (*gff);
  for_const_contents (GFF_index, by_seqname, seqname_gfflist)
    {
      Fold_envelope fold_env;
      fold_env.initialise_from_gff_list (seqname_gfflist->second, seqname_gfflist->first);
      (*this)[seqname_gfflist->first] = fold_env;
    }
}

void Fold_envelope_database::to_gff (GFF_list& gff_list) const
{
  for_const_contents (Fold_envelope_database, *this, seqname_foldenv)
    seqname_foldenv->second.to_gff (gff_list, seqname_foldenv->first);
}

void Fold_envelope_database::load (const char* filename)
{
  ifstream in_file (filename);
  if (!in_file) THROWEXPR ("Fold envelope file not found: " << filename);
  load (in_file);
}

void Fold_envelope_database::load (istream& ins)
{
  sstring s;
  while (ins && !ins.eof())
    {
      s.getline (ins);
      const vector<sstring> f = s.split (" ", TRUE, 2);
      if (f.size() != 2) continue;
      (*this) [f[0]] = Fold_envelope();
      (*this) [f[0]] . initialise_from_envelope_string (f[1]);
    }
}

void Fold_envelope_database::save (ostream& outs) const
{
  for_const_contents (Fold_envelope_database, *this, name_env)
    outs << name_env->first << ' ' << name_env->second.to_envelope_string() << '\n';
}

void Fold_envelope_database::from_gff (const char* gff_filename)
{
  GFF_list gff_list;
  gff_list.load (gff_filename);
  from_gff (gff_list);
}

void Fold_envelope_database::to_gff (const char* gff_filename) const
{
  GFF_list gff_list;
  to_gff (gff_list);
  gff_list.save (gff_filename);
}
