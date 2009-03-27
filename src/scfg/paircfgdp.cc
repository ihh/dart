#include <stack>
#include <deque>
#include "scfg/paircfgdp.h"
#include "util/math_fn.h"
#include "util/vector_output.h"

#define CFGDP_REPORT_INTERVAL 3000000  /* number of bifurcations between logfile messages during DP */

typedef unsigned long long Bifcount;  /* Bifurcation count */

typedef vector<const vector<int>*> Subseq_index_array;

Pair_CFG_filter::Pair_CFG_filter (const Pair_CFG_scores& cfg) :
  allowed_in_states (FlagPairRange),
  allowed_out_states (FlagPairRange),
  allowed_dest_states (FlagPairRange, vector<vector<int> > (cfg.states())),
  allowed_src_states (FlagPairRange, vector<vector<int> > (cfg.states())),
  left_parent (cfg.left_parent()),
  right_parent (cfg.right_parent())
{
  const vector<int> emit_states = cfg.emit_states();
  const vector<int> nonemit_states = cfg.nonemit_states();

  for (int xflag = 0; xflag < FlagRange; ++xflag)
    for (int yflag = 0; yflag < FlagRange; ++yflag)
      {
	const int fp = flag_pair (xflag, yflag);
	for_const_contents (vector<int>, emit_states, s)
	  {
	    bool allowed = 1;
	    const State_type t = cfg.state_type[*s];
	    switch (t & EmitXLR)
	      {
	      case EmitXL:
		if (!(xflag & Subseq::CFLAG_L)) allowed = 0;
		break;
	      case EmitXR:
		if (!(xflag & Subseq::CFLAG_R)) allowed = 0;
		break;
	      case EmitXLR:
		if (!(xflag & Subseq::CFLAG_LR)) allowed = 0;
		break;
	      default:
		break;
	      }
	    switch (t & EmitYLR)
	      {
	      case EmitYL:
		if (!(yflag & Subseq::CFLAG_L)) allowed = 0;
		break;
	      case EmitYR:
		if (!(yflag & Subseq::CFLAG_R)) allowed = 0;
		break;
	      case EmitYLR:
		if (!(yflag & Subseq::CFLAG_LR)) allowed = 0;
		break;
	      default:
		break;
	      }
	    if (allowed)
	      {
		allowed_in_states[fp].push_back(*s);
		allowed_out_states[fp].push_back(*s);
	      }
	  }
	allowed_out_states[fp].insert (allowed_out_states[fp].end(), nonemit_states.begin(), nonemit_states.end());  // outside algorithm fills null states in topological order
	allowed_in_states[fp].insert (allowed_in_states[fp].end(), nonemit_states.rbegin(), nonemit_states.rend());  // inside algorithm fills null states in reverse topological order
      }
  for (int xflag = 0; xflag < FlagRange; ++xflag)
    for (int yflag = 0; yflag < FlagRange; ++yflag)
      {
	const int fp = flag_pair (xflag, yflag);
	allowed_dest_states[fp] = cfg.selected_outgoing_states (allowed_in_states[fp]);
	allowed_src_states[fp] = cfg.selected_incoming_states (allowed_in_states[fp]);
      }
}

Pair_CFG_DP_allocator_base::Pair_CFG_DP_allocator_base (const Named_profile& npx,
							const Named_profile& npy,
							const Fold_envelope& xenv,
							const Fold_envelope& yenv,
							const Pair_envelope& pair_env,
							const Pair_CFG_scores& cfg)
  :  xenv (xenv),
     yenv (yenv),
     pair_env (pair_env),
     npx (npx),
     npy (npy),
     xdsq (npx.dsq),
     ydsq (npy.dsq),
     xmeta (npx.meta_sc),
     ymeta (npy.meta_sc),
     cfg (cfg),
     filter (cfg),
     xsubseqs (xenv.subseq.size()),
     ysubseqs (yenv.subseq.size()),
     states (cfg.states()),
     xlen (npx.size()),
     ylen (npy.size()),
     total_cells (0),
     total_bifs (0),
     bytes_allocated (0)
{ }

void Pair_CFG_DP_allocator_base::count_cells()
{
  // count cells inside pair envelope
  register Bifcount n_cells = 0;  // I don't normally use 'register', but I want this to be invisibly fast
  register Bifcount n_bifs = 0;
  Subseq_index_array valid_ystartpos;
  valid_ystartpos.reserve (ylen + 1);
  for (int x = 0; x <= xlen; ++x)
    {
      // prepare list of valid y start positions for this xstartpos
      valid_ystartpos.clear();
      for (int y = 0; y <= ylen; ++y)
	if (startpos_in_pair_envelope (x, y))
	  valid_ystartpos.push_back (&yenv.by_start[y]);
      // loop through xsubseqs with this startpos
      for_const_contents (vector<int>, xenv.by_start[x], xi)
	{
	  const Subseq& xsubseq = xenv.subseq[*xi];
	  const int xendpos = xsubseq.end();
	  for_const_contents (Subseq_index_array, valid_ystartpos, ysp)
	    for_const_contents (vector<int>, **ysp, yi)
	    {
	      const Subseq& ysubseq = yenv.subseq[*yi];
	      if (endpos_in_pair_envelope (xendpos, ysubseq.end()))
		{
		  ++n_cells;
		  n_bifs += xsubseq.bif_in.size() * ysubseq.bif_in.size();
		}
	    }
	}
    }
  total_cells = n_cells;
  total_bifs = n_bifs;
}

Pair_CFG_DP_sparse_allocator::Pair_CFG_DP_sparse_allocator (const Named_profile& npx,
							    const Named_profile& npy,
							    const Fold_envelope& xenv,
							    const Fold_envelope& yenv,
							    const Pair_envelope& pair_env,
							    const Pair_CFG_scores& cfg)
  : Pair_CFG_DP_allocator_base (npx, npy, xenv, yenv, pair_env, cfg),
    _dummy_inf_score (-InfinityScore),
    xenv_begin (xenv.subseq.begin()),
    yenv_begin (yenv.subseq.begin())
{ }

void Pair_CFG_DP_sparse_allocator::alloc()
{
  // allocate lookup table
  bytes_allocated = sizeof (array2d<By_start_indices>) + (xlen+1) * (ylen+1) * sizeof (By_start_indices);
  CTAG(4,ALLOC) << "Trying to allocate " << xlen+1 << "*" << ylen+1 << " lookup table (approx " << bytes_allocated << " bytes)\n";
  by_start_indices.resize (xlen+1, ylen+1, By_start_indices());

  // figure out how many cells we want to allocate
  unsigned int cells_wanted = 0;
  for (int x = 0; x <= xlen; ++x)
    for (int y = 0; y <= ylen; ++y)
      if (startpos_in_pair_envelope (x, y))
	{
	  const int xs = xenv.by_start[x].size();
	  const int ys = yenv.by_start[y].size();
	  By_start_indices& ind = by_start_indices (x, y);
	  ind.offset = cells_wanted;
	  ind.xsize = xs;
	  cells_wanted += states * xs * ys;
	}

  // print a pre-allocation log message
  unsigned int bytes_wanted = cells_wanted * sizeof(Score);
  CTAG(4,ALLOC) << "Trying to allocate main DP table: " << cells_wanted << " cells (" << bytes_wanted << " bytes)\n";
  
  // allocate the cells
  _cell.resize (cells_wanted, -InfinityScore);
  
  // success
  bytes_allocated += bytes_wanted;
}

Pair_CFG_DP_dense_allocator::Pair_CFG_DP_dense_allocator (const Named_profile& npx,
							  const Named_profile& npy,
							  const Fold_envelope& xenv,
							  const Fold_envelope& yenv,
							  const Pair_envelope& pair_env,
							  const Pair_CFG_scores& cfg)
  : Pair_CFG_DP_allocator_base (npx, npy, xenv, yenv, pair_env, cfg)
{ }

void Pair_CFG_DP_dense_allocator::alloc()
{
  const int sz = states * xsubseqs * ysubseqs;
  const unsigned int bytes_wanted = sz * sizeof(Score);
  CTAG(4,ALLOC) << "Trying to allocate approx " << bytes_wanted << " bytes\n";
  _cell.resize (sz, -InfinityScore);
  bytes_allocated = bytes_wanted;
}

Pair_CFG_DP_matrix::Pair_CFG_DP_matrix (const Named_profile& npx,
					const Named_profile& npy,
					const Fold_envelope& xenv,
					const Fold_envelope& yenv,
					const Pair_CFG_scores& cfg,
					const Pair_envelope& pair_env,
					bool local) :
  CFG_DP_MATRIX_BASE_CLASS (npx, npy, xenv, yenv, pair_env, cfg),
  local (local),
  _show_out (0)
{
  // test validity of Pair SCFG
  if (!cfg.test_valid()) THROWEXPR ("Pair CFG invalid");
  // initialise final score
  if (xlen == 0 && ylen == 0)
    final_score = cfg.start_to_end();
  else
    final_score = -InfinityScore;
  // log the sequences, envelopes & CFG
  if (CTAGGING(-3,CFGDP_GRAMMAR))
    { CL << "Grammar (in CFG matrix constructor):\n"; cfg.show (CL); }
  if (CTAGGING(-3,CFGDP_SEQ))
    { CL << "Sequences (in CFG matrix constructor):\nX: " << npx.seq << "\nY: "<< npy.seq << "\n"; }
  if (CTAGGING(-3,CFGDP_FOLDENV))
    { CL << "X-envelope (in CFG matrix constructor):\n"; xenv.dump(CL); CL << "Y-envelope (in CFG matrix constructor):\n"; yenv.dump(CL); }
  // if requested (via logger), check consistency of fold envelopes, including bifurcations
  if (CTAGGING(1,CFGDP_TEST_FOLDENV))
    {
      CTAG(1,CFGDP_TEST_FOLDENV) << "Testing consistency of fold envelope X\n";
      if (!xenv.test_bifs())
	{
	  xenv.dump(CL);
	  THROWEXPR ("X envelope failed bifurcation consistency test");
	}
      CTAG(1,CFGDP_TEST_FOLDENV) << "Testing consistency of fold envelope Y\n";
      if (!yenv.test_bifs())
	{
	  yenv.dump(CL);
	  THROWEXPR ("Y envelope failed bifurcation consistency test");
	}
      CTAG(1,CFGDP_TEST_FOLDENV) << "Fold envelopes passed consistency tests\n";
    }
  // print ALLOC log messages, allocate & count cells
  if (CTAGGING(5,ALLOC PREALLOC))
    CL << "Allocating Pair SCFG DP matrix:  " << xsubseqs << " * " << ysubseqs << "  cells by  " << states << "  states\n";
  alloc();
  if (CTAGGING(5,ALLOC POSTALLOC))
    CL << "Allocated Pair SCFG DP matrix (approx " << bytes_allocated << " bytes); counting subseqs\n";
  count_cells();
  if (CTAGGING(5,ALLOC CELLCOUNT))
    CL << "Counted " << total_cells << " subseqs and " << total_bifs << " bifurcations in the envelope\n";
}

void Pair_CFG_DP_matrix::show (ostream& o) const
{
  int old_prec = o.precision(3);
  save_flags (o);

  const int page_len = 40;
  sstring page_header;
  page_header << "xi  start+len  xseq XL   XR   XLR  yi  start+len  yseq YL   YR   YLR  ";
  for (int s = 0; s < cfg.states(); ++s)
    {
      sstring state_desc;
      state_desc << s << '(' << state_type_string (cfg.state_type[s]) << ')';
      right_align (page_header);
      page_header.width (15);
      page_header << state_desc;
    }
  page_header << "\n";

  Biosequence xseq;
  Biosequence yseq;
  cfg.alphabet().dsq2seq (xdsq, xseq);
  cfg.alphabet().dsq2seq (ydsq, yseq);

  int lines = 0;
  for (int xsubseq_idx = 0; xsubseq_idx < xsubseqs; ++xsubseq_idx)
    for (int ysubseq_idx = 0; ysubseq_idx < ysubseqs; ++ysubseq_idx)
      {
	if (!in_pair_envelope (xenv.subseq[xsubseq_idx], yenv.subseq[ysubseq_idx]))
	  continue;

	if (lines++ % page_len == 0)
	  {
	    if (lines > 1) o << "\n";
	    o << page_header;
	  }

	o << xenv.subseq[xsubseq_idx].terse_desc (xsubseq_idx, xseq, _show_out);
	o << yenv.subseq[ysubseq_idx].terse_desc (ysubseq_idx, yseq, _show_out);
	right_align (o);
	for (int s = 0; s < cfg.states(); ++s)
	  {
	    o.width (15);
	    ShowScore (read_cell (s, xsubseq_idx, ysubseq_idx), o);
	  }
	o << "\n";
      }
  o << "Final score: ";
  ShowScore (final_score, o);
  o << "\n";
  restore_flags (o);
  o.precision (old_prec);
}

sstring Pair_CFG_DP_matrix::cfg_dump() const
{
  // until sstring::operator<<(int) bug is fixed, we use the following hack
  ostream& dump = CLOGERR;
  //  sstring dump;
  dump << "Pair CFG:\n";
  cfg.show (dump);
  dump << "Envelope X:\n";
  xenv.dump (dump);
  dump << "Envelope Y:\n";
  yenv.dump (dump);
  dump << "Pair CFG DP matrix:\n";
  show (dump);
  //  return dump;
  return sstring();
}

void Pair_CFG_DP_matrix::show_expanded_trace (const Pair_CFG_parse_tree& parse_tree, ostream& out, bool outside) const
{
  typedef pair<int,sstring> Branch_indent;
  stack<Branch_indent> next_branch;
  next_branch.push (Branch_indent (0, sstring("")));
  while (!next_branch.empty())
    {
      const Branch_indent bi = next_branch.top();
      const Pair_CFG_branch& b = parse_tree[bi.first];
      const sstring& indent = bi.second;
      next_branch.pop();
      int xl = b.xl;
      int xr = b.xr;
      int yl = b.yl;
      int yr = b.yr;
      for (int pos = 0; pos < (int) b.path.size(); ++pos)
	{
	  if (indent.size())
	    {
	      for (int i = 0; i < (int) indent.size() - 1; ++i)
		out << indent[i];
	      if (pos == 0)
		out << "+-";
	      else
		out << indent[indent.size()-1] << ' ';
	    }
	  const int s = b.path[pos];
	  const State_type t = s < 0 ? Undefined : cfg.state_type[s];
	  int xsubseq_idx = xenv.find_subseq_idx (xl, xr-xl);
	  int ysubseq_idx = yenv.find_subseq_idx (yl, yr-yl);
	  if (outside && is_emit_type(t))
	    {
	      xsubseq_idx = dest_xsubseq_idx (xsubseq_idx, t);
	      ysubseq_idx = dest_ysubseq_idx (ysubseq_idx, t);
	    }
	  out << "State " << s;
	  out << ", xsubseq " << xenv.subseq[xsubseq_idx].start << "+" << xenv.subseq[xsubseq_idx].len;
	  out << ", ysubseq " << yenv.subseq[ysubseq_idx].start << "+" << yenv.subseq[ysubseq_idx].len;
	  out << ", score ";
	  Score cell_sc;
	  if (s == Grammar_state_enum::Start)
	    cell_sc = final_score;
	  else if (s == Grammar_state_enum::End)
	    cell_sc = 0;
	  else
	    cell_sc = read_cell (s, xsubseq_idx, ysubseq_idx);
	  ShowScore(cell_sc,out);
	  out << ".";
	  if (outside)
	    {
	      if (pos > 0)
		{
		  const int r = b.path[pos-1];
		  const Score trans_sc = cfg.transition (r, s);
		  out << "\tTransition ";
		  ShowScore(trans_sc,out);
		  out << " (" << r << "=>" << s << ")";
		}
	    }
	  else
	    if (pos < (int) b.path.size() - 1)
	      {
		const int d = b.path[pos+1];
		const Score trans_sc = cfg.transition (s, d);
		out << "\tTransition ";
		ShowScore(trans_sc,out);
		out << " (" << s << "=>" << d << ")";
	      }
	  if (s >= 0)
	    {
	      if (is_emit_type(t))
		{
		  int emit_idx = 0;
		  vector<int> emitted (4, -1);
		  if (t & EmitXL) emit_idx += cfg.emit_xl_mul(t) * (emitted[0] = xdsq[xl++]);
		  if (t & EmitXR) emit_idx += cfg.emit_xr_mul(t) * (emitted[1] = xdsq[--xr]);
		  if (t & EmitYL) emit_idx += cfg.emit_yl_mul(t) * (emitted[2] = ydsq[yl++]);
		  if (t & EmitYR) emit_idx += cfg.emit_yr_mul(t) * (emitted[3] = ydsq[--yr]);
		  sstring emit_string ("....");
		  for (int i = 0; i < 4; ++i)
		    if (emitted[i] >= 0)
		      emit_string[i] = cfg.alphabet().int2char_uc (emitted[i]);
		  const Score emit_sc = cfg.emit[s][emit_idx];
		  out << "\tEmit ";
		  ShowScore(emit_sc,out);
		  out << " (" << emit_string << ")";
		}
	    }
	  out << "\n";
	}
      if (b.has_children())
	{
	  out << indent << " |\n";
	  sstring indentl (indent), indentr (indent);
	  indentl << " |";
	  indentr << "  ";
	  next_branch.push (Branch_indent (b.rchild, indentr));
	  next_branch.push (Branch_indent (b.lchild, indentl));
	}
      else
	out << indent << "\n";
    }
}

Pair_CFG_parse_stats::Pair_CFG_parse_stats (const Pair_CFG_DP_matrix& matrix)
  : matrix (matrix)
{ }

void Pair_CFG_parse_stats::count_subseqs (const Pair_CFG_parse_tree& parse_tree)
{
  // loop through subseqs in the path
  stack<int> next_branch;
  next_branch.push (0);
  while (!next_branch.empty())
    {
      const Pair_CFG_branch& b = parse_tree[next_branch.top()];
      next_branch.pop();
      int xi = matrix.xenv.find_subseq_idx (b.xl, b.xr - b.xl);
      int yi = matrix.yenv.find_subseq_idx (b.yl, b.yr - b.yl);
      if (xi < 0) THROWEXPR ("Couldn't find x-subseq " << b.xl << ".." << b.xr);
      if (yi < 0) THROWEXPR ("Couldn't find y-subseq " << b.yl << ".." << b.yr);
      for_const_contents (vector<int>, b.path, s)
	{
	  const State_type t = *s >= 0 ? matrix.cfg.state_type[*s] : Null;
	  if (*s == End || t != Null)
	    ++pair_count [Subseq_index_pair (xi, yi)];  // increment count for this subseq
	  xi = matrix.dest_xsubseq_idx (xi, t);
	  yi = matrix.dest_ysubseq_idx (yi, t);
	}
      if (b.has_children())
	{
	  next_branch.push (b.rchild);
	  next_branch.push (b.lchild);
	}
    }
}

Pair_inside_matrix::Pair_inside_matrix (const Named_profile& npx,
					const Named_profile& npy,
					const Fold_envelope& xenv,
					const Fold_envelope& yenv,
					const Pair_CFG_scores& cfg,
					const Pair_envelope& pair_env,
					bool local,
					bool fill_now) :
  Pair_CFG_DP_matrix (npx, npy, xenv, yenv, cfg, pair_env, local)
{
  if (fill_now)
    fill();
}

Pair_inside_matrix::Pair_inside_matrix (const Named_profile& npx,
					const Named_profile& npy,
					const Fold_envelope& xenv,
					const Fold_envelope& yenv,
					const Pair_CFG_scores& cfg,
					bool local,
					bool fill_now) :
  Pair_CFG_DP_matrix (npx, npy, xenv, yenv, cfg, Pair_envelope (npx.dsq.size(), npy.dsq.size(), 1), local)
{
  if (fill_now)
    fill();
}

void Pair_inside_matrix::fill()
{
  Bifcount last_report = 0;
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;
  // cell (STATE,SUBSEQ) is the score for all parse trees of SUBSEQ rooted at STATE
  Subseq_index_array valid_ystartpos;
  valid_ystartpos.reserve (ylen + 1);
  for (int xstartpos = xlen; xstartpos >= 0; --xstartpos)
    if (xenv.by_start[xstartpos].size())
      {
	// prepare list of valid y start positions for this xstartpos
	valid_ystartpos.clear();
	for (int ystartpos = ylen; ystartpos >= 0; --ystartpos)
	  if (startpos_in_pair_envelope (xstartpos, ystartpos))
	    valid_ystartpos.push_back (&yenv.by_start[ystartpos]);
	// loop through xsubseqs with this startpos
	for_const_contents (vector<int>, xenv.by_start[xstartpos], xsi)
	  {
	    const int xsubseq_idx = *xsi;
	    const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
	    const int xendpos = xsubseq.end();
	    for_const_contents (Subseq_index_array, valid_ystartpos, ysp)
	      for_const_contents (vector<int>, **ysp, ysi)
	      {
		const int ysubseq_idx = *ysi;
		const Subseq& ysubseq = yenv.subseq[ysubseq_idx];
		if (endpos_in_pair_envelope (xendpos, ysubseq.end()))
		  {
		    for_const_contents (vector<int>, allowed_in_states (xsubseq, ysubseq), s)
		      {
			// check state type
			const State_type t = cfg.state_type[*s];
			Score sc;
			if (is_bifurc_type(t))
			  {
			    const int l = cfg.bifurc[*s].l;
			    const int r = cfg.bifurc[*s].r;
			    sc = -InfinityScore;
			    if (t == Bifurc)  // Bifurc
			      for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
				for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
				ScorePSumAcc (sc, ScorePMul (read_cell (l, bx->l, by->l), read_cell (r, bx->r, by->r)));
			    else  // BifurcRevY
			      for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
				for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
				ScorePSumAcc (sc, ScorePMul (read_cell (l, bx->l, by->r), read_cell (r, bx->r, by->l)));
			  }
			else
			  {
			    // get the skinny on this cell
			    const int destx_idx = dest_xsubseq_idx (xsubseq_idx, t);
			    const int desty_idx = dest_ysubseq_idx (ysubseq_idx, t);
			    const Subseq& destx = xenv.subseq[destx_idx];
			    const Subseq& desty = yenv.subseq[desty_idx];
			    // do end transitions
			    if (destx.len == 0 && desty.len == 0)
			      sc = cfg.end[*s];
			    else
			      sc = -InfinityScore;
			    // loop over outgoing transitions
			    for_const_contents (vector<int>, allowed_dest_in_states (*s, destx, desty), d)
			      ScorePSumAcc (sc, ScorePMul (read_cell (*d, destx_idx, desty_idx), cfg.transition (*s, *d)));
			    // add emit score
			    if (t != Null)
			      {
				ScorePMulAcc (sc, cfg.emit[*s] [in_emit_idx (t, xsubseq, ysubseq)]);
				ScorePMulAcc (sc, in_meta_sc (*s, xsubseq, ysubseq));
			      }
			  }
			cell (*s, xsubseq_idx, ysubseq_idx) = sc;
			// do start transitions
			if (local_start_allowed (xsubseq, ysubseq))
			  ScorePSumAcc (final_score, ScorePMul (sc, cfg.start[*s]));
		      }
		    // count the cell
		    ++cells_done;
		    bifs_done += xsubseq.bif_in.size() * ysubseq.bif_in.size();
		  }
	      }
	    if (bifs_done - last_report >= CFGDP_REPORT_INTERVAL || (bifs_done == total_bifs && bifs_done > last_report))
	      {
		last_report = bifs_done;
		CTAG(1,CFGDP) << "Inside: finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)total_bifs))/10. << "%)\n";
	      }
	  }
      }
  if (CTAGGING(-1,INSIDE_MATRIX))
    show(CL);
}

void Pair_inside_matrix::traceback_from (int xsubseq_idx, int ysubseq_idx, int traceback_state, vector<int>& path, bool prepend_Start_state) const
{
  // print log message
  CTAG(4,INSIDE_TRACEBACK) << "Starting sampled traceback at xsubseq #" << xsubseq_idx << ", ysubseq #" << ysubseq_idx << ", state #" << traceback_state << "\n";

  // prepend start state
  if (prepend_Start_state)
    {
      path.clear();
      path.push_back (HMM_state_enum::Start);
    }

  // prepare the stacks
  stack<int> xsubseq_idx_stack;
  stack<int> ysubseq_idx_stack;
  stack<int> state_stack;

  // start the traceback
  while (1)
    {
      path.push_back (traceback_state);

      if (traceback_state == HMM_state_enum::End)
	{
	  if (state_stack.empty()) break;

	  xsubseq_idx = xsubseq_idx_stack.top();
	  ysubseq_idx = ysubseq_idx_stack.top();
	  traceback_state = state_stack.top();

	  xsubseq_idx_stack.pop();
	  ysubseq_idx_stack.pop();
	  state_stack.pop();

	  continue;
	}

      const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
      const Subseq& ysubseq = yenv.subseq[ysubseq_idx];

      const Score traceback_score = read_cell (traceback_state, xsubseq_idx, ysubseq_idx);
      const State_type traceback_type = cfg.state_type [traceback_state];

      const Prob tprob = Score2Prob (traceback_score);
      Prob p = Rnd::prob() * tprob;
      Prob q = 0;
      
      if (is_bifurc_type (traceback_type))
	{
	  const int l = cfg.bifurc[traceback_state].l;
	  const int r = cfg.bifurc[traceback_state].r;
	  while (q < p)
	    {
	      if (traceback_type == Bifurc)  // Bifurc
		{
		  for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
		    {
		      if (q >= p) break;
		      for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
			if ((q += Score2Prob (ScorePMul (read_cell (l, bx->l, by->l), read_cell (r, bx->r, by->r)))) >= p)
			  {
			    xsubseq_idx_stack.push (bx->r);
			    ysubseq_idx_stack.push (by->r);
			    xsubseq_idx = bx->l;
			    ysubseq_idx = by->l;
			    break;
			  }
		    }
		}
	      else  // BifurcRevY
		for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
		{
		  if (q >= p) break;
		  for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
		    if ((q += Score2Prob (ScorePMul (read_cell (l, bx->l, by->r), read_cell (r, bx->r, by->l)))) >= p)
		      {
			xsubseq_idx_stack.push (bx->r);
			ysubseq_idx_stack.push (by->l);
			xsubseq_idx = bx->l;
			ysubseq_idx = by->r;
			break;
		      }
		}
	      if (q < p)
		CLOGERR << "Pair CFG: sampled traceback failed (couldn't find bifurcation) at xi=" << xsubseq_idx << ", yi=" << ysubseq_idx << ", state=" << traceback_state << "; p=" << p << ", q=" << q << ", tprob=" << tprob << ", tscore=" << traceback_score << "; trying again\n";
	    }
	  state_stack.push (r);
	  traceback_state = l;
	}
      else
	{
	  // get the skinny on this cell
	  const int destx_idx = dest_xsubseq_idx (xsubseq_idx, traceback_type);
	  const int desty_idx = dest_ysubseq_idx (ysubseq_idx, traceback_type);
	  const Subseq& destx = xenv.subseq[destx_idx];
	  const Subseq& desty = yenv.subseq[desty_idx];
	  // calculate emit score
	  const Score emit_sc =
	    traceback_type == Null
	    ? 0
	    : ScorePMul (cfg.emit[traceback_state] [in_emit_idx (traceback_type, xsubseq, ysubseq)],
			 in_meta_sc (traceback_state, xsubseq, ysubseq));
	  while (q < p)
	    {
	      // do end transitions
	      if (destx.len == 0 && desty.len == 0 ? (q += Score2Prob (ScorePMul (cfg.end[traceback_state], emit_sc))) >= p : FALSE)
		traceback_state = HMM_state_enum::End;
	      else
		{
		  // loop over outgoing transitions
		  for_const_contents (vector<int>, allowed_dest_in_states (traceback_state, destx, desty), d)
		    if ((q += Score2Prob (ScorePMul3 (read_cell (*d, destx_idx, desty_idx), cfg.transition (traceback_state, *d), emit_sc))) >= p)
		      {
			traceback_state = *d;
			break;
		      }
		}
	      if (q < p)
		CLOGERR << "Pair CFG: sampled traceback failed at xi=" << xsubseq_idx << ", yi=" << ysubseq_idx << ", state=" << traceback_state << "; p=" << p << ", q=" << q << ", tprob=" << tprob << ", tscore=" << traceback_score << "; trying again\n";
	    }
	  xsubseq_idx = destx_idx;
	  ysubseq_idx = desty_idx;
	}
    }
  if (CTAGGING(3,INSIDE_TRACEBACK))
    CL << "Sampled traceback: (" << path << ")\n";
}

Pair_inside_cell_sorter::Pair_inside_cell_sorter (const Pair_inside_matrix& inside)
  : inside (inside)
{
  // print log message
  CTAG(5,INSIDE_TRACEBACK) << "Sorting cells by probability\n";
  // copy some types
  typedef vector<Subseq>::const_iterator Subseq_iterator;
  // copy some variables
  const int xsubseqs = inside.xsubseqs;
  const int ysubseqs = inside.ysubseqs;
  const Subseq_iterator xenv_begin = inside.xenv.subseq.begin();
  const Subseq_iterator yenv_begin = inside.yenv.subseq.begin();
  const Pair_CFG_scores& cfg = inside.cfg;
  // find co-ord range (trivial if nonlocal)
  int minx = 0, miny = 0;
  if (!inside.local)
    {
      minx = inside.xenv.first_start_subseq();
      miny = inside.yenv.first_start_subseq();
    }
  // collect start scores for each cell
  deque<Cell_coords> unsorted_cells;
  deque<Prob> unsorted_cells_prob;
  for (int xsubseq_idx = minx; xsubseq_idx < xsubseqs; ++xsubseq_idx)
    for (int ysubseq_idx = miny; ysubseq_idx < ysubseqs; ++ysubseq_idx)
      {
	// find state probability
	const Subseq& xsubseq = xenv_begin[xsubseq_idx];
	const Subseq& ysubseq = yenv_begin[ysubseq_idx];
	if (inside.in_pair_envelope (xsubseq, ysubseq))
	  for_const_contents (vector<int>, inside.allowed_in_states (xsubseq, ysubseq), s)
	  {
	    const Score sc = ScorePMul (inside.read_cell (*s, xsubseq_idx, ysubseq_idx), cfg.start[*s]);
	    if (sc > -InfinityScore)
	      {
		unsorted_cells.push_back (Cell_coords (*s, xsubseq_idx, ysubseq_idx));
		unsorted_cells_prob.push_back (Score2Prob (sc));
	      }
	  }
      }
  // sort cells
  vector<Prob> unsorted_cells_prob_vec (unsorted_cells_prob.begin(), unsorted_cells_prob.end());
  unsorted_cells_prob = deque<Prob>();  // free up memory
  Reverse_Schwartzian<Prob> by_prob (unsorted_cells_prob_vec);
  const vector<int> cell_order = by_prob.sorted();
  // write sorted_cells and sorted_cells_prob vectors
  const int N = unsorted_cells.size();
  sorted_cells = vector<Cell_coords> (N, Cell_coords(-1,-1,-1));
  sorted_cells_prob = vector<Prob> (N);
  total_prob = 0.;
  for (int i = 0; i < N; ++i)
    {
      sorted_cells[i] = unsorted_cells[cell_order[i]];
      total_prob += (sorted_cells_prob[i] = unsorted_cells_prob_vec[cell_order[i]]);
    }
}

Pair_CFG_DP_matrix::Cell_coords Pair_inside_cell_sorter::sample_coords() const
{
  // set co-ords to null values
  Cell_coords coords (HMM_state_enum::UndefinedState, -1, -1);

  // prepare initial coords
  if (inside.final_score == -InfinityScore)
    THROWEXPR ("Final score is -infinity; traceback likely to break");

  if (inside.xlen == 0 && inside.ylen == 0 && inside.final_score == inside.cfg.start_to_end())
    {
      coords.state = HMM_state_enum::End;
      coords.xsubseq = inside.xsubseqs - 1;
      coords.ysubseq = inside.ysubseqs - 1;
    }
  else
    {
      // sample the final state
      Prob p = Rnd::prob() * total_prob;
      Prob q = 0;
      while (q < p)
	{
	  for (int i = 0; i < (int) sorted_cells.size(); ++i)
	    {
	      if ((q += sorted_cells_prob[i]) >= p)
		{
		  coords = sorted_cells[i];
		  break;
		}
	    }
	  if (q >= p) break;
	  CLOGERR << "Warning -- failed to sample a DP cell; p=" << p << ", q=" << q << ", total_prob=" << total_prob << "; trying again\n";
	}
    }
  return coords;
}

Pair_CFG_local_path Pair_inside_cell_sorter::traceback_with_coords() const
{
  // create local path object
  Pair_CFG_local_path local_path;

  // sample start co-ords
  const Cell_coords coords = sample_coords();
  const int xsubseq_idx = coords.xsubseq;
  const int ysubseq_idx = coords.ysubseq;
  const int traceback_state = coords.state;

  // set local path coords
  local_path.xstart = inside.xenv.subseq[xsubseq_idx].start;
  local_path.xlen   = inside.xenv.subseq[xsubseq_idx].len;
  local_path.ystart = inside.yenv.subseq[ysubseq_idx].start;
  local_path.ylen   = inside.yenv.subseq[ysubseq_idx].len;

  // do traceback, set local path
  inside.traceback_from (xsubseq_idx, ysubseq_idx, traceback_state, local_path.path);

  // return path
  return local_path;
}

Pair_CFG_parse_tree Pair_inside_cell_sorter::parse_tree() const
{
  const Pair_CFG_local_path lp = traceback_with_coords();
  return inside.cfg.parse (lp);
}

Pair_CFG_alignment Pair_inside_cell_sorter::alignment() const
{
  const Pair_CFG_parse_tree parse = parse_tree();
  Pair_CFG_alignment align = parse.alignment (inside.cfg.state_type, inside.npx, inside.npy);
  align.score = inside.cfg.path_score (parse, inside.xdsq, inside.ydsq);
  return align;
}

Pair_outside_matrix::Pair_outside_matrix (const Pair_inside_matrix& inside, bool fill_now) :
  Pair_CFG_DP_matrix (inside.npx, inside.npy, inside.xenv, inside.yenv, inside.cfg, inside.pair_env, inside.local),
  inside (inside),
  count (inside.cfg),
  xmetacount (inside.cfg.max_xmeta_idx() + 1, Metaprob (inside.npx.size(), 0.0)),
  ymetacount (inside.cfg.max_ymeta_idx() + 1, Metaprob (inside.npy.size(), 0.0))
{
  _show_out = 1;
  if (fill_now)
    fill();
}

void Pair_outside_matrix::fill()
{
  // check that final score of inside algorithm was finite
  const Score inside_final_sc = inside.final_score;
  if (inside_final_sc <= -InfinityScore)
    {
      CLOGERR << "WARNING: inside algorithm returned score of -infinity; skipping outside algorithm\n";
      return;
    }
  final_score = inside_final_sc;
  count.log_likelihood = Score2Nats (inside_final_sc);
  // do the DP
  Bifcount last_report = 0;
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;
  // meaning of a cell is different for outside matrix than for inside matrix
  // for outside, cell (STATE,SUBSEQ) is the score for the parse tree landing in STATE
  // having emitted all but SUBSEQ (including the emission from STATE)
  // thus the emissions are outside of SUBSEQ (c.f. inside matrix where emissions are inside)
  if (xlen == 0 && ylen == 0)
    count.start_to_end() = Score2Prob (ScorePMul (cfg.start_to_end(), -inside_final_sc));
  Subseq_index_array valid_ystartpos;
  valid_ystartpos.reserve (ylen + 1);
  for (int xstartpos = 0; xstartpos <= xlen; ++xstartpos)
    if (xenv.by_start[xstartpos].size())
      {
	// prepare list of valid y start positions for this xstartpos
	valid_ystartpos.clear();
	for (int ystartpos = 0; ystartpos <= ylen; ++ystartpos)
	  if (startpos_in_pair_envelope (xstartpos, ystartpos))
	    valid_ystartpos.push_back (&yenv.by_start[ystartpos]);
	// loop through xsubseqs with this startpos
	for_const_reverse_contents (vector<int>, xenv.by_start[xstartpos], xsi)
	  {
	    const int xsubseq_idx = *xsi;
	    const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
	    const int xendpos = xsubseq.end();
	    for_const_contents (Subseq_index_array, valid_ystartpos, ysp)
	      for_const_reverse_contents (vector<int>, **ysp, ysi)
	      {
		const int ysubseq_idx = *ysi;
		const Subseq& ysubseq = yenv.subseq[ysubseq_idx];
		if (endpos_in_pair_envelope (xendpos, ysubseq.end()))
		  {
		    for_const_contents (vector<int>, allowed_out_states (xsubseq, ysubseq), d)
		      {
			// get state type
			const State_type t = cfg.state_type[*d];
			// find source subseqs
			const int srcx_idx = src_xsubseq_idx (xsubseq_idx, t);
			const int srcy_idx = src_ysubseq_idx (ysubseq_idx, t);
			const Subseq& srcx = xenv.subseq[srcx_idx];
			const Subseq& srcy = yenv.subseq[srcy_idx];
			// check that source is in pair envelope
			if (!in_pair_envelope (srcx, srcy))
			  continue;
			// get inside score
			const Score inside_minus_final_sc = ScorePMul (inside.read_cell (*d, srcx_idx, srcy_idx), -inside_final_sc);
			// do start transitions
			Score sc;
			if (local_start_allowed (srcx, srcy))
			  {
			    sc = cfg.start[*d];  // NB must assign sc first, as next line uses new value
			    count.start[*d] += Score2Prob (ScorePMul (sc, inside_minus_final_sc));
			  }
			else
			  sc = -InfinityScore;
			// loop over incoming transitions
			for_const_contents (vector<int>, allowed_src_out_states (*d, srcx, srcy), s)
			  {
			    const Score incoming_sc = ScorePMul (read_cell (*s, srcx_idx, srcy_idx), cfg.transition (*s, *d));
			    count.transition (*s, *d) += Score2Prob (ScorePMul (incoming_sc, inside_minus_final_sc));
			    ScorePSumAcc (sc, incoming_sc);
			  }
			// loop over incoming bifurcations to left & right
			for_const_contents (vector<Bifurcation_left_parent>, filter.left_parent[*d], b)
			  if (cfg.state_type[b->p] == Bifurc)  // Bifurc
			    for_const_contents (Subseq::Bifurc_outl_pseudovec, xsubseq.bif_out_l, bx)
			      for_const_contents (Subseq::Bifurc_outl_pseudovec, ysubseq.bif_out_l, by)
			      ScorePSumAcc (sc, ScorePMul (read_cell (b->p, bx->out, by->out), inside.read_cell (b->l, bx->l, by->l)));
			  else  // BifurcRevY
			    for_const_contents (Subseq::Bifurc_outl_pseudovec, xsubseq.bif_out_l, bx)
			      for_const_contents (Subseq::Bifurc_outr_pseudovec, ysubseq.bif_out_r, by)
			      ScorePSumAcc (sc, ScorePMul (read_cell (b->p, bx->out, by->out), inside.read_cell (b->l, bx->l, by->r)));
			for_const_contents (vector<Bifurcation_right_parent>, filter.right_parent[*d], b)
			  if (cfg.state_type[b->p] == Bifurc)  // Bifurc
			    for_const_contents (Subseq::Bifurc_outr_pseudovec, xsubseq.bif_out_r, bx)
			      for_const_contents (Subseq::Bifurc_outr_pseudovec, ysubseq.bif_out_r, by)
			      ScorePSumAcc (sc, ScorePMul (read_cell (b->p, bx->out, by->out), inside.read_cell (b->r, bx->r, by->r)));
			  else  // BifurcRevY
			    for_const_contents (Subseq::Bifurc_outr_pseudovec, xsubseq.bif_out_r, bx)
			      for_const_contents (Subseq::Bifurc_outl_pseudovec, ysubseq.bif_out_l, by)
			      ScorePSumAcc (sc, ScorePMul (read_cell (b->p, bx->out, by->out), inside.read_cell (b->r, bx->r, by->l)));
			// add emit score
			if (is_emit_type(t))
			  {
			    const int emit_idx = out_emit_idx (t, xsubseq, ysubseq);
			    const Prob post_prob = Score2Prob (ScorePMul (sc, inside_minus_final_sc));
			    count.emit[*d][emit_idx] += post_prob;
			    out_count_meta (post_prob, *d, xsubseq, ysubseq);
			    // NB must update sc after updating counts, as post_prob uses old value
			    ScorePMulAcc (sc, cfg.emit[*d][emit_idx]);
			    ScorePMulAcc (sc, out_meta_sc (*d, xsubseq, ysubseq));
			  }
			// store
			if (sc < -InfinityScore) THROWEXPR("Overflow");
			cell (*d, xsubseq_idx, ysubseq_idx) = sc;
			// do end transitions
			if (xsubseq.len == 0 && ysubseq.len == 0)
			  count.end[*d] += Score2Prob (ScorePMul3 (sc, cfg.end[*d], -inside_final_sc));
		      }
		    // count the cell
		    ++cells_done;
		    bifs_done += xsubseq.bif_in.size() * ysubseq.bif_in.size();
		  }
	      }
	    if (bifs_done - last_report >= CFGDP_REPORT_INTERVAL || (bifs_done == total_bifs && bifs_done > last_report))
	      {
		last_report = bifs_done;
		CTAG(1,CFGDP) << "Outside: finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)total_bifs))/10. << "%)\n";
	      }
	  }
      }
  if (CTAGGING(-1,OUTSIDE_MATRIX))
    show(CL);
  if (CTAGGING(2,CFG_COUNTS))
    count.show (CL);
}

Score Pair_inside_outside_matrix::post_transition_sc (int src_state, int dest_state, int xsubseq_idx, int ysubseq_idx) const
{
  // make copies of some variables
  const Pair_CFG_scores& cfg = inside.cfg;
  const Subseq& xsubseq = inside.xenv.subseq[xsubseq_idx];
  const Subseq& ysubseq = inside.yenv.subseq[ysubseq_idx];
  // calculate P(inside) / P(final)
  Score inside_minus_final_sc = -inside.final_score;
  if (dest_state != HMM_state_enum::End)
    ScorePMulAcc (inside_minus_final_sc, inside.read_cell (dest_state, xsubseq_idx, ysubseq_idx));
  else  // dest_state == End
    if (xsubseq.len > 0 || ysubseq.len > 0)
      return -InfinityScore;
  // calculate P(outside) * P(transition)
  Score incoming_sc = cfg.transition (src_state, dest_state);
  if (src_state != HMM_state_enum::Start)
    ScorePMulAcc (incoming_sc, outside.read_cell (src_state, xsubseq_idx, ysubseq_idx));
  else  // src_state == Start
    if (!inside.local_start_allowed (xsubseq, ysubseq))
      return -InfinityScore;
  // multiply the probabilities and return
  return ScorePMul (incoming_sc, inside_minus_final_sc);
}

Score Pair_inside_outside_matrix::post_state_sc (int dest_state, int xsubseq_idx, int ysubseq_idx) const
{
  // get subseqs
  const Subseq& xsubseq = inside.xenv.subseq[xsubseq_idx];
  const Subseq& ysubseq = inside.yenv.subseq[ysubseq_idx];
  // sum posterior probabilities of incoming transitions
  Score sc = post_transition_sc (HMM_state_enum::Start, dest_state, xsubseq_idx, ysubseq_idx);
  for_const_contents (vector<int>, inside.allowed_src_out_states (dest_state, xsubseq, ysubseq), s)
    ScorePSumAcc (sc, post_transition_sc (*s, dest_state, xsubseq_idx, ysubseq_idx));
  return sc;
}

Pair_inside_outside_matrix::Pair_inside_outside_matrix (const Named_profile& npx,
							const Named_profile& npy,
							const Fold_envelope& xenv,
							const Fold_envelope& yenv,
							const Pair_CFG_scores& cfg,
							bool local,
							bool fill_now) :
  inside (npx, npy, xenv, yenv, cfg, local, fill_now),
  outside (inside, fill_now),
  count (outside.count),
  xmetacount (outside.xmetacount),
  ymetacount (outside.ymetacount)
{ }

Pair_inside_outside_matrix::Pair_inside_outside_matrix (const Named_profile& npx,
							const Named_profile& npy,
							const Fold_envelope& xenv,
							const Fold_envelope& yenv,
							const Pair_CFG_scores& cfg,
							const Pair_envelope& pair_env,
							bool local,
							bool fill_now) :
  inside (npx, npy, xenv, yenv, cfg, pair_env, local, fill_now),
  outside (inside, fill_now),
  count (outside.count),
  xmetacount (outside.xmetacount),
  ymetacount (outside.ymetacount)
{ }

void Pair_inside_outside_matrix::show (ostream& o) const
{
  Biosequence xseq, yseq;
  inside.cfg.alphabet().dsq2seq (inside.xdsq, xseq);
  inside.cfg.alphabet().dsq2seq (inside.ydsq, yseq);
  o << "Sequence X: " << xseq << "\n";
  o << "Sequence Y: " << yseq << "\n";
  o << "Envelope X:\n";
  inside.xenv.dump(o);
  o << "Envelope Y:\n";
  inside.yenv.dump(o);
  o << "Inside matrix:\n";
  inside.show(o);
  o << "Outside matrix:\n";
  outside.show(o);
  o << "Counts:\n";
  count.show(o);
}


Pair_CYK_matrix::Pair_CYK_matrix (const Named_profile& npx,
				  const Named_profile& npy,
				  const Fold_envelope& xenv,
				  const Fold_envelope& yenv,
				  const Pair_CFG_scores& cfg,
				  const Pair_envelope& pair_env,
				  bool local,
				  bool fill_now) :
  Pair_CFG_DP_matrix (npx, npy, xenv, yenv, cfg, pair_env, local),
  final_xsubseq_idx (-1),
  final_ysubseq_idx (-1),
  final_state (HMM_state_enum::UndefinedState)
{
  if (fill_now)
    fill();
}

Pair_CYK_matrix::Pair_CYK_matrix (const Named_profile& npx,
				  const Named_profile& npy,
				  const Fold_envelope& xenv,
				  const Fold_envelope& yenv,
				  const Pair_CFG_scores& cfg,
				  bool local,
				  bool fill_now) :
  Pair_CFG_DP_matrix (npx, npy, xenv, yenv, cfg, Pair_envelope (npx.dsq.size(), npy.dsq.size(), 1), local),
  final_xsubseq_idx (-1),
  final_ysubseq_idx (-1),
  final_state (HMM_state_enum::UndefinedState)
{
  if (fill_now)
    fill();
}

void Pair_CYK_matrix::fill()
{
  Bifcount last_report = 0;
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;

  // handle the special case of zero-length sequences
  if (xlen == 0 && ylen == 0)
    {
      final_state = HMM_state_enum::End;
      final_xsubseq_idx = 0;
      final_ysubseq_idx = 0;
    }

  // cell (STATE,SUBSEQ) is the score for the best parse tree of SUBSEQ rooted at STATE
  Subseq_index_array valid_ystartpos;
  valid_ystartpos.reserve (ylen + 1);
  for (int xstartpos = xlen; xstartpos >= 0; --xstartpos)
    if (xenv.by_start[xstartpos].size())
      {
	// prepare list of valid y start positions for this xstartpos
	valid_ystartpos.clear();
	for (int ystartpos = ylen; ystartpos >= 0; --ystartpos)
	  if (startpos_in_pair_envelope (xstartpos, ystartpos))
	    valid_ystartpos.push_back (&yenv.by_start[ystartpos]);
	// loop through xsubseqs with this startpos
	for_const_contents (vector<int>, xenv.by_start[xstartpos], xsi)
	  {
	    const int xsubseq_idx = *xsi;
	    const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
	    const int xendpos = xsubseq.end();
	    for_const_contents (Subseq_index_array, valid_ystartpos, ysp)
	      for_const_contents (vector<int>, **ysp, ysi)
	      {
		const int ysubseq_idx = *ysi;
		const Subseq& ysubseq = yenv.subseq[ysubseq_idx];
		if (endpos_in_pair_envelope (xendpos, ysubseq.end()))
		  {
		    for_const_contents (vector<int>, allowed_in_states (xsubseq, ysubseq), s)
		      {
			// check state type
			const State_type t = cfg.state_type[*s];
			Score sc;
			if (is_bifurc_type(t))
			  {
			    const int l = cfg.bifurc[*s].l;
			    const int r = cfg.bifurc[*s].r;
			    sc = -InfinityScore;
			    if (t == Bifurc)  // Bifurc
			      for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
				for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
				sc = max (sc, ScorePMul (read_cell (l, bx->l, by->l), read_cell (r, bx->r, by->r)));
			    else  // BifurcRevY
			      for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
				for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
				sc = max (sc, ScorePMul (read_cell (l, bx->l, by->r), read_cell (r, bx->r, by->l)));
			  }
			else
			  {
			    // get the skinny on this cell
			    const int destx_idx = dest_xsubseq_idx (xsubseq_idx, t);
			    const int desty_idx = dest_ysubseq_idx (ysubseq_idx, t);
			    const Subseq& destx = xenv.subseq[destx_idx];
			    const Subseq& desty = yenv.subseq[desty_idx];
			    // do end transitions
			    if (destx.len == 0 && desty.len == 0)
			      sc = cfg.end[*s];
			    else
			      sc = -InfinityScore;
			    // loop over outgoing transitions
			    for_const_contents (vector<int>, allowed_dest_in_states (*s, destx, desty), d)
			      sc = max (sc, ScorePMul (read_cell (*d, destx_idx, desty_idx), cfg.transition (*s, *d)));
			    // add emit score
			    if (t != Null)
			      {
				ScorePMulAcc (sc, cfg.emit[*s] [in_emit_idx (t, xsubseq, ysubseq)]);
				ScorePMulAcc (sc, in_meta_sc (*s, xsubseq, ysubseq));
			      }
			  }
			cell (*s, xsubseq_idx, ysubseq_idx) = sc;
			if (local_start_allowed (xsubseq, ysubseq))
			  {
			    const Score f_sc = ScorePMul (sc, cfg.start[*s]);
			    if (f_sc > final_score)
			      {
				final_score = f_sc;
				final_xsubseq_idx = xsubseq_idx;
				final_ysubseq_idx = ysubseq_idx;
				final_state = *s;
			      }
			  }
		      }
		    // count the cell
		    ++cells_done;
		    bifs_done += xsubseq.bif_in.size() * ysubseq.bif_in.size();
		  }
	      }
	    if (bifs_done - last_report >= CFGDP_REPORT_INTERVAL || (bifs_done == total_bifs && bifs_done > last_report))
	      {
		last_report = bifs_done;
		CTAG(1,CFGDP) << "CYK: finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)total_bifs))/10. << "%)\n";
	      }
	  }
      }
  if (CTAGGING(-1,CYK_MATRIX))
    show(CL);
}

vector<int> Pair_CYK_matrix::traceback() const
{
  vector<int> path;
  traceback_from (final_xsubseq_idx, final_ysubseq_idx, final_state, path, TRUE);
  return path;
}

void Pair_CYK_matrix::traceback_from (int xsubseq_idx, int ysubseq_idx, int traceback_state, vector<int>& path, bool prepend_Start_state) const
{
  if (traceback_state == HMM_state_enum::UndefinedState)
    THROWEXPR (cfg_dump() << "Traceback start co-ords undefined");

  if (final_score == -InfinityScore)
    THROWEXPR ("Final score is -infinity; traceback likely to break");

  // get invariants
  const Subseq& xsubseq_init = xenv.subseq[xsubseq_idx];
  const Subseq& ysubseq_init = yenv.subseq[ysubseq_idx];
  const int xsubseq_init_len = xsubseq_init.len;
  const int ysubseq_init_len = ysubseq_init.len;

  // print log message
  CTAG(4,CYK_TRACEBACK)
    << "Starting CYK traceback at xsubseq #" << xsubseq_idx << ", ysubseq #" << ysubseq_idx << ", state #" << traceback_state << "\n"
    << "xsubseq=(" << xsubseq_init.start << '+' << xsubseq_init.len << ") "
    << "ysubseq=(" << ysubseq_init.start << '+' << ysubseq_init.len << ")\n";

  // prepend start state
  path.clear();
  if (prepend_Start_state)
    path.push_back (HMM_state_enum::Start);

  vector<int> xsubseq_idx_stack;
  vector<int> ysubseq_idx_stack;
  vector<int> state_stack;

  int x_emitted = 0, y_emitted = 0;  // total emitted

  while (1)
    {
      // push this state onto the path
      path.push_back (traceback_state);

      // check if End state reached
      if (traceback_state == HMM_state_enum::End)
	{
	  if (state_stack.size() == 0) break;

	  xsubseq_idx = xsubseq_idx_stack.back();
	  ysubseq_idx = ysubseq_idx_stack.back();
	  traceback_state = state_stack.back();

	  xsubseq_idx_stack.pop_back();
	  ysubseq_idx_stack.pop_back();
	  state_stack.pop_back();

	  continue;
	}

      // get subseqs
      const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
      const Subseq& ysubseq = yenv.subseq[ysubseq_idx];

      // get score & type
      const Score traceback_score = read_cell (traceback_state, xsubseq_idx, ysubseq_idx);
      const State_type traceback_type = cfg.state_type [traceback_state];

      // print absurdly detailed log message
      if (CTAGGING(-1,CYK_TRACEBACK_DETAIL))
	{
	  // Compute total length emitted, stacked & yet-to-emit: should be invariant
	  int xtot = xsubseq.len + x_emitted;
	  int ytot = ysubseq.len + y_emitted;
	  vector<int> xstacked, ystacked;
	  for_const_contents (vector<int>, xsubseq_idx_stack, xsi)
	    {
	      const int xl = xenv.subseq[*xsi].len;
	      xstacked.push_back (xl);
	      xtot += xl;
	    }
	  for_const_contents (vector<int>, ysubseq_idx_stack, ysi)
	    {
	      const int yl = yenv.subseq[*ysi].len;
	      ystacked.push_back (yl);
	      ytot += yl;
	    }
	  // print
	  CL << "In CYK traceback: state=" << traceback_state << " type=" << traceback_type
	     << " xsubseq=(" << xsubseq.start << '+' << xsubseq.len << ") ysubseq=(" << ysubseq.start << '+' << ysubseq.len << ")\n"
	     << "x: emitted (" << x_emitted << ") + stacked (" << xstacked << ") + pending (" << xsubseq.len << ") = " << xtot << "; should be " << xsubseq_init_len << "\n"
	     << "y: emitted (" << y_emitted << ") + stacked (" << ystacked << ") + pending (" << ysubseq.len << ") = " << ytot << "; should be " << ysubseq_init_len << "\n";
	}

      // update x_emitted, y_emitted
      if (traceback_type & EmitXL) ++x_emitted;
      if (traceback_type & EmitXR) ++x_emitted;
      if (traceback_type & EmitYL) ++y_emitted;
      if (traceback_type & EmitYR) ++y_emitted;

      // check for bifurcation
      if (is_bifurc_type (traceback_type))
	{
	  const int l = cfg.bifurc[traceback_state].l;
	  const int r = cfg.bifurc[traceback_state].r;
	  if (traceback_type == Bifurc)  // Bifurc
	    {
	      for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
		for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
		if (traceback_score == ScorePMul (read_cell (l, bx->l, by->l), read_cell (r, bx->r, by->r)))
		  {
		    xsubseq_idx_stack.push_back (bx->r);
		    ysubseq_idx_stack.push_back (by->r);
		    xsubseq_idx = bx->l;
		    ysubseq_idx = by->l;
		    goto CYK_FoundBifurcation;
		  }
	    }
	  else  // BifurcRevY
	    for_const_contents (Subseq::Bifurc_in_pseudovec, xsubseq.bif_in, bx)
	      for_const_contents (Subseq::Bifurc_in_pseudovec, ysubseq.bif_in, by)
	      if (traceback_score == ScorePMul (read_cell (l, bx->l, by->r), read_cell (r, bx->r, by->l)))
		{
		  xsubseq_idx_stack.push_back (bx->r);
		  ysubseq_idx_stack.push_back (by->l);
		  xsubseq_idx = bx->l;
		  ysubseq_idx = by->r;
		  goto CYK_FoundBifurcation;
		}
	  THROWEXPR (cfg_dump() << "Pair CFG: traceback failed (couldn't find bifurcation) at xi=" << xsubseq_idx << ", yi=" << ysubseq_idx << ", state=" << traceback_state);
	CYK_FoundBifurcation:
	  state_stack.push_back (r);
	  traceback_state = l;
	}
      else
	{
	  // get the skinny on this cell
	  const int destx_idx = dest_xsubseq_idx (xsubseq_idx, traceback_type);
	  const int desty_idx = dest_ysubseq_idx (ysubseq_idx, traceback_type);
	  const Subseq& destx = xenv.subseq[destx_idx];
	  const Subseq& desty = yenv.subseq[desty_idx];
	  // calculate emit score
	  const Score emit_sc =
	    traceback_type == Null
	    ? 0
	    : ScorePMul (cfg.emit[traceback_state] [in_emit_idx (traceback_type, xsubseq, ysubseq)],
			 in_meta_sc (traceback_state, xsubseq, ysubseq));
	  // do end transitions
	  if (destx.len == 0 && desty.len == 0 && traceback_score == ScorePMul (cfg.end[traceback_state], emit_sc))
	    {
	      traceback_state = HMM_state_enum::End;
	      goto CYK_FoundEmit;
	    }
	  else
	    {
	      // loop over outgoing transitions
	      for_const_contents (vector<int>, allowed_dest_in_states (traceback_state, destx, desty), d)
		if (traceback_score == ScorePMul3 (read_cell (*d, destx_idx, desty_idx), cfg.transition (traceback_state, *d), emit_sc))
		  {
		    traceback_state = *d;
		    goto CYK_FoundEmit;
		  }
	      THROWEXPR (cfg_dump() << "Pair CFG: traceback failed at xi=" << xsubseq_idx << ", yi=" << ysubseq_idx << ", state=" << traceback_state);
	    }
	CYK_FoundEmit:
	  xsubseq_idx = destx_idx;
	  ysubseq_idx = desty_idx;
	}
    }
  // print log messages
  if (CTAGGING(3,CYK_TRACEBACK))
    CL << "Finished CYK traceback: path=(" << path << ")\n"
       << "x_emitted=" << x_emitted << " y_emitted=" << y_emitted << "\n";
}

Pair_CFG_local_path Pair_CYK_matrix::traceback_with_coords() const
{
  Pair_CFG_local_path local_path;
  local_path.xstart = xenv.subseq[final_xsubseq_idx].start;
  local_path.xlen   = xenv.subseq[final_xsubseq_idx].len;
  local_path.ystart = yenv.subseq[final_ysubseq_idx].start;
  local_path.ylen   = yenv.subseq[final_ysubseq_idx].len;
  local_path.path   = traceback();
  return local_path;
}

Pair_CFG_parse_tree Pair_CYK_matrix::parse_tree() const
{
  const Pair_CFG_local_path lp = traceback_with_coords();
  const Pair_CFG_parse_tree parse = cfg.parse (lp);
  if (CTAGGING(1,CYK_PARSE_TREE))
    {
      CL << "CYK parse tree:\n";
      parse.show (CL, &cfg.state_type);
    }
  if (CTAGGING(-2,CYK_EXPANDED_TRACE))
    {
      CL << "Expanded CYK parse tree:\n";
      show_expanded_trace (parse, CL);
    }
  if (CTAGGING(1,CYK_BREAKDOWN))
    {
      CL << "CYK parse tree score breakdown:\n";
      cfg.show_score_breakdown (parse, xdsq, ydsq, CL);
    }
  return parse;
}

Pair_CFG_alignment Pair_CYK_matrix::alignment() const
{
  const Pair_CFG_parse_tree parse = parse_tree();
  Pair_CFG_alignment align = parse.alignment (cfg.state_type, npx, npy);
  align.score = final_score;
  return align;
}

Pair_KYC_matrix::Pair_KYC_matrix (const Pair_CYK_matrix& cyk, bool fill_now) :
  Pair_CFG_DP_matrix (cyk.npx, cyk.npy, cyk.xenv, cyk.yenv, cyk.cfg, cyk.pair_env, cyk.local),
  cyk (cyk)
{
  _show_out = 1;
  if (fill_now)
    fill();
}

void Pair_KYC_matrix::fill()
{
  // check that final score of inside algorithm was finite
  const Score cyk_final_sc = cyk.final_score;
  if (cyk_final_sc <= -InfinityScore)
    {
      CLOGERR << "WARNING: CYK algorithm returned score of -infinity; skipping KYC algorithm\n";
      return;
    }
  final_score = cyk_final_sc;
  // do the DP
  Bifcount last_report = 0;
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;
  // meaning of a cell is different for KYC matrix than for CYK matrix
  // for KYC, cell (STATE,SUBSEQ) is the max score for any parse tree landing in STATE
  // having emitted all but SUBSEQ (including the emission from STATE)
  // thus the emissions are outside of SUBSEQ (c.f. CYK matrix where emissions are inside)
  Subseq_index_array valid_ystartpos;
  valid_ystartpos.reserve (ylen + 1);
  for (int xstartpos = 0; xstartpos <= xlen; ++xstartpos)
    if (xenv.by_start[xstartpos].size())
      {
	// prepare list of valid y start positions for this xstartpos
	valid_ystartpos.clear();
	for (int ystartpos = 0; ystartpos <= ylen; ++ystartpos)
	  if (startpos_in_pair_envelope (xstartpos, ystartpos))
	    valid_ystartpos.push_back (&yenv.by_start[ystartpos]);
	// loop through xsubseqs with this startpos
	for_const_reverse_contents (vector<int>, xenv.by_start[xstartpos], xsi)
	  {
	    const int xsubseq_idx = *xsi;
	    const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
	    const int xendpos = xsubseq.end();
	    for_const_contents (Subseq_index_array, valid_ystartpos, ysp)
	      for_const_reverse_contents (vector<int>, **ysp, ysi)
	      {
		const int ysubseq_idx = *ysi;
		const Subseq& ysubseq = yenv.subseq[ysubseq_idx];
		if (endpos_in_pair_envelope (xendpos, ysubseq.end()))
		  {
		    for_const_contents (vector<int>, allowed_out_states (xsubseq, ysubseq), d)
		      {
			// get state type
			const State_type t = cfg.state_type[*d];
			// find source subseqs
			const int srcx_idx = src_xsubseq_idx (xsubseq_idx, t);
			const int srcy_idx = src_ysubseq_idx (ysubseq_idx, t);
			const Subseq& srcx = xenv.subseq[srcx_idx];
			const Subseq& srcy = yenv.subseq[srcy_idx];
			// check that source is in pair envelope
			if (!in_pair_envelope (srcx, srcy))
			  continue;
			// do start transitions
			Score sc;
			if (local_start_allowed (srcx, srcy))
			  sc = cfg.start[*d];
			else
			  sc = -InfinityScore;
			// loop over incoming transitions
			for_const_contents (vector<int>, allowed_src_out_states (*d, srcx, srcy), s)
			  {
			    const Score incoming_sc = ScorePMul (read_cell (*s, srcx_idx, srcy_idx), cfg.transition (*s, *d));
			    sc = max (sc, incoming_sc);
			  }
			// loop over incoming bifurcations to left & right
			for_const_contents (vector<Bifurcation_left_parent>, filter.left_parent[*d], b)
			  if (cfg.state_type[b->p] == Bifurc)  // Bifurc
			    for_const_contents (Subseq::Bifurc_outl_pseudovec, xsubseq.bif_out_l, bx)
			      for_const_contents (Subseq::Bifurc_outl_pseudovec, ysubseq.bif_out_l, by)
			      sc = max (sc, ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->l, bx->l, by->l)));
			  else  // BifurcRevY
			    for_const_contents (Subseq::Bifurc_outl_pseudovec, xsubseq.bif_out_l, bx)
			      for_const_contents (Subseq::Bifurc_outr_pseudovec, ysubseq.bif_out_r, by)
			      sc = max (sc, ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->l, bx->l, by->r)));
			for_const_contents (vector<Bifurcation_right_parent>, filter.right_parent[*d], b)
			  if (cfg.state_type[b->p] == Bifurc)  // Bifurc
			    for_const_contents (Subseq::Bifurc_outr_pseudovec, xsubseq.bif_out_r, bx)
			      for_const_contents (Subseq::Bifurc_outr_pseudovec, ysubseq.bif_out_r, by)
			      sc = max (sc, ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->r, bx->r, by->r)));
			  else  // BifurcRevY
			    for_const_contents (Subseq::Bifurc_outr_pseudovec, xsubseq.bif_out_r, bx)
			      for_const_contents (Subseq::Bifurc_outl_pseudovec, ysubseq.bif_out_l, by)
			      sc = max (sc, ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->r, bx->r, by->l)));
			// add emit score
			if (is_emit_type(t))
			  {
			    const int emit_idx = out_emit_idx (t, xsubseq, ysubseq);
			    ScorePMulAcc (sc, cfg.emit[*d][emit_idx]);
			    ScorePMulAcc (sc, out_meta_sc (*d, xsubseq, ysubseq));
			  }
			// store
			if (sc < -InfinityScore) THROWEXPR("Overflow");
			cell (*d, xsubseq_idx, ysubseq_idx) = sc;
		      }
		    // count the cell
		    ++cells_done;
		    bifs_done += xsubseq.bif_in.size() * ysubseq.bif_in.size();
		  }
	      }
	    if (bifs_done - last_report >= CFGDP_REPORT_INTERVAL || (bifs_done == total_bifs && bifs_done > last_report))
	      {
		last_report = bifs_done;
		CTAG(1,CFGDP) << "KYC: finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)total_bifs))/10. << "%)\n";
	      }
	  }
      }
  if (CTAGGING(-1,KYC_MATRIX))
    show(CL);
}

void Pair_KYC_matrix::traceback_from (int xsubseq_idx, int ysubseq_idx, int dest_state, vector<int>& path) const
{
  // this method only works for global DP matrices (since it doesn't return the start coords),
  // so first check that the envelopes & matrix are global, & throw an error if not
  if (local || !xenv.is_global() || !yenv.is_global())
    THROWEXPR ("Attempted a KYC traceback in a local DP matrix without keeping track of traceback start coords");
  // get subseqs
  const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
  const Subseq& ysubseq = yenv.subseq[ysubseq_idx];
  // find max incoming transition
  int best_src = HMM_state_enum::UndefinedState;
  Score best_sc = -InfinityScore;
  if (local_start_allowed (xsubseq, ysubseq))
    {
      best_src = HMM_state_enum::Start;
      best_sc = cfg.start[dest_state];
    }
  for_const_contents (vector<int>, allowed_src_out_states (dest_state, xsubseq, ysubseq), s)
    {
      const Score sc = ScorePMul (read_cell (*s, xsubseq_idx, ysubseq_idx),
				  cfg.transition (*s, dest_state));
      if (sc > best_sc)
	{
	  best_sc = sc;
	  best_src = *s;
	}
    }
  // start traceback from this transition
  traceback_from_transition (xsubseq_idx, ysubseq_idx, dest_state, best_src, path);
}

void Pair_KYC_matrix::traceback_from_transition (int xsubseq_idx, int ysubseq_idx, int dest_state, int src_state, vector<int>& path) const
{
  // print log message
  CTAG(4,KYC_TRACEBACK) << "Starting KYC traceback at xsubseq #" << xsubseq_idx << ", ysubseq #" << ysubseq_idx << ", dest state #" << dest_state << ", source state #" << src_state << "\n";

  // start with a CYK traceback from dest state
  cyk.traceback_from (xsubseq_idx, ysubseq_idx, dest_state, path, FALSE);

  // KYC traceback uses a reversed path; so, start by reversing the CYK traceback
  vector<int> revpath (path.rbegin(), path.rend());

  // now trace outwards through the KYC matrix
  int traceback_state = src_state;
  while (1)
    {
      // push state onto the reverse path
      revpath.push_back (traceback_state);

      // check if Start state reached
      if (traceback_state == HMM_state_enum::Start)
	break;

      // get the skinny on this cell
      const Subseq& xsubseq = xenv.subseq[xsubseq_idx];
      const Subseq& ysubseq = yenv.subseq[ysubseq_idx];
      const Score traceback_sc = read_cell (traceback_state, xsubseq_idx, ysubseq_idx);
      
      // get state type
      const State_type traceback_type = cfg.state_type[traceback_state];
      
      // find source subseqs
      const int srcx_idx = src_xsubseq_idx (xsubseq_idx, traceback_type);
      const int srcy_idx = src_ysubseq_idx (ysubseq_idx, traceback_type);
      const Subseq& srcx = xenv.subseq[srcx_idx];
      const Subseq& srcy = yenv.subseq[srcy_idx];
      
      // calculate emit score
      const Score emit_sc =
	is_emit_type (traceback_type)
	? ScorePMul (cfg.emit[traceback_state][out_emit_idx (traceback_type, xsubseq, ysubseq)],
		     out_meta_sc (traceback_state, xsubseq, ysubseq))
	: 0;

      // set up some variables for recording bifurcation CYK tracebacks
      vector<int> lpath, rpath;

      // do start transitions
      if (local_start_allowed (srcx, srcy) && traceback_sc == ScorePMul (cfg.start[traceback_state], emit_sc))
	{
	  traceback_state = HMM_state_enum::Start;
	  goto KYC_FoundTransition;
	}
      
      // loop over incoming transitions
      for_const_contents (vector<int>, allowed_src_out_states (traceback_state, srcx, srcy), s)
	if (traceback_sc == ScorePMul3 (read_cell (*s, srcx_idx, srcy_idx), cfg.transition (*s, traceback_state), emit_sc))
	  {
	    traceback_state = *s;
	    goto KYC_FoundTransition;
	  }
      
      // loop over incoming bifurcations to left & right
      for_const_contents (vector<Bifurcation_left_parent>, filter.left_parent[traceback_state], b)
	if (cfg.state_type[b->p] == Bifurc)  // Bifurc
	  for_const_contents (Subseq::Bifurc_outl_pseudovec, xsubseq.bif_out_l, bx)
	    for_const_contents (Subseq::Bifurc_outl_pseudovec, ysubseq.bif_out_l, by)
	  {
	    if (traceback_sc == ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->l, bx->l, by->l)))
	      {
		cyk.traceback_from (bx->l, by->l, b->l, lpath, FALSE);
		traceback_state = b->p;
		xsubseq_idx = bx->out;
		ysubseq_idx = by->out;
		goto KYC_FoundBifurcation;
	      }
	  }
	else  // BifurcRevY
	  for_const_contents (Subseq::Bifurc_outl_pseudovec, xsubseq.bif_out_l, bx)
	    for_const_contents (Subseq::Bifurc_outr_pseudovec, ysubseq.bif_out_r, by)
	    if (traceback_sc == ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->l, bx->l, by->r)))
	      {
		cyk.traceback_from (bx->l, by->r, b->l, lpath, FALSE);
		traceback_state = b->p;
		xsubseq_idx = bx->out;
		ysubseq_idx = by->out;
		goto KYC_FoundBifurcation;
	      }

      for_const_contents (vector<Bifurcation_right_parent>, filter.right_parent[traceback_state], b)
	if (cfg.state_type[b->p] == Bifurc)  // Bifurc
	  for_const_contents (Subseq::Bifurc_outr_pseudovec, xsubseq.bif_out_r, bx)
	    for_const_contents (Subseq::Bifurc_outr_pseudovec, ysubseq.bif_out_r, by)
	  {
	    if (traceback_sc == ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->r, bx->r, by->r)))
	      {
		cyk.traceback_from (bx->r, by->r, b->r, rpath, FALSE);
		traceback_state = b->p;
		xsubseq_idx = bx->out;
		ysubseq_idx = by->out;
		goto KYC_FoundBifurcation;
	      }
	  }
	else  // BifurcRevY
	  for_const_contents (Subseq::Bifurc_outr_pseudovec, xsubseq.bif_out_r, bx)
	    for_const_contents (Subseq::Bifurc_outl_pseudovec, ysubseq.bif_out_l, by)
	    if (traceback_sc == ScorePMul (read_cell (b->p, bx->out, by->out), cyk.read_cell (b->r, bx->r, by->l)))
	      {
		cyk.traceback_from (bx->r, by->l, b->r, rpath, FALSE);
		traceback_state = b->p;
		xsubseq_idx = bx->out;
		ysubseq_idx = by->out;
		goto KYC_FoundBifurcation;
	      }

      // if we get here, then traceback failed
      THROWEXPR (cfg_dump() << "Pair CFG: traceback failed at xi=" << xsubseq_idx << ", yi=" << ysubseq_idx << ", state=" << traceback_state);

    KYC_FoundTransition:
      xsubseq_idx = srcx_idx;
      ysubseq_idx = srcy_idx;
      continue;

    KYC_FoundBifurcation:
      // if there was a rightward bifurcation, insert the reversed right-branch traceback at the beginning of revpath
      if (rpath.size())
	revpath.insert (revpath.begin(), rpath.rbegin(), rpath.rend());
      // if there was a leftward bifurcation, insert the reversed left-branch traceback at the end of revpath
      if (lpath.size())
	revpath.insert (revpath.end(), lpath.rbegin(), lpath.rend());
      // continue with the traceback
      continue;
    }

  // path <-- reverse(revpath)
  path.clear();
  path.insert (path.begin(), revpath.rbegin(), revpath.rend());

  // print log messages
  if (CTAGGING(3,KYC_TRACEBACK))
    CL << "KYC traceback: (" << path << ")\n";
}

Score Pair_CYK_KYC_matrix::max_transition_sc (int src_state, int dest_state, int xsubseq_idx, int ysubseq_idx) const
{
  // get subseqs
  const Subseq& xsubseq = cyk.xenv.subseq[xsubseq_idx];
  const Subseq& ysubseq = cyk.yenv.subseq[ysubseq_idx];

  // evaluate
  const Score incoming_sc =
    src_state == HMM_state_enum::Start
    ? (cyk.local_start_allowed (xsubseq, ysubseq) ? 0 : -InfinityScore)
    : kyc.read_cell (src_state, xsubseq_idx, ysubseq_idx);

  const Score transition_sc = cyk.cfg.transition (src_state, dest_state);

  const Score outgoing_sc =
    dest_state == HMM_state_enum::End
    ? (xsubseq.len==0 && ysubseq.len==0 ? 0 : -InfinityScore)
    : cyk.read_cell (dest_state, xsubseq_idx, ysubseq_idx);
  
  return ScorePMul3 (incoming_sc, transition_sc, outgoing_sc);
}

Score Pair_CYK_KYC_matrix::max_state_sc (int dest_state, int xsubseq_idx, int ysubseq_idx) const
{
  // get subseqs
  const Subseq& xsubseq = cyk.xenv.subseq[xsubseq_idx];
  const Subseq& ysubseq = cyk.yenv.subseq[ysubseq_idx];
  // find max over incoming transitions
  Score sc = max_transition_sc (HMM_state_enum::Start, dest_state, xsubseq_idx, ysubseq_idx);
  for_const_contents (vector<int>, cyk.allowed_src_out_states (dest_state, xsubseq, ysubseq), s)
    sc = max (sc, max_transition_sc (*s, dest_state, xsubseq_idx, ysubseq_idx));
  return sc;
}

Pair_CYK_KYC_matrix::Pair_CYK_KYC_matrix (const Named_profile& npx,
					  const Named_profile& npy,
					  const Fold_envelope& xenv,
					  const Fold_envelope& yenv,
					  const Pair_CFG_scores& cfg,
					  bool local,
					  bool fill_now) :
  cyk (npx, npy, xenv, yenv, cfg, local, fill_now),
  kyc (cyk, fill_now)
{ }

Pair_CYK_KYC_matrix::Pair_CYK_KYC_matrix (const Named_profile& npx,
					  const Named_profile& npy,
					  const Fold_envelope& xenv,
					  const Fold_envelope& yenv,
					  const Pair_CFG_scores& cfg,
					  const Pair_envelope& pair_env,
					  bool local,
					  bool fill_now) :
  cyk (npx, npy, xenv, yenv, cfg, pair_env, local, fill_now),
  kyc (cyk, fill_now)
{ }

void Pair_CYK_KYC_matrix::show (ostream& o) const
{
  Biosequence xseq, yseq;
  cyk.cfg.alphabet().dsq2seq (cyk.xdsq, xseq);
  cyk.cfg.alphabet().dsq2seq (cyk.ydsq, yseq);
  o << "Sequence X: " << xseq << "\n";
  o << "Sequence Y: " << yseq << "\n";
  o << "Envelope X:\n";
  cyk.xenv.dump(o);
  o << "Envelope Y:\n";
  cyk.yenv.dump(o);
  o << "CYK matrix:\n";
  cyk.show(o);
  o << "KYC matrix:\n";
  kyc.show(o);
}

Pair_CYK_KYC_cell_sorter::Pair_CYK_KYC_cell_sorter (const Pair_CYK_KYC_matrix& cyk_kyc)
  : cyk_kyc (cyk_kyc)
{
  // print log message
  CTAG(5,CYK_KYC) << "Sorting cells by CYK-KYC score\n";
  // copy some types
  typedef vector<Subseq>::const_iterator Subseq_iterator;
  // copy some variables
  const Pair_CYK_matrix& cyk = cyk_kyc.cyk;
  const int xsubseqs = cyk.xsubseqs;
  const int ysubseqs = cyk.ysubseqs;
  const Subseq_iterator xenv_begin = cyk.xenv.subseq.begin();
  const Subseq_iterator yenv_begin = cyk.yenv.subseq.begin();
  // collect start scores for each cell
  deque<Cell_coords> unsorted_cells;
  deque<Score> unsorted_cells_sc;
  for (int xsubseq_idx = 0; xsubseq_idx < xsubseqs; ++xsubseq_idx)
    for (int ysubseq_idx = 0; ysubseq_idx < ysubseqs; ++ysubseq_idx)
      {
	// find state score
	const Subseq& xsubseq = xenv_begin[xsubseq_idx];
	const Subseq& ysubseq = yenv_begin[ysubseq_idx];
	if (cyk.in_pair_envelope (xsubseq, ysubseq))
	  for_const_contents (vector<int>, cyk.allowed_in_states (xsubseq, ysubseq), s)
	  {
	    const Score sc = cyk_kyc.max_state_sc (*s, xsubseq_idx, ysubseq_idx);
	    if (sc > -InfinityScore)
	      {
		unsorted_cells.push_back (Cell_coords (*s, xsubseq_idx, ysubseq_idx));
		unsorted_cells_sc.push_back (sc);
	      }
	  }
      }
  // sort cells
  vector<Score> unsorted_cells_score_vec (unsorted_cells_sc.begin(), unsorted_cells_sc.end());
  unsorted_cells_sc = deque<Score>();  // free up memory
  Reverse_Schwartzian<Score> by_score (unsorted_cells_score_vec);
  const vector<int> cell_order = by_score.sorted();
  unsorted_cells_score_vec = vector<Score>();  // free up memory
  // write sorted_cells vector
  const int N = unsorted_cells.size();
  sorted_cells = vector<Cell_coords> (N, Cell_coords(-1,-1,-1));
  for (int i = 0; i < N; ++i)
    sorted_cells[i] = unsorted_cells[cell_order[i]];
}
