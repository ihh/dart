#include "scfg/foldenv.h"
#include "indiegram/tripletscfgdp.h"

#define SCFGDP_REPORT_INTERVAL 3000000  /* number of bifurcations between logfile messages during DP */

typedef unsigned long long Bifcount; /// bifurcation count
typedef vector<const vector<int>*> Subseq_index_array;

// constructor
Triplet_SCFG_filter::Triplet_SCFG_filter (const Triplet_SCFG& scfg, const SCFG_state_sorter& sorter)
  : allowed_in_states (FlagTripletRange),
    allowed_out_states (FlagTripletRange),
    allowed_dest_states (FlagTripletRange, vector<vector<int> > (scfg.num_states())),
    allowed_src_states (FlagTripletRange, vector<vector<int> > (scfg.num_states())),
    left_parent (scfg.left_parent()),
    right_parent (scfg.right_parent())
{
  const vector<int> emit_states = sorter.emit_states();
  const vector<int> nonemit_states = sorter.nonemit_states_sorted();

  for (int xflag = 0; xflag < FlagRange; ++xflag)
    for (int yflag = 0; yflag < FlagRange; ++yflag)
      for (int zflag = 0; zflag < FlagRange; ++zflag)
	{
	  const int ft = flag_triplet (xflag, yflag, zflag);
	  for_const_contents (vector<int>, emit_states, s)
	    {
	      bool allowed = true;
	      const State_type t = scfg.state_type[*s];
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
	      switch (t & EmitZLR)
		{
		case EmitZL:
		  if (!(zflag & Subseq::CFLAG_L)) allowed = 0;
		  break;
		case EmitZR:
		  if (!(zflag & Subseq::CFLAG_R)) allowed = 0;
		  break;
		case EmitZLR:
		  if (!(zflag & Subseq::CFLAG_LR)) allowed = 0;
		  break;
		default:
		  break;
		}
	      if (allowed)
		{
		  allowed_in_states[ft].push_back (*s);
		  allowed_out_states[ft].push_back (*s);
		}

	    }
	  // insert Null and Bifurc states in reverse topological (Inside) order:
	  //   NB: void insert( iterator loc, input_iterator start, input_iterator end )
	  //      inserts the elements from start to end before loc
	  //   nonemit_states_sorted() returns states sorted topologically
	  allowed_in_states[ft].insert (allowed_in_states[ft].end(), nonemit_states.rbegin(), nonemit_states.rend());
	  // insert Null and Bifurc states in topological (Outside) order:
	  allowed_out_states[ft].insert (allowed_out_states[ft].end(), nonemit_states.begin(), nonemit_states.end());
	}

  //
  for (int xflag = 0; xflag < FlagRange; ++xflag)
    for (int yflag = 0; yflag < FlagRange; ++yflag)
      for (int zflag = 0; zflag < FlagRange; ++zflag)
	{
	  const int ft = flag_triplet (xflag, yflag, zflag);
	  // scfg.selected_outgoing_states (allowed_in_states[ft]) returns a vector<vector <int> >
	  // mapping each state whose emissions match the Flag_triplet ft to a vector of its possible
	  // destination states
	  allowed_dest_states[ft] = scfg.selected_outgoing_states (allowed_in_states[ft]);
	  allowed_src_states[ft] = scfg.selected_incoming_states (allowed_in_states[ft]);

	  // debugging...
#ifdef DART_DEBUG
	  if (CTAGGING (-1,INDIEGRAM_DP)) {
	    if ((xflag == Subseq::CFLAG_NONE) && (yflag == Subseq::CFLAG_NONE) && (zflag == Subseq::CFLAG_NONE)) {
	      CL << "allowed_in_states[" << ft << "] = ";
	      for (vector<int>::const_iterator s = allowed_in_states[ft].begin(); s != allowed_in_states[ft].end(); ++s)
		CL << *s << "(" << state_type_string (scfg.state_type[*s]) << "), ";
	      CL << endl;
	    }
	  }
#endif /* DART_DEBUG */

	}

}

Triplet_DP_allocator_base::Triplet_DP_allocator_base (const Triplet_SCFG& scfg,
						      const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
						      const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z)
  : scfg (scfg),
    sorter (SCFG_state_sorter (scfg)),
    np_x (np_x), np_y (np_y), np_z (np_z),
    foldenv_x (foldenv_x), foldenv_y (foldenv_y), foldenv_z (foldenv_z),
    dsq_x (np_x.dsq), dsq_y (np_y.dsq), dsq_z (np_z.dsq),
    seqlen_x (np_x.size()), seqlen_y (np_y.size()), seqlen_z (np_z.size()),
    filter (scfg, sorter),
    megabytes_allocated (0),
    num_cells (0),
    num_bifs (0)
{

  // make sure the fold envelopes are global
  if (!(foldenv_x.is_global() && foldenv_y.is_global() && foldenv_z.is_global()))
    THROWEXPR ("Fold envelopes not global");

}

void Triplet_DP_allocator_base::count_cells()
{

  // counts
  register Bifcount n_cells = 0;
  register Bifcount n_bifs = 0;

  // The loop is just the inside ordering used in e.g. CYK.
  for (int startpos_x = seqlen_x; startpos_x >= 0; --startpos_x)
    {
      if (foldenv_x.by_start[startpos_x].size())
	{
	  // loop over X subseqs with start position startpos_x
	  for_const_contents (vector<int>, foldenv_x.by_start[startpos_x], xsi)
	    {
	      const Subseq& subseq_x = foldenv_x.subseq[*xsi];

	      // loop over Y start positions and corresponding subseqs
	      for (int startpos_y = seqlen_y; startpos_y >= 0; --startpos_y)
		{
		  if (foldenv_y.by_start[startpos_y].size())
		    {
		      for_const_contents (vector<int>, foldenv_y.by_start[startpos_y], ysi)
			{
			  const Subseq& subseq_y = foldenv_y.subseq[*ysi];


			  // loop over Z start positions and corresponding subseqs
			  for (int startpos_z = seqlen_z; startpos_z >= 0; --startpos_z)
			    {
			      if (foldenv_z.by_start[startpos_x].size())
				{
				  for_const_contents (vector<int>, foldenv_z.by_start[startpos_z], zsi)
				    {
				      const Subseq& subseq_z = foldenv_z.subseq[*zsi];
				      
				      // do the counting
				      ++n_cells;
				      n_bifs += subseq_x.bif_in.size() * subseq_y.bif_in.size() * subseq_z.bif_in.size();

				    }
				}
			    } // z loop

			}
		    }
		} // y loop

	    }
	}
    } // x loop

  // store the counts
  num_cells = n_cells;
  num_bifs = n_bifs;

}

Triplet_DP_dense_allocator::Triplet_DP_dense_allocator (const Triplet_SCFG& scfg,
							const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
							const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z)
  : Triplet_DP_allocator_base (scfg, np_x, np_y, np_z, foldenv_x, foldenv_y, foldenv_z)
{ }

void Triplet_DP_dense_allocator::alloc()
{

  const int num_states = scfg.num_states(); // doesn't count start and end states (don't get their own DP tables)

  // NB: foldenv_x.subseq.size() is of course equal to cells_needed_x as:
  //  int cells_needed_x = 0;
  //  for (int x = 0; x <= seqlen_x; ++x)
  //    cells_needed_x += foldenv_x.by_start[x].size();

  // pre-allocation log message
  const int cells_needed_x = foldenv_x.subseq.size();
  const int cells_needed_y = foldenv_y.subseq.size();
  const int cells_needed_z = foldenv_z.subseq.size();
  const unsigned long int mb_wanted = (int) floor ((1.0 * cells_needed_x * cells_needed_y * cells_needed_z / 1048576.0) * num_states * sizeof(Score));
  CTAG(7,ALLOC) << "Trying to allocate approx " << num_states << " * (" << cells_needed_x << " * " << cells_needed_y << " * " << cells_needed_z << ") cells * " << sizeof(Score) << " bytes => " << mb_wanted << " Mb.\n";

  // allocate the memory
  this->_cell = XVec (cells_needed_x, YVec (cells_needed_y, ZVec (cells_needed_z, StateVec (num_states, -InfinityScore))));
  megabytes_allocated = mb_wanted;

}



Triplet_DP_sparse_allocator::Triplet_DP_sparse_allocator (const Triplet_SCFG& scfg,
							  const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
							  const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z)
  : Triplet_DP_allocator_base (scfg, np_x, np_y, np_z, foldenv_x, foldenv_y, foldenv_z),
    foldenv_x_begin (foldenv_x.subseq.begin()),
    foldenv_y_begin (foldenv_y.subseq.begin()),
    foldenv_z_begin (foldenv_z.subseq.begin())
{ }

void Triplet_DP_sparse_allocator::alloc()
{
  // lookup table
  megabytes_allocated = sizeof (array2d<By_start_indices>) + (seqlen_x+1) * (seqlen_y+1) * sizeof (By_start_indices);
  CTAG(7,ALLOC) << "Trying to allocate " << seqlen_x+1 << " * " << seqlen_y+1 << " * " << seqlen_z+1 << " lookup table (approx " << megabytes_allocated << " Mb)\n";

  // resize by_start_indices
  by_start_indices.resize (seqlen_z+1);
  for (int z = 0; z <= seqlen_z; ++z)
    by_start_indices[z].resize (seqlen_x+1, seqlen_y+1, By_start_indices());

  // how many cells should we allocate?
  unsigned int cells_wanted = 0;
  for (int x = 0; x <= seqlen_x; ++x)
    for (int y = 0; y <= seqlen_y; ++y)
      for (int z = 0; z <= seqlen_z; ++z)
	if (1) 	//	if (startpos_in_triplet_alignment (x, y, z))
	  {
	    const int xs = foldenv_x.by_start[x].size();
	    const int ys = foldenv_y.by_start[y].size();
	    const int zs = foldenv_z.by_start[z].size();
	    By_start_indices& ind = by_start_indices[z] (x, y);
	    ind.offset = cells_wanted;
	    ind.xsize = xs;
	    ind.zsize = zs;
	    cells_wanted += scfg.num_states() * xs * ys * zs;
	  }
  
  const unsigned long int mb_wanted = (int) floor (1.0 * cells_wanted * sizeof (Score) / 1048576.0);
  CTAG(7,ALLOC) << "Trying to allocate main DP table: " << cells_wanted << " cells (" << mb_wanted << " Mb)\n";

  // allocate
  _cell.resize (cells_wanted, -InfinityScore);

  // we're done
  megabytes_allocated += mb_wanted;

}


// constructor
Triplet_DP_matrix_base::Triplet_DP_matrix_base (const Triplet_SCFG& scfg,
						const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
						const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z)
  : TRIPLET_DP_MATRIX_BASE (scfg, np_x, np_y, np_z, foldenv_x, foldenv_y, foldenv_z),
    _show_out (false)
{
  
  // initialize final score
  if (seqlen_x == 0 && seqlen_y == 0 && seqlen_z == 0)
    final_sc = scfg.transition_scores.start_to_end();
  else
    final_sc = -InfinityScore;

  // count cells and bifurcations and log accordingly
  count_cells();
  if (CTAGGING(5,ALLOC CELLCOUNT))
    CL << "Counted " << num_cells << " subseqs and " << num_bifs << " bifurcations.\n";

  // allocate DP matrix and print log messages
  if (CTAGGING(7,ALLOC PREALLOC))
    CL << "Allocating Triplet SCFG DP matrix: " << foldenv_x.subseq.size() << " * " << foldenv_y.subseq.size() << " * " << foldenv_z.subseq.size() << " cells by " << scfg.num_states() << " states.\n";
  alloc();
  if (CTAGGING(7,ALLOC POSTALLOC))
    CL << "Allocated Triplet SCFG DP matrix (approx " << megabytes_allocated << " Mb).\n";

}

void Triplet_DP_matrix_base::show (ostream& o) const
{
  int old_prec = o.precision (3);
  save_flags (o);

  const int page_len = 40;
  sstring page_header;
  page_header << "xi  start+len  xseq XL   XR   XLR  yi  start+len  yseq YL   YR   YLR  zi  start+len  zseq ZL   ZR   ZLR  ";
  for (int s = 0; s < scfg.num_states(); ++s)
    {
      sstring state_desc; // state description
      state_desc << s << '(' << state_type_string (scfg.state_type[s]) << ')';
      right_align (page_header);
      page_header.width (15);
      page_header << state_desc;
    }
  page_header << "\n";

  // get sequences
  Biosequence xseq;
  Biosequence yseq;
  Biosequence zseq;
  scfg.alphabet().dsq2seq (dsq_x, xseq);
  scfg.alphabet().dsq2seq (dsq_y, yseq);
  scfg.alphabet().dsq2seq (dsq_z, zseq);

  // display the DP matrix
  int lines = 0;
  for (int subseq_idx_x = 0; subseq_idx_x < (int) foldenv_x.subseq.size(); ++subseq_idx_x)
    for (int subseq_idx_y = 0; subseq_idx_y < (int) foldenv_y.subseq.size(); ++subseq_idx_y)
      for (int subseq_idx_z = 0; subseq_idx_z < (int) foldenv_z.subseq.size(); ++subseq_idx_z)
	{
	  // break into pages
	  if (lines++ % page_len == 0)
	    {
	      if (lines > 1) o << "\n";
	      o << page_header;
	    }

	  o << foldenv_x.subseq[subseq_idx_x].terse_desc (subseq_idx_x, xseq, _show_out);
	  o << foldenv_y.subseq[subseq_idx_y].terse_desc (subseq_idx_y, yseq, _show_out);
	  o << foldenv_z.subseq[subseq_idx_z].terse_desc (subseq_idx_z, zseq, _show_out);
	  right_align (o);
	  for (int s = 0; s < scfg.num_states(); ++s)
	    {
	      o.width (15);
	      ShowScore (read_cell (s, subseq_idx_x, subseq_idx_y, subseq_idx_z), o);
	    }
	  o << "\n";
	}

  // add score information
  o << "Final score: ";
  ShowScore (final_sc, o);
  o << "\n";
  restore_flags (o);
  o.precision (old_prec);  

}


void Triplet_DP_matrix_base::show_compact (ostream& o) const
{
  int old_prec = o.precision (3);
  save_flags (o);

  // get sequences
  Biosequence xseq;
  Biosequence yseq;
  Biosequence zseq;
  scfg.alphabet().dsq2seq (dsq_x, xseq);
  scfg.alphabet().dsq2seq (dsq_y, yseq);
  scfg.alphabet().dsq2seq (dsq_z, zseq);

  // display the DP matrix
  for (int subseq_idx_x = 0; subseq_idx_x < (int) foldenv_x.subseq.size(); ++subseq_idx_x)
    for (int subseq_idx_y = 0; subseq_idx_y < (int) foldenv_y.subseq.size(); ++subseq_idx_y)
      for (int subseq_idx_z = 0; subseq_idx_z < (int) foldenv_z.subseq.size(); ++subseq_idx_z)
	{

	  for (int s = 0; s < scfg.num_states(); ++s)
	    {

	      Score sc = read_cell (s, subseq_idx_x, subseq_idx_y, subseq_idx_z);
	      if (sc > -InfinityScore)
		{
		  sstring state_desc; // state description
		  state_desc << s << '(' << state_type_string (scfg.state_type[s]) << ')';
		  o << "(" << state_desc << ", "
		    << foldenv_x.subseq[subseq_idx_x].terser_desc() << ", "
		    << foldenv_y.subseq[subseq_idx_y].terser_desc() << ", "
		    << foldenv_z.subseq[subseq_idx_z].terser_desc()
		    << ") => ";
		  ShowScore (sc, o);
		  o << "\n";
		}
	    }

	  o << "\n";

	}

  // add score information
  o << "Final score: ";
  ShowScore (final_sc, o);
  o << "\n";
  restore_flags (o);
  o.precision (old_prec);  

}


sstring Triplet_DP_matrix_base::scfg_dump() const
{
  ostream& dump = CLOGERR;
  dump << "Triplet SCFG DP matrix:\n";
  //  show (dump);
  show_compact (dump);
  dump << "Triplet SCFG:\n";
  scfg.show (dump);
  dump << "Envelope X:\n";
  foldenv_x.dump (dump);
  dump << "Envelope Y:\n";
  foldenv_y.dump (dump);
  dump << "Envelope Z:\n";
  foldenv_z.dump (dump);
  
  return sstring();
}





Triplet_inside_matrix::Triplet_inside_matrix (const Triplet_SCFG& scfg,
					      const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
					      const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z,
					      bool fill_now /* = true */)
  : Triplet_DP_matrix_base (scfg,
			    np_x, np_y, np_z,
			    foldenv_x, foldenv_y, foldenv_z)
{
  if (fill_now)
    fill();
}

// main loop of Inside algorithm
// The interpretation of Inside scores:
// cell (a, subseq_idx_x, subseq_idx_y, subseq_idx_z) is the
// summed score that subseqs (subseq_x, subseq_y, subseq_z) were generated
// by a parse tree rooted at state a
// => P (subseq_x, subseq_y, subseq_z | subtree generating subseq_x, subseq_y, subseq_z is rooted at a).
void Triplet_inside_matrix::fill()
{

  // print log message: starting Inside fill
  CTAG(7,INDIEGRAM_DP) << "Starting Inside fill\n";

  // for logging progress of DP
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;
  Bifcount last_report = 0;

  // loop over X start positions
  for (int startpos_x = seqlen_x; startpos_x >= 0; --startpos_x)
    {
      if (foldenv_x.by_start[startpos_x].size())
	{
	  // loop over X subseqs with start position startpos_x
	  // to do: Ian notes in foldenv.h that .by_start[START][] isn't in a valid
	  // inside-outside ordering.  Is this true, and if so then why is the below loop ok?
	  for_const_contents (vector<int>, foldenv_x.by_start[startpos_x], xsi)
	    {
	      const int subseq_idx_x = *xsi;
	      const Subseq& subseq_x = foldenv_x.subseq[subseq_idx_x];

	      // loop over Y start positions and corresponding subseqs
	      for (int startpos_y = seqlen_y; startpos_y >= 0; --startpos_y)
		{
		  if (foldenv_y.by_start[startpos_y].size())
		    {
		      for_const_contents (vector<int>, foldenv_y.by_start[startpos_y], ysi)
			{
			  const int subseq_idx_y = *ysi;
			  const Subseq& subseq_y = foldenv_y.subseq[subseq_idx_y];

			  // loop over Z start positions and corresponding subseqs
			  for (int startpos_z = seqlen_z; startpos_z >= 0; --startpos_z)
			    {
			      if (foldenv_z.by_start[startpos_z].size())
				{
				  for_const_contents (vector<int>, foldenv_z.by_start[startpos_z], zsi)
				    {
				      const int subseq_idx_z = *zsi;
				      const Subseq& subseq_z = foldenv_z.subseq[subseq_idx_z];
				      
				      // Fill in all cells with the subsequence indices (subseq_x, subseq_y, subseq_z)
				      //  for all states whose production rules match the connection flags for this subsequence triplet.
				      //  NB: this is the "source" state (left-hand side of production rules)
				      for_const_contents (vector<int>, allowed_in_states (subseq_x, subseq_y, subseq_z), s)
					{

					  const State_type t = scfg.state_type[*s];
					  Score score;

					  // bifurcation
					  if (is_bifurc_type (t))
					    {
					      // get left and right child states
					      const int l = scfg.bifurc[*s].l;
					      const int r = scfg.bifurc[*s].r;

					      score = -InfinityScore;
					      // loop over all "bifurcation connections" and sum over all inside subseq pairs
					      //  (possible inside subsequence pairs corresponding to the left and right child states of the current bifurcation state)
					      // (*bx).l is the subsequence index for the left inside subseq (corresponding to the left child state)
					      for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_x.bif_in, bx)
						{
						  for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_y.bif_in, by)
						    {
						      for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_z.bif_in, bz)
							{
							  // no transition score associated with bifurcations in our framework
							  ScorePSumAcc (score, ScorePMul (read_cell (l, bx->l, by->l, bz->l), read_cell (r, bx->r, by->r, bz->r)));
							}
						    }
						}
					    }
					  // emit or otherwise
					  else
					    {
					      // now it's time to perform the summation over possible destination states corresponding 
					      // to the inside subsequence (contained in this one) reached by an emission from the current "source" state (*s)

					      // find the subsequence triplet reached by an emission from the current state (*s)
					      // dest_idx_x is the X subseq index of reachable cells
					      const int dest_idx_x = dest_subseq_idx_x (subseq_idx_x, t);
					      const int dest_idx_y = dest_subseq_idx_y (subseq_idx_y, t);
					      const int dest_idx_z = dest_subseq_idx_z (subseq_idx_z, t);
					      // dest_x is the corresponding X subseq
					      const Subseq& dest_x = foldenv_x.subseq[dest_idx_x];
					      const Subseq& dest_y = foldenv_y.subseq[dest_idx_y];
					      const Subseq& dest_z = foldenv_z.subseq[dest_idx_z];
				      
					      // if the subseqs are empty, account for transitions to end
					      if (dest_x.len == 0 && dest_y.len == 0 && dest_z.len == 0)
						score = scfg.transition_scores.end[*s];
					      else
						score = -InfinityScore;

					      // Sum over the scores corresponding to parse trees rooted at the 
					      // possible destination state as well as the associated transition scores:
					      //  Sum over all possible destination states
					      //  corresponding to the inside subseq pair (dest_x, dest_y, dest_z)
					      //  reached by an emission from the current "source" state (*s).
					      // allowed_dest_in_states (state, dest_x, dest_y, dest_z) returns a vector
					      //  of all possible destination states from 'state', given 
					      //  the connection flags associated with the destination subsequence triplet (dest_x, dest_y, dest_z).
					      // Ie, we're considering all cells cell (s', dest_idx_x, dest_idx_y, dest_idx_z) reachable 
					      // from the current cell cell (*s, subseq_idx_x, subseq_idx_y, subseq_idx_z) and summing over s'.
					      for_const_contents (vector<int>, allowed_dest_in_states (*s, dest_x, dest_y, dest_z), d)
						ScorePSumAcc (score, ScorePMul (read_cell (*d, dest_idx_x, dest_idx_y, dest_idx_z), scfg.transition_scores.transition (*s, *d)));
					      
					      // Multiply by the emit score for the current state.
					      // This factor must come here (rather than at the start of the block).
					      // NB: remember that we're using a formulation where the source state (*s) does the emitting
					      // Also note that we're only looping over states (*s) whose production rules match the connection flags 
					      // for the subsequence triplet (subseq_x, subseq_y, subseq_z), so the emit score emit[*s][...]
					      // is necessarily defined for this state (*s) (it won't in general be defined for an arbitrary emit state).
					      if (is_emit_type (t))
						ScorePMulAcc (score, scfg.emit[*s][in_emit_idx (t, subseq_x, subseq_y, subseq_z)]);

					    }

					  // now store the score
					  cell (*s, subseq_idx_x, subseq_idx_y, subseq_idx_z) = score;

					  // is a start transition allowed for the current subsequence triplet?
					  // because we're doing global alignment, this is only allowed for the maximal subseqs
					  if (local_start_allowed (subseq_x, subseq_y, subseq_z))
					    ScorePSumAcc (final_sc, ScorePMul (score, scfg.transition_scores.start[*s])); // scfg.start[*s] is the score of the transition Start -> (*s)

					} // loop over states

				      // count the cell
				      ++cells_done;
				      bifs_done += subseq_x.bif_in.size() * subseq_y.bif_in.size() * subseq_z.bif_in.size();
				    }
				}
			    } // z loop

			  // log progress
			  if (bifs_done - last_report >= SCFGDP_REPORT_INTERVAL || (bifs_done == num_bifs && bifs_done > last_report))
			    {
			      last_report = bifs_done;
			      CTAG(7,INDIEGRAM_DP) << "Inside: Finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)num_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)num_bifs))/10. << "%)\n";
			    }
			}
		    }
		} // y loop

	    }
	}
    } // x loop

  // all done!  : )

}


Triplet_outside_matrix::Triplet_outside_matrix (const Triplet_inside_matrix& inside, bool fill_now /* = true */)
  : Triplet_DP_matrix_base (inside.scfg,
			    inside.np_x, inside.np_y, inside.np_z,
			    inside.foldenv_x, inside.foldenv_y, inside.foldenv_z),
    inside (inside)
{
  _show_out = 1;
  if (fill_now)
    fill();
}


// The interpretation of Outside scores:
// cell (a, subseq_idx_x, subseq_idx_y, subseq_idx_z) is the
// summed score of a (incomplete) parse tree rooted at the Start state
// which generates (X \ {subseq_x}, Y \ {subseq_y}, Z \ {subseq_z}) and terminates
// at state a for the subseqs (subseq_x, subseq_y, subseq_z)
// (so the subseqs (subseq_x, subseq_y, subseq_z) are generated by a parse tree rooted at a).
// => P (SUBSEQ, subtree generating SUBSEQ is rooted at a)

/* If state a is an Emit state, then we associate its emission (e.g. (x_{i-1} x_{j+1})) to /outside/ of SUBSEQ,
 * so a is the root of the subtree generating SUBSEQ in the sense of
 *          a
 *   /      |					\
 * (x_{i-1}) (subtree) (x_{j+1})
 *         /					\
 *      (SUBSEQ)
 */

// This formulation of the Outside probabilities (compared to the more standard version in Durbin)
// makes calculation of posterior transition probabilities easier but posterior state probabilities 
// more difficult (vice versa for the Durbin version).
// See my notes for more details.  Either interpretation is valid; I think that the dart one is nicer,
// at least when dealing with fold envelopes.
// 
// Ian says:
// meaning of a cell is different for outside matrix than for inside matrix
// for outside, cell (STATE,SUBSEQ) is the score for the parse tree landing in STATE
// having emitted all but SUBSEQ (including the emission from STATE)
// thus the emissions are outside of SUBSEQ (c.f. inside matrix where emissions are inside).
void Triplet_outside_matrix::fill()
{

  // print log message: starting Outside fill
  CTAG(7,INDIEGRAM_DP) << "Starting Outside fill\n";

  // check that final score of inside algorithm was finite
  const Score inside_final_sc = inside.final_sc;
  if (inside_final_sc <= -InfinityScore)
    {
      CLOGERR << "WARNING: inside algorithm returned score of -infinity; skipping outside algorithm\n";
      return;
    }
  final_sc = inside_final_sc;

  // do the DP
  Bifcount last_report = 0;
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;

  // loop over X start positions
  for (int startpos_x = 0; startpos_x <= seqlen_x; ++startpos_x)
    if (foldenv_x.by_start[startpos_x].size())
      {
	// loop over X subseqs with start position startpos_x
	for_const_reverse_contents (vector<int>, foldenv_x.by_start[startpos_x], xsi)
	  {
	    const int subseq_idx_x = *xsi;
	    const Subseq& subseq_x = foldenv_x.subseq[subseq_idx_x];

	    // loop over Y start positions and corresponding subseqs
	    for (int startpos_y = 0; startpos_y <= seqlen_y; ++startpos_y)
	      if (foldenv_y.by_start[startpos_y].size())
		{
		  for_const_reverse_contents (vector<int>, foldenv_y.by_start[startpos_y], ysi)
		    {
		      const int subseq_idx_y = *ysi;
		      const Subseq& subseq_y = foldenv_y.subseq[subseq_idx_y];

		      // loop over Z start positions and corresponding subseqs
		      for (int startpos_z = 0; startpos_z <= seqlen_z; ++startpos_z)
			if (foldenv_z.by_start[startpos_z].size())
			  {
			    for_const_reverse_contents (vector<int>, foldenv_z.by_start[startpos_z], zsi)
			      {
				const int subseq_idx_z = *zsi;
				const Subseq& subseq_z = foldenv_z.subseq[subseq_idx_z];

				for_const_contents (vector<int>, allowed_out_states (subseq_x, subseq_y, subseq_z), d)
				  {
				    const State_type t = scfg.state_type[*d];
				    Score score;

				    // get source subseqs
				    const int src_idx_x = src_subseq_idx_x (subseq_idx_x, t);
				    const int src_idx_y = src_subseq_idx_y (subseq_idx_y, t);
				    const int src_idx_z = src_subseq_idx_z (subseq_idx_z, t);

				    const Subseq& src_x = foldenv_x.subseq[src_idx_x];
				    const Subseq& src_y = foldenv_y.subseq[src_idx_y];
				    const Subseq& src_z = foldenv_z.subseq[src_idx_z];

				    // start transitions
				    if (local_start_allowed (src_x, src_y, src_z))
				      score = scfg.transition_scores.start[*d];
				    else
				      score = -InfinityScore;

				    // loop over incoming transitions
				    for_const_contents (vector<int>, allowed_src_out_states (*d, src_x, src_y, src_z), s)
				      {
					const Score incoming_sc = ScorePMul (read_cell (*s, src_idx_x, src_idx_y, src_idx_z), scfg.transition_scores.transition (*s, *d));
					ScorePSumAcc (score, incoming_sc);
				      }

				    // loop over incoming bifurcations to left & right
				    for_const_contents (vector<Bifurcation_left_parent>, filter.left_parent[*d], b)
				      if (scfg.state_type[b->p] == Bifurc) // bifurc
					for_const_contents (Subseq::Bifurc_outl_pseudovec, subseq_x.bif_out_l, bx)
					  {
					    for_const_contents (Subseq::Bifurc_outl_pseudovec, subseq_y.bif_out_l, by)
					      {
						for_const_contents (Subseq::Bifurc_outl_pseudovec, subseq_z.bif_out_l, bz)
						  ScorePSumAcc (score, ScorePMul (read_cell (b->p, bx->out, by->out, bz->out), inside.read_cell (b->l, bx->l, by->l, bz->l)));
					      }
					  }
				    for_const_contents (vector<Bifurcation_right_parent>, filter.right_parent[*d], b)
				      if (scfg.state_type[b->p] == Bifurc) // bifurc
					for_const_contents (Subseq::Bifurc_outr_pseudovec, subseq_x.bif_out_r, bx)
					  {
					    for_const_contents (Subseq::Bifurc_outr_pseudovec, subseq_y.bif_out_r, by)
					      {
						for_const_contents (Subseq::Bifurc_outr_pseudovec, subseq_z.bif_out_r, bz)
						  ScorePSumAcc (score, ScorePMul (read_cell (b->p, bx->out, by->out, bz->out), inside.read_cell (b->r, bx->r, by->r, bz->r)));
					      }
					  }

				    // add emit score
				    if (is_emit_type(t))
				      {
					const int emit_idx = out_emit_idx (t, subseq_x, subseq_y, subseq_z);
					ScorePMulAcc (score, scfg.emit[*d][emit_idx]);
				      }

				    // store
				    if (score < -InfinityScore) THROWEXPR ("Overflow");
				    cell (*d, subseq_idx_x, subseq_idx_y, subseq_idx_z) = score;

				  } // loop over states

				++cells_done;
				bifs_done += subseq_x.bif_in.size() * subseq_y.bif_in.size() * subseq_z.bif_in.size();

			      }
			  } // z loop

		      // log progress
		      if (bifs_done - last_report >= SCFGDP_REPORT_INTERVAL || (bifs_done == num_bifs && bifs_done > last_report))
			{
			  last_report = bifs_done;
			  CTAG(7,INDIEGRAM_DP) << "Outside: Finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)num_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)num_bifs))/10. << "%)\n";
			}

		    }
		} // y loop

	  }
      } // x loop

  if (CTAGGING(-1,OUTSIDE_MATRIX))
    show(CL);

}


// CYK matrix constructor
Triplet_CYK_matrix::Triplet_CYK_matrix (const Triplet_SCFG& scfg,
					const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
					const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z,
					bool fill_now /* = true */)
  : Triplet_DP_matrix_base (scfg,
			    np_x, np_y, np_z,
			    foldenv_x, foldenv_y, foldenv_z),
    final_subseq_idx_x (-1),
    final_subseq_idx_y (-1),
    final_subseq_idx_z (-1),
    final_state (Grammar_state_enum::UndefinedState)
{
  if (fill_now)
    fill();
}


// main loop of CYK algorithm
// Compare with Pair_CYK_matrix::fill(): We currently don't allow alignment envelope constraints,
// so there's no need to fool with the 'startpos_in_pair_envelope' stuff.
// Note also that there are two constructors for Pair_CYK_matrix, one which takes
// a Pair_envelope as a parameter and one which doesn't.  The one which doesn't
// creates a full (unconstrained) alignment envelope in the constructor body.
// This is equivalent to what we're doing.
// cell (state, subseq_idx_x, subseq_idx_y, subseq_idx_z) is the
// score for the best parse tree of (subseq_idx_x, subseq_idx_y, subseq_idx_z)
// rooted at state.
void Triplet_CYK_matrix::fill()
{

  // print log message: starting CYK fill
  CTAG(7,INDIEGRAM_DP) << "Starting CYK fill\n";

  // for logging progress of DP
  Bifcount cells_done = 0;
  Bifcount bifs_done = 0;
  Bifcount last_report = 0;

  // loop over X start positions
  for (int startpos_x = seqlen_x; startpos_x >= 0; --startpos_x)
    {
      if (foldenv_x.by_start[startpos_x].size())
	{
	  // loop over X subseqs with start position startpos_x
	  // to do: Ian notes in foldenv.h that .by_start[START][] isn't in a valid
	  // inside-outside ordering.  Is this true, and if so then why is the below loop ok?
	  for_const_contents (vector<int>, foldenv_x.by_start[startpos_x], xsi)
	    {
	      const int subseq_idx_x = *xsi;
	      const Subseq& subseq_x = foldenv_x.subseq[subseq_idx_x];

	      // loop over Y start positions and corresponding subseqs
	      for (int startpos_y = seqlen_y; startpos_y >= 0; --startpos_y)
		{
		  if (foldenv_y.by_start[startpos_y].size())
		    {
		      for_const_contents (vector<int>, foldenv_y.by_start[startpos_y], ysi)
			{
			  const int subseq_idx_y = *ysi;
			  const Subseq& subseq_y = foldenv_y.subseq[subseq_idx_y];

			  // loop over Z start positions and corresponding subseqs
			  for (int startpos_z = seqlen_z; startpos_z >= 0; --startpos_z)
			    {
			      if (foldenv_z.by_start[startpos_z].size())
				{
				  for_const_contents (vector<int>, foldenv_z.by_start[startpos_z], zsi)
				    {
				      const int subseq_idx_z = *zsi;
				      const Subseq& subseq_z = foldenv_z.subseq[subseq_idx_z];
				      

				      // debugging
#ifdef DART_DEBUG
				      if (((subseq_x.len == 0) && (subseq_y.len == 0) && (subseq_z.len == 0)) || ((subseq_x.len == seqlen_x) && (subseq_y.len == seqlen_y) && (subseq_z.len == seqlen_z)))
					{
					  CL << "(" << subseq_x.terser_desc() << ", "<< subseq_y.terser_desc() << ", "<< subseq_z.terser_desc() << "):\n";
					  CL << "  allowed_in_states  = ";
					  for_const_contents (vector<int>, allowed_in_states (subseq_x, subseq_y, subseq_z), s)
					    {
					      sstring state_desc; // state description
					      state_desc << *s << '(' << state_type_string (scfg.state_type[*s]) << ')';
					      CL << state_desc << ", ";
					    }
					  CL << "\n";
					}
#endif /* DART_DEBUG */

				      // Fill in all cells with the subsequence indices (subseq_x, subseq_y, subseq_z)
				      //  for all states whose production rules match the connection flags for this subsequence triplet.
				      //  NB: this is the "source" state (left-hand side of production rules)
				      for_const_contents (vector<int>, allowed_in_states (subseq_x, subseq_y, subseq_z), s)
					{

					  const State_type t = scfg.state_type[*s];
					  Score score;

					  // bifurcation
					  if (is_bifurc_type (t))
					    {
					      // get left and right child states
					      const int l = scfg.bifurc[*s].l;
					      const int r = scfg.bifurc[*s].r;

					      score = -InfinityScore;
					      // loop over all "bifurcation connections" and take the highest-scoring inside subseq pair
					      //  (possible inside subsequence pairs corresponding to the left and right child states of the current bifurcation state)
					      // (*bx).l is the subsequence index for the left inside subseq (corresponding to the left child state)
					      for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_x.bif_in, bx)
						{
						  for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_y.bif_in, by)
						    {
						      for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_z.bif_in, bz)
							{
							  // no transition score associated with bifurcations in our framework
							  // NB: they may be set to -InfinityScore anyways when e.g. 
							  // Triplet_SCFG::selected_outgoing_states() is called by the
							  // Triplet_SCFG_filter constructor
							  score = max (score, ScorePMul (read_cell (l, bx->l, by->l, bz->l), read_cell (r, bx->r, by->r, bz->r)));
							}
						    }
						}
					    }
					  // emit or otherwise
					  else
					    {
					      // now it's time to perform the summation (max, actually) over possible destination states corresponding 
					      // to the inside subsequence (contained in this one) reached by an emission from the current "source" state (*s)

					      // find the subsequence triplet reached by an emission from the current state (*s)
					      // dest_idx_x is the X subseq index of reachable cells
					      const int dest_idx_x = dest_subseq_idx_x (subseq_idx_x, t);
					      const int dest_idx_y = dest_subseq_idx_y (subseq_idx_y, t);
					      const int dest_idx_z = dest_subseq_idx_z (subseq_idx_z, t);
					      // dest_x is the corresponding X subseq
					      const Subseq& dest_x = foldenv_x.subseq[dest_idx_x];
					      const Subseq& dest_y = foldenv_y.subseq[dest_idx_y];
					      const Subseq& dest_z = foldenv_z.subseq[dest_idx_z];
				      
					      // if the subseqs are empty, account for transitions to end
					      if (dest_x.len == 0 && dest_y.len == 0 && dest_z.len == 0)
						score = scfg.transition_scores.end[*s];
					      else
						score = -InfinityScore;

					      // Multiply by the score corresponding to the parse tree rooted at the 
					      // most likely destination state as well as the associated transition score:
					      // Perform a max operation over all possible destination states
					      // corresponding to the inside subseq pair (dest_x, dest_y, dest_z)
					      // reached by an emission from the current "source" state (*s).
					      // allowed_dest_in_states (state, dest_x, dest_y, dest_z) returns a vector
					      //  of all possible destination states from 'state', given 
					      //  the connection flags associated with the destination subsequence triplet (dest_x, dest_y, dest_z).
					      // Ie, we're considering all cells cell (s', dest_idx_x, dest_idx_y, dest_idx_z) reachable 
					      // from the current cell cell (*s, subseq_idx_x, subseq_idx_y, subseq_idx_z) and maxing over s'.
					      for_const_contents (vector<int>, allowed_dest_in_states (*s, dest_x, dest_y, dest_z), d)
						score = max (score, ScorePMul (read_cell (*d, dest_idx_x, dest_idx_y, dest_idx_z), scfg.transition_scores.transition (*s, *d)));
					      
					      // Multiply by the emit score for the current state.
					      // This factor must come here (rather than at the start of the block).
					      // NB: remember that we're using a formulation where the source state (*s) does the emitting
					      // Also note that we're only looping over states (*s) whose production rules match the connection flags 
					      // for the subsequence triplet (subseq_x, subseq_y, subseq_z), so the emit score emit[*s][...]
					      // is necessarily defined for this state (*s) (it won't in general be defined for an arbitrary emit state).
					      if (is_emit_type (t))
						ScorePMulAcc (score, scfg.emit[*s][in_emit_idx (t, subseq_x, subseq_y, subseq_z)]);

					    }

					  // debugging
#ifdef DART_DEBUG
					  if (((subseq_x.len == 0) && (subseq_y.len == 0) && (subseq_z.len == 0)) || ((subseq_x.len == seqlen_x) && (subseq_y.len == seqlen_y) && (subseq_z.len == seqlen_z)))
					    CL << "     " << *s << " => " << score << "\n";
#endif /* DART_DEBUG */

					  // now store the score
					  cell (*s, subseq_idx_x, subseq_idx_y, subseq_idx_z) = score;

					  // is a start transition allowed for the current subsequence triplet?
					  // because we're doing global alignment, this is only allowed for the maximal subseqs
					  if (local_start_allowed (subseq_x, subseq_y, subseq_z))
					    {
					      const Score f_sc = ScorePMul (score, scfg.transition_scores.start[*s]); // scfg.start[*s] is the score of the transition Start -> (*s)
					      // if this is a new best score, store the relevant information
					      if (f_sc > final_sc)
						{
						  final_sc = f_sc;
						  final_subseq_idx_x = subseq_idx_x;
						  final_subseq_idx_y = subseq_idx_y;
						  final_subseq_idx_z = subseq_idx_z;
						  final_state = *s;
						}
					    }

					} // loop over states

				      // count the cell
				      ++cells_done;
				      bifs_done += subseq_x.bif_in.size() * subseq_y.bif_in.size() * subseq_z.bif_in.size();
				    }
				}
			    } // z loop

			  // log progress
			  if (bifs_done - last_report >= SCFGDP_REPORT_INTERVAL || (bifs_done == num_bifs && bifs_done > last_report))
			    {
			      last_report = bifs_done;
			      CTAG(7,INDIEGRAM_DP) << "CYK: Finished " << cells_done << " subseqs (" << ((int)(1000.*(double)cells_done/(double)num_cells))/10. << "%) and " << bifs_done << " bifurcations (" << ((int)(1000.*(double)bifs_done/(double)num_bifs))/10. << "%)\n";
			    }
			}
		    }
		} // y loop

	    }
	}
    } // x loop

  // all done!  : )

}

// Traceback debugging note: I've checked that this pretty much mirrors Ian's code.

// this performs the same function as Pair_CYK_matrix::traceback()
// I integrate the helper function traceback_from() into this
// because I'm only interested in global tracebacks (for now, at least).

// Remember that CYK fills the DP matrix in an Inside-Outside order.  Traceback
// therefore proceeds in an Outside-Inside order, which means that 
// a qualitative sketch of the algorithm is:
// -Look at the state type of the state of the current cell.
// -If it's a bifurcation, push the child states onto the stack.
// -Otherwise, loop over the possible destination states (child states of the current state)
//  to find the best one (the one whose score).
// We're building the parse tree top-down.  This is nice.
vector<int> Triplet_CYK_matrix::traceback (bool prepend_Start_state /* = true */) const
{
  // this will hold the traceback state path
  vector<int> path;

  // if no successful parse, dump DP matrix and exit
  if (final_sc == -InfinityScore)
    THROWEXPR (scfg_dump() << "Final score is -InfinityScore; traceback is likely to break.\n");

  // set initial subsequences and lengths (corresponding to last cell filled by CYK)
  // NB: We're now doing only global tracebacks, so these are the outermost subsequences.
  const Subseq& subseq_init_x = foldenv_x.subseq[final_subseq_idx_x];
  const Subseq& subseq_init_y = foldenv_y.subseq[final_subseq_idx_y];
  const Subseq& subseq_init_z = foldenv_z.subseq[final_subseq_idx_z];

  // print log message: starting CYK traceback
  CTAG(4,CYK_TRACEBACK)
    << "Starting CYK traceback at X subseq #" << final_subseq_idx_x << ", Y subseq #" << final_subseq_idx_y << ", Z subseq #" << final_subseq_idx_z << ", state #" << final_state << "\n"
    << "X subseq = (" << subseq_init_x.start << '+' << subseq_init_x.len << ") "
    << "Y subseq = (" << subseq_init_y.start << '+' << subseq_init_y.len << ") "
    << "Z subseq = (" << subseq_init_z.start << '+' << subseq_init_z.len << ")\n";

  // prepend start state if requested
  path.clear();
  if (prepend_Start_state)
    path.push_back (Grammar_state_enum::Start);

  // stacks for subseqs and states
  vector<int> subseq_idx_x_stack;
  vector<int> subseq_idx_y_stack;
  vector<int> subseq_idx_z_stack;
  vector<int> state_stack;
  
  // total # of symbols emitted to sequences
  int emitted_x = 0, emitted_y = 0, emitted_z = 0;

  // placekeepers for state and subseqs of current cell in traceback
  int traceback_state;
  int subseq_idx_x, subseq_idx_y, subseq_idx_z;

  // initialize placekeepers to the first cell in traceback (last cell filled in CYK)
  traceback_state = final_state;
  subseq_idx_x = final_subseq_idx_x;
  subseq_idx_y = final_subseq_idx_y;
  subseq_idx_z = final_subseq_idx_z;

  // now do the traceback: for all states in path
  while (1)
    {
      // push current state onto path
      path.push_back (traceback_state);
      
      // check if End state
      if (traceback_state == Grammar_state_enum::End)
	{
	  // if no more states left in the state path, then we're done with the traceback
	  if (state_stack.size() == 0) break;

	  // otherwise, we've just finished a branch, so
	  // get the next subseq triplet and state in the traceback,
	  subseq_idx_x = subseq_idx_x_stack.back();
	  subseq_idx_y = subseq_idx_y_stack.back();
	  subseq_idx_z = subseq_idx_z_stack.back();
	  traceback_state = state_stack.back();
	  // remove them from the stacks,
	  subseq_idx_x_stack.pop_back();
	  subseq_idx_y_stack.pop_back();
	  subseq_idx_z_stack.pop_back();
	  state_stack.pop_back();
	  // and then continue to deal with this next subseq triplet and state
	  continue;
	}

      // get subseqs for the current cell
      const Subseq& subseq_x = foldenv_x.subseq[subseq_idx_x];
      const Subseq& subseq_y = foldenv_y.subseq[subseq_idx_y];
      const Subseq& subseq_z = foldenv_z.subseq[subseq_idx_z];

      // get score and state type for the current cell
      const Score traceback_sc = read_cell (traceback_state, subseq_idx_x, subseq_idx_y, subseq_idx_z);
      const State_type traceback_type = scfg.state_type[traceback_state];

      // update emitted_x, etc.
      if (traceback_type & EmitXL) ++emitted_x;
      if (traceback_type & EmitXR) ++emitted_x;
      if (traceback_type & EmitYL) ++emitted_y;
      if (traceback_type & EmitYR) ++emitted_y;
      if (traceback_type & EmitZL) ++emitted_z;
      if (traceback_type & EmitZR) ++emitted_z;

      // check for bifurcation
      if (is_bifurc_type (traceback_type))
	{
	  const int l = scfg.bifurc[traceback_state].l;
	  const int r = scfg.bifurc[traceback_state].r;
	  // loop over all "bifurcation connections" to look for the correct (best-scoring) subseq pair
	  for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_x.bif_in, bx)
	    {
	      for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_y.bif_in, by)
		{
		  for_const_contents (Subseq::Bifurc_in_pseudovec, subseq_z.bif_in, bz)
		    {
		      if (traceback_sc == ScorePMul (read_cell (l, bx->l, by->l, bz->l), read_cell (r, bx->r, by->r, bz->r)))
			{
			  // we want path to correspond to a pre-order traversal of the parse tree,
			  // so push the subseq info for the right child onto the subseq stacks
			  subseq_idx_x_stack.push_back (bx->r);
			  subseq_idx_y_stack.push_back (by->r);
			  subseq_idx_z_stack.push_back (bz->r);
			  // and then store the left subseq info to start constructing the left subtree
			  subseq_idx_x = bx->l;
			  subseq_idx_y = by->l;
			  subseq_idx_z = bz->l;
			  goto CYK_FoundBifurcation;
			}
		    }
		}
	    }
	  // sanity check: if we're here, then we couldn't find the correct subseq pair
	  THROWEXPR (scfg_dump() << "Triplet SCFG: traceback failed (couldn't find bifurcation) at xi = " << subseq_idx_x << ", yi = " << subseq_idx_y << ", zi = " << subseq_idx_z << ", state = " << traceback_state << " (" << state_type_string (traceback_type) << ")\n");
	CYK_FoundBifurcation:
	  // push the right child onto the state stack
	  state_stack.push_back (r);
	  // and store the left child state to start constructing the left subtree
	  traceback_state = l;
	}
      // if emit or null state
      else
	{
	  // find the destination subsequence triplet reached by an emission from the current state traceback_state
	  // dest_idx_x is the X subseq index of reachable cells
	  const int dest_idx_x = dest_subseq_idx_x (subseq_idx_x, traceback_type);
	  const int dest_idx_y = dest_subseq_idx_y (subseq_idx_y, traceback_type);
	  const int dest_idx_z = dest_subseq_idx_z (subseq_idx_z, traceback_type);
	  // dest_x is the corresponding X subseq
	  const Subseq& dest_x = foldenv_x.subseq[dest_idx_x];
	  const Subseq& dest_y = foldenv_y.subseq[dest_idx_y];
	  const Subseq& dest_z = foldenv_z.subseq[dest_idx_z];

	  // calculate emit score
	  const Score emit_sc =
	    traceback_type == Null
	    ? 0
	    : scfg.emit[traceback_state][in_emit_idx (traceback_type, subseq_x, subseq_y, subseq_z)];

	  // do transitions to End if necesssary (if we've finished a branch)
	  if (dest_x.len == 0 && dest_y.len == 0 && dest_z.len == 0 && traceback_sc == ScorePMul (scfg.transition_scores.end[traceback_state], emit_sc))
	    {
	      traceback_state = Grammar_state_enum::End;
	      goto CYK_FoundEmit;
	    }
	  // otherwise find the next state on the branch
	  else
	    {
	      // loop over outgoing transitions (destination states)
	      for_const_contents (vector<int>, allowed_dest_in_states (traceback_state, dest_x, dest_y, dest_z), d)
		// if we've found the correct destination state
		if (traceback_sc == ScorePMul3 (read_cell (*d, dest_idx_x, dest_idx_y, dest_idx_z), scfg.transition_scores.transition (traceback_state, *d), emit_sc))
		  {
		    traceback_state = *d; // store this as the next state to add to the state path
		    goto CYK_FoundEmit;
		  }
	      // sanity check: if we're here, then we couldn't find the correct destination state
	      THROWEXPR (scfg_dump() << "Triplet SCFG: traceback failed at xi = " << subseq_idx_x << ", yi = " << subseq_idx_y << ", zi = " << subseq_idx_z << ", state = " << traceback_state << " (" << state_type_string (traceback_type) << ")");
	    }
	CYK_FoundEmit:
	  // update the placekeepers
	  subseq_idx_x = dest_idx_x;
	  subseq_idx_y = dest_idx_y;
	  subseq_idx_z = dest_idx_z;
	}
    }
  
  // log messages for finished traceback
  if (CTAGGING(3,CYK_TRACEBACK))
    CL << "Finished CYK traceback: path = (" << path << ")\n"
       << "emitted_x = " << emitted_x << " emitted_y = " << emitted_y << " emitted_z = " << emitted_z << "\n";

  return path;

}


Triplet_SCFG_parse_tree Triplet_CYK_matrix::parse_tree() const
{
  const vector<int> path = traceback (true);
  const Triplet_SCFG_parse_tree parse = scfg.parse (path);
  return parse;
}

Triplet_SCFG_alignment Triplet_CYK_matrix::alignment() const
{
  const Triplet_SCFG_parse_tree parse = parse_tree();
  Triplet_SCFG_alignment align = parse.alignment (scfg.state_type, scfg.state_type_ancestral_map, np_x, np_y, np_z);
  return align;
}


