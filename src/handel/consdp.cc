#include "handel/consdp.h"

// Transducer_constrained_DP_matrix
void Transducer_constrained_DP_matrix::alloc()
{
  Transducer_DP_base::alloc();

  /* commented out log message because we don't have a tree any more and I might deprecate this source file anyway
  if (CTAGGING(-1,TRANSDUCER_CONSTRAINED_DP TRANSDUCER_DP_CONSTRAINT))
    {
      CL << "Alignment path constraint for transducer DP (before removing unobserved sequences):\n";
      path->show (CL, tree->node_name);
    }
  */

  collapsed_path = Alignment_path (observed_seqs.size(), path->columns());
  for (int observed_row = 0; observed_row < (int) observed_seqs.size(); ++observed_row)
    collapsed_path[observed_row] = (*path)[observed_seqs[observed_row]];
  collapsed_path.erase_empty_columns();

  cell_sc.resize (trans_sc->states(), collapsed_path.columns(), -InfinityScore);
  end_sc = -InfinityScore;

  Sequence_coords coords = collapsed_path.create_seq_coords();
  seq_coords.clear();
  for (int col = 0; col < collapsed_path.columns(); ++col)
    {
      collapsed_path.inc_seq_coords (coords, col);
      seq_coords.push_back (coords);
    }
}


// Transducer_constrained_forward_matrix
void Transducer_constrained_forward_matrix::fill()
{
  /* commented out log message because we don't have a tree any more and I might deprecate this source file anyway
  if (CTAGGING(2,TRANSDUCER_CONSTRAINED_DP TRANSDUCER_DP_CONSTRAINT))
    {
      CL << "Alignment path constraint for transducer DP:\n";
      vector<sstring> row_name (collapsed_path.rows());
      for (int i = 0; i < (int) observed_seqs.size(); ++i)
	row_name[i] = tree->node_specifier (observed_seqs[i]);
      collapsed_path.show (CL, row_name);
      CL << "Tree for transducer DP:\n";
      tree->write (CL);
    }
  */

  if (CTAGGING(2,TRANSDUCER_CONSTRAINED_DP))
    {
      CL << "Transducer for constrained DP:\n";
      trans_sc->show (CL);
    }

  // loop over columns of constraint alignment
#define TRANSDPTAGS TRANSDUCER_CONSTRAINED_DP TRANSDUCER_CONSTRAINED_DP_MATRIX
  const bool logging = CTAGGING(2,TRANSDPTAGS);
  if (logging)
    CL << "Constrained DP matrix:\n";
  for (int col = 0; col < columns(); ++col)
    {
      // create a string to hold DP matrix log message
      sstring log_string;

      // loop over states, testing for consistency with observed column
      for (int dest = 0; dest < states(); ++dest)
	{
	  const State_type dest_type = trans_sc->state_type[dest];
	  bool consistent = true;
	  for (int row = 0; row < sequences(); ++row)
	    if (collapsed_path(row,col) != type_node_emit (dest_type, observed_seqs[row]))
	      {
		consistent = false;
		break;
	      }
	  if (consistent)
	    {
	      // get address of cell
	      Score& dest_cell = cell_sc (dest, col);
	      // handle the first column separately
	      if (col == 0)
		dest_cell = trans_sc->transition (Grammar_state_enum::Start, dest);
	      else
		for_const_contents (vector<int>, incoming[dest], src)
		  ScorePSumAcc (dest_cell, ScorePMul (cell_sc (*src, col - 1), trans_sc->transition (*src, dest)));
	      // add emit score
	      const Score emit_sc = calc_cell_emit (dest_type, seq_coords[col]);
	      ScorePMulAcc (dest_cell, emit_sc);
	      // log
	      if (logging && dest_cell > -InfinityScore)
		log_string << ' ' << dest << '{' << trans_sc->state_name[dest] << "}[" << dest_cell << ']';
	    }
	}
      if (logging)
	CTAG(2,TRANSDPTAGS) << "Column #" << col
			    << ", seq_coords (" << seq_coords[col]
			    << "):" << log_string
			    << '\n';
    }

  // do end cell
  if (columns() == 0)
    {
      end_sc = trans_sc->start_to_end();
      for (int s = 0; s < states(); ++s)
	if (trans_sc->state_type[s] == 0)  // loop over wait states
	  ScorePSumAcc (end_sc, ScorePMul (trans_sc->transition (Grammar_state_enum::Start, s),
					   trans_sc->transition (s, Grammar_state_enum::End)));
    }
  else
    for (int src = 0; src < states(); ++src)
      ScorePSumAcc (end_sc, ScorePMul (cell_sc (src, columns() - 1), trans_sc->transition (src, Grammar_state_enum::End)));
  if (logging)
    CTAG(2,TRANSDPTAGS) << "End: " << end_sc << '\n';
}

vector<int> Transducer_constrained_forward_matrix::sample_traceback()
{
  if (end_sc <= -InfinityScore)
    THROWEXPR ("Tried to sample traceback when Forward likelihood is zero");

  vector<int> trace;
  int dest_state = Grammar_state_enum::End;
  for (int src_col = columns() - 1; src_col >= 0; --src_col)
    {
      // add current state to traceback
      trace.push_back (dest_state);

      // calculate emit likelihood
      Score emit_sc = 0;
      if (dest_state != Grammar_state_enum::End)
	{
	  const int dest_col = src_col + 1;
	  emit_sc = calc_cell_emit (trans_sc->state_type[dest_state], seq_coords[dest_col]);
	}

      // calculate incoming log-likelihoods
      vector<Score> sc;
      for (int src_state = 0; src_state < states(); ++src_state)
	sc.push_back (ScorePMul3 (cell_sc (src_state, src_col),
				  trans_sc->transition (src_state, dest_state),
				  emit_sc));

      // normalize incoming log-likelihoods
      Score total_sc = -InfinityScore;
      for_const_contents (vector<Score>, sc, s)
	ScorePSumAcc (total_sc, *s);

      // sample next state
      Prob p = Rnd::prob();
      for (dest_state = 0; dest_state < states(); ++dest_state)
	if ((p -= Score2Prob (ScorePMul (sc[dest_state], -total_sc))) <= 0)
	  break;
    }

  // add final(first) & Start states; reverse trace; return
  trace.push_back (dest_state);
  trace.push_back (Grammar_state_enum::Start);
  reverse (trace.begin(), trace.end());
  return trace;
}
