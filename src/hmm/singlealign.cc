#include "hmm/singlealign.h"

void Single_HMM_alignment::emit_row (const sstring& name, Sequence_database& db)
{
  Named_profile tmp_np;
  db.push_back(tmp_np);
  Named_profile& np = db.back();

  row_name.push_back (np.name = name);
  prof.push_back (&np);
  state_path.push_back (hmm.sample_state_path());

  hmm.sample_sequence (state_path.back(), np.dsq);
  np.dsq2score (hmm.alphabet());

  prof_null_sc.push_back (np.prof_sc.null_score (null_emit_sc));

  total_counts.add_counts_from_state_path (hmm, *prof.back(), state_path.back());
  ++aligned_rows;
}

void Single_HMM_alignment::add_row (const sstring& name, const Named_profile* profile, const vector<int>* path)
{
  row_name.push_back (name);
  prof.push_back (profile);
  prof_null_sc.push_back (profile->prof_sc.null_score (null_emit_sc));
  if (path)
    {
      state_path.push_back (*path);
      total_counts.add_counts_from_state_path (hmm, *profile, *path);
      ++aligned_rows;
    }
  else
    state_path.push_back (vector<int>());
}

void Single_HMM_alignment::set_null_model (const vector<Score>& null_emit, Score null_prior)
{
  null_emit_sc = null_emit;
  null_prior_sc = null_prior;
  model_prior_sc = Prob2Score (1 - Score2Prob (null_prior));
  for (int row = 0; row < rows(); ++row) prof_null_sc[row] = prof[row]->prof_sc.null_score (null_emit);
}

bool Single_HMM_alignment::optimise_row (int row)
{
  Single_Viterbi_matrix matrix (hmm, *prof[row]);
  // compare to null model score
  vector<int> optimal_path;
  if (matrix.final_score + model_prior_sc > prof_null_sc[row] + null_prior_sc)
    optimal_path = matrix.optimal_state_path();
  // update alignment
  if (state_path[row] == optimal_path) return 0;
  if (state_path[row].size())
    {
      Single_HMM_counts row_counts (hmm);
      row_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
      total_counts.subtract_counts (row_counts);
      --aligned_rows;
    }
  state_path[row] = optimal_path;
  total_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
  ++aligned_rows;
  return 1;
}

bool Single_HMM_alignment::sample_row (int row, double kT)
{
  Single_HMM_scores hmm_kT = hmm;
  hmm_kT.scale_all_scores (1/kT);
  Single_forward_matrix matrix (hmm_kT, *prof[row]);
  // choose between this and null model
  vector<int> sampled_path;
  const Score model_evidence_sc = matrix.final_score + model_prior_sc;
  const Score null_evidence_sc = prof_null_sc[row] + null_prior_sc;
  if (Rnd::decide (Score2Prob (model_evidence_sc - (ScorePSum (model_evidence_sc, null_evidence_sc)))))
    sampled_path = matrix.sample_state_path();
  // update counts
  if (state_path[row] == sampled_path) return 0;
  if (state_path[row].size())
    {
      Single_HMM_counts row_counts (hmm);
      row_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
      total_counts.subtract_counts (row_counts);
      --aligned_rows;
    }
  state_path[row] = sampled_path;
  total_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
  ++aligned_rows;
  return 1;
}

bool Single_HMM_alignment::Gibbs_optimise_row (int row, const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts)
{
  if (got_path(row))
    {
      Single_HMM_counts row_counts (hmm);
      row_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
      total_counts.subtract_counts (row_counts);
      --aligned_rows;
    }
  optimise_HMM (mask, pseudocounts);
  Single_Viterbi_matrix matrix (hmm, *prof[row]);
  // compare to null model score
  vector<int> optimal_path;
  if (matrix.final_score + model_prior_sc > prof_null_sc[row] + null_prior_sc)
    optimal_path = matrix.optimal_state_path();
  const bool changed = state_path[row] == optimal_path;
  state_path[row] = optimal_path;
  // update alignment
  total_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
  ++aligned_rows;
  return changed;
}

bool Single_HMM_alignment::Gibbs_sample_row (int row, const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts)
{
  if (got_path(row))
    {
      Single_HMM_counts row_counts (hmm);
      row_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
      total_counts.subtract_counts (row_counts);
      --aligned_rows;
    }
  optimise_HMM (mask, pseudocounts);
  Single_forward_matrix matrix (hmm, *prof[row]);
  // choose between this and null model
  vector<int> sampled_path;
  const Score model_evidence_sc = matrix.final_score + model_prior_sc;
  const Score null_evidence_sc = prof_null_sc[row] + null_prior_sc;
  if (Rnd::decide (Score2Prob (model_evidence_sc - (ScorePSum (model_evidence_sc, null_evidence_sc)))))
    sampled_path = matrix.sample_state_path();
  // update counts
  const bool changed = state_path[row] == sampled_path;
  state_path[row] = sampled_path;
  total_counts.add_counts_from_state_path (hmm, *prof[row], state_path[row]);
  return changed;
}

void Single_HMM_alignment::optimise_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts) const
{
  Single_HMM_counts counts = total_counts;
  counts.add_counts (pseudocounts);
  counts.update_HMM_scores (hmm, mask, 0);
}

void Single_HMM_alignment::sample_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts, double kT) const
{
  Single_HMM_counts counts = total_counts;
  counts.add_counts (pseudocounts);
  counts.update_HMM_scores (hmm, mask, 1, kT);
}


double Single_HMM_alignment::accuracy (const Single_HMM_alignment& reference_alignment) const
{
  if (rows() != reference_alignment.rows()) THROW Standard_exception ("Can't compare two alignments that don't have the same number of rows");
  vector<bool> row_done (rows(), (bool) 0);
  int matches = 0;
  int total = 0;
  for (int ref_row = 0; ref_row < reference_alignment.rows(); ++ref_row)
    {
      int row = find_row_by_name (reference_alignment.row_name[ref_row].c_str());
      if (row_done[row]) THROW Standard_exception ("Two rows in reference alignment mapped onto the same row in test alignment");
      for (int i = 0; i < (int) min (state_path[row].size(), reference_alignment.state_path[ref_row].size()); ++i)
	if (state_path[row][i] == reference_alignment.state_path[row][i])
	  ++matches;
      total += reference_alignment.state_path[ref_row].size();
      row_done[row] = 1;
    }
  return total == 0 ? 0 : ((double) matches) / ((double) total);
}

void Single_HMM_alignment::make_local_alignment (const vector<int>& state_sequence, const sstring& label_sequence, Local_alignment& ret_alignment, sstring& ret_column_labels) const
{
  ret_alignment.reset (rows());
  ret_alignment.row_name = row_name;

  ret_alignment.prof.clear();
  for (int i = 0; i < rows(); ++i) ret_alignment.prof.push_back (&prof[i]->prof_sc);

  ret_column_labels.clear();
  vector<int> sp_cursor (rows(), -1);
  if (state_sequence.size() > 0)
    {
      for (int r = 0; r < rows(); ++r)
	{
	  while ((ret_alignment.start_coord[r] = ++sp_cursor[r]) < (int) state_path[r].size())
	    if (state_path[r][sp_cursor[r]] == state_sequence.front())
	      break;
	}
      for (int i = 0; i < (int) state_sequence.size(); ++i)
	{
	  vector<int>  match_rows (rows());
	  vector<int>  match_rows_in_next_column;
	  for (int r = 0; r < rows(); ++r) match_rows[r] = r;
	  while (match_rows.size())
	    {
	      vector<bool> column (rows(), (bool) 0);
	      bool found_match = 0;
	      for_const_contents (vector<int>, match_rows, r)
		if (sp_cursor[*r] < (int) state_path[*r].size())
		  if (state_path[*r][sp_cursor[*r]] == state_sequence[i])
		    {
		      found_match = column[*r] = 1;
		      if (++sp_cursor[*r] < (int) state_path[*r].size())
			if (state_path[*r][sp_cursor[*r]] == state_sequence[i])
			  match_rows_in_next_column.push_back (*r);
		    }
	      if (found_match)
		{
		  ret_alignment.path.append_column (column);
		  ret_column_labels.push_back (label_sequence[i]);
		}
	      match_rows.swap (match_rows_in_next_column);
	      match_rows_in_next_column.clear();
	    }
	}
    }
}

void Single_HMM_alignment_counts::add_row (const sstring& name, const Named_profile* profile)
{
  row_name.push_back (name);
  prof.push_back (profile);
  counts.push_back (Single_HMM_counts (hmm));
}

void Single_HMM_alignment_counts::reestimate_row_counts (int row, double kT)
{
  total_counts.subtract_counts (counts[row]);
  Single_HMM_scores hmm_kT = hmm;
  hmm_kT.scale_all_scores (1/kT);
  Single_forward_backward_matrix matrix (hmm_kT, *prof[row]);
  counts[row] = matrix.counts;
  total_counts.add_counts (counts[row]);
}

void Single_HMM_alignment_counts::optimise_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts) const
{
  Single_HMM_counts counts = total_counts;
  counts.add_counts (pseudocounts);
  counts.update_HMM_scores (hmm, mask, 0);
}

void Single_HMM_alignment_counts::sample_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts, double kT) const
{
  Single_HMM_counts counts = total_counts;
  counts.add_counts (pseudocounts);
  counts.update_HMM_scores (hmm, mask, 1, kT);
}

Single_HMM_alignment_counts& Single_HMM_alignment_counts::operator= (const Single_HMM_alignment_counts& c)
{
  prof.clear();
  counts.clear();
  prof.reserve (c.rows());
  counts.reserve (c.rows());
  for (int r = 0; r < c.rows(); ++r)
  {
    prof[r] = ((Single_HMM_alignment_counts&) c) . prof[r];
    counts[r] = ((Single_HMM_alignment_counts&) c) . counts[r];
  }
  total_counts = c.total_counts;
  return *this;
}
