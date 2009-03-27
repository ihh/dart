#include <iomanip>

#include "kimono/gibbs.h"
#include "util/math_fn.h"

Gibbs_dataset::Gibbs_dataset (const Sequence_database_index& db_index, const Sequence_database_index& revcomp_index, const vector<sstring>& row_name, const Alphabet& alphabet)
  : db_index (db_index),
    revcomp_index (revcomp_index),
    profile (row_name.size(), (const Digitized_biosequence*) 0),
    revcomp(),
    profile_seq (row_name.size(), (const Biosequence*) 0),
    revcomp_seq(),
    row_name (row_name)
{
  CLOG(6) << "Cross-referencing sequences by name\n";
  
  if (revcomp_index.size() > 0)
    {
      revcomp = vector<const Digitized_biosequence*> (row_name.size(), (const Digitized_biosequence*) 0);
      revcomp_seq = vector<const Biosequence*> (row_name.size(), (const Biosequence*) 0);
    }

  for (int i = 0; i < (int) row_name.size(); ++i)
    {
      int j = db_index.name2profile_index.lookup (row_name[i]);
      if (j < 0)
	CLOG(8) << "Warning - sequence '" << row_name[i] << "' not found in sequence database\n";
      else
	{
	  profile[i] = &db_index.profile[j]->dsq;
	  profile_seq[i] = &db_index.profile[j]->seq;
	  if (revcomp_index.size() > 0)
	    {
	      revcomp[i] = &revcomp_index.profile[j]->dsq;
	      revcomp_seq[i] = &revcomp_index.profile[j]->seq;
	    }
	}
    }
}

int Gibbs_dataset::row_index (const sstring& name) const
{
  vector<sstring>::const_iterator f = find (row_name.begin(), row_name.end(), name);
  return f == row_name.end() ? -1 : f - row_name.begin();
}

Gibbs_null_model::Gibbs_null_model (const Gibbs_dataset& dataset, const Alphabet& alphabet, double composition_pseudocount, double motif_del_prob, int context_order)
  : dataset (dataset),
    context_order (context_order),
    context_string ((int) pow ((double) alphabet.size(), (double) (context_order + 1))),
    profile_context (dataset.profile.size()),
    profile_rev_context (dataset.profile.size()),
    composition (alphabet.size(), composition_pseudocount),
    composition_context ((int) pow ((double) alphabet.size(), (double) (context_order + 1))),
    context_null_sc ((int) pow ((double) alphabet.size(), (double) (context_order + 1))),
    profile_null_score (dataset.profile.size()),
    motif_delete_probability (motif_del_prob)
{
  CLOG(6) << "Making null model for sequences\n";

  if (motif_delete_probability == 0) CLOGERR << "Warning! setting motif_delete_probability==0 will probably screw things up; try e.g. 1e-100 instead\n";
  
  double total_seq_len = 0;
  for (int i = 0; i < (int) dataset.profile.size(); ++i)
    if (dataset.includes_row (i))
      {
	const Digitized_biosequence& dsq = *dataset.profile[i];
	dsq.count_symbol_frequencies (composition);

	profile_context[i] = dsq.make_context_dsq (alphabet.size(), context_order);
	profile_rev_context[i] = dsq.make_reverse_context_dsq (alphabet.size(), context_order);
	
	profile_context[i].count_symbol_frequencies (composition_context);
	profile_rev_context[i].count_symbol_frequencies (composition_context);

	if (dataset.revcomp_enabled())
	  {
	    const Digitized_biosequence& rev_dsq = *dataset.revcomp[i];
	    rev_dsq.count_symbol_frequencies (composition);
	    Digitized_biosequence tmp_dsq = rev_dsq.make_context_dsq (alphabet.size(), context_order);
	    Digitized_biosequence tmp_rev_dsq = rev_dsq.make_reverse_context_dsq (alphabet.size(), context_order);
	    tmp_dsq.count_symbol_frequencies (composition_context);
	    tmp_rev_dsq.count_symbol_frequencies (composition_context);
	  }

	if (CLOGGING(-2)) {
	  CL << "dsq: (" << dsq << ")\n";
	  CL << "context: (" << profile_context[i] << ")\n";
	  CL << "rev_context: (" << profile_rev_context[i] << ")\n";
	}

	total_seq_len += dsq.size();
      }

  NormalisePr (composition);

  // normalise composition_context cliques
  for (int prev_context = 0; prev_context < (int) composition_context.size(); prev_context += alphabet.size()) {
    double norm = 0;
    for (int c = 0; c < alphabet.size(); ++c) norm += composition_context[prev_context + c];
    if (norm == 0) norm = 1;
    for (int c = 0; c < alphabet.size(); ++c) composition_context[prev_context + c] /= norm;
  }
  context_null_sc = Prob2ScoreVec (composition_context);

  mean_sequence_length = total_seq_len / (double) dataset.profile.size();

  for (int context = 0; context < (int) composition_context.size(); ++context) {
    Digitized_biosequence nmer = Digitized_biosequence::context_to_nmer (context, alphabet.size(), context_order);
    Biosequence nmer_seq;
    alphabet.dsq2seq (nmer, nmer_seq);
    context_string[context] << nmer_seq;
  }

  if (CLOGGING(5)) {
    CL << "Composition = (" << composition << ")\n";
    CL << "Order " << context_order <<  " Markov chain = (";
    vector<sstring> cc_string = context_string;
    for (int context = 0; context < (int) composition_context.size(); ++context)
      cc_string[context] << ":" << composition_context[context];
    CL << cc_string << ")\n";
    CL << "Mean sequence length = " << mean_sequence_length << "\n";
  }
  
  for (int i = 0; i < (int) dataset.profile.size(); ++i)
    {
      if (dataset.includes_row (i))
	profile_null_score[i] = (profile_context[i].null_score (context_null_sc) + profile_rev_context[i].null_score (context_null_sc)) / 2;
      else
	profile_null_score[i] = 0;
      if (CLOGGING(4)) CL << "Row '" << dataset.row_name[i] << "' null score = " << Score2Bits (profile_null_score[i]) << " bits\n";
    }
}

int Gibbs_alignment::left_col() const
{
  int left_col = 0;
  bool first = 1;
  for (int row = 0; row < rows(); ++row)
    if (row_aligned[row])
      if (first || start_col(row) < left_col)
	{
	  left_col = start_col(row);
	  first = 0;
	}
  return left_col;
}

int Gibbs_alignment::right_col() const
{
  int right_col = 0;
  bool first = 1;
  for (int row = 0; row < rows(); ++row)
    if (row_aligned[row])
      if (first || end_col(row) > right_col)
	right_col = end_col(row);
  return right_col;
}

void Gibbs_alignment::optimise_motif()
{
  CLOG(3) << "Optimising motif\n";

  if (motif_length > columns())
    {
      CLOG(3) << "Warning - alignment too narrow for motif\n";
      motif_start = left_col();
      motif_end = right_col();
      recalculate_counts();
      return;
    }

  if (aligned_row_count < 2)
    {
      motif_start = left_col();
      motif_end = motif_start + motif_length - 1;
      recalculate_counts();
      return;
    }

  vector<int> end_sc_matrix;
  calculate_motif_matrix (end_sc_matrix);

  motif_start = left_col() + max_element (end_sc_matrix.begin(), end_sc_matrix.end()) - end_sc_matrix.begin();
  motif_end = motif_start - motif_length + 1;
  recalculate_counts();

  if (CLOGGING(1)) display_motif_matrix (CL, end_sc_matrix);
}

void Gibbs_alignment::sample_motif (double kT)
{
  CLOG(3) << "Sampling motif\n";
  if (CLOGGING(-1)) { CL << "Old motif:\n"; display(CL); CL << "left_col=" << left_col() << " motif_start=" << motif_start << "\n"; }

  if (motif_length > columns())
    {
      CLOG(3) << "Warning - alignment too narrow for motif\n";
      motif_start = left_col();
      motif_end = right_col();
      recalculate_counts();
      return;
    }

  if (aligned_row_count < 2)
    {
      CLOG(3) << "Choosing random motif position\n";
      motif_start = Rnd::rnd_int (columns() - motif_length + 1);
      motif_end = motif_start + motif_length - 1;
      recalculate_counts();
      return;
    }
  
  vector<int> end_sc_matrix;
  calculate_motif_matrix (end_sc_matrix);
  
  motif_start = left_col() + Rnd::choose (Score2BoltzmannVec (end_sc_matrix, kT));
  motif_end = motif_start + motif_length - 1;
  recalculate_counts();

  if (CLOGGING(1)) display_motif_matrix (CL, end_sc_matrix);
}

void Gibbs_alignment::calculate_motif_matrix (vector<int>& matrix) const
{
  const int cols = columns() - motif_length + 1;
  const int lc = left_col();
  matrix = vector<int> (cols);
  vector<int> recent_match_sc (motif_length);
  int end_sc = 0;
  vector<double> sym_freq;
  vector<double> context_freq;
  double total_freq;
  for (int m = 0; m < columns(); ++m)
    {
      const int col = m + lc;
      recalculate_column_counts (col, sym_freq, context_freq, total_freq);
      const double gap_count = total_aligned_in_prob - total_freq;
      const double total_gap_score = max ((double) -InfinityScore, gap_count * (double) delete_extend_score);
      const int match_sc = ScorePMul (column_score (sym_freq, context_freq), (int) total_gap_score);
      ScorePMulAcc (end_sc, match_sc);
      if (m >= motif_length) ScorePMulAcc (end_sc, -recent_match_sc[m % motif_length]);
      if (m >= motif_length - 1) matrix [m - motif_length + 1] = ScorePMul (end_sc, 2 * delete_end_score);
      recent_match_sc[m % motif_length] = match_sc;
    }
}

void Gibbs_alignment::display_motif_matrix (ostream& o, const vector<int>& matrix, double display_threshold) const
{
  o << "Start  End    L/bits Alignment\n";
  int old_prec = o.precision (2);
  save_flags (o);
  right_align (o);
  vector<double> sym_freq;
  vector<double> context_freq;
  double total_freq;
  for (int m = 0; m < columns(); ++m)
    {
      int col = m + left_col();

      if (m < motif_length - 1) o << "    "; else o << setw(4) << col - motif_length + 1;
      o << (col == motif_end ? '*' : ' ');
      o << setw(4) << col;
      o << (col == motif_end ? '*' : ' ');
      
      o << setw(9);
      if (m < motif_length - 1) o << ""; else o << Score2Bits (matrix [m - motif_length + 1]);
      o << " ";
      
      for (int row = 0; row < rows(); ++row)
	if (row_aligned[row] && in_prob[row] >= display_threshold)
	  o << (not_gap(row,col) ? residue_char(row,col) : '.');

      recalculate_column_counts (col, sym_freq, context_freq, total_freq);
      o << (col == motif_end ? '*' : ' ');
      const double gap_count = total_aligned_in_prob - total_freq;
      const double total_gap_score = max ((double) -InfinityScore, gap_count * (double) delete_extend_score);
      
      o << " column_score = " << column_score (sym_freq, context_freq);
      o << " gap_score = " << total_gap_score;

      vector<sstring> cfreq_str = null_model.context_string;
      for (int context = 0; context < (int) cfreq_str.size(); ++context)
	cfreq_str[context] << ":" << context_freq[context];
      o << " sym_freq=(" << sym_freq << ") context_freq=(" << cfreq_str << ")";

      o << "\n";
    }
  o.precision (old_prec);
  restore_flags (o);
}

void Gibbs_alignment::unalign_row (int row)
{
  if (!row_aligned[row]) return;     // bail now if already unaligned
  --aligned_row_count;

  row_aligned[row] = 0;

  const double row_weight = in_prob[row];
  for (int col = max (start_col(row), motif_start); col <= min (end_col(row), motif_end); ++col)
    add_counts (residue(row,col), residue_context(row,col), residue_rev_context(row,col), col, -row_weight);
  total_aligned_in_prob -= row_weight;

  if (CLOGGING(1))
    {
      CL << "Removed row #" << row << " (" << dataset.row_name[row] << "), weight " << in_prob[row] << "; new counts:\n";
      row_aligned[row] = 1;          // hack so that sequence will display properly
      display_counts (CL, row);
      row_aligned[row] = 0;          // restore the proper state of affairs
      if (CLOGGING(0))  display_column_log_odds_scores (CL);
      if (CLOGGING(-1)) verify_counts (CL);
    }
}

void Gibbs_alignment::align_row (int row, int offset, bool reversed)
{
  if (!dataset.includes_row (row)) THROW Standard_exception ("Attempt to align a sequence that doesn't exist");

  if (!row_aligned[row]) ++aligned_row_count;

  row_offset[row] = offset;
  row_reversed[row] = reversed;
  row_aligned[row] = 1;

  // add counts for this row to col_sym_freq[][]
  //
  double row_weight = in_prob[row];
  for (int col = max (start_col(row), motif_start); col <= min (end_col(row), motif_end); ++col)
    add_counts (residue(row,col), residue_context(row,col), residue_rev_context(row,col), col, row_weight);
  total_aligned_in_prob += row_weight;

  if (CLOGGING(1))
    {
      CL << "Added row #" << row << " (" << dataset.row_name[row] << "), weight " << row_weight << "; new counts:\n";
      display_counts (CL, row);
      if (CLOGGING(0))  display_column_log_odds_scores (CL);
      if (CLOGGING(-1)) verify_counts (CL);
    }
}

void Gibbs_alignment::calculate_DP_matrix (const Digitized_biosequence& dsq,
					   const Digitized_biosequence& fwd_context,
					   const Digitized_biosequence& rev_context,
					   vector<int>::iterator matrix_begin,
					   int matrix_size,
					   double weight) const
{
  const int sequence_length = dsq.size();

  if (CLOGGING(-1)) display_motif_scores (CL);
  if (CLOGGING(-2)) display_motif_counts (CL);

  // row_offset can run from (motif_start - sequence_length) to (motif_end + 1)
  //
  for (int matrix_col = 0; matrix_col < matrix_size; ++matrix_col)
    {
      int& score = *(matrix_begin + matrix_col);

      const int first_seq_col = matrix_col + (motif_start - (sequence_length - 1));           // column index of the first residue in the sequence
      const int last_seq_col  = first_seq_col + sequence_length - 1;                          // column index of the last residue in the sequence
      
      const int ldel = max (min (motif_end+1, first_seq_col) - motif_start, 0);               // deleted residues on LHS of motif
      const int rdel = max (motif_end+1 - max (last_seq_col+1, motif_start), 0);              // deleted residues on RHS of motif
      
      score = (ldel + rdel) * delete_extend_score;
      
      const int first_overlapping_motif_pos = max (motif_start, first_seq_col) - motif_start; // first motif index that overlaps with the sequence
      const int last_overlapping_motif_pos  = min (motif_end, last_seq_col) - motif_start;    // last motif index that overlaps
      const int motif_start_seq_pos         = motif_start - first_seq_col;                    // the residue in the sequence that's aligned with the first position in the motif
      
      for (int motif_pos = first_overlapping_motif_pos; motif_pos <= last_overlapping_motif_pos; ++motif_pos)
	{
	  const int dsq_pos = motif_pos + motif_start_seq_pos;
	  ScorePMulAcc (score, motif_score (dsq[dsq_pos],
								 fwd_context[dsq_pos],
								 rev_context[dsq_pos],
								 motif_pos + motif_start,
								 weight));
	}
    }
}

void Gibbs_alignment::display_DP_matrix (ostream& o, const vector<int>& matrix, int row, int highlight_col) const
{
  o << "Alignment log-likelihoods for '" << dataset.row_name[row] << "':\n";
  int old_prec = o.precision (2);
  save_flags (o);
  left_align (o);

  Digitized_biosequence cons_dsq = alphabet_consensus();
  Biosequence cons_seq;
  alphabet.dsq2seq (cons_dsq, cons_seq);

  int forward_size = DP_matrix_size (row);

  sstring banner;
  banner << "+/- Offset Log-like/bits " << cons_seq << "\n";
  for (int i = 0; i < (int) matrix.size(); ++i)
    {
      const int matrix_col = i < forward_size ? forward_size-1 - i : matrix.size()-1 - i + forward_size;     // going backwards looks better (as though a sliding window is moving right)
      if (i % 20 == 0) { if (i > 0) o << "\n"; o << banner; }
      char c = (matrix_col == highlight_col ? '*' : ' ');
      bool rev = matrix_col >= forward_size;
      o << c << (rev ? '-' : '+') << c << " ";
      right_align (o);
      const int motif_start_seq_pos = motif_start - convert_matrix_column_to_row_offset (matrix_col % forward_size, row);
      const Biosequence& row_seq = *(rev ? dataset.revcomp_seq[row] : dataset.profile_seq[row]);
      o << setw(6) << motif_start_seq_pos << " " << setw(12) << Score2Bits (matrix[matrix_col]) << "  ";
      for (int j = motif_start_seq_pos; j < motif_start_seq_pos + (int) cons_dsq.size(); ++j)
	o << (j < 0 || j >= (int) row_seq.size() ? '.' : row_seq[j]);
      o << "\n";
    }
  o << "\n";
  
  restore_flags (o);
  o.precision (old_prec);
}

Digitized_biosequence Gibbs_alignment::alphabet_consensus() const
{
  Digitized_biosequence dsq (motif.size());
  for (int i = 0; i < (int) motif.size(); ++i)
    dsq[i] = max_element (motif[i].begin(), motif[i].end()) - motif[i].begin();
  return dsq;
}

sstring Gibbs_alignment::mixture_consensus() const
{
  sstring consensus;

  vector<int> optimal_sc;
  vector<double> log_dirichlet_evidence;
  double log_total_dirichlet_evidence;
  
  for (int m = 0; m < motif_length; ++m) {
    calculate_optimal_scores (col_sym_freq[m], optimal_sc, log_dirichlet_evidence, log_total_dirichlet_evidence);
    int max_component = max_element (log_dirichlet_evidence.begin(), log_dirichlet_evidence.end()) - log_dirichlet_evidence.begin();
    const sstring& max_cpt_name = dirichlet_component_name[max_component];
    if (max_cpt_name.size() > 1) consensus << "[" << max_cpt_name << "]";
    else consensus << max_cpt_name;
  }
  return consensus;
}

bool Gibbs_alignment::sample_row (int row, double kT)
{
  if (!dataset.includes_row (row))
    {
      if (CLOGGING(3)) CL << "No data for sequence '" << dataset.row_name[row] << "' - skipping alignment sampling step\n";
      return 0;
    }
  else
    if (CLOGGING(3)) CL << "Sampling alignment of sequence '" << dataset.row_name[row] << "'\n";

  if (aligned_row_count == 0)
    {
      align_row (row, 0, 0);
      sample_motif();
      return 1;
    }
  
  bool row_was_aligned = row_aligned[row];
  int  old_offset      = row_offset[row];
  bool old_reversed    = row_reversed[row];

  unalign_row (row);

  int forward_size = DP_matrix_size (row);
  vector<int> matrix ((revcomp_enabled() ? 2 : 1) * forward_size);
  calculate_DP_matrix (*dataset.profile[row],
		       null_model.profile_context[row],
		       null_model.profile_rev_context[row],
		       matrix.begin(),
		       forward_size,
		       in_prob[row]);
  if (revcomp_enabled())
    calculate_DP_matrix (*dataset.revcomp[row],
			 null_model.profile_rev_context[row],
			 null_model.profile_context[row],
			 matrix.begin() + forward_size,
			 forward_size,
			 in_prob[row]);

  int sample_col = Rnd::choose (Score2BoltzmannVec (matrix, kT));
  if (CLOGGING(2)) display_DP_matrix (CL, matrix, row, sample_col);
  
  bool reversed = sample_col >= forward_size;
  sample_col = sample_col % forward_size;

  int new_offset = convert_matrix_column_to_row_offset (sample_col, row);
  align_row (row, new_offset, reversed);
  return !(row_was_aligned && new_offset == old_offset && reversed == old_reversed);
}

bool Gibbs_alignment::optimise_row (int row)
{
  if (!dataset.includes_row (row))
    {
      if (CLOGGING(3)) CL << "No data for sequence '" << dataset.row_name[row] << "' - skipping alignment optimisation step\n";
      return 0;
    }
  else
    if (CLOGGING(3)) CL << "Optimising alignment of sequence '" << dataset.row_name[row] << "'\n";

  if (aligned_row_count == 0)
    {
      align_row (row, 0, 0);
      optimise_motif();
      return 1;
    }

  bool row_was_aligned = row_aligned[row];
  int  old_offset      = row_offset[row];
  bool old_reversed    = row_reversed[row];

  unalign_row (row);

  int forward_size = DP_matrix_size (row);
  vector<int> matrix ((revcomp_enabled() ? 2 : 1) * forward_size);
  calculate_DP_matrix (*dataset.profile[row],
		       null_model.profile_context[row],
		       null_model.profile_rev_context[row],
		       matrix.begin(),
		       forward_size,
		       in_prob[row]);
  if (revcomp_enabled())
    calculate_DP_matrix (*dataset.revcomp[row],
			 null_model.profile_rev_context[row],
			 null_model.profile_context[row],
			 matrix.begin() + forward_size,
			 forward_size,
			 in_prob[row]);
  
  int max_col = max_element (matrix.begin(), matrix.end()) - matrix.begin();
  if (CLOGGING(2)) display_DP_matrix (CL, matrix, row, max_col);

  bool reversed = max_col >= forward_size;
  max_col = max_col % forward_size;

  int new_offset = convert_matrix_column_to_row_offset (max_col, row);
  align_row (row, new_offset, reversed);
  return !(row_was_aligned && new_offset == old_offset && reversed == old_reversed);
}

Gibbs_alignment::Gibbs_alignment (const Gibbs_dataset& dataset,
				  const Alphabet& alphabet,
				  const Gibbs_null_model& null_model,
				  const vector<double>& log_dirichlet_component_weight,
				  const vector<vector<double> >& dirichlet_mixture,
				  const vector<sstring>& dirichlet_component_name,
				  int motif_length,
				  double motif_delete_probability)
  : dataset (dataset),
    alphabet (alphabet),
    null_model (null_model),
    context_null_sc (null_model.context_null_sc),
    profile_null_score (null_model.profile_null_score),
    log_dirichlet_component_weight (log_dirichlet_component_weight),
    dirichlet_mixture (dirichlet_mixture),
    dirichlet_component_name (dirichlet_component_name),
    motif_length (motif_length),
    delete_extend_score (Prob2Score (motif_delete_probability)),
    delete_end_score (Prob2Score (1.0 - motif_delete_probability)),
    row_aligned (rows(), 0),
    row_offset (rows()),
    row_reversed (rows()),
    aligned_row_count (0),
    motif_start (0),
    motif_end (motif_length - 1),
    col_sym_freq (motif_length, vector<double> (alphabet.size(), (double) 0)),
    col_context_freq (motif_length, vector<double> (context_null_sc.size(), (double) 0)),
    total_col_freq (motif_length, (double) 0),
    motif (motif_length, vector<int> (alphabet.size(), (int) 0)),
    in_sc (rows(), (int) -InfinityScore),
    in_prob (rows(), (double) 0),
    total_aligned_in_prob (0)
{
  return;
}

void Gibbs_alignment::set_membership_probability (int row, int score)
{
  // change row weight
  //
  double old_weight = in_prob[row];
  double new_weight = Score2Prob (score);
  in_sc[row] = score;
  in_prob[row] = new_weight;

  // change col_sym_freq[][] appropriately
  //
  if (row_aligned[row])
    {
      for (int col = max (start_col(row), motif_start); col <= min (end_col(row), motif_end); ++col)
	add_counts (residue(row,col), residue_context(row,col), residue_rev_context(row,col), col, new_weight - old_weight);
      total_aligned_in_prob += new_weight - old_weight;
    }

  if (CLOGGING(1))
    {
      CL << "Changed weight of row #" << row << " (" << dataset.row_name[row] << ") from " << old_weight << " to " << new_weight << "; new counts:\n";
      display_counts (CL, row);
      if (CLOGGING(0))  display_column_log_odds_scores (CL);
      if (CLOGGING(-1)) verify_counts (CL);
    }
}

int Gibbs_alignment::log_likelihood_sc (int row) const
{
  if (!row_aligned[row]) return -InfinityScore;
  int score = 0;
  score += delete_extend_score * left_delete (row) + delete_end_score;
  score += delete_extend_score * right_delete (row) + delete_end_score;
  for (int col = motif_start; col <= motif_end; ++col)
    if (not_gap (row, col))
      ScorePMulAcc (score, motif_score (residue(row,col),
							     residue_context(row,col),
							     residue_rev_context(row,col),
							     col));
  return ScorePMul (score, profile_null_score[row]);
}

int Gibbs_alignment::log_parameter_prior_sc() const
{
  int sc = -InfinityScore;
  for (int col = motif_start; col <= motif_end; ++col)
    {
      vector<int> motif_sc = motif [col-motif_start];
      vector<double> emit_prob = Score2ProbVecNorm (motif_sc);
      for (int i = 0; i < (int) dirichlet_mixture.size(); ++i)
	ScorePSumAcc (sc, Nats2Score (log_dirichlet_component_weight[i] + Math_fn::log_dirichlet (emit_prob, dirichlet_mixture[i])));
    }
  return sc;
}

void Gibbs_alignment::display (ostream& o, double threshold, int columns)
{
  int old_prec = o.precision(3);
  save_flags (o);
  left_align (o);

  int align_cols = columns - 28;
  if (align_cols <= 0) THROW Standard_exception ("Display too narrow");
  for (int col = motif_start; col <= motif_end; col += align_cols)
    {
      int e = min (col + align_cols, motif_end + 1);
      for (int row = 0; row < rows(); ++row)
	if (row_aligned[row] && in_prob[row] >= threshold)
	  {
	    sstring n = dataset.row_name[row];
	    if (n.size() > 10) n.erase (n.begin() + 10, n.end());
	    o << setw(10) << n << " ";
	    if (revcomp_enabled()) o << (row_reversed[row] ? "(-) " : "(+) ");
	    o.width(4);
	    right_align (o);
	    if (e < start_col(row) || col > end_col(row)) o << ""; else o << max (col - row_offset[row] + 1, 1);
	    left_align (o);
	    o << " ";
	    for (int c = col; c < e; ++c) o << (not_gap (row, c) ? residue_char (row, c) : '.');
	    o << " ";
	    o.width(4);
	    if (e < start_col(row) || col > end_col(row)) o << ""; else o << min (e - row_offset[row], residues(row));
	    right_align (o);
	    o << " " << setw(3) << (int) (100.0 * in_prob[row]) << "%\n";
	    left_align (o);
	  }
      o << "\n";
    }

  restore_flags (o);
  o.precision (old_prec);
}

void Gibbs_alignment::display_scores (ostream& o, double threshold)
{
  for (int row = 0; row < rows(); ++row)
    if (row_aligned[row] && in_prob[row] >= threshold)
      {
	o << "Sequence '" << dataset.row_name[row] << "' log-likelihood = " << Score2Bits (log_likelihood_sc (row)) << " bits";
	o << "  (null = " << Score2Bits (null_model.profile_null_score[row]) << " bits)\n";
      }
  o << "Log-parameter prior = " << Score2Bits(log_parameter_prior_sc()) << "\n";
}

void Gibbs_alignment::display_counts (ostream& o, int row)
{
  int old_prec = o.precision(1);
  int old_fill = o.fill();
  save_flags (o);
  right_align (o);
  
  o << "col:";
  for (int col = motif_start; col <= motif_end; ++col) o << " " << setw(4) << col;
  o << "\n";
  for (int sym = 0; sym < alphabet.size(); ++sym)
    {
      o << "  " << alphabet.int2char(sym) << ":";
      for (int col = motif_start; col <= motif_end; ++col)
	{
	  if (col == motif_start) o << "<" << setfill('<');
	  else if (col == motif_end+1) o << ">" << setfill('>');
	  else if (col > motif_start && col <= motif_end) o << "^" << setfill('^');
	  else o << " " << setfill(' ');
	  o << setw(4) << col_sym_freq[col - motif_start][sym];
	}
      if (motif_end == right_col()) o << ">>";
      o << setfill(' ') << "\n";
    }
  if (row < 0 ? 0 : row_aligned[row])
    {
      o << "seq:";
      for (int col = motif_start; col < start_col(row); ++col) o << setw(5) << "";
      o << "    " << residue_char (row, start_col(row));
      for (int col = start_col(row)+1; col <= min (end_col(row), motif_end); ++col) o << "----" << residue_char(row,col);
      o << "  (" << dataset.row_name[row] << ")\n";
    }
  
  restore_flags (o);
  o.fill (old_fill);
  o.precision (old_prec);
}

void Gibbs_alignment::display_column_log_odds_scores (ostream& o)
{
  int old_prec = o.precision(2);
  save_flags (o);
  right_align (o);
  
  /*
  o << "bits";
  for (int col = motif_start; col <= motif_end; ++col)
    {
      const int m = col - motif_start;
      int sc = 0;
      for (int sym = 0; sym < alphabet.size(); ++sym)
	{
	  ScorePMulAcc (sc, (int) (col_sym_freq[m][sym] / (total_col_freq[m] == 0 ? 1 : total_col_freq[m]) * (double) motif_score (sym, -1, col, 0)));
	  ScorePMulAcc (sc, (int) (stut_col_sym_freq[m][sym] / (total_col_freq[m] == 0 ? 1 : total_col_freq[m]) * (double) motif_score (sym, sym, col, 0)));
	}
      o << " " << setw(4) << Score2Bits (sc);
    }
  o << "\n";
  */
  o << "WARNING -- display_column_log_odds_scores skipped\n";
  
  restore_flags (o);
  o.precision (old_prec);
}

void Gibbs_alignment::recalculate_counts()
{
  col_sym_freq = vector<vector<double> > (motif_length);
  col_context_freq = vector<vector<double> > (motif_length);
  total_col_freq = vector<double> (motif_length);
  motif = vector<vector<int> > (motif_length);
  for (int m = 0; m < motif_end + 1 - motif_start; ++m) {
    recalculate_column_counts (m + motif_start, col_sym_freq[m], col_context_freq[m], total_col_freq[m]);
    recalculate_motif_column (m + motif_start);
  }
}

void Gibbs_alignment::verify_counts (ostream& o, double tol)
{
  /*
  int old_prec = o.precision(1);
  save_flags (o);
  right_align (o);
  
  o << "vrfy";
  for (int col = motif_start; col <= motif_end; ++col)
    {
      vector<double> tally = pseudocount;
      for (int row = 0; row < rows(); ++row)
	if (not_gap (row, col))
	  {
	    tally[residue(row,col)] += in_prob[row];
	  }
      bool ok = 1;
      for (int sym = 0; sym < alphabet.size() && ok; ++sym)
	if (abs (col_sym_freq[sym][col-motif_start] + stut_col_sym_freq[sym][col-motif_start] - tally[sym]) >= tol)
	  ok = 0;
      o << " " << setw(4) << (ok ? "OK" : "BAD");
    }
  o << "\n";
  
  restore_flags (o);
  o.precision (old_prec);
  */
  o << "WARNING -- verify_counts skipped\n";
}

void Gibbs_alignment::display_all_scores (ostream& o, const Gibbs_null_model& null, const vector<Gibbs_alignment*>& model)
{
  Stream_saver ss;
  
  int old_prec = o.precision(3);
  ss.save_flags (o);
  ss.left_align (o);

  o << "Gibbs alignment log-likelihoods (in bits):\n";
  o << "[gibbs] Sequence   Null      ";
  for (int m = 0; m < (int) model.size(); ++m) o << " Model #" << setw(3) << m+1;
  o << "\n";
  
  for (int g = 0; g < (int) null.dataset.profile.size(); ++g)
    {
      o << "[gibbs] ";
      sstring n = null.dataset.row_name[g];
      if (n.size() > 10) n.erase (n.begin() + 10, n.end());
      o << setw(11) << n;
      ss.right_align (o);
      o << setw(10) << Score2Bits (null.profile_null_score[g]);
      for (int m = 0; m < (int) model.size(); ++m)
	o << " " << setw(10) << Score2Bits (model[m]->log_likelihood_sc (g));
      ss.left_align (o);
      o << "\n";
    }

  ss.restore_flags (o);
  o.precision (old_prec);
}

void Gibbs_alignment::display_all_mixture_consensi (ostream& o, const vector<Gibbs_alignment*>& model)
{
  o << "Consensi:";
  for_const_contents (vector<Gibbs_alignment*>, model, g) {
    sstring c = (*g)->mixture_consensus();
    o << " " << c;
  }
  o << "\n";
}

void Gibbs_alignment::display_all_alphabet_consensi (ostream& o, const vector<Gibbs_alignment*>& model)
{
  o << "Consensi:";
  for_const_contents (vector<Gibbs_alignment*>, model, g) {
    Digitized_biosequence dsq = (*g)->alphabet_consensus();
    Biosequence seq;
    (*g)->alphabet.dsq2seq (dsq, seq);
    o << " " << seq;
  }
  o << "\n";
}

void Gibbs_alignment::display_motif_scores (ostream& o) const
{
  o << "Motif scores:\n";
  for (int col = motif_start; col <= motif_end; ++col) {
    o << "Column #" << col << ":";
    for (int context = 0; context < (int) context_null_sc.size(); ++context) {
      Digitized_biosequence nmer = Digitized_biosequence::context_to_nmer (context, alphabet.size(), null_model.context_order);
      Biosequence nmer_seq;
      alphabet.dsq2seq (nmer, nmer_seq);
      o << " " << nmer_seq << "=" << motif_score (nmer.back(), context, context, col, 0);
    }
    o << "\n";
  }
}

void Gibbs_alignment::display_motif_counts (ostream& o) const
{
  o << "Motif counts:\n";
  for (int col = motif_start; col <= motif_end; ++col) {
    o << "Column #" << col << ": col_sym_freq=(" << col_sym_freq[col-motif_start];
    o << ") col_context_freq=(" << col_context_freq[col-motif_start];
    o << ") total_col_freq=" << total_col_freq[col-motif_start] << "\n";
  }
}

Gibbs_alignment_factory::Gibbs_alignment_factory (const Gibbs_dataset& dataset,
						  const Alphabet& alphabet,
						  const Gibbs_null_model& null_model,
						  const vector<double>& log_dirichlet_component_weight,
						  const vector<vector<double> >& dirichlet_mixture,
						  const vector<sstring>& dirichlet_component_name,
						  int motif_length)
  : dataset (dataset),
    alphabet (alphabet),
    null_model (null_model),
    context_null_sc (null_model.context_null_sc),
    profile_null_score (null_model.profile_null_score),
    log_dirichlet_component_weight (log_dirichlet_component_weight),
    dirichlet_mixture (dirichlet_mixture),
    dirichlet_component_name (dirichlet_component_name),
    motif_length (motif_length),
    motif_delete_probability (null_model.motif_delete_probability)
{
  if (CLOGGING(5)) {
    for (int i = 0; i < (int) log_dirichlet_component_weight.size(); ++i) {
      CL << "Dirichlet mixture component #" << i << ": log(weight) = " << log_dirichlet_component_weight[i] << ", pseudocounts = (" << dirichlet_mixture[i] << ")\n";
    }
  }
  CLOG(7) << "Gibbs_alignment_factory constructed\n";
}

Gibbs_alignment* Gibbs_alignment_factory::new_model() const
{
  return new Gibbs_alignment (dataset, alphabet, null_model, log_dirichlet_component_weight, dirichlet_mixture, dirichlet_component_name, motif_length, motif_delete_probability);
}

