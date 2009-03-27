#ifndef GIBBS_INCLUDED
#define GIBBS_INCLUDED

#include "kimono/expr.h"
#include "util/math_fn.h"
#include "util/vector_output.h"
#include "util/strsaver.h"

#include <deque>

// Gibbs_dataset contains the sequences and row names for a Gibbs_alignment
//
struct Gibbs_dataset
{
  const Sequence_database_index& db_index;
  const Sequence_database_index& revcomp_index;

  vector<const Digitized_biosequence*>  profile;               // our sequences (subset of all sequences in db_index)
  vector<const Digitized_biosequence*>  revcomp;               // reverse complemented sequences

  vector<const Biosequence*>  profile_seq;                     // text representations
  vector<const Biosequence*>  revcomp_seq;

  const vector<sstring>& row_name;

  int row_index (const sstring& name) const;

  bool includes_row (int row) const { return profile[row] == 0 ? 0 : profile[row]->size() > 0; }
  bool revcomp_enabled() const { return revcomp.size() > 0; }

  Gibbs_dataset (const Sequence_database_index& db_index, const Sequence_database_index& revcomp_index, const vector<sstring>& row_name, const Alphabet& alphabet);
};

struct Gibbs_null_model
{
  const Gibbs_dataset& dataset;

  int context_order;  /* the 'n' in the 'nth-order' phrases that follow... */
  vector<sstring> context_string;  /* e.g. (for 2nd-order) aaa, aac, aag, aat, aca... */
  
  vector<Digitized_biosequence> profile_context;      /* nth-order description of profile */
  vector<Digitized_biosequence> profile_rev_context;  /* nth-order description of profile reversed (but not complemented!) */

  vector<double> composition;  /* 1th-order composition */
  vector<double> composition_context;  /* nth-order composition */
  vector<int>    context_null_sc;        /* = log compositition_context */
  double         mean_sequence_length;

  vector<int>    profile_null_score;         // null emission score for each sequence (averaged over forward & backward)

  double         motif_delete_probability;
  
  Gibbs_null_model (const Gibbs_dataset& dataset, const Alphabet& alphabet, double composition_pseudocount, double motif_del_prob, int context_order = 0);
};

struct Gibbs_alignment : Stream_saver
{
  // variables that aren't ours to change

  const Gibbs_dataset&          dataset;
  const Alphabet&               alphabet;
  const Gibbs_null_model&       null_model;

  const vector<int>&            context_null_sc;                  // Markov chain null model
  const vector<int>&            profile_null_score;       // scores of each sequence under the null model

  const vector<double>& log_dirichlet_component_weight;
  const vector<vector<double> >& dirichlet_mixture;
  const vector<sstring>& dirichlet_component_name;
  
  int                           motif_length;
  int                           delete_extend_score;
  int                           delete_end_score;

  // our variables
  
  vector<int>                   row_aligned;        // flag indicating whether each row is aligned
  vector<int>                   row_offset;         // column of alignment in which first residue of each sequence is located
  vector<int>                   row_reversed;       // flag indicating whether, for each row, the sequence or its reverse complement is aligned
  int                           aligned_row_count;  // number of aligned rows

  // motif data

  int                           motif_start;        // start column for motif
  int                           motif_end;          // end column for motif
  vector<vector<double> >       col_sym_freq;       // col_sym_freq [col] [sym] = expected count for symbol <sym> in motif column <col>
  vector<vector<double> >       col_context_freq;   // col_context_freq [col] [context] = expected count for Markov-chain context <context> in motif column <col>
  vector<double>                total_col_freq;     // total_col_freq [col] = sum_{sym} col_sym_freq[col][sym]
  vector<vector<int> >          motif;              // motif[C][R] = log likelihood score for finding symbol R in column C of motif
  
  // the following variables are set by set_membership_probability()
  
  vector<int>                   in_sc;              // in_sc[R] = score for row R to be in the alignment
  vector<double>                in_prob;            // same, in probability space

  double                        total_aligned_in_prob;   // sum of in_prob's for all aligned rows (i.e. rows for which row_aligned[row] == true)

  // methods

  int            left_col() const;
  int            right_col() const;

  int flush_left()   // shift whole alignment over so that left_col() == 0; return alignment offset
  {
    const int lc = left_col();
    for (int row = 0; row < rows(); ++row) row_offset[row] -= lc;
    motif_start -= lc;
    motif_end -= lc;
    return -lc;
  }

  int            rows() const { return dataset.profile.size(); }
  int            columns() const { return right_col() + 1 - left_col(); };
  bool           revcomp_enabled() const { return dataset.revcomp_enabled(); }

  void           unalign_row (int row);
  void           align_row (int row, int offset, bool reversed);

  int            residues (int row) const { return dataset.profile[row]->size(); }
  int            start_col (int row) const { return row_offset[row]; }                                                   // first column with a residue in
  int            end_col (int row) const { return row_offset[row] + residues (row) - 1; }                                // last column with a residue in

  bool           not_gap (int row, int col) const { const int offset = row_offset[row]; return row_aligned[row] && col >= offset && col < offset + residues (row); }

  int            left_pad (int row) const { return max (min (end_col(row)+1, motif_start) - start_col(row), 0); }   // padding residues emitted before entering motif on LHS
  int            left_delete (int row) const { return max (min (motif_end+1, start_col(row)) - motif_start, 0); }   // deleted residues on LHS of motif
  int            right_pad (int row) const { return max (end_col(row)+1 - max (start_col(row), motif_end+1), 0); }  // padding residues emitted after leaving motif on RHS
  int            right_delete (int row) const { return max (motif_end+1 - max (end_col(row)+1, motif_start), 0); }  // deleted residues on RHS of motif

  const Digitized_biosequence& row_data (int row) const { return *(row_reversed[row] ? dataset.revcomp[row] : dataset.profile[row]); }
  const Biosequence&           row_data_seq (int row) const { return *(row_reversed[row] ? dataset.revcomp_seq[row] : dataset.profile_seq[row]); }
  int            residue (int row, int col) const { return row_data (row) [col - row_offset[row]]; }
  char           residue_char (int row, int col) const { return row_data_seq (row) [col - row_offset[row]]; }

  const Digitized_biosequence& context_row_data (int row) const { return row_reversed[row] ? null_model.profile_rev_context[row] : null_model.profile_context[row]; }
  const Digitized_biosequence& rev_context_row_data (int row) const { return row_reversed[row] ? null_model.profile_context[row] : null_model.profile_rev_context[row]; }
  int            residue_context (int row, int col) const { return context_row_data(row) [col - row_offset[row]]; }
  int            residue_rev_context (int row, int col) const {
    const Digitized_biosequence& revcon = rev_context_row_data(row);
    const int m = col - row_offset[row];
    return revcon [revcon.size() - 1 - m];
  }

  // CHANGED
  void           add_counts (int sym, int fwd_context, int rev_context, int col, double weight)
    {
      if (sym < 0) THROW Standard_exception ("Can't handle ambiguous symbols");
      const int m = col - motif_start;
      col_sym_freq[m][sym] += weight;
      col_context_freq[m][fwd_context] += weight / 2;
      col_context_freq[m][rev_context] += weight / 2;
      total_col_freq[m] += weight;
      recalculate_motif_column (col);
    }

  // CHANGED
  int            motif_score (int sym, int fwd_context, int rev_context, int col, double weight = 0) const             // works out the optimal column score for the given residue (weight is currently unused)
    {
      if (sym < 0) THROW Standard_exception ("Can't handle ambiguous symbols");
      return motif [col - motif_start] [sym] - (context_null_sc [fwd_context] + context_null_sc [rev_context]) / 2;
      /*
       * Old code:
       *
      if (sym < 0 || prev_sym < 0) THROW Standard_exception ("Can't handle ambiguous symbols");
      const int m = col - motif_start;

      return ScorePMul (Prob2Score (col_sym_freq[m][sym] / total_col_freq[m]), 
					  -(sym==prev_sym ? stutter_sc : null_sc[sym]));
      */

      /*
       * The following, commented-out definition adds the weight for this residue onto the current column counts.
       * I took this out because it was destabilising fragile motifs.
       *
      return ScorePMul (Prob2Score ((col_sym_freq[m][sym] + weight) / (total_col_freq[m] + weight)), 
					  -(sym==prev_sym ? stutter_sc : null_sc[sym]));
      */
    }
  
  // CHANGED
  int            column_score (const vector<double>& sym_freq, const vector<double>& context_freq) const
    {
      vector<int> optimal_sc;
      vector<double> log_dirichlet_evidence;
      double log_total_dirichlet_evidence;
      calculate_optimal_scores (sym_freq, optimal_sc, log_dirichlet_evidence, log_total_dirichlet_evidence);
      int sc = 0;
      for (int sym = 0; sym < alphabet.size(); ++sym)
	ScorePMulAcc (sc, (int) (sym_freq[sym] * (double) optimal_sc[sym]));
      for (int context = 0; context < (int) context_freq.size(); ++context)
	ScorePMulAcc (sc, (int) (-context_freq[context] * (double) context_null_sc[context]));
      return sc;
    }
  
  void           optimise_motif();
  void           sample_motif (double kT = 1);       // NB only samples the motif *position*, not the probabilities
  void           calculate_motif_matrix (vector<int>& matrix) const;  // calculates DP matrix for finding motif position
  void           display_motif_matrix (ostream& o, const vector<int>& matrix, double display_threshold = .01) const;
  
  int            DP_matrix_size (int row) const { return residues (row) + motif_end - motif_start; }
  int            convert_matrix_column_to_row_offset (int matrix_col, int row) const { return matrix_col + (motif_start - (residues(row) - 1)); }
  void           calculate_DP_matrix (const Digitized_biosequence& dsq,
				      const Digitized_biosequence& fwd_context,
				      const Digitized_biosequence& rev_context,
				      vector<int>::iterator matrix_begin,
				      int matrix_size,
				      double weight) const;   // assumes that motif_start & motif_end are valid
  void           display_DP_matrix (ostream& o, const vector<int>& matrix, int row, int highlight_col = -1) const;

  void           recalculate_counts();  // calls recalculate_column_counts() and recalculate_motif_column(); call after motif_start and motif_end changed
  // CHANGED
  void           recalculate_column_counts (int col, vector<double>& sym_freq, vector<double>& context_freq, double& total_freq) const
  {
    sym_freq = vector<double> (alphabet.size(), (double) 0);
    context_freq = vector<double> (context_null_sc.size(), (double) 0);
    total_freq = 0;
    for (int row = 0; row < rows(); ++row) {
      if (not_gap (row, col)) {
	const double w = in_prob[row];
	sym_freq [residue (row, col)] += w;
	context_freq [residue_context (row, col)] += w / 2;
	context_freq [residue_rev_context (row, col)] += w / 2;
	total_freq += w;
      }
    }
  }
  void recalculate_motif_column (int col) {
    const int m = col - motif_start;
    vector<double> log_dirichlet_evidence;
    double log_total_dirichlet_evidence;
    calculate_optimal_scores (col_sym_freq[m], motif[m], log_dirichlet_evidence, log_total_dirichlet_evidence);
  }
  void calculate_optimal_scores (const vector<double>& sym_freq, vector<int>& motif_col, vector<double>& log_dirichlet_evidence, double& log_total_dirichlet_evidence) const {

    vector<vector<double> > freq_with_pseudocounts (dirichlet_mixture.size(), sym_freq);
    
    log_dirichlet_evidence = vector<double> (dirichlet_mixture.size());
    log_total_dirichlet_evidence = -InfinityScore;

    for (int i = 0; i < (int) dirichlet_mixture.size(); ++i) {
      for (int s = 0; s < alphabet.size(); ++s) freq_with_pseudocounts[i][s] += dirichlet_mixture[i][s];
      log_dirichlet_evidence[i] = log_dirichlet_component_weight[i] + Math_fn::log_dirichlet_normaliser (freq_with_pseudocounts[i]) - Math_fn::log_dirichlet_normaliser (dirichlet_mixture[i]);
      NatsPSumAcc (log_total_dirichlet_evidence, log_dirichlet_evidence[i]);
    }

    motif_col = vector<int> (alphabet.size(), -InfinityScore);

    for (int i = 0; i < (int) log_dirichlet_evidence.size(); ++i) {
      double norm = 0;
      for (int s = 0; s < alphabet.size(); ++s) norm += freq_with_pseudocounts[i][s];
      if (norm == 0) norm = 1;
      int posterior_component_weight_sc = Nats2Score (log_dirichlet_evidence[i] - log_total_dirichlet_evidence);
      for (int s = 0; s < alphabet.size(); ++s)
	ScorePSumAcc (motif_col[s],
				       ScorePMul (posterior_component_weight_sc,
								    Prob2Score (freq_with_pseudocounts[i][s] / norm)));
    }
    /* logging code commented out because it disrupted display_motif_matrix logging status
     * (wishlist: push & pop logging status)
     *
    if (CLOGGING(-2)) {
      vector<sstring> ldp_string (dirichlet_component_name);
      for (int i = 0; i < log_dirichlet_evidence.size(); ++i) ldp_string[i] << ":" << log_dirichlet_evidence[i] - log_total_dirichlet_evidence;
      CL << "sym_freq = (" << sym_freq << ") log_dirichlet_posterior = (" << ldp_string << ") motif_col=(" << motif_col << ")\n";
    }
    */
  }
  
  bool           sample_row (int row, double kT = 1);       // returns TRUE if row changed
  bool           optimise_row (int row);                    // returns TRUE if row changed

  void           set_membership_probability (int row, int score);    // modifies col_sym_freq[][] accordingly
  int            log_likelihood_sc (int row) const;                  // returns the row log-likelihood score, conditional on the row being in the alignment
  int            log_parameter_prior_sc() const;                     // returns the log of the prior probability of the model parameters

  void           display (ostream& o, double threshold = .01, int columns = 60);   // sequences whose membership probability are less than (threshold) will not be displayed
  void           display_scores (ostream& o, double threshold = .01);
  void           display_counts (ostream& o, int row = -1);
  void           display_column_log_odds_scores (ostream& o);
  void           verify_counts (ostream& o, double tol = .01);
  static void    display_all_scores (ostream& o, const Gibbs_null_model& null, const vector<Gibbs_alignment*>& model);
  static void    display_all_mixture_consensi (ostream& o, const vector<Gibbs_alignment*>& model);
  static void    display_all_alphabet_consensi (ostream& o, const vector<Gibbs_alignment*>& model);

  void display_motif_scores (ostream& o) const;
  void display_motif_counts (ostream& o) const;

  Digitized_biosequence alphabet_consensus() const;
  sstring mixture_consensus() const;
  
  Gibbs_alignment (const Gibbs_dataset& dataset,
		   const Alphabet& alphabet,
		   const Gibbs_null_model& null_model,
		   const vector<double>& log_dirichlet_component_weight,
		   const vector<vector<double> >& dirichlet_mixture,
		   const vector<sstring>& dirichlet_component_name,
 		   int motif_length,
		   double motif_delete_probability);

};

struct Gibbs_alignment_factory
{
  const Gibbs_dataset&    dataset;                  // sequences
  const Alphabet&         alphabet;                 // alphabet
  const Gibbs_null_model& null_model;               // null model
  const vector<int>&      context_null_sc;          // Markov chain from null model
  const vector<int>&      profile_null_score;       // null emit scores by sequence

  const vector<double>& log_dirichlet_component_weight;
  const vector<vector<double> >& dirichlet_mixture;
  const vector<sstring>& dirichlet_component_name;

  int                     motif_length;
  double                  motif_delete_probability;

  Gibbs_alignment_factory (const Gibbs_dataset& dataset,
			   const Alphabet& alphabet,
			   const Gibbs_null_model& null_model,
			   const vector<double>& log_dirichlet_component_weight,
			   const vector<vector<double> >& dirichlet_mixture,
			   const vector<sstring>& dirichlet_component_name,
			   int motif_length);

  Gibbs_alignment* new_model() const;
};

#endif
