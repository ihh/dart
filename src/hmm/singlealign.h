#ifndef SINGLE_ALIGN_INCLUDED
#define SINGLE_ALIGN_INCLUDED

#include "hmm/singledp.h"

class Single_HMM_alignment : public Named_rows
{
public:
  Single_HMM_scores&           hmm;
  vector<const Named_profile*> prof;
  vector<vector<int> >         state_path;
  Single_HMM_counts            total_counts;
  int                          aligned_rows;

protected:
  vector<Score>                null_emit_sc;    // null model; -infinity by default
  Score                        null_prior_sc;   // prior score for null model; -infinity by default
  Score                        model_prior_sc;  // prior score for Single_HMM  (= prob2score(1-score2prob(null_prior_sc)))
  vector<Score>                prof_null_sc;    // null scores for each sequence

public:
  Single_HMM_alignment (Single_HMM_scores& hmm) :
    hmm (hmm),
    total_counts (hmm),
    aligned_rows (0),
    null_emit_sc (hmm.alphabet().size(), -Prob2Score (1.0 / (double) hmm.alphabet().size())),
    null_prior_sc (-InfinityScore),
    model_prior_sc (0)
  { }

  void              emit_row (const sstring& name, Sequence_database& db);  // generates a new, aligned sequence
  void              add_row (const sstring& name, const Named_profile* profile, const vector<int>* path = 0);               // adds a row (with no alignment, unless path != 0)
  bool              got_path (int row) const { return ((Single_HMM_alignment&)*this) . state_path[row].size() != 0; }

  void              set_null_model (const vector<Score>& null_emit_sc, Score null_prior_sc);   // recalculates all the null scores

  void              optimise_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts) const;
  void              sample_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts, double kT = 1) const;   // samples from Dirichlet approximation to posterior

  bool              optimise_row (int row);                                                                // returns TRUE if row changed
  bool              sample_row (int row, double kT = 1);

  bool              Gibbs_optimise_row (int row, const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts);   // optimises HMM for all other rows, then aligns this row
  bool              Gibbs_sample_row (int row, const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts);     // optimises HMM for all other rows, then samples this row

  int               path_transition_score (int row) const { return hmm.path_transition_score (((Single_HMM_alignment&)*this).state_path[row]); }
  int               path_emit_score (int row) const { return hmm.path_emit_score (((Single_HMM_alignment&)*this).state_path[row], *((Single_HMM_alignment&)*this).prof[row]); }
  int               path_score (int row) const  { return hmm.path_score (((Single_HMM_alignment&)*this).state_path[row], *((Single_HMM_alignment&)*this).prof[row]); } 
  int               total_path_score() const { int sc = 0; for (int r = 0; r < rows(); ++r) if (got_path(r)) ScorePMulAcc (sc, path_score(r)); return sc; }

  double            accuracy (const Single_HMM_alignment& reference_alignment) const;      // returns proportion of residues aligned to the same states as in the reference alignment
  void              make_local_alignment (const vector<int>& state_sequence, const sstring& label_sequence, Local_alignment& ret_alignment, sstring& ret_column_labels) const;

  bool                  operator== (const Single_HMM_alignment& a) const { return state_path == a.state_path; }
  Single_HMM_alignment& operator= (const Single_HMM_alignment& a) { row_name = a.row_name; prof = a.prof; state_path = a.state_path; total_counts = a.total_counts; return *this; }
};

struct Single_HMM_alignment_counts : Named_rows
{
  Single_HMM_scores&           hmm;
  vector<const Named_profile*> prof;
  vector<Single_HMM_counts>    counts;
  Single_HMM_counts            total_counts;
  
  Single_HMM_alignment_counts (Single_HMM_scores& hmm) : hmm(hmm), total_counts(hmm) { }

  void              add_row (const sstring& name, const Named_profile* profile);                                            // adds a row (with no counts)
  void              reestimate_row_counts (int row, double kT = 1);
  void              optimise_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts) const;
  void              sample_HMM (const Single_HMM_mask& mask, const Single_HMM_counts& pseudocounts, double kT = 1) const;   // samples from Dirichlet approximation to posterior

  Single_HMM_alignment_counts& operator= (const Single_HMM_alignment_counts& c);
};

#endif
