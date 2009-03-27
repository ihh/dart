#ifndef MOTIF_INCLUDED
#define MOTIF_INCLUDED

#include "hmm/singlehmm.h"
#include "seq/dirichlet.h"

// WAY messy. Convert to trainer.h in NEAR future
class Motif : public HMM_state_enum
{
public:
  const int len;
  const double pseud_mul;

private:
  const Alphabet&        alphabet;
  Single_PHMM            hmm;
  Single_HMM_scores      hmm_scores;
  PScores                em_scores;
  Dirichlet_prior        prior;
  Single_matrix_factory& matrix_factory;

  Alphabet_group         null_emit;
  vector<Alphabet_group> motif_emit;
  Boolean_group          null_extend;
  Boolean_group          null_decide;
  
  const int START; // start state
  const int END;   // end state
  const int NPAD;  // null pad state
  const int LPAD;  // left pad state
  const int RPAD;  // right pad state
  const int MBEG;  // model begin (first) state
  const int MFIN;  // model end (last) state

  // this next method appears to be misnamed. it doesn't increment the motif_emit counts, it optimises them.
  // also, why do we need a model_count here?
  // mmmm. bitrot.
  void inc_emit_counts (const Single_HMM_counts& hmm_counts, Prob model_count = 1.0);

  vector<Single_forward_backward_interface*> fb;    // forward-backward matrices
  void clear_fb();

public:

  Motif (int len, double pseud_mul, Single_matrix_factory& matrix_factory, const Alphabet& alphabet);
  ~Motif();

  void set_null_emit_prob (const vector<Prob>& null_emit_prob);
  void set_null_length (const double length);
  void set_model_prior (const Prob p);

  void seed (const Named_profile& np, int pos);
  Loge do_EM (const Sequence_database_index& index);  // does one round of EM (calculates FB matrices & counts); returns log-likelihood
  void mask (const Sequence_database_index& index) const;   // call AFTER a call to do_EM() (assumes FB matrices pre-calculated)

  Score_profile prof_sc() const;

  void display (ostream& o) const;
};

#endif
