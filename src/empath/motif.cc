#include "empath/motif.h"
#include "seq/seqlogo.h"
#include "util/vector_output.h"
#include "util/math_fn.h"

Motif::Motif (int len, double pseud_mul, Single_matrix_factory& matrix_factory, const Alphabet& alphabet)
  : len (len),
    pseud_mul (pseud_mul),
    alphabet (alphabet),
    hmm (len+3, alphabet),
    hmm_scores(),
    em_scores(),
    prior (em_scores),
    matrix_factory (matrix_factory),
    motif_emit (len),
    START(Start), END(End), NPAD(0), LPAD(1), RPAD(2), MBEG(3), MFIN(3+len-1)
{
  // set up HMM dimensions

  for (int i = 0; i < hmm.states(); ++i) hmm.emit[i] = vector<PFunc> (alphabet.size());

  // make null emit group

  null_emit = em_scores.new_alphabet_group (alphabet, "Null emit");

  // set padding state emit scores to 0

  for (int sym = 0; sym < alphabet.size(); ++sym)
    {
      hmm.emit[NPAD][sym] = 1.0;   // explicit null model state; will replace this with separate null model
      hmm.emit[LPAD][sym] = 1.0;   // left padding state
      hmm.emit[RPAD][sym] = 1.0;   // right padding state
    }

  // make motif emit groups and assign to match states
  
  for (int m = 0; m < len; ++m)
    {
      const int match_state = MBEG + m;
      sstring match_state_name;
      match_state_name << "Match state #" << m + 1;
      motif_emit[m] = em_scores.new_alphabet_group (alphabet, match_state_name.c_str());

      for (int sym = 0; sym < alphabet.size(); ++sym)
	hmm.emit[match_state][sym] = motif_emit[m][sym] / null_emit[sym];
    }

  // make transition groups

  // null model extend probability
  null_extend = em_scores.new_boolean_group ("Null extend");

  // null/prior model decision probability (this will end up elsewhere)
  null_decide = em_scores.new_boolean_group ("Null begin");

  // do transition assignments

  hmm.transition (START, END)  = null_decide.YES * null_extend.NO;

  hmm.transition (START, NPAD) = null_decide.YES * null_extend.YES;
  hmm.transition (NPAD,  NPAD) = null_extend.YES;
  hmm.transition (NPAD,  END)  = null_extend.NO;

  hmm.transition (START, MBEG) = null_decide.NO * null_extend.NO;

  hmm.transition (START, LPAD) = null_decide.NO * null_extend.YES;
  hmm.transition (LPAD,  LPAD) = null_extend.YES;
  hmm.transition (LPAD,  MBEG) = null_extend.NO;

  for (int m = MBEG; m < MFIN; ++m)
    hmm.transition (m, m + 1)  = 1.0;

  hmm.transition (MFIN,  END)  = null_extend.NO;

  hmm.transition (MFIN,  RPAD) = null_extend.YES;
  hmm.transition (RPAD,  RPAD) = null_extend.YES;
  hmm.transition (RPAD,  END)  = null_extend.NO;

  // set up the metascore indices for masking

  const int mask_idx = 0;   // index of mask metascore vector
  for (int m = MBEG; m <= MFIN; ++m)
    hmm.metascore_idx[m].push_back (mask_idx);
}

Motif::~Motif() { clear_fb(); }

void Motif::clear_fb()
{
  for_contents (vector<Single_forward_backward_interface*>, fb, fb_iter) if (*fb_iter) delete *fb_iter;
  fb.clear();
}

void Motif::set_null_emit_prob (const vector<Prob>& null_emit_prob)
{
  em_scores[null_emit] = Prob2ScoreVec (null_emit_prob);

  vector<Prob> pseud = null_emit_prob;
  for_contents (vector<Prob>, pseud, pc) *pc *= pseud_mul;
  Dirichlet_mixture emit_mixture (pseud);

  for_const_contents (vector<Alphabet_group>, motif_emit, g)
    prior.assign (*g, emit_mixture);
}

void Motif::set_null_length (const double length)
{
  const Prob null_extend_prob = Math_fn::p_extend (length);

  const Score null_extend_sc = Prob2Score (null_extend_prob);
  const Score null_end_sc = Prob2Score (1.0 - null_extend_prob);
  
  em_scores[null_extend][0] = null_end_sc;
  em_scores[null_extend][1] = null_extend_sc;
}

void Motif::set_model_prior (const Prob p)
{
  const Score null_prior_sc = Prob2Score (1.0 - p);
  const Score model_prior_sc = Prob2Score (p);

  em_scores[null_decide][0] = model_prior_sc;
  em_scores[null_decide][1] = null_prior_sc;
}

void Motif::seed (const Named_profile& np, int pos)
{
  if (CLOGGING(5)) CL << "Seeding at sequence '" << np.name << "' position " << pos << "\n";
  prior.initialise();
  hmm_scores = hmm.eval_hmm_scores (em_scores);

  vector<int> state_path (np.prof_sc.size());
  for (int i = 0; i < pos; ++i) state_path[i] = LPAD;
  for (int i = 0; i < len; ++i) state_path[pos+i] = MBEG + i;
  for (int i = pos+len; i < (int) np.prof_sc.size(); ++i) state_path[i] = RPAD;
  
  Single_HMM_counts hmm_counts (hmm_scores);
  hmm_counts.add_counts_from_state_path (hmm_scores, np, state_path);
  inc_emit_counts (hmm_counts, 1.0);
}

Loge Motif::do_EM (const Sequence_database_index& index)
{
  clear_fb();
  fb.reserve (index.size());

  hmm_scores = hmm.eval_hmm_scores (em_scores);

  Single_HMM_counts total_hmm_counts (hmm_scores);
  for (int i = 0; i < index.size(); ++i)
    {
      if (CLOGGING(4)) CL << "Getting EM counts for sequence '" << index.name[i] << "'\n";
      fb.push_back (matrix_factory.new_forward_backward (hmm_scores, *index.profile[i]));
      Single_HMM_counts dsq_hmm_counts = fb.back()->get_expected_counts();
      total_hmm_counts.add_counts (dsq_hmm_counts);
    }

  // update EM counts for emit state only
  inc_emit_counts (total_hmm_counts, 1.0);

  return total_hmm_counts.log_likelihood;
}

void Motif::mask (const Sequence_database_index& index) const
{
  if ((int) fb.size() != index.size())
    THROW Standard_exception ("Forward-backward matrices don't match sequence database index");
  Metascore mask;
  for (int i = 0; i < (int) fb.size(); ++i)
    {
      const Metaprob& in_motif_prob = fb[i]->get_expected_metacounts()[0];
      if ((int) in_motif_prob.size() != index.profile[i]->size())
	THROW Standard_exception ("Forward-backward matrix doesn't match sequence database index");

      mask.clear();
      mask.reserve (in_motif_prob.size());

      for (int pos = 0; pos < (int) in_motif_prob.size(); ++pos)
	mask.push_back (Prob2Score (1.0 - in_motif_prob[pos]));

      if (CLOGGING(4)) CL << "Masking '" << index.name[i] << "' with (" << mask << ")\n";
      
      index.profile[i]->add_metascores (0, mask);
    }
}

void Motif::inc_emit_counts (const Single_HMM_counts& hmm_counts, Prob model_count)
{
  // create a counts object for the EM parameter set
  PCounts em_counts (em_scores);

  // map the counts for emit states back onto the EM parameter set
  for (int m = MBEG; m <= MFIN; ++m)
    for (int sym = 0; sym < alphabet.size(); ++sym)
      hmm.emit[m][sym].inc_var_counts (em_counts, em_scores, hmm_counts.emit[m][sym] * model_count);

  // debug output
  if (CTAGGING(2,PCOUNTS)) em_counts.show(CL);

  // ask the Dirichlet prior object to optimise the EM parameters, given the counts
  prior.optimise (em_counts);

  // more debug output
  if (CTAGGING(3,PSCORES)) em_scores.show(CL);
}

Score_profile Motif::prof_sc() const
{
  Score_profile prof_sc;
  for (int pos = 0; pos < len; ++pos)
    {
      Symbol_score_map ssm;
      for (int sym = 0; sym < alphabet.size(); ++sym)
	ssm[sym] = em_scores[motif_emit[pos]][sym];
      prof_sc.push_back (ssm);
    }
  return prof_sc;
}

void Motif::display (ostream& o) const
{
  em_scores.show(o);

  Score_profile prof = prof_sc();
  Sequence_logo logo (prof, alphabet);
  logo.show (o);
}

