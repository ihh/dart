#include "empath/trainer.h"
#include "seq/seqlogo.h"
#include "util/math_fn.h"
#include "util/vector_output.h"

Trainable::Trainable (int states, const Alphabet& alphabet, PScores& pscore) :
  emit_group(),
  hmm (states, alphabet),
  pscore (pscore),
  prior (pscore),
  mask_metascore_idx (-1),
  null_emit (pscore.new_alphabet_group (alphabet, "Null emit")),
  null_extend (pscore.new_boolean_group ("Null extend")),
  seed_path()
{
  for (int i = 0; i < StateOffset; ++i)
    hmm.state_type[i] = Single_state_typing::Null;
  hmm.transition (Start, FwdStart) = 1;
  hmm.transition (Start, RevStart) = 1;
  hmm.transition (FwdEnd, End) = 1;
  hmm.transition (RevEnd, End) = 1;
}

Trainable::~Trainable()
{ }

void Trainable::instruct_turtle (ostream& turtle_stream) const
{
  int old_prec = turtle_stream.precision (10);
  save_flags (turtle_stream);
  right_align (turtle_stream);
  turtle_stream.fill (' ');

  for (int pos = -1; pos < (int) emit_group.size(); ++pos)
    {
      turtle_stream.width (4);
      if (pos < 0)
	turtle_stream << "null";
      else
	turtle_stream << pos;
      const vector<Score>& score_vec = pscore [pos < 0 ? null_emit : emit_group[pos]];
      for (int sym = 0; sym < hmm.alphabet().size(); ++sym)
	{
	  turtle_stream.width (16);
	  turtle_stream << Score2Prob (score_vec[sym]);
	}
      turtle_stream << '\n';
    }

  restore_flags (turtle_stream);
  turtle_stream.precision (old_prec);
}

Alphabet_group Trainable::new_emit_group (const char* name)
{
  Alphabet_group g = pscore.new_alphabet_group (hmm.alphabet(), name);
  emit_group.push_back (g);
  return g;
}

void Trainable::optimise_pscore (const Single_HMM_counts& hmm_counts)
{
  // get counts for PVar's
  PCounts pcount (pscore);
  hmm.inc_var_counts (pcount, pscore, hmm_counts, 1.0);
  if (CTAGGING(2,PCOUNTS)) pcount.show(CL);

  // update PVar scores
  prior.optimise (pcount);
  if (CTAGGING(3,PSCORES)) pscore.show(CL);
}

int Trainable::seed_path_residues() const
{
  int res = 0;
  for_const_contents (vector<int>, seed_path, s)
    if (hmm.state_type[*s] == Single_PHMM::Emit)
      ++res;
  return res;
}

Score_profile Trainable::prof_sc() const
{
  Single_HMM_scores hmm_sc = hmm.eval_hmm_scores (pscore);
  vector<int> path = hmm_sc.consensus_state_path();
  Score_profile prof_sc;
  hmm_sc.make_profile (path, prof_sc);
  return prof_sc;
}

void Trainable::set_null_emit_prob (const vector<Prob>& null_emit_prob, double pseudocount_multiplier)
{
  pscore[null_emit] = Prob2ScoreVec (null_emit_prob);

  vector<Prob> pseud = null_emit_prob;
  for_contents (vector<Prob>, pseud, pc) *pc *= pseudocount_multiplier;
  Dirichlet_mixture emit_mixture (pseud);

  for_const_contents (vector<Alphabet_group>, emit_group, g)
    prior.assign (*g, emit_mixture);
}

void Trainable::set_null_length (const double length)
{
  const Prob null_extend_prob = Math_fn::p_extend (length);

  const Score null_extend_sc = Prob2Score (null_extend_prob);
  const Score null_end_sc = Prob2Score (1.0 - null_extend_prob);
  
  pscore[null_extend.NO] = null_end_sc;
  pscore[null_extend.YES] = null_extend_sc;
}

void Trainable::optimise_null_model (const Sequence_database& db, double pseudocount_multiplier)
{
  const vector<Prob> null_emit_prob = db.get_null_model (hmm.alphabet().size());
  const double null_length = db.mean_length();
  set_null_emit_prob (null_emit_prob, pseudocount_multiplier);
  set_null_length (null_length);
}

void Trainable::reset_to_prior()
{
  CLOG(4) << "Resetting model probabilities to prior mean\n";
  prior.initialise();
}

void Trainable::seed (const Named_profile& np)
{
  reset_to_prior();
  if (np.seq.size()) CLOG(4) << "Seeding model on sequence '" << np.seq << "'";
  // get HMM scores
  Single_HMM_scores hmm_scores = hmm.eval_hmm_scores (pscore);
  // get counts for seed path
  Single_HMM_counts hmm_counts (hmm_scores);
  hmm_counts.add_counts_from_state_path (hmm_scores, np, seed_path);
  // blank out transition counts (vvv hacky! this whole seed routine should be rewritten)
  for (int i = 0; i < hmm.states(); ++i)
    {
      for (int j = 0; j < hmm.states(); ++j)
	hmm_counts.transition(i,j) = 0;
      hmm_counts.transition(Start,i) = 0;
      hmm_counts.transition(i,End) = 0;
    }
  // update pscore
  optimise_pscore (hmm_counts);
}

void Trainable::train (const Sequence_database& db, Single_matrix_factory& dp_factory)
{
  reset_to_prior();
  CLOG(6) << "Training model\n";
  Loge loglike = -InfinityLoge;
  Loge prev_loglike;
  do
    {
      prev_loglike = loglike;
      loglike = do_EM (db, dp_factory);
    }
  while (loglike > prev_loglike);
}

Loge Trainable::do_EM (const Sequence_database& db, Single_matrix_factory& dp_factory)
{
  // get HMM scores
  Single_HMM_scores hmm_scores = hmm.eval_hmm_scores (pscore);
  Single_HMM_counts total_hmm_counts (hmm_scores);
  // tally counts for each sequence
  for_const_contents (Sequence_database, db, np)
    {
      if (CLOGGING(4)) CL << "Getting EM counts for sequence '" << (*np).name << "'\n";
      const Single_fast_forward_backward_matrix dsq_fb (hmm_scores, *np);
      const Single_HMM_counts dsq_hmm_counts = dsq_fb.get_expected_counts();
      total_hmm_counts.add_counts (dsq_hmm_counts);
    }
  // update pscore & return
  optimise_pscore (total_hmm_counts);
  return total_hmm_counts.log_likelihood;
}

Local_trainer::Local_trainer (Trainable& model, Single_matrix_factory& dp_factory) :
  model (model),
  dp_factory (dp_factory),
  local_hmm (model.hmm, model.null_emit, model.null_extend, model.pscore, model.mask_metascore_idx),
  prior (model.prior)
{ }

void Local_trainer::global_train (const Sequence_database& db)
{
  prior.initialise();
  model.train (db, dp_factory);
}

void Local_trainer::local_search (const Sequence_database& db, GFF_list& results, const char* source, const char* feature)
{
  // create a consensus string
  Score_profile prof_sc = model.prof_sc();
  Biosequence seq;
  model.hmm.alphabet().score2seq (prof_sc, seq);
  // leave out the regexp version of the model, as most just end up being /[acgt][acgt][acgt].../
  sstring pattern;
  // old, regexp version:
  //  sstring pattern = model.hmm.alphabet().degenerate_seq_pattern (seq);
  // do the search
  local_hmm.search (model.pscore, db, results, model.reverse_states, dp_factory, source, feature, pattern.c_str());
}

void Local_trainer::local_seed (const Named_profile& np, int pos)
{
  if (CTAGGING(5,LOCAL_SEED)) CL << "Seeding at sequence '" << np.name << "' position " << pos << "\n";
  Named_profile np_subseq = np.subseq (pos, model.seed_path_residues());
  prior.initialise();
  model.seed (np_subseq);
}

Loge Local_trainer::local_EM (const Sequence_database& db)
{
  // get HMM scores
  Single_HMM_scores hmm_scores = local_hmm.eval_hmm_scores (model.pscore);
  Single_HMM_counts total_hmm_counts (hmm_scores);
  // do forward-backward on each sequence
  for_const_contents (Sequence_database, db, np)
    {
      if (CLOGGING(4)) CL << "Getting EM counts for sequence '" << np->name << "'\n";
      const Single_forward_backward_interface* fb = dp_factory.new_forward_backward (hmm_scores, *np);
      const Single_HMM_counts hmm_counts = fb->get_expected_counts();
      total_hmm_counts.add_counts (hmm_counts);
      delete fb;
    }
  // update and return
  optimise_pscore (total_hmm_counts);
  return total_hmm_counts.log_likelihood;
}

Loge Local_trainer::local_EM (const Named_profile& np)
{
  // get HMM scores
  Single_HMM_scores hmm_scores = local_hmm.eval_hmm_scores (model.pscore);
  // do forward-backward on single sequence
  if (CLOGGING(4)) CL << "Getting EM counts for sequence '" << np.name << "'\n";
  const Single_forward_backward_interface* fb = dp_factory.new_forward_backward (hmm_scores, np);
  const Single_HMM_counts hmm_counts = fb->get_expected_counts();
  delete fb;
  // update and return
  optimise_pscore (hmm_counts);
  return hmm_counts.log_likelihood;
}

Loge Local_trainer::local_forward (const Named_profile& np)
{
  // get HMM scores
  Single_HMM_scores hmm_scores = local_hmm.eval_hmm_scores (model.pscore);
  // do forward on sequence
  if (CLOGGING(4)) CL << "Getting EM counts for sequence '" << np.name << "'\n";
  const Single_forward_interface* fwd = dp_factory.new_forward (hmm_scores, np);
  const Loge loglike = Score2Nats (fwd->get_forward_score());
  delete fwd;
  return loglike;
}

void Local_trainer::local_mask (Sequence_database& db) const
{
  CLOG(6) << "Masking sequences\n";
  if (!local_hmm.uses_metascores()) { CLOG(6) << "HMM has no metascores; skipping masking step\n"; return; }
  // get HMM scores
  Single_HMM_scores hmm_scores = local_hmm.eval_hmm_scores (model.pscore);
  // do forward-backward on each sequence
  for_contents (Sequence_database, db, np)
    {
      const Single_forward_backward_interface* fb = dp_factory.new_forward_backward (hmm_scores, *np);
      const vector<Metaprob>& expected_metacounts = fb->get_expected_metacounts();
      if (CTAGGING(3,LOCAL_MASK)) { CL << "Expected metacounts for sequence '" << np->name << "':\n"; fb->show_metacounts(CL); }
      local_hmm.mask (*np, expected_metacounts);
      delete fb;
    }
}

void Local_trainer::get_local_mask (const Sequence_database& db, GFF_list& gff_list, Prob min_prob, const char* source, const char* feature) const
{
  CLOG(6) << "Getting sequence masks\n";
  if (!local_hmm.uses_metascores()) { CLOG(6) << "HMM has no metascores; skipping mask-finding step\n"; return; }
  // get HMM scores
  Single_HMM_scores hmm_scores = local_hmm.eval_hmm_scores (model.pscore);
  // do forward-backward on each sequence
  for_const_contents (Sequence_database, db, np)
    {
      const Single_forward_backward_interface* fb = dp_factory.new_forward_backward (hmm_scores, *np);
      const vector<Metaprob>& expected_metacounts = fb->get_expected_metacounts();
      if (CTAGGING(3,LOCAL_MASK)) { CL << "Expected metacounts for sequence '" << np->name << "':\n"; fb->show_metacounts(CL); }
      local_hmm.get_mask_gff (gff_list, expected_metacounts, min_prob, np->name.c_str(), np->seq.c_str(), source, feature);
      delete fb;
    }
}

void Local_trainer::optimise_pscore (const Single_HMM_counts& local_hmm_counts)
{
  // get counts for PVar's
  PCounts pcount (model.pscore);
  local_hmm.inc_var_counts (pcount, model.pscore, local_hmm_counts, 1.0);
  if (CTAGGING(2,PCOUNTS)) pcount.show(CL);

  // update PVar scores
  prior.optimise (pcount);
  if (CTAGGING(3,PSCORES)) model.pscore.show(CL);
}

void Local_trainer::display_vars (ostream& o) const
{
  model.pscore.show(o);
}

sstring Local_trainer::seqlogo() const
{
  sstring logo_str;
  Score_profile prof = model.prof_sc();
  Sequence_logo logo (prof, model.hmm.alphabet());
  logo.show (logo_str);
  return logo_str;
}
