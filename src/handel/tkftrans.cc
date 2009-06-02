#include "handel/tkftrans.h"

TKF91_transducer_factory::TKF91_transducer_factory (double lambda, double mu)
  : Transducer_alignment_with_subst_model(),
    lambda (lambda),
    mu (mu)
{
  sort_indels = false;  // TKF91 transducer has something to say about indel ordering
}

Pair_transducer_scores TKF91_transducer_factory::prior_pair_trans_sc()
{
  Pair_transducer_scores prior (2);

  prior.state_type[0] = Transducer_state_type_enum::TransducerInsertType;
  prior.state_type[1] = Transducer_state_type_enum::TransducerWaitType;

  prior.alphabet_size = subst_model.m();
  prior.alloc_pair_emit (Prob2Score(1.));

  prior.state_name[0] = "I";
  prior.state_name[1] = "W";

  for (int s = 0; s < prior.states(); ++s)
    if (prior.state_type[s] != TransducerWaitType)
      prior.emit_label[s] = prior.state_name[s];

  const vector<Prob> ins = subst_model.create_prior();
  for (int i = 0; i < subst_model.m(); ++i)
    prior.insert_val(0,i) = Prob2Score(ins[i]);

  prior.transition (Start, 0) = Prob2Score (lambda / mu);
  prior.transition (Start, 1) = Prob2Score (1. - lambda / mu);

  prior.transition (0, 0) = Prob2Score (lambda / mu);
  prior.transition (0, 1) = Prob2Score (1. - lambda / mu);

  prior.transition (1, End) = Prob2Score (1.);

  return prior;
}

Pair_transducer_scores TKF91_transducer_factory::branch_pair_trans_sc (double time)
{
  Pair_transducer_scores branch (4);

  branch.state_type[0] = Transducer_state_type_enum::TransducerInsertType;
  branch.state_type[1] = Transducer_state_type_enum::TransducerWaitType;
  branch.state_type[2] = Transducer_state_type_enum::TransducerMatchType;
  branch.state_type[3] = Transducer_state_type_enum::TransducerDeleteType;

  branch.alphabet_size = subst_model.m();
  branch.alloc_pair_emit (Prob2Score(1.));

  branch.state_name[0] = "I";
  branch.state_name[1] = "W";
  branch.state_name[2] = "M";
  branch.state_name[3] = "D";

  for (int s = 0; s < branch.states(); ++s)
    if (branch.state_type[s] != TransducerWaitType)
      branch.emit_label[s] = branch.state_name[s];

  const vector<Prob> ins = subst_model.create_prior();
  const array2d<Prob> mat = subst_model.create_conditional_substitution_matrix (time);

  for (int i = 0; i < subst_model.m(); ++i)
    {
      branch.insert_val(0,i) = Prob2Score(ins[i]);
      branch.delete_val(3,i) = Prob2Score(1.);
      for (int j = 0; j < subst_model.m(); ++j)
	branch.match_val(2,i,j) = Prob2Score (mat(i,j));
    }

  const double alpha = exp (-mu * time);
  const double beta = (lambda * (1 - exp ((lambda - mu) * time))) / (mu - lambda * exp ((lambda - mu) * time));
  const double gamma = 1. + ((1 - exp (-mu * time)) * (mu - lambda * exp ((lambda - mu) * time))) / (-mu * (1 - exp ((lambda - mu) * time)));

  branch.transition (Start, 0) = Prob2Score (beta);
  branch.transition (Start, 1) = Prob2Score (1. - beta);

  branch.transition (0, 0) = Prob2Score (beta);
  branch.transition (0, 1) = Prob2Score (1. - beta);

  branch.transition (1, 2) = Prob2Score (alpha);
  branch.transition (1, 3) = Prob2Score (1. - alpha);
  branch.transition (1, End) = Prob2Score (1.);

  branch.transition (2, 0) = Prob2Score (beta);
  branch.transition (2, 1) = Prob2Score (1. - beta);

  branch.transition (3, 0) = Prob2Score (gamma);
  branch.transition (3, 1) = Prob2Score (1. - gamma);

  return branch;
}

double TKF91_transducer_factory::gap_rate() { return mu; }
double TKF91_transducer_factory::mean_gap_size() { return 1.; }

TKF91_transducer_factory* TKF91_transducer_factory::clone()
{
  TKF91_transducer_factory* ttf = new TKF91_transducer_factory (lambda, mu);
  ttf->subst_model = subst_model;
  // clone Handel_base
  (Handel_base&) *ttf = (Handel_base&) *this;
  // clone selected Tree_alignment members (this should be in Tree_alignment, seriously...)
  ttf->tree = tree;
  ttf->align = align;
  ttf->node2row = node2row;
  ttf->row2node = row2node;
  ttf->handel_seq_scores = handel_seq_scores;
  ttf->handel_branch_scores = handel_branch_scores;
  // return
  return ttf;
}
