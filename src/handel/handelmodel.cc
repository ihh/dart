#include "handelmodel.h"

Pair_HMM_scores Handel_model::conditional_pair_HMM (double t)
{
  Wildcard_transducer trans = transducer (t);
  Pair_HMM_scores hmm (trans.states(), &alphabet());
  hmm.assign_transition_matrix (trans);
  hmm.state_type = trans.state_type;

  const array2d<Score> condmat = Prob2ScoreArray2d (create_conditional_substitution_matrix (t));
  const vector<Score> prior = Prob2ScoreVec (create_prior());
  const vector<Score> one (prior.size(), (Score) 0);

  for (int s = 0; s < hmm.states(); ++s)
    switch (hmm.state_type[s])
      {
      case Pair_HMM_state_type_enum::EmitXY:
	hmm.pair_emit[s] = condmat;
	break;
      case Pair_HMM_state_type_enum::EmitX:
	hmm.single_emit[s] = one;
	break;
      case Pair_HMM_state_type_enum::EmitY:
	hmm.single_emit[s] = prior;
	break;
      case Pair_HMM_state_type_enum::Null:
      default:
	break;
      }

  return hmm;
}

Pair_HMM_scores Handel_model::joint_pair_HMM (double t)
{
  Pair_HMM_scores hmm = conditional_pair_HMM (t);

  const array2d<Score> jointmat = Prob2ScoreArray2d (create_joint_substitution_matrix (t));
  const vector<Score> prior = Prob2ScoreVec (create_prior());

  for (int s = 0; s < hmm.states(); ++s)
    switch (hmm.state_type[s])
      {
      case Pair_HMM_state_type_enum::EmitXY:
	hmm.pair_emit[s] = jointmat;
	break;
      case Pair_HMM_state_type_enum::EmitX:
	hmm.single_emit[s] = prior;
	break;
      default:
	break;
      }

  return hmm;
}

