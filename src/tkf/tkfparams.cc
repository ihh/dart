#include "tkf/tkfparams.h"

double TKF_branch_transition_scores::effective_zero_time = .00001;

TKF_branch_transition_scores::TKF_branch_transition_scores (const TKF_params& params, double time) : TKF_seq_scores (params)
{
  time = max (time, effective_zero_time);
  
  ancest = Prob2Score (params.ancestral_survival_prob(time));
  not_ancest = Prob2Score (1 - params.ancestral_survival_prob(time));
    
  descen = Prob2Score (params.more_descendants_prob(time));
  not_descen = Prob2Score (1 - params.more_descendants_prob(time));
  
  orphan = Prob2Score (params.orphan_survival_prob(time));
  not_orphan = Prob2Score (1 - params.orphan_survival_prob(time));
}

int TKF_branch_transition_scores::pair_HMM_end_state = 4;

int TKF_branch_transition_scores::conditional_path_score (const Pairwise_path& path) const
{
  vector<bool>& parent = ((Pairwise_path&)path).parent();
  vector<bool>& child  = ((Pairwise_path&)path).child();

  if (path.columns() == 0) return 0;

  int score = 0;
  int i = 0;
  while (i < path.columns())
    {
      if (parent[i]) { score += not_descen; break; }
      if (child[i]) score += descen;
      i++;
    }

  int d;
  int not_d = not_descen;
  while (i < path.columns())
    {
      if (child[i]) { score += ancest;     d = descen; not_d = not_descen; }
      else          { score += not_ancest; d = orphan; not_d = not_orphan; }
      while (++i < path.columns())
	{
	  if (parent[i]) { score += not_d; break; }
	  if (child[i]) { score += d; d = descen; not_d = not_descen; }
	}
    }
  score += not_d;

  return score;
}

int TKF_branch_transition_scores::joint_pair_HMM_transition_score (int src, int dest) const
{
  // States: 0=start, 1=-C, 2=P-, 3=PC, 4=end   (P=parent, C=child)

  int t = -InfinityScore;
  switch (src)
    {
    case 0:                         // start
    case 1:                         // -C
    case 3: switch (dest)           // PC
      {
      case 1: t = descen; break;                             // -C
      case 2: t = eqmlen + not_descen + not_ancest; break;   // P-
      case 3: t = eqmlen + not_descen + ancest; break;       // PC
      case 4: t = not_eqmlen + not_descen; break;            // end
      default: break;
      }
    break;
    case 2: switch (dest)           // P-
      {
      case 1: t = orphan; break;                             // -C
      case 2: t = eqmlen + not_orphan + not_ancest; break;   // P-
      case 3: t = eqmlen + not_orphan + ancest; break;       // PC
      case 4: t = not_eqmlen + not_orphan; break;            // end
      default: break;
      }
    break;
    default: break;
    }

  return max(-InfinityScore,t);
}

int TKF_branch_transition_scores::conditional_pair_HMM_transition_score (int src, int dest) const
{
  // States: 0=start, 1=-C, 2=P-, 3=PC, 4=end   (P=parent, C=child)

  int t = -InfinityScore;
  switch (src)
    {
    case 0:                         // start
    case 1:                         // -C
    case 3: switch (dest)           // PC
      {
      case 1: t = descen; break;                    // -C
      case 2: t = not_descen + not_ancest; break;   // P-
      case 3: t = not_descen + ancest; break;       // PC
      case 4: t = not_descen; break;                // end
      default: break;
      }
    break;
    case 2: switch (dest)           // P-
      {
      case 1: t = orphan; break;                    // -C
      case 2: t = not_orphan + not_ancest; break;   // P-
      case 3: t = not_orphan + ancest; break;       // PC
      case 4: t = not_orphan; break;                // end
      default: break;
      }
    break;
    default: break;
    }

  return max(-InfinityScore,t);
}

void TKF_branch_transition_scores::explain_labels (ostream& o)
{
  o << "TKF log-likelihoods (in bits): a = ancestral survival, o = first orphan, d = more descendants, e = equilibrium length\n";
}

void TKF_branch_transition_scores::show (ostream& o) const
{
  o << "a = " << Score2Bits(ancest) << ", !a = " << Score2Bits(not_ancest) << ", ";
  o << "o = " << Score2Bits(orphan) << ", !o = " << Score2Bits(not_orphan) << ", ";
  o << "d = " << Score2Bits(descen) << ", !d = " << Score2Bits(not_descen) << ", ";
  o << "e = " << Score2Bits(eqmlen) << ", !e = " << Score2Bits(not_eqmlen) << "\n";
}

TKF_branch_scores::TKF_branch_scores (const TKF_params& params, double time) :
  TKF_branch_transition_scores (params, time),
  cond_submat (Prob2ScoreArray2d (params.submat_factory.create_conditional_substitution_matrix(time)))
{ }
