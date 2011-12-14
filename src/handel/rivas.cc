#include "handel/rivas.h"

// comment the following #define to enable D->I transitions
// #define ALLOW_DELETE_TO_INSERT_TRANSITIONS

// number of times to sample parameters during parameter-sampling step
#define RIVAS_PRIOR_SAMPLES 100

// static inits
// May want to allow these to be configurable from cmdline?
// IH (3/30/09): configurable from cmdline would be good...

// alpha, beta are hyperparameters for the gamma 
// distribution from which delete rates are sampled
// default rate = alpha/beta, weight = 1/beta

// TODO: reestimate convex hyperparams!
double Affine_transducer_factory::alpha = 195.437;
double Affine_transducer_factory::beta = 3.861;
double Convex_transducer_factory::alpha = 195.437;
double Convex_transducer_factory::beta = 3.861;

// _alpha and _beta are the hyperparameters
// for the beta distribution from which gamma and 
// delete extend probabilities are sampled
double Affine_transducer_factory::gamma_alpha = 8.597;
double Affine_transducer_factory::gamma_beta = 0.488;
double Convex_transducer_factory::gamma_alpha = 8.597;
double Convex_transducer_factory::gamma_beta = 0.488;
double Affine_transducer_factory::delete_extend_alpha = 10.313;
double Affine_transducer_factory::delete_extend_beta = 1.639;

// Component_prior singleton constructor

Convex_transducer_factory::Component_prior::Component_prior()
  :  weight_count (2),  // 2-component mixture
     alpha (2),
     beta (2)
{
  // mild (2:1) prior on mixture component weights
  weight_count[0] = 2.;
  weight_count[1] = 1.;

  // mild del_extend prior; hand-tweaked hyperparams
  alpha[0] = 19.72;
  beta[0] = 21.94;
  alpha[1] = 7.83;
  beta[1] = 1.64;
}


Convex_transducer_factory::Component_prior Convex_transducer_factory::cpt_prior = Convex_transducer_factory::Component_prior();

Rivas_transducer_factory::Rivas_transducer_factory (int branch_states, double gamma_val) :
  prior (2),
  branch (branch_states),
  del_rate (0.),
  mean_del_size (0.),
  gamma_val (gamma_val)
{
  // initialise prior transducer
  init_prior_ins (0, "Emit");   // Ins
  init_prior_wait (1, "Wait");  // Wait

  prior.transition (Start, 0) = Prob2Score (gamma_val);
  prior.transition (Start, 1) = Prob2Score (1. - gamma_val);

  prior.transition (0, 0) = Prob2Score (gamma_val);
  prior.transition (0, 1) = Prob2Score (1. - gamma_val);

  prior.transition (1, End) = 0;
}

Rivas_transducer_factory::~Rivas_transducer_factory()
{
  for_contents (Trans_model_map, trans_model, tm)
    if (tm->second)
      delete tm->second;
}

Rivas_transducer_factory* Rivas_transducer_factory::clone()
{
  Rivas_transducer_factory* rtf = new Rivas_transducer_factory (prior.states(), branch.states());
  rtf->prior = prior;
  rtf->branch = branch;
  rtf->assign_subst_model (*this);
  rtf->del_rate = del_rate;
  rtf->mean_del_size = mean_del_size;
  // clone Handel_base
  (Handel_base&) *rtf = (Handel_base&) *this;
  // clone selected Tree_alignment members (this should be in Tree_alignment, seriously...)
  rtf->tree = tree;
  rtf->align = align;
  rtf->node2row = node2row;
  rtf->row2node = row2node;
  rtf->handel_seq_scores = handel_seq_scores;
  rtf->handel_branch_scores = handel_branch_scores;
  // return
  return rtf;
}

void Rivas_transducer_factory::init_trans_model (int src_state, int n_dest_states)
{
  Irrev_EM_matrix*& irrev_ptr = trans_model[src_state];
  if (irrev_ptr)
    delete irrev_ptr;
  irrev_ptr = new Irrev_EM_matrix (1, n_dest_states);
  dest_states[src_state] = vector<int>();
}

int Rivas_transducer_factory::get_dest_state_index (int src_state, int dest_state)
{
  vector<int>& dest_vec = dest_states[src_state];
  if (find (dest_vec.begin(), dest_vec.end(), dest_state) == dest_vec.end())
    dest_vec.push_back (dest_state);
  return find (dest_vec.begin(), dest_vec.end(), dest_state) - dest_vec.begin();
}

void Rivas_transducer_factory::set_trans_prob_init (int src_state, int dest_state, Prob prob)
{
  const int dest_n = get_dest_state_index (src_state, dest_state);

  Irrev_EM_matrix* irrev = trans_model[src_state];
  irrev->pi[dest_n] = prob;
}

void Rivas_transducer_factory::set_trans_prob_rate (int src_state, int old_dest_state, int new_dest_state, double rate)
{
  const int old_dest_n = get_dest_state_index (src_state, old_dest_state);
  const int new_dest_n = get_dest_state_index (src_state, new_dest_state);

  Irrev_EM_matrix* irrev = trans_model[src_state];
  irrev->X[0] (old_dest_n, new_dest_n) += rate;
  irrev->X[0] (old_dest_n, old_dest_n) -= rate;
}

Pair_transducer_scores Rivas_transducer_factory::prior_pair_trans_sc()
{
  prior.alphabet_size = subst_model().m();
  prior.alloc_pair_emit (0);  // prior.state_type already initialized by calls to init_prior_* from constructor

  const vector<Prob> ins = subst_model().create_prior();
  for (int i = 0; i < subst_model().m(); ++i)
    prior.insert_val(0,i) = Prob2Score(ins[i]);

  for (int s = 0; s < prior.states(); ++s)
    if (prior.state_type[s] != TransducerWaitType)
      prior.emit_label[s] = prior.state_name[s];

  return prior;
}

void Rivas_transducer_factory::update_trans_probs()
{
  for_contents (Trans_model_map, trans_model, tm)
    if (tm->second)
      tm->second->update();
}

Pair_transducer_scores Rivas_transducer_factory::branch_pair_trans_sc (double time)
{
  branch.alphabet_size = subst_model().m();
  branch.alloc_pair_emit (0);  // branch.state_type already initialized by calls to init_branch_* from set_transitions_method (called by constructor)

  const vector<Prob> ins = subst_model().create_prior();
  const array2d<Prob> match = subst_model().create_conditional_substitution_matrix (time);

  for (int i = 0; i < subst_model().m(); ++i)
    for (int s = 0; s < branch.states(); ++s)
      switch (branch.state_type[s])
	{
	case TransducerMatchType:
	  for (int j = 0; j < subst_model().m(); ++j)
	    branch.match_val(s,i,j) = Prob2Score(match(i,j));
	  break;

	case TransducerInsertType:
	  branch.insert_val(s,i) = Prob2Score(ins[i]);
	  break;

	case TransducerDeleteType:
	  branch.delete_val(s,i) = 0;
	  break;

	case TransducerWaitType:
	  break;

	default:
	  THROWEXPR ("Unexpected state type");
	  break;
	}

  for_contents (Trans_model_map, trans_model, s_tm)
    {
      const int src_state = s_tm->first;
      Irrev_EM_matrix* tm = s_tm->second;
      const vector<int>& dest_vec = dest_states[src_state];

      const vector<Prob> trans_prob = tm->create_state_vector (time);
      for (int dest_n = 0; dest_n < (int) dest_vec.size(); ++dest_n)
	branch.transition (src_state, dest_vec[dest_n]) = Prob2Score (trans_prob[dest_n]);
    }

  for (int s = 0; s < branch.states(); ++s)
    if (branch.state_type[s] != TransducerWaitType)
      branch.emit_label[s] = branch.state_name[s];

  if (CTAGGING(3,RIVAS))
    {
      CL << "Branch transducer for time " << time << ":\n";
      branch.show (CL);
    }
  return branch;
}

void Affine_transducer_factory::set_transitions(double gamma_val, double delete_rate, double del_extend_prob)
{
  gamma_param = gamma_val;
  delete_rate_param = delete_rate;
  delete_extend_prob_param = del_extend_prob;


  // initialise branch transducer
  init_branch_wait (0, "MWait");   // match-Wait
  init_branch_match (1, "Match");  // Match
  init_branch_wait (2, "DWait");   // del-Wait
  init_branch_del (3, "Delete");    // Del
  init_branch_ins (4, "Insert");    // Ins

  // Detailed balance:       Rate(insertion of size n) = Rate(deletion of size n) * gamma^n
  // Geometric distribution: Rate(deletion of size n) = total_delete_rate * del_extend_prob^{n-1} * (1 - del_extend_prob)
  // Therefore:              total_insert_rate = sum_{n=1}^infty (total_delete_rate * del_extend_prob^{n-1} * (1 - del_extend_prob) * gamma_val^n)
  //                                           = total_delete_rate * gamma_val * (1 - del_extend_prob) / (1 - del_extend_prob * gamma)
  double ins_extend_prob = del_extend_prob * gamma_val;
  double insert_rate = delete_rate * gamma_val * (1 - del_extend_prob) / (1 - ins_extend_prob);

  init_trans_model (Start, 2);
  set_trans_prob_init (Start, 0, 1.);
  set_trans_prob_rate (Start, 0, 4, insert_rate);

  init_trans_model (0, 2);
  set_trans_prob_init (0, 1, 1.);
  set_trans_prob_rate (0, 1, 3, delete_rate);

  init_trans_model (1, 2);
  set_trans_prob_init (1, 0, 1.);
  set_trans_prob_rate (1, 0, 4, insert_rate);

  branch.transition (2, 3) = Prob2Score (1.);

#ifdef ALLOW_DELETE_TO_INSERT_TRANSITIONS
  init_trans_model (3, 3);
  set_trans_prob_init (3, 0, 1. - del_extend_prob);
  set_trans_prob_init (3, 2, del_extend_prob);
  set_trans_prob_rate (3, 0, 4, insert_rate);
#else /* ALLOW_DELETE_TO_INSERT_TRANSITIONS */
  branch.transition (3, 0) = Prob2Score (1. - del_extend_prob);
  branch.transition (3, 2) = Prob2Score (del_extend_prob);
#endif /* ALLOW_DELETE_TO_INSERT_TRANSITIONS */

  branch.transition (4, 0) = Prob2Score (1. - ins_extend_prob);
  branch.transition (4, 4) = Prob2Score (ins_extend_prob);

  update_trans_probs();

  // set up banding coefficients
  del_rate = delete_rate;
  mean_del_size = 1. / (1. - del_extend_prob);
}

Affine_transducer_factory_param_container Affine_transducer_factory::propose_indel_params()
{
  Affine_transducer_factory_param_container container;
  // Sample a normal distribution with variance=indel_proposal_variance
  // and exp it to generate from log-normal
  multiplier = exp(gennor(0, indel_proposal_variance)); 
  //  CTAG(5, PARAM_SAMPLE) << "\nvariance, multiplier, new rate: " << indel_proposal_variance << " " << multiplier << " " << delete_rate_param*multiplier << endl; 

  // Currently holding  fixed all parameters except delete rate param
  // They seem to slow the mixing, and are generally not investigated downstream
  // These could be given their own multipliers as above, but it would
  // complicate things (very) slightly. 
  container.gamma_param = gamma_param; //Rnd::sample_beta (gamma_alpha, gamma_beta);
  container.delete_extend_prob_param = delete_extend_prob_param; //Rnd::sample_beta (delete_extend_alpha, delete_extend_beta);

  // Old sampling from prior:
  //container.delete_rate_param = Rnd::sample_gamma (alpha, beta);

  // New - scale rate parameter by log-normal multiplier
  container.delete_rate_param = delete_rate_param*multiplier; 
  return container;
} 

Convex_transducer_factory_param_container Convex_transducer_factory::propose_indel_params()
{
  Convex_transducer_factory_param_container container;
  container.gamma_param = Rnd::sample_beta (gamma_alpha, gamma_beta);
  
  container.cpt_weight_param = Rnd::sample_dirichlet (cpt_prior.weight_count);
  for (unsigned int cpt = 0; cpt < cpt_weight_param.size(); ++cpt)
    container.cpt_delete_extend.push_back (Rnd::sample_beta (cpt_prior.alpha.at(cpt), cpt_prior.beta.at(cpt)));

  container.delete_rate_param = Rnd::sample_gamma(alpha, beta);
  return container;
}

sstring Affine_transducer_factory::indel_parameter_string() const
{
  sstring param_string;
  param_string << "Eqm_extend: " << gamma_param << ", Delete_rate: " << delete_rate_param << ", Del_extend: " << delete_extend_prob_param;
  return param_string;
}

sstring Convex_transducer_factory::indel_parameter_string() const
{
  sstring param_string;
  param_string << "Eqm_extend: " << gamma_param << ", Delete_rate: " << delete_rate_param
	       << ", Del_size_class:( " << cpt_weight_param << "), Del_extend:( " << cpt_delete_extend << ')';
  return param_string;
}

Score  Affine_transducer_factory::alignment_path_score_with_prior()
{
  Score path_score = alignment_path_score(); 
  path_score += Nats2Score(Math_fn::gamma_density(delete_rate_param, alpha, beta));
  path_score += Nats2Score( Math_fn::beta_density(delete_extend_prob_param, delete_extend_alpha, delete_extend_beta));
  path_score += Nats2Score( Math_fn::beta_density(gamma_param, gamma_alpha, gamma_beta));
  return path_score; 
}

double Affine_transducer_factory::proposal_prob(double old_param, double new_param)
{
  // This assumes a log-normal proposal density
  // and returns the probability of proposing new_param given the chain
  // is at old_param.  
  return Math_fn::normal_density(log(new_param/old_param), 0, indel_proposal_variance); 
}


void Affine_transducer_factory::sample_indel_params()
{
  //  Score old_sc = alignment_path_score();
  // Use score with prior on indel params
  Score old_sc = alignment_path_score_with_prior();

  for (int sample = 0; sample < RIVAS_PRIOR_SAMPLES; ++sample)
    {
      Affine_transducer_factory_param_container
	old_params = *this,
	// Note the variable 'multiplier' gets set within this function:
	new_params = propose_indel_params();

      set_transitions (new_params.gamma_param, new_params.delete_rate_param, new_params.delete_extend_prob_param);
      //      const Score new_sc = alignment_path_score();
      const Score new_sc = alignment_path_score_with_prior();
      
      // Proposal probabilities for hastings ratio.  
      // P(propose(new|old))
      double propose_new_prob = proposal_prob(new_params.delete_rate_param,
					      old_params.delete_rate_param); 
      // P(propose(old|new))
      double propose_old_prob = proposal_prob(old_params.delete_rate_param,
					      new_params.delete_rate_param);
      double hastingsRatio = 
	Score2Prob (ScorePMul (new_sc, -old_sc)) * 
	(propose_old_prob / propose_new_prob) * (1.0/multiplier); // added jacobian multipler as per Lakner et al 2008 Sys Bio
      // Bumped  this down to 4, since 5 is now the default log level.  -OW
      CTAG(4, PARAM_SAMPLE) << "Sampled parameters ("
			    << sample+1 << " of " << RIVAS_PRIOR_SAMPLES
			    << "): gamma=" << gamma_param
			    << " delete_rate=" << delete_rate_param
			    << " delete_extend=" << delete_extend_prob_param
			    << " new_bitscore=" << Score2Bits(new_sc)
			    << " old_bitscore=" << Score2Bits(old_sc)
			    << " old |new proposal: " << propose_old_prob 
			    << " new|old proposal: " << propose_new_prob 
			    << " Jacobian multiplier " << multiplier 
			    << " ; \n\tHastings ratio:"  << hastingsRatio
			    << ";  move ";
      
  
      if (hastingsRatio > 1.0 ? true : Rnd::decide(hastingsRatio))
	{
	  CL << "accepted\n";
	  old_sc = new_sc;
	}
      else
	{
	  CL << "rejected\n";
	  set_transitions (old_params.gamma_param, old_params.delete_rate_param, old_params.delete_extend_prob_param);
	}
    }
}

void Convex_transducer_factory::sample_indel_params()
{
  for (int sample = 0; sample < RIVAS_PRIOR_SAMPLES; ++sample)
    {
      Score old_sc = alignment_path_score();

      Convex_transducer_factory_param_container
	old_params = *this,
	new_params = propose_indel_params();

      set_transitions (new_params.gamma_param, new_params.delete_rate_param, new_params.cpt_weight_param, new_params.cpt_delete_extend);
      const Score new_sc = alignment_path_score();

      CTAG(5, PARAM_SAMPLE) << "Sampled parameters ("
			    << sample+1 << " of " << RIVAS_PRIOR_SAMPLES
			    << "): gamma=" << gamma_param
			    << " delete_rate=" << delete_rate_param
			    << " cpt_weight=(" << cpt_weight_param
			    << ") delete_extend=(" << cpt_delete_extend
			    << ") new_bitscore=" << Score2Bits(new_sc)
			    << " old_bitscore=" << Score2Bits(old_sc)
			    << "; move ";
  
      if (new_sc > old_sc ? true : Rnd::decide (Score2Prob (ScorePMul (new_sc, -old_sc))))
	{
	  CL << "accepted\n";
	  old_sc = new_sc;
	}
      else
	{
	  CL << "rejected\n";
	  set_transitions (old_params.gamma_param, old_params.delete_rate_param, old_params.cpt_weight_param, old_params.cpt_delete_extend);
	}
    }
}

Affine_transducer_factory::Affine_transducer_factory (double gamma, double delete_rate, double del_extend_prob)
  : Rivas_transducer_factory (5, gamma)
{
  set_transitions(gamma, delete_rate, del_extend_prob);
  sort_indels = true;
}

void Convex_transducer_factory::set_transitions (double gamma_val, double delete_rate, const vector<double>& cpt_weight, const vector<double>& cpt_del_extend)
{
  if (cpt_weight.size() != cpt_del_extend.size())
    THROWEXPR ("Bad parameters in Convex_transducer_factory set_transitions");
  const int cpts = cpt_weight.size();
  
  gamma_param = gamma_val;
  delete_rate_param = delete_rate;
  cpt_weight_param = cpt_weight;
  cpt_delete_extend = cpt_del_extend;

  const int mwait = 0;
  const int match = 1;
  init_branch_wait (mwait, "MWait");   // match-Wait
  init_branch_match (match, "Match");  // Match
  vector<int> dwait, del, ins;
  for (int c = 0; c < cpts; ++c)
    {
      const int s = 2 + 3*c;
      sstring dwait_str, del_str, ins_str;
      dwait_str << "DWait" << c+1;
      del_str << "Delete" << c+1;
      ins_str << "Insert" << c+1;

      dwait.push_back (s);
      del.push_back (s+1);
      ins.push_back (s+2);

      init_branch_wait (dwait[c], dwait_str.c_str());   // del-Wait
      init_branch_del (del[c], del_str.c_str());    // Del
      init_branch_ins (ins[c], ins_str.c_str());    // Ins
    }

  vector<double> cpt_ins_extend, cpt_ins_rate, cpt_del_rate;
  for (int c = 0; c < cpts; ++c)
    {
      // Detailed balance:       Rate(insertion of size n) = Rate(deletion of size n) * gamma^n
      // Geometric distribution: Rate(deletion of size n) = total_delete_rate * del_extend_prob^{n-1} * (1 - del_extend_prob)
      // Therefore:              total_insert_rate = sum_{n=1}^infty (total_delete_rate * del_extend_prob^{n-1} * (1 - del_extend_prob) * gamma^n)
      //                                           = total_delete_rate * gamma * (1 - del_extend_prob) / (1 - del_extend_prob * gamma)
      const double c_ins_ext = gamma_val * cpt_del_extend[c];
      const double c_del_rate = delete_rate * cpt_weight[c];
      const double c_ins_rate = c_del_rate * gamma_val * (1 - cpt_del_extend[c]) / (1 - c_ins_ext);
      cpt_ins_extend.push_back (c_ins_ext);
      cpt_del_rate.push_back (c_del_rate);
      cpt_ins_rate.push_back (c_ins_rate);
    }

  init_trans_model (Start, 1 + cpts);
  set_trans_prob_init (Start, mwait, 1.);
  for (int c = 0; c < cpts; ++c)
    set_trans_prob_rate (Start, mwait, ins[c], cpt_ins_rate[c]);

  init_trans_model (mwait, 1 + cpts);
  set_trans_prob_init (mwait, match, 1.);
  for (int c = 0; c < cpts; ++c)
    set_trans_prob_rate (mwait, match, del[c], cpt_del_rate[c]);

  init_trans_model (match, 1 + cpts);
  set_trans_prob_init (match, mwait, 1.);
  for (int c = 0; c < cpts; ++c)
    set_trans_prob_rate (match, mwait, ins[c], cpt_ins_rate[c]);

  for (int c = 0; c < cpts; ++c)
    branch.transition (dwait[c], del[c]) = Prob2Score (1.);

#ifdef ALLOW_DELETE_TO_INSERT_TRANSITIONS
  for (int c = 0; c < cpts; ++c)
    {
      init_trans_model (del[c], 2 + cpts);

      set_trans_prob_init (del[c], mwait, 1. - cpt_del_extend[c]);
      set_trans_prob_init (del[c], dwait[c], cpt_del_extend[c]);
      for (int c2 = 0; c2 < cpts; ++c2)
	set_trans_prob_rate (del[c], mwait, ins[c2], cpt_ins_rate[c2]);
    }
#else /* ALLOW_DELETE_TO_INSERT_TRANSITIONS */
  for (int c = 0; c < cpts; ++c)
    {
      branch.transition (del[c], mwait) = Prob2Score (1. - cpt_del_extend[c]);
      branch.transition (del[c], dwait[c]) = Prob2Score (cpt_del_extend[c]);
    }
#endif /* ALLOW_DELETE_TO_INSERT_TRANSITIONS */

  for (int c = 0; c < cpts; ++c)
    {
      branch.transition (ins[c], mwait) = Prob2Score (1. - cpt_ins_extend[c]);
      branch.transition (ins[c], ins[c]) = Prob2Score (cpt_ins_extend[c]);
    }

  update_trans_probs();

  // set up banding coefficients
  del_rate = delete_rate;
  mean_del_size = 0.;
  for (int c = 0; c < cpts; ++c)
    mean_del_size += cpt_weight[c] / (1. - cpt_del_extend[c]);
}

Convex_transducer_factory::Convex_transducer_factory (double gamma, double delete_rate, const vector<double>& cpt_weight, const vector<double>& cpt_del_extend)
  : Rivas_transducer_factory (2 + 3 * cpt_weight.size(), gamma)
{
  set_transitions (gamma, delete_rate, cpt_weight, cpt_del_extend);
  sort_indels = true;
}
