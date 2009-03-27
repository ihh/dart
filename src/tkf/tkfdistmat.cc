#include "tkf/tkfdistmat.h"
#include "tkf/tkfhmm.h"
#include "util/maximise.h"

TKF_dist_func_factory::TKF_dist_func* TKF_dist_func_factory::create_dist_func (const Alignment& align)
{
  return new TKF_dist_func (*this, align);
}

TKF_dist_func_factory::TKF_dist_func_factory (const TKF_params& params,
					      bool use_indels, bool use_subst,
					      double tres, double tmax,
					      int n_realign, double default_time)
  : params (params),
    use_indels (use_indels), use_subst (use_subst),
    tres (tres), tmax (tmax),
    n_realign (n_realign), default_time (default_time)
{ }

double TKF_dist_func_factory::TKF_dist_func::operator() (int i, int j)
{
  // get i & j Score_profiles
  const Score_profile& iprof = *align.prof[i];
  const Score_profile& jprof = *align.prof[j];

  // create the distance function object
  TKF_counts_function* f = 0;
  if (factory.n_realign < 0)
    f = new TKF_unaligned_counts_function (iprof, jprof, &factory.params.submat_factory.alphabet());
  else if (factory.n_realign == 0)
    f = new TKF_aligned_counts_function (align, i, j, &factory.params.submat_factory.alphabet());
  else if (factory.n_realign > 0)
    {
      // estimate a pairwise alignment
      TKF_joint_pair_HMM_scores tkf_pair_hmm (factory.params, factory.default_time);
      Pair_Viterbi_DP_matrix tkf_pair_matrix (tkf_pair_hmm, iprof, jprof);
      vector<int> vtraceback = tkf_pair_matrix.optimal_state_path();
      Pairwise_path vpath = tkf_pair_hmm.convert_state_path_to_alignment (vtraceback);
      // create distance function from alignment
      f = new TKF_aligned_counts_function (iprof, jprof, vpath, &factory.params.submat_factory.alphabet());
    }

  // create TKF_functions object
  TKF_functions funcs (factory.params, *f, factory.use_indels, factory.use_subst, factory.tres);

  // do Brent maximisation
  double t1, t2, t3, tbest, fbest;
  bracket_maximum (funcs.log_like, t1, t2, t3, 0., factory.tmax);
  brent_deriv (funcs.log_like, funcs.log_like_dt, t1, t2, t3, factory.tres, tbest, fbest);

  // delete the distance function object
  delete f;

  // return the ML distance estimate
  return tbest;
}
