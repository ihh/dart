#ifndef TKFHMM_INCLUDED
#define TKFHMM_INCLUDED

#include <math.h>
#include <functional>
#include "tkf/tkfparams.h"
#include "util/maximise.h"

// Standalone pair HMM structure
//
struct TKF_joint_pair_HMM_scores : Pair_HMM_scores
{
  TKF_joint_pair_HMM_scores (const TKF_params& params, double time, bool use_indels = TRUE, bool use_subst = TRUE);
};

// Derivatives of pair HMM transition probabilities w.r.t. time
//
struct TKF_joint_pair_HMM_dt : Pair_HMM_derivatives
{
  TKF_joint_pair_HMM_dt (const TKF_params& params, double time, bool use_indels = TRUE, bool use_subst = TRUE);
};

// Counts function objects
//
struct TKF_counts_function : unary_function<const Pair_HMM_scores&,Pair_HMM_counts>
{
  const Score_profile& xseq;
  const Score_profile& yseq;
  const Alphabet* alphabet;
  TKF_counts_function (const Score_profile& xseq, const Score_profile& yseq, const Alphabet* alphabet = 0)
    : xseq(xseq), yseq(yseq), alphabet(alphabet) { }
  virtual ~TKF_counts_function();
  virtual Pair_HMM_counts operator() (const Pair_HMM_scores& scores) = 0;
};


struct TKF_unaligned_counts_function : TKF_counts_function
{
  TKF_unaligned_counts_function (const Score_profile& xseq, const Score_profile& yseq, const Alphabet* alphabet = 0)
    : TKF_counts_function(xseq,yseq,alphabet) { }
  Pair_HMM_counts operator() (const Pair_HMM_scores& scores);
};


struct TKF_aligned_counts_function : TKF_counts_function
{
  Pairwise_path path;

  TKF_aligned_counts_function (const Score_profile& xseq, const Score_profile& yseq,
			       const Pairwise_path& path, const Alphabet* alphabet = 0) :
    TKF_counts_function (xseq, yseq, alphabet),
    path (path)
    { }

  TKF_aligned_counts_function (const Alignment& align, int xrow, int yrow, const Alphabet* alphabet = 0) :
    TKF_counts_function (*align.prof[xrow], *align.prof[yrow], alphabet),
    path (align.path, xrow, yrow, 1)
    { }

  Pair_HMM_counts operator() (const Pair_HMM_scores& scores);
};


// Log-likelihood, derivative & related function objects with encapsulating counts-caching class
//
struct TKF_functions
{
  bool use_indels;
  bool use_subst;

  struct Pair_scores_function : unary_function<double,TKF_joint_pair_HMM_scores>
  {
    TKF_functions& tf;
    Pair_scores_function (TKF_functions& tf) : tf(tf) { }
    result_type operator() (argument_type time) { return TKF_joint_pair_HMM_scores (tf.params, time, tf.use_indels, tf.use_subst); }
  } pair_scores_function;

  struct Pair_dt_function : unary_function<double,TKF_joint_pair_HMM_dt>
  {
    TKF_functions& tf;
    Pair_dt_function (TKF_functions& tf) : tf(tf) { }
    result_type operator() (argument_type time) { return TKF_joint_pair_HMM_dt (tf.params, time, tf.use_indels, tf.use_subst); }
  } pair_dt_function;

  struct Pair_counts_function : unary_function<double,Pair_HMM_counts>
  {
    TKF_functions& tf;
    Pair_counts_function (TKF_functions& tf) : tf(tf) { }
    result_type operator() (argument_type time) { return tf.counts_function (tf.pair_scores(time)); }
  } pair_counts_function;

  struct Log_like_function : unary_function<double,double>
  {
    TKF_functions& tf;
    Log_like_function (TKF_functions& tf) : tf(tf) { }
    result_type operator() (argument_type time)
    {
      const Pair_HMM_counts hmm_counts = tf.pair_counts(time);
      const Loge loglike = hmm_counts.log_likelihood;
      if (CTAGGING(3,TKFHMM TKFHMM_LOGLIKE))
	CL << "At time " << time << ": loglike = " << loglike << "\n";
      return loglike;
    }
  } log_like_function;
  
  struct Log_like_dt_function : unary_function<double,double>
  {
    TKF_functions& tf;
    Log_like_dt_function (TKF_functions& tf) : tf(tf) { }
    result_type operator() (argument_type time)
    {
      const Loge dlog_dt = tf.pair_counts(time).dloglike_dx(tf.pair_scores(time), tf.pair_dt(time));
      if (CTAGGING(3,TKFHMM TKFHMM_LOGLIKE))
	{ CL << "At time " << time << ": d(loglike)/dt = " << dlog_dt << "\n"; }
      return dlog_dt;
    }
  } log_like_dt_function;
  
  const TKF_params&     params;
  TKF_counts_function&  counts_function;

  Function_cache<Pair_scores_function>  pair_scores;
  Function_cache<Pair_dt_function>      pair_dt;
  Function_cache<Pair_counts_function>  pair_counts;
  Function_cache<Log_like_function>     log_like;
  Function_cache<Log_like_dt_function>  log_like_dt;

  TKF_functions (const TKF_params& params, TKF_counts_function& counts_function, bool use_indels, bool use_subst, double resolution) :
    use_indels (use_indels),
    use_subst (use_subst),
    
    pair_scores_function (*this),
    pair_dt_function     (*this),
    pair_counts_function (*this),
    log_like_function    (*this),
    log_like_dt_function (*this),
    
    params (params),
    counts_function (counts_function),

    pair_scores (pair_scores_function, resolution),
    pair_dt     (pair_dt_function, resolution),
    pair_counts (pair_counts_function, resolution),
    log_like    (log_like_function, resolution),
    log_like_dt (log_like_dt_function, resolution)

    { }

};

#endif
