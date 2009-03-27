#ifndef SCORE_INCLUDED
#define SCORE_INCLUDED

#include <math.h>
#include <algorithm>
#include <vector>
#include <set>
#include "util/array2d.h"
#include "util/logfile.h"
#include "util/rnd.h"
#include "util/math_fn.h"

// typedefs

typedef int    Score; // Scores really HAVE to be integers, because DP tracebacks do equality tests which may fail on floating points
typedef double Prob;
typedef long double Rate;  /* changed from "double" to "long double" to fix xrate large-matrix precision problems at suggestion of Peter Schattner -- 12/15/2008, IH   */
typedef double Loge;
typedef double Log2;

// public global constant macros
// the only guaranteed relationship between the infinity's is:    Score2Nats(InfinityScore) <= InfinityLoge
// i.e. they are NOT necessarily equal!

#define Loge2         0.693147180559945
#define InfinityScore 987654321
#define InfinityLoge  9.87654321e+99
#define InfinityProb  1e+300
#define ZeroProb      1e-300

// private global constant macros
// DartScore2BitsRatio is assumed to be a power of 10 by string display routines

#define DartScore2BitsRatio   1000
#define Log10DartScore2BitsRatio 3
#define DartScorePSumLookup   20000
#define DartScorePSumOverflow 40000
#define DartScore2ProbLookup  20000
#define DartScore2ProbTiny    9.5367431640625e-07  /* = Nats2Prob (Score2Nats (-DartScore2ProbLookup)) */
#define DartScore2NatsRatio   (DartScore2BitsRatio / Loge2)
#define DartLogeHuge          700         /* roughly Prob2Nats(InfinityProb) */ 
#define DartScoreHuge         1010000     /* roughly Prob2Score(InfinityProb) */ 

// basic score conversion macros

#define Prob2Score(P) Score_fns::prob2score(P)
#define Score2Prob(S) Score_fns::score2prob(S)

#define Nats2Score(N) Score_fns::loglike2score(N)
#define Score2Nats(S) Score_fns::score2loglike(S)

#define Score2Bits(S) Score_fns::score2bits(S)
#define Bits2Score(B) Score_fns::bits2score(B)

#define Nats2Bits(N)  Score_fns::loglike2bits(N)
#define Bits2Nats(B)  Score_fns::bits2loglike(B)

#define Prob2Nats(P)  Score_fns::prob2loglike(P)
#define Nats2Prob(N)  Score_fns::loglike2prob(N)

#define Bits2Prob(B)  Score_fns::bits2prob(B)
#define Prob2Bits(P)  Score_fns::prob2bits(P)

// slightly more exotic conversion macros

#define ScaleScore(SC,BETA) Score_fns::scale_score(SC,BETA)

#define Prob2ScoreVec(PV)       Score_fns::prob2score_vector(PV)
#define Score2ProbVecNorm(SV)   Score_fns::score2prob_vector_normalised(SV)
#define Score2ProbVecUnnorm(SV) Score_fns::score2prob_vector_unnormalised(SV)

#define Prob2ScoreArray2d(PA) Score_fns::prob2score_array2d(PA)
#define Score2ProbArray2d(SA) Score_fns::score2prob_array2d(SA)

#define Prob2NatsArray2d(PA) Score_fns::prob2loglike_array2d(PA)
#define Nats2ProbArray2d(NA) Score_fns::loglike2prob_array2d(NA)

#define Score2Boltzmann(S,KT)     Score_fns::score2boltzmann(S,KT)
#define Score2BoltzmannVec(SV,KT) Score_fns::score2boltzmann_vector(SV,KT)

#define NormaliseSc(SVEC) Score_fns::normalise_sc(SVEC)
#define NormalisePr(PVEC) Score_fns::normalise_pr(PVEC)

// probability-space arithmetic macros

#define ScorePMul(A,B)    Score_fns::score_pr_product(A,B)
#define ScorePMul3(A,B,C) Score_fns::score_pr_product(A,B,C)
#define ScorePMulAcc(A,B) Score_fns::score_pr_product_accum(A,B)

#define ScorePSum(A,B)    Score_fns::score_pr_sum(A,B)
#define ScorePSum3(A,B,C) Score_fns::score_pr_sum(A,B,C)
#define ScorePSumAcc(A,B) Score_fns::score_pr_sum_accum(A,B)
#define ScorePSumSet(SET) Score_fns::score_pr_sum(SET)

#define ScorePDiff(A,B)   Score_fns::score_pr_diff_slow(A,B)

#define NatsPMul(A,B)     Score_fns::loge_pr_product(A,B)
#define NatsPMul3(A,B,C)  Score_fns::loge_pr_product(A,B,C)
#define NatsPMulAcc(A,B)  Score_fns::loge_pr_product_accum(A,B)

#define NatsPSum(A,B)     Score_fns::loge_pr_sum(A,B)
#define NatsPSum3(A,B,C)  Score_fns::loge_pr_sum(A,B,C)
#define NatsPSumAcc(A,B)  Score_fns::loge_pr_sum_accum(A,B)
#define NatsPSumSet(SET)  Score_fns::loge_pr_sum(SET)

#define NatsPSumSlow(A,B)    Score_fns::loge_pr_sum_slow(A,B)
#define NatsPSumAccSlow(A,B) Score_fns::loge_pr_sum_accum_slow(A,B)
#define NatsPSumSetSlow(SET) Score_fns::loge_pr_sum_slow(SET)

#define NatsPDiff(A,B)    Score_fns::loge_pr_diff_slow(A,B)

// score display macro; prints "Infinity" or "OVERFLOW" where appropriate

#define ShowScore(SC,OUT) Score_fns::show_score(SC,OUT)

// score<->log2 string conversion macros

#define Score2BitsString(SC) Score_fns::score2bits_string(SC)
#define Bits2ScoreString(SC) Score_fns::bits2score_string(SC)

class Score_fns
{
private:

  // nifty counter (ensures singleton initialisation)
  static int nifty_counter;

  // lookup tables
  static Prob*  score2prob_lookup;
  static Score* score_pr_sum_lookup;

  // overflow counters
  static unsigned int overflows;
  static unsigned int overflow_notify;

  static inline void inc_overflows (int lookup)
  {
    if (lookup < DartScorePSumOverflow)
      if (++overflows == overflow_notify)
	{
	  CLOGERR << "warning -- hit " << overflows << " log(1+e^-x) lookup overflows; consider using a multiset\n";
	  overflow_notify *= 100;
	}
  }

  static void cancelled_infinities_error();

public:

  static inline Score prob2score (const Prob p)
  { return p >= InfinityProb ? InfinityScore : (p <= ZeroProb ? -InfinityScore : (Score) (DartScore2NatsRatio * log(p))); }
  static inline Prob  score2prob (const Score sc)  // exp(...) calls in this function were formerly DartScore2ProbTiny
    {
      if (sc <= 0)
	{
	  if (sc <= -DartScore2ProbLookup)
	    return sc <= -DartScoreHuge ? 0. : exp ((double) sc / DartScore2NatsRatio);
	  return score2prob_lookup[-sc];
	}
      else
	return sc >= DartScore2ProbLookup ? (sc >= DartScoreHuge ? InfinityProb : exp ((double) sc / DartScore2NatsRatio)) : 1. / score2prob_lookup[sc];
    }

  static inline Prob  score2boltzmann (const Score sc, const double kT) { return score2prob ((Score) (((double)sc) / kT)); }
  static vector<Prob> score2boltzmann_vector (const vector<Score>& sc, const double kT);

  static inline Score loglike2score (const Loge loglike)
    {
      return (Score) minmax (loglike * DartScore2NatsRatio, (Loge) -InfinityScore, (Loge) InfinityScore);
    }
  static inline Loge  score2loglike (const Score sc)
    {
      if (sc >= InfinityScore) return InfinityLoge;
      if (sc <= -InfinityScore) return -InfinityLoge;
      return ((Loge) sc) / DartScore2NatsRatio;
    }
  
  static inline Log2  score2bits (const Score score) { return ((double) score) / ((double) DartScore2BitsRatio); }
  static inline Score bits2score (const Log2 bits) { return (Score) round (bits * ((double) DartScore2BitsRatio)); }

  static inline Log2 loglike2bits (const Loge loglike) { return loglike / Loge2; }
  static inline Loge bits2loglike (const Log2 bits) { return bits * Loge2; }

  static inline Prob loglike2prob (const Loge loglike) { return loglike >= DartLogeHuge ? InfinityProb : (loglike <= -DartLogeHuge ? 0. : exp(loglike)); }
  static inline Loge prob2loglike (const Prob prob) { return prob <= ZeroProb ? -InfinityLoge : log(prob); }
  
  static inline Prob bits2prob (const Log2 bits) { return loglike2prob (bits2loglike (bits)); }
  static inline Log2 prob2bits (const Prob prob) { return loglike2bits (prob2loglike (prob)); }

  static inline void scale_score (Score& score, const double beta) { if (score > -InfinityScore && score < InfinityScore) score = (Score) (beta * (double) score); }
  
  static inline vector<Score> prob2score_vector (const vector<Prob>& p_vector)
    {
      vector<Score> v (p_vector.size());
      for (int i = 0; i < (int) p_vector.size(); ++i) v[i] = prob2score (p_vector[i]);
      return v;
    }

  static vector<Prob> score2prob_vector_normalised (const vector<Score>& sc_vector);    // normalises probabilities
  static vector<Prob> score2prob_vector_unnormalised (const vector<Score>& sc_vector);  // doesn't normalise probabilities

  static void normalise_pr (vector<Prob>& weight_vector);
  static void normalise_sc (vector<Score>& score_vector);

  static inline array2d<Score> prob2score_array2d (const array2d<Prob>& p_array)
    {
      array2d<Score> a (p_array.xsize(), p_array.ysize());
      for (int x = 0; x < p_array.xsize(); ++x)
	for (int y = 0; y < p_array.ysize(); ++y)
	  a(x,y) = prob2score (p_array(x,y));
      return a;
    }

  static inline array2d<Prob> loglike2prob_array2d (const array2d<Loge>& ll_array)
    {
      array2d<Prob> p (ll_array.xsize(), ll_array.ysize());
      for (int x = 0; x < ll_array.xsize(); ++x)
	for (int y = 0; y < ll_array.ysize(); ++y)
	  p(x,y) = loglike2prob (ll_array(x,y));
      return p;
    }

  // probability-space arithmetic for Scores

  static inline Score score_pr_product (const Score a, const Score b)
  {
#ifdef DART_DEBUG
    if (a <= -InfinityScore)
      {
	if (b >= InfinityScore) cancelled_infinities_error();
	return -InfinityScore;
      }
    else if (b <= -InfinityScore)
      {
	if (a >= InfinityScore) cancelled_infinities_error();
	return -InfinityScore;
      }
    else if (a >= InfinityScore || b >= InfinityScore)
      return InfinityScore;
    return a + b;
#else /* DART_DEBUG */
    if (a >= InfinityScore || b >= InfinityScore)
      return InfinityScore;
    if (a <= -InfinityScore || b <= -InfinityScore)
      return -InfinityScore;
    return a + b;
#endif /* DART_DEBUG */
  }
  static inline Score score_pr_product (const Score a, const Score b, const Score c) { return score_pr_product (score_pr_product (a, b), c); }    // just for convenience
  static inline void  score_pr_product_accum (Score&a, const Score b) { a = score_pr_product(a,b); }

  static Score score_pr_sum (const Score a, const Score b)
  {
    if (a <= -InfinityScore) return b;
    else if (b <= -InfinityScore) return a;
    else if (a >= InfinityScore || b >= InfinityScore) return InfinityScore;
    Score diff = a - b;
    if (diff > 0) {   //  a > b
      if (diff > DartScorePSumLookup) {
#if 0  /* commented out, because the warnings started to get on my nerves */
#ifdef DART_DEBUG
	inc_overflows (diff);
#endif /* DART_DEBUG */
#endif /* 0 */
	return a;
      } else {
	return a + score_pr_sum_lookup[diff];
      }
    } else {          //  b > a
      if (diff < -DartScorePSumLookup) {
#if 0  /* commented out, because the warnings started to get on my nerves */
#ifdef DART_DEBUG
	inc_overflows (-diff);
#endif /* DART_DEBUG */
#endif /* 0 */
	return b;
      } else {
	return b + score_pr_sum_lookup[-diff];
      }
    }
    return -InfinityScore - 1; // should never get here
  }
  static Score score_pr_sum (const Score a, const Score b, const Score c) { return score_pr_sum (score_pr_sum (a, b), c); }
  static inline void score_pr_sum_accum (Score& a, const Score b) { a = score_pr_sum (a, b); }

  static Score score_pr_sum (multiset<Score>& x)     // accurate, efficient method for summing long series; corrupts x
    {
      if (x.empty()) return -InfinityScore;
      while (x.size() > 1)
	{
	  multiset<Score>::iterator i = x.begin();
	  multiset<Score>::iterator j = i;
	  ++j;
	  const Score ij_sum = score_pr_sum (*i, *j);
	  x.erase (i, ++j);
	  x.insert (ij_sum);
	}
      return *(x.begin());
    }

  static Score score_pr_diff_slow (const Score a, const Score b)
  { return max (a, b) + prob2score (1. - score2prob (-abs (b - a))); }

  // the same probability-space arithmetic, but for Loge's

  static inline Loge loge_pr_product (const Loge a, const Loge b) { return a + b; }  // no real danger of overflow
  static inline Loge loge_pr_product (const Loge a, const Loge b, Loge c) { return a + b + c; }
  static inline void loge_pr_product_accum (Loge& a, const Loge b) { a = loge_pr_product(a,b); }

  static Loge loge_pr_sum (const Loge a, const Loge b)
    {
      Score diff = loglike2score (a - b);
      if (diff > 0) {
	if (diff > DartScorePSumLookup) {
#if 0  /* commented out, because the warnings started to get on my nerves */
#ifdef DART_DEBUG
	  inc_overflows (diff);
#endif /* DART_DEBUG */
#endif /* 0 */
	  return a;
	} else {
	  return a + score2loglike (score_pr_sum_lookup[diff]);
	}
      } else {
	if (diff < -DartScorePSumLookup) {
#if 0  /* commented out, because the warnings started to get on my nerves */
#ifdef DART_DEBUG
	  inc_overflows (-diff);
#endif /* DART_DEBUG */
#endif /* 0 */
	  return b;
	} else {
	  return b + score2loglike (score_pr_sum_lookup[-diff]);
	}
      }
    }
  static inline void loge_pr_sum_accum (Loge& a, const Loge b) { a = loge_pr_sum (a, b); }

  static Loge loge_pr_sum (multiset<Loge>& x)    // accurate(ish), efficient method for summing long series; corrupts x
    {
      if (x.empty()) return -InfinityLoge;
      while (x.size() > 1)
	{
	  multiset<Loge>::iterator i = x.begin();
	  multiset<Loge>::iterator j = i;
	  ++j;
	  const Loge ij_sum = loge_pr_sum (*i, *j);
	  x.erase (i, ++j);
	  x.insert (ij_sum);
	}
      return *(x.begin());
    }

  // slow (but hopefully more accurate) probability-sum methods for Loge's
  
  static inline Loge loge_pr_sum_slow (const Loge a, const Loge b) { return a > b ? a + log(1+exp(b-a)) : b + log(1+exp(a-b)); }
  static inline void loge_pr_sum_accum_slow (Loge& a, const Loge b) { a = loge_pr_sum_slow (a, b); }

  static Loge loge_pr_sum_slow (multiset<Loge>& x)
    {
      if (x.empty()) return -InfinityLoge;
      while (x.size() > 1)
	{
	  multiset<Loge>::iterator i = x.begin();
	  multiset<Loge>::iterator j = i;
	  ++j;
	  const Loge ij_sum = loge_pr_sum_slow (*i, *j);
	  x.erase (i, ++j);
	  x.insert (ij_sum);
	}
      return *(x.begin());
    }

  static Loge loge_pr_diff_slow (const Loge a, const Loge b)
  { return max (a, b) + prob2loglike (1. - loglike2prob (-abs (b - a))); }

  // display methods

  static void describe_scoring_scheme (ostream& o);      // reports the internal format to the logfile
  static void show_score (Score sc, ostream& o);         // displays the score, or "infinity", or "OVERFLOW"

  // string conversion methods
  static sstring score2bits_string (Score sc);
  static Score bits2score_string (const sstring& str);

  // constructor & destructor
  Score_fns();
  ~Score_fns();
};

// dummy initialiser, so that Score_fns constructor is called
static Score_fns _dummy_score_initialiser;

// Kronecker delta class for Scores

class Kronecker_score : public vector<Score>
{
public:
  Kronecker_score (int size, int val) : vector<Score> (size, (Score) -InfinityScore) { (*this)[val] = 0; }
};

#endif
