#ifndef LOGPROB_INCLUDED
#define LOGPROB_INCLUDED

#include "util/score.h"

// LogProb:
// Template wrapper allowing logspace types (i.e. Loge & Score) to be handled transparently as though they were Prob's.
// The "Algebra" argument is a class with "value_type" typedef'd to Loge or Score,
// that supports the following static methods:
//   to_prob(), from_prob(), pmul(), psum(), pdiff(), pmulacc(), psumacc(), pdiffacc()
// See below for Algebra examples.
//
template <class Algebra>
class LogProb
{
public:
  // typedefs
  typedef typename Algebra::value_type value_type;

private:
  // value
  value_type val;

public:
  // constructors
  LogProb() { }  // no initialisation, for speed
  LogProb (Prob px) : val(from_prob(px)) { }
  LogProb (const LogProb& lx) : val(lx.val) { }

  // assignment operators
  inline LogProb& operator= (const LogProb& lx) { val = lx.val; return *this; }
  inline LogProb& operator= (Prob px) { val = from_prob(px); return *this; }

  // arithmetic operators; all combinations of LogProb and Prob are covered
  inline friend LogProb operator+ (const LogProb& lx, const LogProb& ly) { return from_log (psum (lx.val, ly.val)); }
  inline friend LogProb operator+ (const LogProb& lx, Prob py) { return from_log (psum (lx.val, from_prob(py))); }
  inline friend LogProb operator+ (Prob px, const LogProb& ly) { return from_log (psum (from_prob(px), ly.val)); }
  inline LogProb& operator+= (const LogProb& lx) { psumacc (val, lx.val); return *this; }
  inline LogProb& operator+= (Prob px) { psumacc (val, from_prob(px)); return *this; }

  inline friend LogProb operator- (const LogProb& lx, const LogProb& ly) { return from_log (pdiff (lx.val, ly.val)); }
  inline friend LogProb operator- (const LogProb& lx, Prob py) { return from_log (pdiff (lx.val, from_prob(py))); }
  inline friend LogProb operator- (Prob px, const LogProb& ly) { return from_log (pdiff (from_prob(px), ly.val)); }
  inline LogProb& operator-= (const LogProb& lx) { pdiffacc (val, lx.val); return *this; }
  inline LogProb& operator-= (Prob px) { pdiffacc (val, from_prob(px)); return *this; }

  inline friend LogProb operator* (const LogProb& lx, const LogProb& ly) { return from_log (pmul (lx.val, ly.val)); }
  inline friend LogProb operator* (const LogProb& lx, Prob py) { return from_log (pmul (lx.val, from_prob(py))); }
  inline friend LogProb operator* (Prob px, const LogProb& ly) { return from_log (pmul (from_prob(px), ly.val)); }
  inline LogProb& operator*= (const LogProb& lx) { pmulacc (val, lx.val); return *this; }
  inline LogProb& operator*= (Prob px) { pmulacc (val, from_prob(px)); return *this; }

  inline friend LogProb operator/ (const LogProb& lx, const LogProb& ly) { return from_log (pmul (lx.val, -ly.val)); }
  inline friend LogProb operator/ (const LogProb& lx, Prob py) { return from_log (pmul (lx.val, -from_prob(py))); }
  inline friend LogProb operator/ (Prob px, const LogProb& ly) { return from_log (pmul (from_prob(px), -ly.val)); }
  inline LogProb& operator/= (const LogProb& lx) { pmulacc (val, -lx.val); return *this; }
  inline LogProb& operator/= (Prob px) { pmulacc (val, -from_prob(px)); return *this; }

  // increment & decremement
  LogProb& operator++() { *this += 1.; return *this; }
  LogProb operator++(int) { LogProb tmp (*this); ++(*this); return tmp; }

  LogProb& operator--() { *this -= 1.; return *this; }
  LogProb operator--(int) { LogProb tmp (*this); --(*this); return tmp; }

  // relational operators
  inline friend int operator== (const LogProb& lx, const LogProb& ly) { return lx.val == ly.val; }
  inline friend int operator== (const LogProb& lx, const Prob py) { return lx.val == from_prob(py); }
  inline friend int operator== (const Prob px, const LogProb& ly) { return from_prob(px) == ly.val; }

  inline friend int operator!= (const LogProb& lx, const LogProb& ly) { return lx.val != ly.val; }
  inline friend int operator!= (const LogProb& lx, const Prob py) { return lx.val != from_prob(py); }
  inline friend int operator!= (const Prob px, const LogProb& ly) { return from_prob(px) != ly.val; }

  inline friend int operator< (const LogProb& lx, const LogProb& ly) { return lx.val < ly.val; }
  inline friend int operator< (const LogProb& lx, const Prob py) { return lx.val < from_prob(py); }
  inline friend int operator< (const Prob px, const LogProb& ly) { return from_prob(px) < ly.val; }

  inline friend int operator> (const LogProb& lx, const LogProb& ly) { return lx.val > ly.val; }
  inline friend int operator> (const LogProb& lx, const Prob py) { return lx.val > from_prob(py); }
  inline friend int operator> (const Prob px, const LogProb& ly) { return from_prob(px) > ly.val; }

  inline friend int operator<= (const LogProb& lx, const LogProb& ly) { return lx.val <= ly.val; }
  inline friend int operator<= (const LogProb& lx, const Prob py) { return lx.val <= from_prob(py); }
  inline friend int operator<= (const Prob px, const LogProb& ly) { return from_prob(px) <= ly.val; }

  inline friend int operator>= (const LogProb& lx, const LogProb& ly) { return lx.val >= ly.val; }
  inline friend int operator>= (const LogProb& lx, const Prob py) { return lx.val >= from_prob(py); }
  inline friend int operator>= (const Prob px, const LogProb& ly) { return from_prob(px) >= ly.val; }

  // stream operators
  inline friend ostream& operator<< (ostream& out, const LogProb& lx) { out << to_prob(lx.val); return out; }
  inline friend istream& operator>> (istream& in, const LogProb& lx) { Prob px; in >> px; lx.val = px; return in; }

  // cast operators
  inline Prob prob() const { return to_prob (val); }
  inline operator double() const { return to_prob (val); }

private:
  // private Algebra method wrappers
  static inline Prob to_prob (value_type X) { return Algebra::to_prob (X); }
  static inline value_type from_prob (Prob P) { return Algebra::from_prob (P); }
  static inline value_type pmul (value_type X, value_type Y) { return Algebra::pmul (X, Y); }
  static inline value_type psum (value_type X, value_type Y) { return Algebra::psum (X, Y); }
  static inline value_type pdiff (value_type X, value_type Y) { return Algebra::pdiff (X, Y); }
  static inline void pmulacc (value_type& X, value_type Y) { Algebra::pmulacc (X, Y); }
  static inline void psumacc (value_type& X, value_type Y) { Algebra::psumacc (X, Y); }
  static inline void pdiffacc (value_type& X, value_type Y) { Algebra::pdiffacc (X, Y); }

  // static constructor from logspace value
  static inline LogProb from_log (value_type X) { LogProb lx; lx.val = X; return lx; }
};

// LogReal: a LogProb with a sign.
//
template <class Algebra>
class LogReal
{
public:
  // typedefs
  typedef typename Algebra::value_type value_type;

private:
  // value
  LogProb<Algebra> val;
  signed short sign;

public:
  // constructors
  LogReal() { }  // no initialisation, for speed
  LogReal (double px) : val(abs(px)), sign(ssgn(px)) { }
  LogReal (const LogReal& lx) : val(lx.val), sign(lx.sign) { }

  // assignment operators
  inline LogReal& operator= (const LogReal& lx) { val = lx.val; sign = lx.sign; return *this; }
  inline LogReal& operator= (double px) { val = abs(px); sign = ssgn(px); return *this; }

  // arithmetic operators; all combinations of LogReal and double are covered
  inline friend LogReal operator+ (const LogReal& lx, const LogReal& ly)
  { return lx.sign == ly.sign ? from_log (lx.val + ly.val, lx.sign) : from_log (lx.val - ly.val, lx.val > ly.val ? lx.sign : ly.sign); }
  inline friend LogReal operator+ (const LogReal& lx, double py) { LogReal ly (py); return lx + ly; }
  inline friend LogReal operator+ (double px, const LogReal& ly) { LogReal lx (px); return lx + ly; }
  inline LogReal& operator+= (const LogReal& lx) { return *this = *this + lx; }
  inline LogReal& operator+= (double px) { LogReal lx (px); return *this = *this + lx; }

  inline friend LogReal operator- (const LogReal& lx, const LogReal& ly)
  { return lx.sign == ly.sign ? from_log (lx.val - ly.val, lx.val > ly.val ? lx.sign : -lx.sign) : from_log (lx.val + ly.val, lx.sign); }
  inline friend LogReal operator- (const LogReal& lx, double py) { LogReal ly (py); return lx - ly; }
  inline friend LogReal operator- (double px, const LogReal& ly) { LogReal lx (px); return lx - ly; }
  inline LogReal& operator-= (const LogReal& lx) { return *this = *this - lx; }
  inline LogReal& operator-= (double px) { LogReal lx (px); return *this = *this - lx; }

  inline friend LogReal operator* (const LogReal& lx, const LogReal& ly) { return from_log (lx.val * ly.val, lx.sign * ly.sign); }
  inline friend LogReal operator* (const LogReal& lx, double py) { return from_log (lx.val * abs(py), lx.sign * ssgn(py)); }
  inline friend LogReal operator* (double px, const LogReal& ly) { return from_log (abs(px) * ly.val, ssgn(px) * ly.sign); }
  inline LogReal& operator*= (const LogReal& lx) { val *= lx.val; sign *= lx.sign; return *this; }
  inline LogReal& operator*= (double px) { val *= abs(px); sign *= ssgn(px); return *this; }

  inline friend LogReal operator/ (const LogReal& lx, const LogReal& ly) { return from_log (lx.val / ly.val, lx.sign * ly.sign); }
  inline friend LogReal operator/ (const LogReal& lx, double py) { return from_log (lx.val / abs(py), lx.sign * ssgn(py)); }
  inline friend LogReal operator/ (double px, const LogReal& ly) { return from_log (abs(px) / ly.val, ssgn(px) * ly.sign); }
  inline LogReal& operator/= (const LogReal& lx) { val /= lx.val; sign *= lx.sign; return *this; }
  inline LogReal& operator/= (double px) { val /= abs(px); sign *= ssgn(px); return *this; }

  // increment & decremement
  LogReal& operator++() { *this += 1.; return *this; }
  LogReal operator++(int) { LogReal tmp (*this); ++(*this); return tmp; }

  LogReal& operator--() { *this -= 1.; return *this; }
  LogReal operator--(int) { LogReal tmp (*this); --(*this); return tmp; }

  // relational operators
  inline friend int operator== (const LogReal& lx, const LogReal& ly) { return lx.val == ly.val && lx.sign == ly.sign; }
  inline friend int operator== (const LogReal& lx, const Prob py) { return lx.val == abs(py) && lx.sign == ssgn(py); }
  inline friend int operator== (const Prob px, const LogReal& ly) { return abs(px) == ly.val && ssgn(px) == ly.sign; }

  inline friend int operator!= (const LogReal& lx, const LogReal& ly) { return !(lx == ly); }
  inline friend int operator!= (const LogReal& lx, const Prob py) { return !(lx == py); }
  inline friend int operator!= (const Prob px, const LogReal& ly) { return !(px == ly); }

  inline friend int operator< (const LogReal& lx, const LogReal& ly)
  { return lx.sign > 0 ? (ly.sign > 0 && lx.val < ly.val) : (ly.sign > 0 || lx.val > ly.val); }
  inline friend int operator< (const LogReal& lx, const Prob py) { LogReal ly (py); return lx < ly; }
  inline friend int operator< (const Prob px, const LogReal& ly) { LogReal lx (px); return lx < ly; }

  inline friend int operator> (const LogReal& lx, const LogReal& ly) { return ly < lx; }
  inline friend int operator> (const LogReal& lx, const Prob py) { return py < lx; }
  inline friend int operator> (const Prob px, const LogReal& ly) { return ly < px; }

  inline friend int operator<= (const LogReal& lx, const LogReal& ly)
  { return lx.sign > 0 ? (ly.sign > 0 && lx.val <= ly.val) : (ly.sign > 0 || lx.val >= ly.val); }
  inline friend int operator<= (const LogReal& lx, const Prob py) { LogReal ly (py); return lx <= ly; }
  inline friend int operator<= (const Prob px, const LogReal& ly) { LogReal lx (px); return lx <= ly; }

  inline friend int operator>= (const LogReal& lx, const LogReal& ly) { return ly <= lx; }
  inline friend int operator>= (const LogReal& lx, const Prob py) { return py <= lx; }
  inline friend int operator>= (const Prob px, const LogReal& ly) { return ly <= px; }

  // stream operators
  inline friend ostream& operator<< (ostream& out, const LogReal& lx) { out << (lx.sign > 0 ? '+' : '-') << lx.val; return out; }
  inline friend istream& operator>> (istream& in, const LogReal& lx) { double px; in >> px; lx = px; return in; }

  // cast operators
  inline double prob() const { return (double) val * (double) sign; }
  inline operator double() const { return (double) val * (double) sign; }

private:
  // static constructor from logspace value and sign
  static inline LogReal from_log (const LogProb<Algebra>& val, signed short sign) { LogReal lx; lx.val = val; lx.sign = sign; return lx; }
  // cast sgn()
  static inline signed short ssgn (double p) { return p >= 0 ? +1 : -1; }
};

// Algebra classes
struct Fast_loge_algebra
{
  typedef Loge value_type;
  static inline Prob to_prob (Loge X) { return Score_fns::loglike2prob(X); }
  static inline Loge from_prob (Prob P) { return Score_fns::prob2loglike(P); }
  static inline Loge pmul (Loge X, Loge Y) { return Score_fns::loge_pr_product(X,Y); }
  static inline Loge psum (Loge X, Loge Y) { return Score_fns::loge_pr_sum(X,Y); }
  static inline Loge pdiff (Loge X, Loge Y) { return Score_fns::loge_pr_diff_slow(X,Y); }
  static inline void pmulacc (Loge& X, Loge Y) { Score_fns::loge_pr_product_accum(X,Y); }
  static inline void psumacc (Loge& X, Loge Y) { Score_fns::loge_pr_sum_accum(X,Y); }
  static inline void pdiffacc (Loge& X, Loge Y) { X = Score_fns::loge_pr_diff_slow(X,Y); }
};
typedef LogProb<Fast_loge_algebra> Loge_prob;
typedef LogReal<Fast_loge_algebra> Loge_real;

struct Slow_loge_algebra
{
  typedef Loge value_type;
  static inline Prob to_prob (Loge X) { return Score_fns::loglike2prob(X); }
  static inline Loge from_prob (Prob P) { return Score_fns::prob2loglike(P); }
  static inline Loge pmul (Loge X, Loge Y) { return Score_fns::loge_pr_product(X,Y); }
  static inline Loge psum (Loge X, Loge Y) { return Score_fns::loge_pr_sum_slow(X,Y); }
  static inline Loge pdiff (Loge X, Loge Y) { return Score_fns::loge_pr_diff_slow(X,Y); }
  static inline void pmulacc (Loge& X, Loge Y) { Score_fns::loge_pr_product_accum(X,Y); }
  static inline void psumacc (Loge& X, Loge Y) { Score_fns::loge_pr_sum_accum_slow(X,Y); }
  static inline void pdiffacc (Loge& X, Loge Y) { X = Score_fns::loge_pr_diff_slow(X,Y); }
};
typedef LogProb<Slow_loge_algebra> Slow_loge_prob;
typedef LogReal<Slow_loge_algebra> Slow_loge_real;

struct Fast_score_algebra
{
  typedef Score value_type;
  static inline Prob to_prob (Score X) { return Score_fns::score2prob(X); }
  static inline Score from_prob (Prob P) { return Score_fns::prob2score(P); }
  static inline Score pmul (Score X, Score Y) { return Score_fns::score_pr_product(X,Y); }
  static inline Score psum (Score X, Score Y) { return Score_fns::score_pr_sum(X,Y); }
  static inline Score pdiff (Score X, Score Y) { return Score_fns::score_pr_diff_slow(X,Y); }
  static inline void pmulacc (Score& X, Score Y) { Score_fns::score_pr_product_accum(X,Y); }
  static inline void psumacc (Score& X, Score Y) { Score_fns::score_pr_sum_accum(X,Y); }
  static inline void pdiffacc (Score& X, Score Y) { X = Score_fns::score_pr_diff_slow(X,Y); }
};
typedef LogProb<Fast_score_algebra> Score_prob;
typedef LogReal<Fast_score_algebra> Score_real;

#endif /* LOGPROB_INCLUDED */
