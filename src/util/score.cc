#include "util/score.h"

// nifty_counter
int Score_fns::nifty_counter = 0;

// lookup tables (null for now)
Prob* Score_fns::score2prob_lookup = 0;
Score* Score_fns::score_pr_sum_lookup = 0;

// overflow counts
unsigned int Score_fns::overflows = 0;
unsigned int Score_fns::overflow_notify = 10000;

// and now, the constructor
Score_fns::Score_fns()
{
  // only initialise once
  if (nifty_counter++ == 0) {
    // check #define's
    if (Log10DartScore2BitsRatio != (int) Log10DartScore2BitsRatio)
      THROWEXPR ("DartScore2BitsRatio expected to be a power of 10");

    // initialise lookup tables
    if (score2prob_lookup != 0) delete[] score2prob_lookup;

    score2prob_lookup = new Prob[DartScore2ProbLookup];
    for (int x = 0; x < DartScore2ProbLookup; ++x)
      score2prob_lookup[x] = loglike2prob (score2loglike ((Score) -x));

    if (score_pr_sum_lookup != 0) delete[] score_pr_sum_lookup;

    score_pr_sum_lookup = new Score[DartScorePSumLookup];
    for (int x = 0; x < DartScorePSumLookup; ++x)
      score_pr_sum_lookup[x] = prob2score(1 + score2prob((Score) -x));
  }
}

void Score_fns::describe_scoring_scheme (ostream& o)
{  
  o << "Calculated score lookup tables (scores stored in 1/" << DartScore2BitsRatio << " bit units, inf=" << InfinityScore << ")\n";
}

void Score_fns::show_score (Score sc, ostream& o)
{
  if (sc <= -InfinityScore) {
    if (sc == -InfinityScore) o << "-oo";
    else o << "-OVERFLOW";
  } else if (sc >= InfinityScore) {
    if (sc == InfinityScore) o << "+oo";
    else o << "+OVERFLOW";
  } else o << sc;
}

Score_fns::~Score_fns()
{
  // only cleanup once
  if (--nifty_counter == 0) {
    // following warning message is commented out cos basically it's just a pain in the arse, especially when using multiple processes
    // if (overflows > 10000) cerr << "warning -- " << overflows << " log(1+e^-x) lookup overflows\n";
    delete[] score_pr_sum_lookup;
    delete[] score2prob_lookup;
  }
}

vector<Prob> Score_fns::score2boltzmann_vector (const vector<Score>& sc, const double kT)
{
  Score max_sc = *max_element (sc.begin(), sc.end());
  vector<Prob> b (sc.size());
  for (int i = 0; i < (int) sc.size(); i++)
    b[i] = sc[i] > -InfinityScore ? score2boltzmann (sc[i] - max_sc, kT) : 0;
  normalise_pr (b);
  return b;
}

vector<Prob> Score_fns::score2prob_vector_normalised (const vector<Score>& sc_vector)
{
  Score max_sc = *max_element (sc_vector.begin(), sc_vector.end());
  vector<Prob> v (sc_vector.size());
  for (int i = 0; i < (int) sc_vector.size(); i++)
    v[i] = score2prob (sc_vector[i] - max_sc);
  normalise_pr (v);
  return v;
}

vector<Prob> Score_fns::score2prob_vector_unnormalised (const vector<Score>& sc_vector)
{
  vector<Prob> v (sc_vector.size());
  for (int i = 0; i < (int) sc_vector.size(); i++)
    v[i] = score2prob (sc_vector[i]);
  return v;
}

void Score_fns::normalise_pr (vector<Prob>& weight_vector)
{
  double n = accumulate (weight_vector.begin(), weight_vector.end(), 0.0);
  if (n != 0) for_contents (vector<Prob>, weight_vector, w) *w = *w / n;
}

void Score_fns::normalise_sc (vector<Score>& score_vector)
{
  multiset<Score> sorted_scores (score_vector.begin(), score_vector.end());
  Score n = score_pr_sum (sorted_scores);
  if (n > -InfinityScore) for_contents (vector<Score>, score_vector, s) score_pr_product_accum (*s, -n);
}

void Score_fns::cancelled_infinities_error()
{
  CLOGERR << "Attempt to cancel infinities!\n";
}

sstring Score_fns::score2bits_string (Score sc)
{
  sstring s;
  if (sc < 0) { sc = -sc; s << '-'; }
  const int div = sc / DartScore2BitsRatio;
  const int mod = sc % DartScore2BitsRatio;
  if (div != 0 || mod == 0) s << div;
  if (mod != 0) s << '.' << mod;
  return s;
}

Regexp minus_bits_re ("^-");
Regexp integer_bits_re ("^(-?)([0123456789]*)$");
Regexp fractional_bits_re ("^(-?)([0123456789]*)\\.([0123456789]*)$");
Score Score_fns::bits2score_string (const sstring& str)
{
  Score sc;
  if (integer_bits_re.Match (str.c_str()))
    sc = integer_bits_re[1].to_int();
  else if (fractional_bits_re.Match (str.c_str()))
    {
      const int div = fractional_bits_re[1].to_int();
      const int modlen = fractional_bits_re[2].size();
      int mod = fractional_bits_re[2].to_int();
      for (int i = modlen; i > Log10DartScore2BitsRatio; --i)
	mod = mod / DartScore2BitsRatio;
      for (int i = modlen; i < Log10DartScore2BitsRatio; ++i)
	mod = mod * DartScore2BitsRatio;
      sc = div * DartScore2BitsRatio + mod;
    }
  else
    THROWEXPR ("Malformed bitscore: " << str);
  if (minus_bits_re.Match (str.c_str()))
    sc = -sc;
  return sc;
}
