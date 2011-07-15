#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "seq/biosequence.h"
#include "util/istream_get_iterator.h"
#include "util/logfile.h"
#include "util/Regexp.h"
#include "util/vector_output.h"
#include "util/nstring.h"

#define CHAR2INT_TABLE_SIZE 256

int Digitized_biosequence::null_score (const vector<int>& null_model) const
{
  int sc = 0;
  for_const_contents (Digitized_biosequence, *this, i)
    if (*i < 0)
      THROW Standard_exception ("Can't handle ambiguous symbols yet");
    else
      ScorePMulAcc (sc, null_model[*i]);
  return sc;
}

int Digitized_biosequence::stuttered_null_score (const vector<int>& null_model, const int stutter_sc) const
{
  int sc = 0;
  int prev_sym = back();
  for_const_contents (Digitized_biosequence, *this, i)
    {
      if (*i < 0)
	THROW Standard_exception ("Can't handle ambiguous symbols yet");
      else
	ScorePMulAcc (sc, *i == prev_sym ? stutter_sc : null_model[*i]);
      prev_sym = *i;
    }
  return sc;
}

void Digitized_biosequence::count_symbol_frequencies (vector<double>& composition) const
{
  for_const_contents (Digitized_biosequence, *this, i)
    if (*i < 0)
      THROW Standard_exception ("Can't handle ambiguous symbols yet");
    else
      composition[*i]++;
}

Digitized_biosequence Digitized_biosequence::make_context_dsq (int alphabet_size, int context_order) const
{
  Digitized_biosequence dsq (size());
  for (int i = 0; i < (int) size(); ++i) {
    if ((*this)[i] < 0 || (*this)[i] >= alphabet_size) THROW Standard_exception ("Bad symbol in sequence");
    int context = 0;
    int mul = 1;
    for (int j = 0; j <= context_order; ++j) {
      int pos = i - j;
      while (pos < 0) pos += size();
      context += mul * (*this)[pos];
      mul *= alphabet_size;
    }
    dsq[i] = context;
  }
  return dsq;
}

Digitized_biosequence Digitized_biosequence::make_reverse_context_dsq (int alphabet_size, int context_order) const
{
  Digitized_biosequence dsq (size());
  for (int i = 0; i < (int) size(); ++i) {
    if ((*this)[i] < 0 || (*this)[i] >= alphabet_size) THROW Standard_exception ("Bad symbol in sequence");
    int context = 0;
    int mul = 1;
    for (int j = 0; j <= context_order; ++j) {
      int pos = i + j;
      while (pos >= (int) size()) pos -= size();
      context += mul * (*this)[pos];
      mul *= alphabet_size;
    }
    dsq [size() - 1 - i] = context;
  }
  return dsq;
}

Digitized_biosequence Digitized_biosequence::context_to_nmer (int context, int alphabet_size, int context_order)
{
  Digitized_biosequence dsq;
  for (int b = context_order; b >= 0; --b)
    dsq.push_back ((context / (int) pow((double)alphabet_size,(double)b)) % alphabet_size);
  return dsq;
}

Weight_profile::Weight_profile (const Digitized_biosequence& dsq)
{
  clear();
  reserve (dsq.size());
  for_const_contents (Digitized_biosequence, dsq, i)
    {
      Symbol_weight_map m;
      m[*i] = 1;
      push_back (m);
    }
}

int Weight_profile::weight_map_consensus (const Symbol_weight_map& w)
{
  if (w.size() == 0) THROW Standard_exception ("Attempt to take consensus of a zero-information weight profile");
  Symbol_weight best (-1, 0);
  for_const_contents (Symbol_weight_map, w, sw)
    if ((*sw).second > best.second) best = *sw;
  return best.first;
}

Digitized_biosequence Weight_profile::consensus_dsq() const
{
  Digitized_biosequence dsq (size(), -1);
  for (int pos = 0; pos < (int) size(); pos++) dsq[pos] = consensus(pos);
  return dsq;
}

void Weight_profile::normalise()
{
  for_contents (Weight_profile, *this, swm)
    {
      Prob norm = 0;
      for_contents (Symbol_weight_map, *swm, sw) norm += (*sw).second;
      if (norm > 0) for_contents (Symbol_weight_map, *swm, sw) (*sw).second /= norm;
    }
}

void Weight_profile::count_symbol_frequencies (vector<double>& composition) const
{
  for_const_contents (Weight_profile, *this, ssm)
    for_const_contents (Symbol_weight_map, *ssm, ss)
    composition [(*ss).first] += (*ss).second;
}

Score_profile::Score_profile (const Digitized_biosequence& dsq)
{
  clear();
  reserve (dsq.size());
  for_const_contents (Digitized_biosequence, dsq, i)
    {
      Symbol_score_map m;
      m[*i] = 0;
      push_back (m);
    }
}

Score_profile::Score_profile (const Weight_profile& profile)
{
  clear();
  reserve (profile.size());
  for_const_contents (Weight_profile, profile, swm)
    {
      Symbol_score_map ssm;
      for_const_contents (Symbol_weight_map, *swm, sw)
	ssm.insert (Symbol_score ((*sw).first, Prob2Score ((*sw).second)));
      push_back (ssm);
    }
}

int Score_profile::score_map_consensus (const Symbol_score_map& s)
{
  if (s.size() == 0) THROW Standard_exception ("Attempt to take consensus of a zero-information profile entry");
  Symbol_score best (-1, 0);
  bool first = 1;
  for_const_contents (Symbol_score_map, s, ss)
    if (first || (*ss).second > best.second) { best = *ss; first = 0; }
  return best.first;
}

int Score_profile::score_map_sample (const Symbol_score_map& s, double kT)
{
  if (s.size() == 0) THROW Standard_exception ("Attempt to sample from a zero-information profile entry");
  vector<int> res;
  vector<int> sc;
  for_const_contents (Symbol_score_map, s, ss)
    {
      res.push_back ((*ss).first);
      sc.push_back ((*ss).second);
    }
  return res [Rnd::choose (Score2BoltzmannVec (sc, kT))];
}

int Score_profile::score_map_score (const Symbol_score_map& ssm, const vector<int>& null_model)
{
  multiset<int> series;
  for_const_contents (Symbol_score_map, ssm, ss)
    series.insert (null_model[(*ss).first] + (*ss).second);
  return ScorePSumSet (series);
}

int Score_profile::null_score (const vector<int>& null_model) const
{
  int sc = 0;
  for_const_contents (Score_profile, *this, ssm)
    ScorePMulAcc (sc, score_map_score (*ssm, null_model));
  return sc;
}

void Score_profile::normalise()
{
  for_contents (Score_profile, *this, ssm)
    {
      int norm = -InfinityScore;
      for_contents (Symbol_score_map, *ssm, ss) ScorePSumAcc (norm, (*ss).second);
      for_contents (Symbol_score_map, *ssm, ss) ScorePMulAcc ((*ss).second, -norm);
    }
}

Digitized_biosequence Score_profile::consensus_dsq() const
{
  Digitized_biosequence dsq (size(), -1);
  for (int pos = 0; pos < (int) size(); pos++) dsq[pos] = consensus(pos);
  return dsq;
}

Digitized_biosequence Score_profile::sample_dsq (double kT) const
{
  Digitized_biosequence dsq (size(), -1);
  for (int pos = 0; pos < (int) size(); pos++) dsq[pos] = sample (pos, kT);
  return dsq;
}

Score_profile Score_profile::evolve (const array2d<int>& substitution_matrix) const
{
  Score_profile result (size(), Symbol_score_map());
  for (int pos = 0; pos < (int) size(); ++pos)
    for (int j = 0; j < substitution_matrix.ysize(); ++j)         // j is target residue
      {
	Score j_score = -InfinityScore;
	for_const_contents (Symbol_score_map, (*this)[pos], ss)
	  ScorePSumAcc (j_score, ScorePMul ((*ss).second, substitution_matrix ((*ss).first, j)));
	if (j_score > -InfinityScore)
	  result[pos].insert (result[pos].end(), Symbol_score (j, j_score));
      }
  return result;
}

vector<int> Score_profile::null_scores_by_residue (const vector<int>& null_model) const
{
  vector<int> scores (size());
  for (int i = 0; i < (int) size(); ++i) scores[i] = score_map_score (((Score_profile&) *this) [i], null_model);
  return scores;
}

int Score_profile::inner_product (const Score_profile& prof) const
{
  if (size() != prof.size()) THROW Standard_exception ("Attempt to take the inner product of Score_profiles of different lengths");
  int score = 0;
  for (int i = 0; i < (int) size(); ++i)
    {
      int residue_score = -InfinityScore;
      for_const_contents (Symbol_score_map, (*this)[i], ss)
	{
	  Symbol_score_map::const_iterator f = prof[i].find ((*ss).first);
	  if (f != prof[i].end()) ScorePSumAcc (residue_score, (*f).second + (*ss).second);
	}
      ScorePMulAcc (score, residue_score);
    }
  return score;
}

Weight_profile Score_profile::score2weight() const
{
  Weight_profile profile;
  profile.reserve (size());
  for_const_contents (Score_profile, *this, ssm)
    {
      Symbol_weight_map swm;
      for_const_contents (Symbol_score_map, *ssm, ss)
	swm.insert (Symbol_weight ((*ss).first, Score2Prob ((*ss).second)));
      profile.push_back (swm);
    }
  return profile;
}


void Score_profile::count_symbol_frequencies (vector<double>& composition) const
{
  for_const_contents (Score_profile, *this, ssm)
    for_const_contents (Symbol_score_map, *ssm, ss)
    composition [(*ss).first] += Score2Prob ((*ss).second);
}

Alphabet::Alphabet (const char* name, int size)
{
  reset (name, size);
}

Alphabet::Alphabet (void)
{
  reset("uninitialized", 1);
}

void Alphabet::reset (const char* new_name, int new_size)
{
  name = new_name;
  alphabet_size = new_size;

  char2int_table = vector<int> (CHAR2INT_TABLE_SIZE, (int) alphabet_size);
  char_aff_table = vector<double> (CHAR2INT_TABLE_SIZE, 0.);
  nondegen_comp = vector<int>();
  degen_comp = vector<int>();

  tok.clear();
  chars_lc.clear();
  chars_uc.clear();
  degen_lc.clear();
  display_lc.clear();

  degen_weight.clear();
  degen_score.clear();
  nondegen_weight.clear();
  nondegen_score.clear();
  for (int c = 0; c < alphabet_size; ++c)
    {
      Symbol_weight_map w;
      w.insert (Symbol_weight (c, 1.));
      nondegen_weight.push_back (w);

      Symbol_score_map sc;
      sc.insert (Symbol_score (c, 0));
      nondegen_score.push_back (sc);
    }
  // alphabet returns (alphabet_size) if symbol unknown,
  // so place a degenerate symbol at the end of nondegen_score & nondegen_weight
  const double w = 1. / (double) size();
  nondegen_weight.push_back (flat_weight_map (w));
  nondegen_score.push_back (flat_score_map (Prob2Score (w)));

  wild_ssm.clear();
  wild_swm.clear();
  for (int c = 0; c < alphabet_size; ++c)
    {
      wild_swm.insert (Symbol_weight (c, 1.));
      wild_ssm.insert (Symbol_score (c, 0));
    }

  case_sensitive = false;
  unknown_int = alphabet_size;
  degenerate_lcd = alphabet_size;
  warn_randomise = 1;
  warn_unknown = 1;
}

void Alphabet::init_chars (const char* chars, const char* comp)
{
  const int len = strlen (chars);
  if (len != 0 && len != alphabet_size) THROWEXPR ("Bad chars string: '" << chars << "'");
  // figure out if alphabet is case-sensitive
  for (int c = 0; c < CHAR2INT_TABLE_SIZE; ++c)
    char2int_table[c] = alphabet_size;
  case_sensitive = false;
  for (int i = 0; i < len; ++i)
    {
      const char c = tolower (chars[i]);
      const int j = char2int_table[c];
      if (j != alphabet_size)
	{
	  if (chars[j] == chars[i])
	    THROWEXPR ("Duplicate character '" << chars[i] << "' in alphabet");
	  case_sensitive = true;
	}
      char2int_table[c] = i;

      sstring c_tok;
      c_tok << c;
      tok.push_back (c_tok);
    }
  // reset char2int_table
  if (case_sensitive)
    for (int i = 0; i < CHAR2INT_TABLE_SIZE; ++i)
      char2int_table[i] = alphabet_size;
  // initialise the alphabet
  chars_lc = chars_uc = chars;
  if (!case_sensitive)
    {
      chars_lc.to_lower();
      chars_uc.to_upper();
    }
  for (int i = 0; i < len; ++i)
    {
      char2int_table[chars_lc[i]] = char2int_table[chars_uc[i]] = i;
      char_aff_table[chars_lc[i]] = char_aff_table[chars_uc[i]] = 1.0;
    }
  nondegen_comp.clear();
  const int comp_len = strlen (comp);
  if (comp_len != 0 && comp_len != alphabet_size) THROWEXPR ("Bad complement string: '" << comp << "'");
  for (int i = 0; i < comp_len; ++i)
    nondegen_comp.push_back (char2int_strict (comp[i]));
}

void Alphabet::init_degen_chars (const char* chars)
{
  const int len = strlen (chars);
  degen_lc = chars;
  if (case_sensitive)
    for (int i = 0; i < len; ++i)
      char2int_table[degen_lc[i]] = -i-1;
  else
    {
      degen_lc.to_lower();
      for (int i = 0; i < len; ++i)
	{
	  char2int_table[degen_lc[i]] = char2int_table[toupper(degen_lc[i])] = -i-1;
	  char_aff_table[degen_lc[i]] = char_aff_table[toupper(degen_lc[i])] = 0.5;
	}
    }
  Symbol_score_map wildcard_score_map = flat_score_map (0);
  Symbol_weight_map wildcard_weight_map = flat_weight_map (1.0);
  degen_weight = vector<Symbol_weight_map> (len, wildcard_weight_map);
  degen_score = vector<Symbol_score_map> (len, wildcard_score_map);
}

void Alphabet::init_degen_complement()
{
  degen_comp.clear();
  if (has_complement())
    for (int i = 0; i < (int) degen_lc.size(); ++i)
      degen_comp.push_back (complement (-i-1));
}

void Alphabet::init_display (const sstring& chars)
{
  if ((int) chars.size() != 0 && (int) chars.size() != alphabet_size)
    THROWEXPR ("Bad display string: '" << chars << "'");
  display_lc = chars;
}

void Alphabet::set_degen_char (char degen_char, const Symbol_weight_map& p)
{
  const int i = char2int_deg (degen_char);
  if (i >= 0) THROWEXPR ("Bad degenerate char: '" << degen_char << "'");
  degen_weight[-i-1] = p;
  Symbol_score_map sc;
  for_const_contents (Symbol_weight_map, p, sw)
    sc[sw->first] = Prob2Score (sw->second);
  degen_score[-i-1] = sc;
  // check to see if this is the unknown char
  if (unknown_int >= alphabet_size && sc == flat_score_map (-Prob2Score(alphabet_size)))
    unknown_int = i;
}

void Alphabet::set_degen_char (char degen_char, const char* possibilities)
{
  Symbol_weight_map p;
  const int n_poss = strlen (possibilities);
  for (int i = 0; i < n_poss; ++i) p[char2int_strict(possibilities[i])] = 1.0 / (double) n_poss;
  set_degen_char (degen_char, p);
}

void Alphabet::set_wildcard (char wildcard_char)
{
  Symbol_weight_map p = flat_weight_map (1.0);
  set_degen_char (wildcard_char, p);
}

void Alphabet::set_unknown (char unknown_char)
{
  Symbol_weight_map p = flat_weight_map (1.0 / (double) alphabet_size);
  set_degen_char (unknown_char, p);
}

bool Alphabet::all_wild (const Score_profile& prof) const
{
  Symbol_score_map wildcard_score_map = flat_score_map (0);
  for (int i = 0; i < (int) prof.size(); ++i) if (prof[i] != wildcard_score_map) return 0;
  return 1;
}

void Alphabet::init_hidden (const Alphabet& base_alphabet, int hidden_classes)
{
  vector<sstring> hidden_class_labels;
  for (int c = 0; c < hidden_classes; ++c)
    {
      sstring s;
      s << c + 1;
      hidden_class_labels.push_back (s);
    }
  init_hidden (base_alphabet, hidden_class_labels);
}

void Alphabet::init_hidden (const Alphabet& base_alphabet, const vector<sstring>& hidden_class_labels)
{
  const int A = base_alphabet.size();
  const int C = hidden_class_labels.size();

  if (C == 0)
    THROWEXPR ("Alphabet::init_hidden called with zero hidden classes (at least one is needed). Complain loudly to your nearest dartcoder!");

  if (C == 1)
    *this = base_alphabet;
  else  // C > 1
    {
      reset (base_alphabet.name.c_str(), C*A);
      case_sensitive = base_alphabet.case_sensitive;

      int_string C_string (C);
      name << C_string;

      const sstring nd_chars (base_alphabet.nondegenerate_chars());
      const sstring d_chars (base_alphabet.degenerate_chars());

      sstring chars;
      chars << nd_chars << d_chars;
      init_degen_chars (chars.c_str());

      for (int a = 0; a < A; ++a)
	{
	  Symbol_weight_map p;
	  for (int c = 0; c < C; ++c)
	    p.insert (Symbol_weight (c*A + a, 1));
	  if ((int) nd_chars.size() == A)
	    set_degen_char (nd_chars[a], p);
	}

      for (int i = 0; i < (int) d_chars.size(); ++i)
	{
	  const Symbol_weight_map deg_w = base_alphabet.char2weight (d_chars[i]);
	  Symbol_weight_map p;
	  for_const_contents (Symbol_weight_map, deg_w, sw)
	    for (int c = 0; c < C; ++c)
	      p.insert (Symbol_weight (c*A + sw->first, sw->second));
	  set_degen_char (d_chars[i], p);
	}

      sstring base_display_chars;
      if ((int) base_alphabet.display_lc.size() == A)
	base_display_chars = base_alphabet.display_lc;
      else if ((int) nd_chars.size() == A)
	base_display_chars = nd_chars;

      sstring display_chars;
      for (int c = 0; c < C; ++c)
	display_chars << base_display_chars;

      init_display (display_chars);

      for (int c = 0; c < C; ++c)
	for (int a = 0; a < A; ++a)
	  {
	    sstring ca_tok;
	    if ((int) base_alphabet.tok.size() > a && base_alphabet.tok[a].size())
	      ca_tok << base_alphabet.tok[a];
	    else if ((int) nd_chars.size() > a)
	      ca_tok << nd_chars[a];
	    else
	      ca_tok << (a + 1) << '_';
	    ca_tok << hidden_class_labels[c];
	    tok.push_back (ca_tok);
	  }
    }
}


bool Alphabet::all_wild (const Weight_profile& prof) const
{
  Symbol_weight_map wildcard_weight_map = flat_weight_map (1.0);
  for (int i = 0; i < (int) prof.size(); ++i) if (prof[i] != wildcard_weight_map) return 0;
  return 1;
}

int Alphabet::weight2int (const Symbol_weight_map& w) const
{
  if (w.size() == 1) return w.begin()->first;
  // minimise sum-of-squares error
  int best = -1;
  double best_err = 0;
  bool first = 1;
  for (int i = 0; i < (int) (size() + degen_weight.size()); ++i)
    {
      const int sym = i < size() ? i : -(i-size())-1;
      Symbol_weight_map x = int2weight(sym);
      double err = 0;
      for (int j = 0; j < size(); ++j)
	{
	  Symbol_weight_map::const_iterator pw = w.find(j);
	  Symbol_weight_map::const_iterator px = x.find(j);
	  const double w1 = pw == w.end() ? 0.0 : pw->second;
	  const double w2 = px == x.end() ? 0.0 : px->second;
	  err += (w1 - w2) * (w1 - w2);
	}
      if (first || err < best_err*.99) { best = sym; best_err = err; first = 0; }  // .99 term favours earlier symbols (hacky)
    }
  return best;
}

char Alphabet::weight2char (const Symbol_weight_map& w) const
{
  return int2char (weight2int (w));
}

int Alphabet::score2int (const Symbol_score_map& sc) const
{
  if (sc.size() == 1) return sc.begin()->first;
  Symbol_weight_map w;
  for_const_contents (Symbol_score_map, sc, ss) w[ss->first] = Score2Prob(ss->second);
  return weight2int (w);
}

char Alphabet::score2char (const Symbol_score_map& sc) const
{
  return int2char (score2int (sc));
}

Symbol_score_map Alphabet::flat_score_map (int score) const
{
  Symbol_score_map m;
  for (int s = 0; s < size(); s++) m.insert (Symbol_score (s, score));
  return m;
}

Symbol_weight_map Alphabet::flat_weight_map (double weight) const
{
  Symbol_weight_map m;
  for (int s = 0; s < size(); s++) m.insert (Symbol_weight (s, weight));
  return m;
}

Digitized_biosequence& Alphabet::seq2dsq(const Biosequence& seq, Digitized_biosequence& dsq) const
{
  dsq.clear();
  dsq.reserve(seq.size());
  for (int i = 0; i < (int) seq.size(); ++i) dsq.push_back (char2int_strict (seq[i]));
  return dsq;
}

Biosequence& Alphabet::dsq2seq(const Digitized_biosequence& dsq, Biosequence& seq) const
{
  seq.clear();
  seq.reserve(dsq.size());
  for (int i = 0; i < (int) dsq.size(); ++i) seq.push_back (int2char (dsq[i]));
  return seq;
}


Weight_profile& Alphabet::seq2weight (const Biosequence& seq, Weight_profile& prof) const
{
  prof.clear();
  prof.reserve (seq.size());
  for_const_contents (Biosequence, seq, c) prof.push_back (char2weight(*c));
  return prof;
}

Biosequence& Alphabet::weight2seq (const Weight_profile& prof, Biosequence& seq) const
{
  seq.clear();
  seq.reserve (prof.size());
  for_const_contents (Weight_profile, prof, w) seq.push_back (weight2char(*w));
  return seq;
}

Score_profile& Alphabet::seq2score (const Biosequence& seq, Score_profile& prof) const
{
  Weight_profile wprof;
  return prof = Score_profile (seq2weight (seq, wprof));
}

Biosequence& Alphabet::score2seq (const Score_profile& prof, Biosequence& seq) const
{
  Weight_profile wprof = prof.score2weight();
  return weight2seq (wprof, seq);
}

Score_profile& Alphabet::dsq2score (const Digitized_biosequence& dsq, Score_profile& prof) const
{
  Biosequence seq;
  dsq2seq (dsq, seq);
  return seq2score (seq, prof);
}

Weight_profile& Alphabet::dsq2weight (const Digitized_biosequence& dsq, Weight_profile& prof) const
{
  Biosequence seq;
  dsq2seq (dsq, seq);
  return seq2weight (seq, prof);
}

Digitized_biosequence Alphabet::new_seq2dsq (const Biosequence& seq) const
{ Digitized_biosequence dsq; seq2dsq (seq, dsq); return dsq; }

Biosequence Alphabet::new_dsq2seq (const Digitized_biosequence& dsq) const
{ Biosequence seq; dsq2seq (dsq, seq); return seq; }

Weight_profile Alphabet::new_seq2weight (const Biosequence& seq) const
{ Weight_profile prof; seq2weight (seq, prof); return prof; }

Biosequence Alphabet::new_weight2seq (const Weight_profile& prof) const
{ Biosequence seq; weight2seq (prof, seq); return seq; }

Score_profile Alphabet::new_seq2score (const Biosequence& seq) const
{ Score_profile prof; seq2score (seq, prof); return prof; }

Biosequence Alphabet::new_score2seq (const Score_profile& prof) const
{ Biosequence seq; score2seq (prof, seq); return seq; }

Weight_profile Alphabet::new_dsq2weight (const Digitized_biosequence& dsq) const
{ Weight_profile prof; dsq2weight (dsq, prof); return prof; }

Score_profile Alphabet::new_dsq2score (const Digitized_biosequence& dsq) const
{ Score_profile prof; dsq2score (dsq, prof); return prof; }

Symbol_weight_map Alphabet::complement_weight (const Symbol_weight_map& m) const
{
  Symbol_weight_map result;
  for_const_contents (Symbol_weight_map, m, sw) result [complement ((*sw).first)] = (*sw).second;
  return result;
}

Symbol_score_map Alphabet::complement_score (const Symbol_score_map& m) const
{
  Symbol_score_map result;
  for_const_contents (Symbol_score_map, m, ss) result [nondegen_comp[ss->first]] = ss->second;
  return result;
}

int Alphabet::complement (int sym) const
{
  if (sym >= 0 && sym < alphabet_size) return nondegen_comp[sym];
  else if (sym < 0 && -sym-1 < (int) degen_comp.size()) return degen_comp[-sym-1];
  // try to guess complement
  Symbol_score_map sc = int2score (sym);
  Symbol_score_map sc_rev = complement_score (sc);
  return score2int (sc_rev);
}

char Alphabet::complement_char (char c) const
{
  char c_rev = int2char (complement (char2int_deg (c)));
  return case_sensitive ? c_rev : (isupper(c) ? toupper(c_rev) : tolower(c_rev));
}

Weight_profile Alphabet::revcomp_weight_profile (const Weight_profile& p) const
{
  Weight_profile result (p.size(), Symbol_weight_map());
  for (int i = 0; i < (int) p.size(); ++i)
    for_const_contents (Symbol_weight_map, p[p.size()-1-i], sw)
      result[i] [complement ((*sw).first)] = (*sw).second;
  return result;
}

Score_profile Alphabet::revcomp_score_profile (const Score_profile& p) const
{
  Score_profile result (p.size(), Symbol_score_map());
  for (int i = 0; i < (int) p.size(); ++i)
    for_const_contents (Symbol_score_map, p[p.size()-1-i], ss)
      result[i] [complement ((*ss).first)] = (*ss).second;
  return result;
}

Digitized_biosequence Alphabet::revcomp_digitized_sequence (const Digitized_biosequence& dsq) const
{
  Digitized_biosequence result (dsq.size());
  for (int i = 0; i < (int) dsq.size(); ++i) result[i] = complement (dsq[dsq.size()-1-i]);
  return result;
}

Biosequence Alphabet::revcomp_sequence (const Biosequence& s) const
{
  Biosequence r (s.size());
  for (int i = 0; i < (int) s.size(); ++i)
    r[i] = complement_char (s[s.size() - 1 - i]);
  return r;
}

int Alphabet::nmer_token (int word_len, const Digitized_biosequence::const_iterator iter)
{
  int nmer = 0;
  const int s = size();
  for (int mul = 1, j = word_len - 1; j >= 0; --j, mul *= s)
    nmer += mul * iter[j];
  return nmer;
}

int Alphabet::nmer_token (const char* nmer)
{
  const Biosequence seq (nmer);
  const Digitized_biosequence dsq = new_seq2dsq (seq);
  return nmer_token (dsq.size(), dsq.begin());
}

sstring Alphabet::degenerate_seq_pattern (const Biosequence& seq) const
{
  sstring pattern;
  for_const_contents (Biosequence, seq, c)
    {
      const int sym = char2int_deg(*c);
      if (sym >= 0) pattern << *c;
      else if (sym == unknown_int) pattern << '.';
      else
	{
	  pattern << '[';
	  for_const_contents (Symbol_weight_map, degen_weight[-sym-1], sw)
	    pattern << int2char (sw->first);
	  pattern << ']';
	}
    }
  return pattern;
}

Simple_alphabet::Simple_alphabet (const char* name, const char* chars)
  : Alphabet (name, strlen (chars))
{
  init_chars (chars);
}

DNA_alphabet_type::DNA_alphabet_type() : Alphabet ("DNA", 4)
{
  init_chars ("acgt", "tgca");
  init_degen_chars ("nx*urymkswhbvd");
  set_unknown ('n');
  set_unknown ('x');
  set_wildcard ('*');
  set_degen_char ('u', "t");
  set_degen_char ('r', "ag");
  set_degen_char ('y', "ct");
  set_degen_char ('m', "ac");
  set_degen_char ('k', "gt");
  set_degen_char ('s', "cg");
  set_degen_char ('w', "at");
  set_degen_char ('h', "act");
  set_degen_char ('b', "cgt");
  set_degen_char ('v', "acg");
  set_degen_char ('d', "agt");
  degenerate_lcd = 12;
  init_degen_complement();
  char_aff_table['n'] = char_aff_table['N'] = 1.;
}


RNA_alphabet_type::RNA_alphabet_type() : Alphabet ("RNA", 4)
{
  init_chars ("acgu", "ugca");
  init_degen_chars ("nx*trymkswhbvd");
  set_unknown ('n');
  set_unknown ('x');
  set_wildcard ('*');
  set_degen_char ('t', "u");
  set_degen_char ('r', "ag");
  set_degen_char ('y', "ct");
  set_degen_char ('m', "ac");
  set_degen_char ('k', "gt");
  set_degen_char ('s', "cg");
  set_degen_char ('w', "at");
  set_degen_char ('h', "act");
  set_degen_char ('b', "cgt");
  set_degen_char ('v', "acg");
  set_degen_char ('d', "agt");
  degenerate_lcd = 12;
  init_degen_complement();
  char_aff_table['n'] = char_aff_table['N'] = 1.;
}

Protein_alphabet_type::Protein_alphabet_type() : Alphabet ("Protein", 20)
{
  // Protein alphabets are arranged in increasing alphabetic order of the initial letter of the amino acid.
  // The old, IUPAC-alphabetical order is "acdefghiklmnpqrstvwy".
  init_chars ("arndcqeghilkmfpstwyv", "");
  init_degen_chars ("x*bz");
  set_unknown ('x');
  set_wildcard ('*');
  set_degen_char ('b', "nd");
  set_degen_char ('z', "qe");
}

Roman_alphabet_type::Roman_alphabet_type() : Alphabet ("RomanText", 27)
{
  init_chars ("abcdefghijklmnopqrstuvwxyz_", "");
  init_degen_chars ("~*");
  set_unknown ('~');
  set_wildcard ('*');
}

Text_alphabet_type::Text_alphabet_type() : Alphabet ("PunctuatedText", 52)
{
  init_chars ("abcdefghijklmnopqrstuvwxyz,.:;+-/'!?_$&()@0123456789", "");
  init_degen_chars ("~*\"");
  set_unknown ('~');
  set_wildcard ('*');
  set_degen_char ('\"', "'");
}

DNA_alphabet_type DNA_alphabet;
RNA_alphabet_type RNA_alphabet;
Protein_alphabet_type Protein_alphabet;
Roman_alphabet_type Roman_alphabet;
Text_alphabet_type Text_alphabet;
Simple_alphabet Dummy_alphabet ("Dummy", "x");

void Named_profile::write_FASTA (ostream& out,
				 int column_width) const
{
  assert_profiles_consistent (SEQ);  // check that everything is consistent with the ASCII sequence
  out << ">" << name.c_str() << " " << cruft.c_str() << "\n";
  int s = seq.size();
  for (int i = 0; i < s; i++)
    {
      if (i > 0 && (i % column_width) == 0) out << "\n";
      out << seq[i];
    }
  if (s > 0) out << "\n";
}


void Named_profile::clear()
{
  name.clear();
  cruft.clear();
  seq.clear();
  prof_sc.clear();
  prof_w.clear();
}

Regexp Named_profile::name_regexp ("^[ \t]*>([^ \t]*)[ \t]*(.*)");
Regexp Named_profile::quick_name_regexp ("^[ \t]*>");
Regexp Named_profile::nse_regexp ("^(.*)/([0-9])+-([0-9])+$");

void Named_profile::read_FASTA (istream& in)
{
  clear();

  sstring s;
  bool got_name;
  do {
    s.getline(in);
    s.chomp();
    got_name = name_regexp.Match (s.c_str());
  } while (!got_name && in && !in.eof());

  if (got_name)
    {
      if (name_regexp.SubStrings() >= 1) name = name_regexp[1];
      if (name_regexp.SubStrings() >= 2) cruft = name_regexp[2];

      while (1)
	{
	  streampos line_start = in.tellg();
	  s.getline(in);
	  if (quick_name_regexp.Match (s.c_str())) { in.seekg (line_start); break; }
	  for_const_contents (sstring, s, i) if (!isspace (*i) && *i != '\0') seq.push_back (*i);
	  if (!in || in.eof()) break;
	}
      CTAG(3,BIOSEQ_IO BIOSEQ) << "Read sequence '" << name << "' (" << seq.size() << " residues)\n";
    }
}

bool Named_profile::detect_FASTA (istream& in)
{
  streampos in_pos = in.tellg();
  istream_get_iterator i (in);
  istream_get_iterator end;
  while (isspace(*i) && i != end)
    i++;
  bool is_FASTA = i != end && *i == '>';
  in.seekg(in_pos);
  CLOG(7) << "Autodetect: " << (is_FASTA ? "looks like FASTA format to me (but what do I know)\n" : "doesn't look like FASTA format to me (but I'm no expert system)\n");
  return is_FASTA;
}

void Named_profile::write_metaprob (int meta_idx, ostream& out, int column_width) const
{
  if (meta_idx >= (int) meta_sc.size()) THROWEXPR ("Metascore index out of range");
  if ((int) meta_sc[meta_idx].size() != size()) THROWEXPR ("Metascore size mismatch");
  out << '>' << name << ' ' << cruft << '\n';
  for (int pos = 0; pos < size(); pos += column_width)
    {
      for (int i = pos; i < min (pos + column_width, (int) size()); ++i)
	{
	  int percent = (int) (100.0 * Score2Prob (meta_sc[meta_idx][i]) + .5);
	  if (percent > 100) percent = 100;
	  unsigned char c = percent + '@';
	  if (percent >= 50) c += 50;
	  out << c;
	}
      out << '\n';
    }
}

void Named_profile::read_metaprob (int meta_idx, istream& in, bool expect_name)
{
  if (meta_idx >= (int) meta_sc.size()) THROWEXPR ("Metascore index out of range");
  meta_sc[meta_idx].clear();


  sstring s;
  if (expect_name)
    {
      bool got_name = 0;
      while (!got_name && in && !in.eof())
	{
	  s.getline(in);
	  s.chomp();
	  got_name = name_regexp.Match (s.c_str());
	}
      if (!got_name) THROWEXPR ("While reading metascore: missing name");
      if (name_regexp[1] != name) THROWEXPR ("While reading metascore: name mismatch");
    }
  
  while (1)
    {
      streampos line_start = in.tellg();
      s.getline(in);
      if (quick_name_regexp.Match (s.c_str())) { in.seekg (line_start); break; }
      for_const_contents (sstring, s, i)
	if (!isspace (*i) && *i != '\0')
	  {
	    int percent = (int) (((unsigned char) *i) - '@');
	    if (percent >= 50) percent -= 50;
	    if (percent < 0 || percent > 100) THROWEXPR ("While reading metascore: bad character '" << *i << "'");
	    meta_sc[meta_idx].push_back (Prob2Score (100.0 / (double) percent));
	  }
    }
  
  if ((int) meta_sc[meta_idx].size() != size())
    THROWEXPR ("While reading metascore: metascore length doesn't match sequence length");
  
  CTAG(3,BIOSEQ_IO BIOSEQ) << "Read metascore '" << name << "'\n";
}

void Named_profile::seq_update (const Alphabet& a, Profile_flags convert)
{
  if (convert & SCORE) seq2score(a);
  if (convert & WEIGHT) seq2weight(a);
  if (convert & DSQ) seq2dsq(a);
}

Named_profile Named_profile::revcomp (const Alphabet& a, const char* name_suffix, const char* cruft_suffix) const
{
  Named_profile np;

  np.name = name;
  np.cruft = cruft;
  np.name << name_suffix;
  np.cruft << cruft_suffix;

  np.seq     = a.revcomp_sequence (seq);
  np.prof_sc = a.revcomp_score_profile (prof_sc);
  np.prof_w  = a.revcomp_weight_profile (prof_w);
  np.dsq     = a.revcomp_digitized_sequence (dsq);

  return np;
}


void Sequence_database::write_FASTA (ostream& out, int column_width) const
{
  for_iterator (const_iterator, i, begin(), end()) (*i).write_FASTA (out, column_width);
}

void Sequence_database::to_lower()
{
  for_iterator (iterator, i, begin(), end()) (*i).seq.to_lower();
}

void Sequence_database::to_upper()
{
  for_iterator (iterator, i, begin(), end()) (*i).seq.to_upper();
}

void Sequence_database::read_FASTA (istream& in)
{
  int count = 0;
  Named_profile profile;
  // currently, the following eof() test misses empty sequences at the end of a file. THIS IS A BUG.
  while (!in.eof())
    {
      profile.read_FASTA (in);
      if (profile.name.size() > 0 || !in.eof())
	{
	  push_back (profile);
	  ++count;
	  if (count % 500 == 0) CTAG(4,BIOSEQ_IO BIOSEQ) << "Read " << count << " sequences\n";
	}
    }
  if (count % 500 != 0) CTAG(4,BIOSEQ_IO BIOSEQ) << "Read " << count << " sequences\n";
}

void Sequence_database::write_metaprob (int meta_idx, ostream& out, int column_width) const
{
  for_iterator (const_iterator, np, begin(), end()) np->write_metaprob (meta_idx, out, column_width);
}

void Sequence_database::read_metaprob (int meta_idx, istream& in)
{
  for_iterator (iterator, np, begin(), end()) np->read_metaprob (meta_idx, in, 1);
}

int Sequence_database::total_residues() const
{
  int sz = 0;
  for_const_contents (list<Named_profile>, *this, prof) sz += prof->size();
  return sz;
}

void Sequence_database::seqs2scores (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting ASCII sequences to score profiles\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).seq2score (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::seqs2weights (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting ASCII sequences to probability profiles\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).seq2weight (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::seqs2dsqs (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting ASCII sequences to integer arrays\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).seq2dsq (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::seqs_update (const Alphabet& a, Profile_flags convert)
{
  if (convert & SCORE) seqs2scores(a);
  if (convert & WEIGHT) seqs2weights(a);
  if (convert & DSQ) seqs2dsqs(a);
}

void Sequence_database::scores2seqs (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting score profiles to ASCII sequences\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).score2seq (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::scores2weights()
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting score profiles to probability profiles\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).score2weight();
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::weights2seqs (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting probability profiles to ASCII sequences\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).weight2seq (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::weights2scores()
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting probability profiles to score profiles\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).weight2score();
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::consensus_scores2seqs (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting score profiles to consensus ASCII sequences\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).consensus_score2seq (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Sequence_database::consensus_weights2seqs (const Alphabet& a)
{
  CTAG(6,BIOSEQ) << "Sequence_database: converting probability profiles to consensus ASCII sequences\n";
  int count = 0;
  for_contents (list<Named_profile>, *this, prof)
    {
      (*prof).consensus_weight2seq (a);
      ++count;
      if (count % 500 == 0) CLOG(4) << "Processed " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Processed " << count << " sequences\n";
}

void Named_profile::reset_metascores (int n_metascores, Score default_sc)
{ meta_sc = vector<Metascore> (n_metascores, Metascore (size(), (Score) default_sc)); }

void Named_profile::add_metascores (int metascore_idx, const Metascore& delta)
{ for (int pos = 0; pos < (int) delta.size(); ++pos) ScorePMulAcc (meta_sc[metascore_idx][pos], delta[pos]); }

int Named_profile::size() const
{
  int seq_sz = seq.size(), dsq_sz = dsq.size(), sc_sz = prof_sc.size(), w_sz = prof_w.size();
  int sz = max (max (max (seq_sz, sc_sz), w_sz), dsq_sz);
  // do a simple consistency test between non-empty profiles
  if ((seq_sz != 0 && seq_sz != sz) || (sc_sz != 0 && sc_sz != sz) || (w_sz != 0 && w_sz != sz) || (dsq_sz != 0 && dsq_sz != sz))
    THROWEXPR ("Named_profile looks internally inconsistent (lengths: seq " << seq_sz << ", dsq " << dsq_sz << ", score " << sc_sz << ", weight " << w_sz << ").\nYou probably need to call a method like Named_profile::seq_update() or Sequence_database::seqs_update() after updating your sequences.");
  return sz;
}

Named_profile Named_profile::subseq (int start, int len) const
{
  if (start < 0 || len <= 0 || start + len > size())
    THROWEXPR ("Tried to take subseq " << start << "+" << len << " of sequence length " << size());
  Named_profile result;
  result.name << name << "/" << start+1 << "-" << start+len;
  if (seq.size()) result.seq = sstring (seq.begin() + start, seq.begin() + start + len);
  if (dsq.size()) result.dsq.insert (result.dsq.end(), dsq.begin() + start, dsq.begin() + start + len);
  if (prof_w.size()) result.prof_w.insert (result.prof_w.end(), prof_w.begin() + start, prof_w.begin() + start + len);
  if (prof_sc.size()) result.prof_sc.insert (result.prof_sc.end(), prof_sc.begin() + start, prof_sc.begin() + start + len);
  result.meta_sc = vector<Metascore> (meta_sc.size(), Metascore());
  for (int i = 0; i < (int) meta_sc.size(); ++i)
    if ((int) meta_sc[i].size() >= start + len)
      result.meta_sc[i].insert (result.meta_sc[i].end(), meta_sc[i].begin() + start, meta_sc[i].begin() + start + len);
    else
      THROWEXPR ("Dodgy Metascore");
  return result;
}

bool Named_profile::name_is_NSE() const
{
  return nse_regexp.Match (name.c_str());
}

NSE Named_profile::parse_name() const
{
  if (nse_regexp.Match (name.c_str()))
    return NSE (nse_regexp[1], nse_regexp[2].to_int(), nse_regexp[3].to_int());
  return NSE();
}

void Named_profile::assert_profiles_consistent (Profile_flags test) const
{
  // check for profile synchronisation by looking at sizes (not perfect but close enough)
  int seq_sz = seq.size(), dsq_sz = dsq.size(), sc_sz = prof_sc.size(), w_sz = prof_w.size();
  int sz;
  if (test & SEQ) sz = seq_sz;
  else if (test & SCORE) sz = sc_sz;
  else if (test & WEIGHT) sz = w_sz;
  else if (test & DSQ) sz = dsq_sz;
  else THROW Standard_exception ("assert_profiles_consistent() called with no profile flags set");
  bool ok = 1;
  if (((test & SEQ) || seq_sz != 0) && sz != seq_sz) ok = 0;
  if (((test & SCORE) || sc_sz != 0) && sz != sc_sz) ok = 0;
  if (((test & WEIGHT) || w_sz != 0) && sz != w_sz) ok = 0;
  if (((test & DSQ) || dsq_sz != 0) && sz != dsq_sz) ok = 0;
  for_const_contents (vector<Metascore>, meta_sc, msc)
    if (seq_sz != (int) (*msc).size()) ok = 0;
  if (!ok)
    THROWEXPR ("Named_profile looks internally inconsistent (lengths: seq " << seq_sz << ", dsq " << dsq_sz << ", score " << sc_sz << ", weight " << w_sz << ").\nYou probably need to call a method like Named_profile::seq_update() or Sequence_database::seqs_update() after reading in or updating your sequences.");
}

Sequence_database Sequence_database::revcomp (const Alphabet& a, const char* name_suffix, const char* cruft_suffix) const
{
  CTAG(6,BIOSEQ) << "Reverse complementing sequence database\n";
  int count = 0;
  Sequence_database rev;
  for_const_contents (list<Named_profile>, *this, prof)
    {
      rev.push_back ((*prof).revcomp (a, name_suffix, cruft_suffix));
      ++count;
      if (count % 500 == 0) CLOG(4) << "Reverse complemented " << count << " sequences\n";
    }
  if (count % 500 != 0) CLOG(4) << "Reverse complemented " << count << " sequences\n";
  return rev;
}

bool Sequence_database::matches_alphabet (const Alphabet& alphabet, double threshold) const
{
  int count = 0;
  int total = 0;
  for_const_contents (Sequence_database, *this, np)
    {
      total += (*np).seq.size();
      for_const_contents (Biosequence, (*np).seq, c)
	if (alphabet.contains_deg (*c)) ++count;
    }
  if (total == 0) { CLOGERR << "Warning: Sequence database is empty & therefore trivially consistent with " << alphabet.name << " alphabet\n"; return 1; }
  double frac = ((double) count) / ((double) total);
  bool is_alphabet = frac >= threshold;
  int perc = (int) (100.0*frac);
  if (is_alphabet)
    CLOG(7) << perc << "% of residues matched " << alphabet.name << " alphabet; looks like " << alphabet.name << " to me\n";
  else
    CLOG(7) << perc << "% of residues matched " << alphabet.name << " alphabet; doesn't look like " << alphabet.name << " to me\n";
  return is_alphabet;
}

double Sequence_database::alphabet_similarity (const Alphabet& alphabet) const
{
  double score = 0;
  for_const_contents (Sequence_database, *this, np)
    for_const_contents (Biosequence, (*np).seq, c)
    score += alphabet.char_affinity (*c);
  return score;
}

const Alphabet& Sequence_database::detect_alphabet() const
{
  vector<Alphabet*> alph;
  alph.push_back (&DNA_alphabet);
  alph.push_back (&RNA_alphabet);
  alph.push_back (&Protein_alphabet);
  alph.push_back (&Text_alphabet);

  vector<double> alph_weight (alph.size(), 1.0);
  alph_weight[2] = .99;    // give Protein_alphabet 99% weight
  alph_weight[3] = .6;  // give Text_alphabet 60% weight
  
  CTAG(2,BIOSEQ) << "Trying to figure out the best alphabet...\n";
  int best = -1;
  double best_sim = 0;
  for (int i = 0; i < (int) alph.size(); ++i)
    {
      const double sim = alph_weight[i] * alphabet_similarity (*alph[i]);
      CTAG(2,BIOSEQ) << "Score for " << alph[i]->name << " alphabet is " << sim << "\n";
      if (i == 0 || sim > best_sim) { best = i; best_sim = sim; }
    }
  const Alphabet& best_alph = *alph[best];
  CTAG(3,BIOSEQ) << "Best alphabet is " << best_alph.name << "\n";
  return best_alph;
}

void Sequence_database::count_symbol_frequencies (vector<double>& composition) const
{
  CTAG(6,BIOSEQ) << "Counting symbol frequencies.\n";
  for_const_contents (Sequence_database, *this, np)
    {
      (*np).assert_profiles_consistent (DSQ);    // this call protects against the caller not having called seqs2dsqs() or whatever
      (*np).dsq.count_symbol_frequencies (composition);
    }
}

vector<Prob> Sequence_database::get_null_model (int alphabet_size, double pseudocount) const
{
  vector<Prob> null_model (alphabet_size, pseudocount);
  count_symbol_frequencies (null_model);
  NormalisePr (null_model);
  if (CLOGGING(5)) CL << "Residue frequencies: (" << null_model << ")\n";
  return null_model;
}

Prob Sequence_database::get_null_extend (double pseudocount) const
{
  const double res = pseudocount + (double) total_residues();
  const double seqs = pseudocount + (double) size();
  return res / (res + seqs);
}

double Sequence_database::mean_length() const
{
  int total_len = 0;
  for_const_contents (Sequence_database, *this, np)
    total_len += (*np).size();
  return ((double) total_len) / ((double) size());
}

Sequence_database::iterator Sequence_database::append_sequence_database (const Sequence_database& seq_db)
{
  insert (end(), seq_db.begin(), seq_db.end());
  // this is a hacky way of finding the first new element, but it works (ug)
  Sequence_database::iterator iter = end();
  for_const_contents (Sequence_database, seq_db, dummy_iter)
    --iter;
  return iter;
}

FASTA_sequence_database::FASTA_sequence_database (const char* filename,
						  const Alphabet* alphabet,
						  Profile_flags convert)
  : my_alphabet (alphabet),
    index()
{
  if (filename)
    read_from_file (filename, alphabet, convert);
}

void FASTA_sequence_database::read_FASTA (istream& in)
{
  ((Sequence_database*) this)->read_FASTA (in);
  index.update (*this);
}

void FASTA_sequence_database::read_from_file (const char* filename,
					      const Alphabet* alphabet,
					      Profile_flags convert)
{
  CLOG(7) << "FASTA_sequence_database: reading from file '" << filename << "'\n";
  ifstream seq_stream (filename);
  if (!seq_stream) { sstring e; e << "Sequence file '" << filename << "' not found"; THROW Standard_exception (e); }
  read_FASTA (seq_stream);
  my_alphabet = alphabet==0 ? &detect_alphabet() : alphabet;
  seqs_update (*my_alphabet, convert);
  index.update (*this);
}

void FASTA_sequence_database::add_seq (const char* name, const char* seq, Profile_flags convert)
{
  Named_profile np;
  np.name = sstring (name);
  np.seq = sstring (seq);
  push_back (np);
  index.update (*this);
  if (my_alphabet == 0) my_alphabet = &detect_alphabet();
  back().seq_update (*my_alphabet, convert);
}

FASTA_sequence_database::iterator FASTA_sequence_database::append_sequence_database (const Sequence_database& seq_db)
{
  FASTA_sequence_database::iterator iter = Sequence_database::append_sequence_database (seq_db);
  update_index();
  return iter;
}

void FASTA_sequence_database::update_alphabet (const Alphabet* alphabet)
{
  my_alphabet = alphabet ? alphabet : &detect_alphabet();
}

void Sequence_database_index::update (Sequence_database& seq_db)
{
  CTAG(6,BIOSEQ) << "Sequence_database_index: building name index\n";

  name.clear();
  profile.clear();
  name2profile_index.clear();

  profile.reserve (seq_db.size());

  for_contents (Sequence_database, seq_db, np) {
    name2profile_index [(*np).name] = profile.size();
    name.push_back ((*np).name);
    profile.push_back (&*np);
  }
}

Named_profile* Sequence_database_index::name2profile (const sstring& name) const
{
  int i = name2profile_index.lookup (name);
  if (i < 0) { sstring e; e << "Sequence " << name << " not found"; THROW Standard_exception (e); }
  return profile[i];
}

void Sequence_database_index::assert_names_unique (bool test_NSE_conflict) const
{
  // barf if any names are duplicated
  // also check that no names are of form "X/Y-Z" where X is another sequence name
  for (int i = 0; i < size(); ++i)
    {
      const Named_profile& np (*profile[i]);
      const sstring& name = np.name;
      if (name2profile_index.lookup (name) != i)
	THROWEXPR ("Duplicate sequence name '" << name << "' in sequence database");
      if (test_NSE_conflict && np.name_is_NSE())
	{
	  const NSE nse = np.parse_name();
	  const int j = name2profile_index.lookup (nse.seqname);
	  if (j >= 0)
	    THROWEXPR ("Sequence names '" << name << "' and '" << nse.seqname << "' could potentially conflict");
	}
    }
}
