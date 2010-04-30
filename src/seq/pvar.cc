#include "seq/pvar.h"
#include "util/dexception.h"

void PVar::show (ostream& o, const vector<vector<sstring> >* group_suffix, bool no_group_prefix, bool quote_special_chars) const
{
  sstring text;
  if (no_group_prefix)
    text << (*group_suffix)[group_idx][var_idx];
  else
    {
      if (group_idx >= 26)
	text << (char) ('a' + group_idx / 26);
      text << (char) ('a' + (group_idx % 26));
      if (group_suffix && (int) group_suffix->size() > group_idx && (int) (*group_suffix)[group_idx].size() > var_idx)
	text << '[' << (*group_suffix)[group_idx][var_idx] << ']';
      else
	text << '[' << var_idx << ']';
    }
  if (quote_special_chars)
    {
      const SExpr_atom quoted_text (text);
      o << quoted_text;
    }
  else
    o << text;
}

Kronecker_score boolean_delta_no  (2, Boolean_group::No_index);
Kronecker_score boolean_delta_yes (2, Boolean_group::Yes_index);

void PScores::normalise()
{
  for (int g = 0; g < groups(); ++g)
    if (group_size (g) > 1)
      {
	// too lazy to write a NormaliseFSc right now - 4/29/2010
	vector<Score> sc (group_size(g));
	for (int v = 0; v < group_size(g); ++v)
	  sc[v] = FScore2Score (group[g][v]);
	NormaliseSc (sc);
	for (int v = 0; v < group_size(g); ++v)
	  group[g][v] = Score2FScore (sc[v]);
      }
}

void PScores::set_null_model (const Sequence_database& seq_db, const Alphabet_group& null_emit, const Boolean_group& null_extend)
{
  const int A = null_emit.group_size;
  const vector<Prob> null_emit_prob = seq_db.get_null_model (A, 0.);
  const Prob null_extend_prob = seq_db.get_null_extend (0.);
  for (int a = 0; a < A; ++a)
    (*this)[null_emit[a]] = Prob2Score (null_emit_prob[a]);
  (*this)[null_extend.n] = Prob2Score (1. - null_extend_prob);
  (*this)[null_extend.y] = Prob2Score (null_extend_prob);
}

void PScores::write (ostream& out) const
{
  out << groups() << "\n";
  for (int g = 0; g < groups(); ++g)
    {
      out << group_size(g);
      for (int p = 0; p < group_size(g); ++p)
	out << " " << group[g][p];
      out << "\n";
    }
}

void PScores::read (istream& in)
{
  int i;
  in >> i;
  if (i != groups()) THROWEXPR ("PScores: read " << i << " groups, expected " << groups());
  for (int g = 0; g < groups(); ++g)
    {
      in >> i;
      if (i != group_size(g)) THROWEXPR ("PScores: read size " << i << " for group " << g << ", expected " << group_size(g));
      for (int p = 0; p < group_size(g); ++p)
	in >> group[g][p];
    }
}

PCounts::PCounts() { }

PCounts::PCounts (const PScores& var_scores)
{
  restructure (var_scores);
}

void PCounts::restructure (const PScores& var_scores)
{
  group.clear();
  group_name.clear();
  group_suffix.clear();
  wait.clear();

  for (int g = 0; g < var_scores.groups(); ++g)
    new_group (var_scores.group[g].size(), var_scores.group_name[g].c_str(), var_scores.group_suffix[g]);
}

void PCounts::clear()
{
  for (int g = 0; g < groups(); ++g)
    {
      if (group_size(g) == 1)
	wait[g] = 0.;
      for (int p = 0; p < group_size(g); ++p)
	group[g][p] = 0.0;
    }
}

void PCounts::show (ostream&o) const
{
  PVar_container<Prob>::show(o,"counts");

  bool first_wait = true;
  for (int g = 0; g < groups(); ++g)
    if (group[g].size() == 1)
      {
	if (first_wait)
	  o << "Wait times:\n";
	first_wait = false;

	PVar id (g, 0);
	id.show (o, &group_suffix);
	o << ' ' << wait[g] << "\n";
      }
}

void PCounts::optimise_group (int group_idx, PScores& var_scores) const
{
  assert_same_dimensions (var_scores);
  const int sz = group_size (group_idx);
  if (sz > 1)
    {
      double norm = 0.;
      for (int p = 0; p < sz; ++p)
	norm += group[group_idx][p];
      if (norm > 0.)
	for (int p = 0; p < sz; ++p)
	  var_scores.group[group_idx][p] = Prob2Score (group[group_idx][p] / norm);
      else
	CLOGERR << "Warning: insufficient data to optimize probability variable group #" << group_idx << "\n";
    }
  else
    if (wait[group_idx] > 0.)
      var_scores.group[group_idx][0] = Prob2Score (group[group_idx][0] / wait[group_idx]);
    else
      CLOGERR << "Warning: insufficient data to optimize rate variable #" << group_idx << "\n";
}

void PCounts_like::clear()
{
  log_likelihood = 0.0;
  PCounts::clear();
}

int PScores::pgroups() const
{
  return groups();
}

int PScores::pgroup_size (int pgroup_idx) const
{
  return group_size (pgroup_idx);
}

PGroup PScores::new_group (unsigned int size, const char* name, const char* suffix)
{
  return PGroup (_new_group_index (size, -InfinityScore, name, suffix), size);
}

PGroup PScores::new_group (unsigned int size, const char* name, const vector<sstring>& suffix)
{
  return PGroup (_new_group_index (size, -InfinityScore, name, suffix), size);
}

PGroup PScores::new_group (const char* name, const char* suffix)
{
  const int g = _new_group_index (-InfinityScore, name, suffix);
  return PGroup (g, group_size(g));
}

Boolean_group PScores::new_boolean_group (const char* name)
{
  sstring name_str;
  if (name) name_str << name;
  return Boolean_group (_new_group_index (2, -InfinityScore, name_str.c_str(), "NY"));
}

Alphabet_group PScores::new_alphabet_group (const Alphabet& a, const char* name, bool bigend)
{
  sstring name_str;
  if (name) name_str << name;
  return Alphabet_group (_new_group_index (a.size(), -InfinityScore, name_str.c_str(), a.nondegenerate_chars_toupper()), a, 1, bigend);
}

Alphabet_group PScores::new_alphabet_group (const Alphabet& a, int word_len, const char* name, bool bigend)
{
  Alphabet_group pg = new_alphabet_group (a.size(), word_len, name, bigend);
  pg._alphabet = &a;
  vector<sstring> suffix (pg.group_size);
  for (int w = 0; w < pg.group_size; ++w)
    suffix[w] = pg.index2word (w);
  group_suffix[pg.group_idx] = suffix;
  return pg;
}

Alphabet_group PScores::new_alphabet_group (const vector<sstring>& alphabet_tokens, int word_len, const char* name, bool bigend)
{
  Alphabet_group pg = new_alphabet_group (alphabet_tokens.size(), word_len, name, bigend);
  vector<sstring> suffix (pg.group_size);
  for (int w = 0; w < pg.group_size; ++w)
    {
      const vector<int> intvec = pg.index2intvec (w);
      for (int i = 0; i < (int) intvec.size(); ++i)
	suffix[w] << alphabet_tokens[intvec[i]] << ' ';
      suffix[w].chomp (' ');
    }
  group_suffix[pg.group_idx] = suffix;
  return pg;
}

Alphabet_group PScores::new_alphabet_group (int alphabet_size, int word_len, const char* name, bool bigend)
{
  sstring name_str;
  if (name) name_str << name;
  const int n_words = (int) pow ((double) alphabet_size, (double) word_len);
  vector<sstring> suffix (n_words);
  const int g = _new_group_index (n_words, -InfinityScore, name_str.c_str(), suffix);
  const Alphabet_group pg = Alphabet_group (g, alphabet_size, word_len, bigend);
  for (int w = 0; w < n_words; ++w)
    {
      const vector<int> intvec = pg.index2intvec (w);
      for (int i = 0; i < (int) intvec.size(); ++i)
	suffix[w] << intvec[i] << ' ';
      suffix[w].chomp (' ');
    }
  group_suffix[g] = suffix;
  return pg;
}
