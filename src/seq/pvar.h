// this module does EM more nicely than in the singlehmm & pairhmm modules
// allowing for redundant dependencies of HMM parameters on a restricted parameter set

#ifndef PVAR_INCLUDED
#define PVAR_INCLUDED

#include "util/score.h"
#include "seq/biosequence.h"
#include "util/strsaver.h"
#include "util/sexpr.h"

// A PVar struct is a tag identifying a parameter by group & index-within-group
struct PVar
{
  int group_idx;
  int var_idx;

  PVar () : group_idx(-1), var_idx(-1) { }
  PVar (unsigned int g, unsigned int p) : group_idx(g), var_idx(p) { }

  bool operator== (const PVar& pv) const
  { return group_idx == pv.group_idx && var_idx == pv.var_idx; }

  bool operator< (const PVar& pv) const
  { return group_idx < pv.group_idx ? 1 : (group_idx == pv.group_idx ? (var_idx < pv.var_idx) : 0); }

  bool is_null() const { return group_idx < 0; }

  // show method maps:
  //   group_idx --> 'a', 'b', 'c' ... 'aa', 'ab' etc
  //   var_idx --> '[0]', '[1]' etc
  //          or --> '[A]', '[C]' etc if (e.g.) (*group_suffix)[group_idx] == { "A", "C", "G", "T" }

  void show (ostream& o, const vector<vector<sstring> >* group_suffix = 0, bool no_group_prefix = false, bool quote_special_chars = true) const;
};

// PGroup has a group index and a size parameter, and can generate PVar's using the [] operator
// the size parameter is basically decoration, providing range checking for the [] operator
struct PGroup
{
  int group_idx;
  int group_size;

  PGroup () : group_idx(-1), group_size(-1) { }
  PGroup (int g_idx) : group_idx(g_idx), group_size(-1) { }
  PGroup (int g_idx, int g_sz) : group_idx(g_idx), group_size(g_sz) { }

  PVar operator[] (int var_idx) const
  {
    if (var_idx < 0 || (group_size >= 0 && var_idx >= group_size))
      THROW Standard_exception ("PVar index out of range");
    return PVar (group_idx, var_idx);
  }

  bool operator== (const PGroup& pg) const
  { return group_idx == pg.group_idx && group_size == pg.group_size; }

  bool operator< (const PGroup& pg) const
  { return group_idx < pg.group_idx ? 1 : (group_idx == pg.group_idx ? (group_size < pg.group_size) : 0); }

  bool is_null() const { return group_idx < 0; }
};

// Boolean_group is a subclass of PGroup for groups of size 2

struct Boolean_group : PGroup
{
  enum { No_index = 0, Yes_index = 1 };
  PVar NO, n;
  PVar YES, y;
  Boolean_group() : PGroup(), NO(), n(), YES(), y() { }
  Boolean_group (int g_idx) : PGroup(g_idx,2),
			      NO(g_idx,No_index), n(g_idx,No_index),
			      YES(g_idx,Yes_index), y(g_idx,Yes_index) { }
  PVar operator[] (bool var_idx) const { return PGroup::operator[] (var_idx ? 1 : 0); }
};

// Kronecker delta Score vectors for Boolean_group's
extern Kronecker_score boolean_delta_no;
extern Kronecker_score boolean_delta_yes;

// Alphabet_group is a subclass of PGroup for Alphabet-ic groups

struct Alphabet_group : PGroup
{
  // data
  const Alphabet* _alphabet;
  int word_len;
  int alphabet_size;
  bool big_endian;  // if true, then Hubertus' leftmost digit ascends first

  // constructors
  Alphabet_group (bool bigend = true)
    : PGroup(),
      _alphabet (0),
      word_len(0),
      alphabet_size(0),
      big_endian (bigend)
  { }

  Alphabet_group (int g_idx, const Alphabet& a, int _word_len = 1, bool bigend = true)
    : PGroup (g_idx, (int) pow ((double) a.size(), (double) _word_len)),
      _alphabet (&a),
      word_len (_word_len),
      alphabet_size (a.size()),
      big_endian (bigend)
  { }

  Alphabet_group (int g_idx, int alph_sz, int _word_len = 1, bool bigend = true)
    : PGroup (g_idx, (int) pow ((double) alph_sz, (double) _word_len)),
      _alphabet (0),
      word_len (_word_len),
      alphabet_size (alph_sz),
      big_endian (bigend)
  { }

  // alphabet accessors
  bool has_alphabet() const { return _alphabet != 0; }
  const Alphabet& alphabet() const { return *_alphabet; }

  // operator[int] accessor, so that operator[char] doesn't override superclass method (why the f**k does this happen? goddamn C++)
  PVar operator[] (int var_idx) const
  {
    return PGroup::operator[] (var_idx);
  }

  // operator[char] accessor for words of length 1
  PVar operator[] (char var_char) const
  {
    if (var_char >= 0 && var_char < group_size) return PGroup::operator[] (var_char);  // doubles up as operator[](int)
    if (word_len != 1) THROWEXPR ("Tried to use a [char] accessor on a multi-char Alphabet PGroup");
    if (!_alphabet->contains_strict (var_char))
      THROWEXPR ("Tried to get unknown PVar '" << var_char << "' from Alphabet PGroup");
    return PGroup::operator[] (_alphabet->char2int_strict (var_char));
  }
  // operator[char*] accessor for words of any length
  PVar operator[] (const char* var_str) const { return PGroup::operator[] (word2index (var_str)); }
  // operator[vector<int>] accessor for words that have already been tokenized (allows us to access even without an Alphabet object)
  PVar operator[] (const vector<int>& intvec) const { return PGroup::operator[] (intvec2index (intvec)); }

  // helper functions to convert words <-> PVar indices / integer vectors
  // ordering is little-endian (i.e. first char has multiplier 1)
  int word2index (const char* var_str) const { return intvec2index (word2intvec (var_str)); }
  vector<int> word2intvec (const char* var_str) const
  {
    if ((int) strlen(var_str) != word_len)
      THROWEXPR ("Tried to get PVar '" << var_str << "' from wordlength-" << word_len << " Alphabet PGroup");
    vector<int> word (word_len);
    for (int i = 0; i < word_len; ++i)
      if (!_alphabet->contains_strict (var_str[i]))
	THROWEXPR ("PVar '" << var_str << "' contains unknown chars")
	  else
	    word[i] = _alphabet->char2int_strict (var_str[i]);
    return word;
  }
  int intvec2index (const vector<int>& int_word) const
  {
    int var_idx = 0;
    int var_mul = 1;
    for (int i = (big_endian ? 0 : word_len - 1);
	 big_endian ? (i < word_len) : (i >= 0);
	 (big_endian ? ++i : --i),
	   var_mul *= alphabet_size)
      var_idx += var_mul * int_word[i];
    return var_idx;
  }
  sstring index2word (int var_idx) const { return intvec2word (index2intvec (var_idx)); }
  vector<int> index2intvec (int var_idx) const
  {
    vector<int> word (word_len);
    int var_mul = 1;
    for (int i = (big_endian ? 0 : word_len - 1);
	 big_endian ? (i < word_len) : (i >= 0);
	 (big_endian ? ++i : --i),
	   var_mul *= alphabet_size)
      word[i] = (var_idx / var_mul) % alphabet_size;
    return word;
  }
  sstring intvec2word (const vector<int>& int_word) const
  {
    sstring word (word_len);
    for (int i = 0; i < word_len; ++i)
      word[i] = _alphabet->int2char_uc (int_word[i]);
    return word;
  }
};

// template for parameter space classes
template<class T>
class PVar_container : Stream_saver
{
private:
  // private data
  set<sstring> group_name_set;  // fast indexing for group_name
protected:
  // protected method _new_group_index allocates a new group index & returns it
  int _new_group_index (unsigned int size, T default_val, const char* name, const vector<sstring>& suffix)
  {
    // get group index
    int group_idx = group.size();
    group.push_back (vector<T> (size, default_val));

    // make group name
    sstring gname;
    for (int gnum = 1;; gnum++)
      {
	gname.clear();
	if (name) gname << name; else gname << "Group" << group_idx;
	if (gnum > 1) gname << gnum;
	if (group_name_set.find (gname) == group_name_set.end())
	  break;
      }

    // store and return
    group_name.push_back (gname);
    group_name_set.insert (gname);
    group_suffix.push_back (suffix);

    // group index
    return group_idx;
  }
  // overloaded version of same method, passing suffix as a const char* for backward compatibility
  int _new_group_index (unsigned int size, T default_val, const char* name = 0, const char* suffix = "")
  {
    vector<sstring> suffix_vec;
    if (suffix)
      for (int i = 0; i < (int) strlen(suffix); ++i)
	{
	  suffix_vec.push_back (sstring());
	  suffix_vec.back().push_back ((char) suffix[i]);
	}
    return _new_group_index (size, default_val, name, suffix_vec);
  }
  // yet another overloaded version, that figures out the group size by splitting suffix on whitespace
  int _new_group_index (T default_val, const char* name, const char* suffix)
  {
    const sstring suffix_string (suffix);
    const vector<sstring> suffix_vec = suffix_string.split();
    return _new_group_index (suffix_vec.size(), default_val, name, suffix_vec);
  }
public:
  // data
  vector<vector<T> >       group;    // access order is group[group_idx][var_idx]
  vector<sstring>          group_name;
  vector<vector<sstring> > group_suffix;

  // accessors
  int groups() const { return group.size(); }
  int group_size (int group_idx) const { return group[group_idx].size(); }
  int total_pvars() const
  {
    int tot = 0;
    for (int g = 0; g < groups(); ++g)
      tot += group_size (g);
    return tot;
  }

  void delete_group (int group_idx)
  {
    const sstring gname = group_name[group_idx];
    group.erase (group.begin() + group_idx);
    group_name.erase (group_name.begin() + group_idx);
    group_name_set.erase (gname);
    group_suffix.erase (group_suffix.begin() + group_idx);
  }

  bool contains_suffix (const sstring& suffix) const
  {
    for (int g = 0; g < groups(); ++g)
      for (int v = 0; v < (int) group[g].size(); ++v)
	if (group_suffix[g][v] == suffix)
	  return true;
    return false;
  }

  T&       operator[] (const PVar& pv) { return group[pv.group_idx][pv.var_idx]; }
  const T& operator[] (const PVar& pv) const { return group[pv.group_idx][pv.var_idx]; }

  T&       operator[] (const sstring& suffix)
  {
    for (int g = 0; g < groups(); ++g)
      for (int v = 0; v < (int) group[g].size(); ++v)
	if (group_suffix[g][v] == suffix)
	  return group[g][v];
    THROWEXPR ("PVar '" << suffix << "' not found");
  }

  const T& operator[] (const sstring& suffix) const
  {
    for (int g = 0; g < groups(); ++g)
      for (int v = 0; v < (int) group[g].size(); ++v)
	if (group_suffix[g][v] == suffix)
	  return group[g][v];
    THROWEXPR ("PVar '" << suffix << "' not found");
  }

  vector<T>&       operator[] (const PGroup& pg) { return group[pg.group_idx]; }
  const vector<T>& operator[] (const PGroup& pg) const { return group[pg.group_idx]; }

  // display method
  void     show (ostream& o, const char* entries_descriptor = "entries") const;

  // dimension equality comparison
  template<class S>
  bool     same_dimensions (const PVar_container<S>& s) const
  {
    if (groups() != s.groups()) return 0;
    for (int g = 0; g < groups(); ++g)
      if (group[g].size() != s.group[g].size()) return 0;
    return 1;
  }

  // dimension equality assertion
  template<class S>
  void     assert_same_dimensions (const PVar_container<S>& s) const
  {
    if (!same_dimensions(s))
      THROW Standard_exception ("Attempt to use parameter spaces of mismatched size");
  }
};

// PScope is a virtual class for allocating PGroups & finding dimensions of a PGroup space
// i.e. an abstract wrapper for the PVar_container template
struct PScope
{
  // query methods
  virtual int pgroups() const = 0;
  virtual int pgroup_size (int pgroup_idx) const = 0;
  // maybe put pgroup_name() & pgroup_suffix() methods here as well?
  // should definitely have separate Telegraph name methods; for now, use PVar defaults ("a", "b" etc)

  // alloc methods
  // it's kind of crazy having all these overloaded virtual new_group() methods that do almost the same thing here...
  // should probably replace most of them with concrete methods that call a single virtual, instead. oh well.
  virtual PGroup         new_group (const char* name, const char* suffix) = 0;  // number of whitespace-separated suffix fields determines size
  virtual PGroup         new_group (unsigned int size, const char* name = 0, const char* suffix = 0) = 0;
  virtual PGroup         new_group (unsigned int size, const char* name, const vector<sstring>& suffix) = 0;
  virtual Boolean_group  new_boolean_group (const char* name = 0) = 0;
  virtual Alphabet_group new_alphabet_group (const Alphabet& a, const char* name = 0, bool bigend = true) = 0;
  virtual Alphabet_group new_alphabet_group (const Alphabet& a, int word_len, const char* name = 0, bool bigend = true) = 0;
  virtual Alphabet_group new_alphabet_group (int alphabet_size, int word_len, const char* name = 0, bool bigend = true) = 0;
  virtual Alphabet_group new_alphabet_group (const vector<sstring>& alphabet_tokens, int word_len, const char* name = 0, bool bigend = true) = 0;

  // virtual destructor
  virtual ~PScope() { }
};

// PScores class
// this is the main class for specifying both the dimensions of a parameter space & the scores for each parameter.

class PScores : public PScope, public PVar_container<FScore>
{
public:
  PScores() { }

  // PScope methods

  int pgroups() const;
  int pgroup_size (int pgroup_idx) const;

  PGroup         new_group (const char* name, const char* suffix);
  PGroup         new_group (unsigned int size, const char* name = 0, const char* suffix = 0);
  PGroup         new_group (unsigned int size, const char* name, const vector<sstring>& suffix);
  Boolean_group  new_boolean_group (const char* name = 0);
  Alphabet_group new_alphabet_group (const Alphabet& a, const char* name = 0, bool bigend = true);
  Alphabet_group new_alphabet_group (const Alphabet& a, int word_len, const char* name = 0, bool bigend = true);
  Alphabet_group new_alphabet_group (int alphabet_size, int word_len, const char* name = 0, bool bigend = true);
  Alphabet_group new_alphabet_group (const vector<sstring>& alphabet_tokens, int word_len, const char* name = 0, bool bigend = true);

  // helper methods
  void normalise();
  void set_null_model (const Sequence_database& seq_db, const Alphabet_group& null_emit, const Boolean_group& null_extend);

  // debugging output method
  void show (ostream&o) const { PVar_container<FScore>::show(o,"scores"); }
  
  // crude I/O
  void write (ostream& out) const;
  void read (istream& in);  // dies if dimensions don't match input
};

// PCounts class
// specifies counts for EM updates
class PCounts : public PVar_container<Prob>
{
public:
  // data
  vector<double> wait;  // wait times for one-variable PGroups

  // constructors
  PCounts();
  PCounts (const PScores& var_scores);

  // virtual destructor
  virtual ~PCounts() { }

  // builders
  void restructure (const PScores& var_scores);
  virtual void clear();

  // helpers
  void optimise_group (int group_idx, PScores& scores) const;    // uses these counts to optimise scores for the specified parameter group
  void optimise (PScores& scores, const set<int>& mutable_groups) const;    // uses these counts to optimise scores for the specified set of parameter groups
  void optimise (PScores& scores) const;    // uses these counts to optimise scores for all parameter groups

  void randomize (PScores& scores, const set<int>& pgroups_to_randomize);
  Loge log_prior (const PScores& scores, const set<int>& pgroups_to_include);

  PCounts& operator+= (const PCounts& counts);

  // output
  void show (ostream&o) const;

private:
  // private builders
  void new_group (unsigned int size, const char* name, const char* suffix)
  {
    _new_group_index (size, 0., name, suffix);
    wait.push_back (0.);
  }

  void new_group (unsigned int size, const char* name, const vector<sstring>& suffix)
  {
    _new_group_index (size, 0.0, name, suffix);
    wait.push_back (0.);
  }
};

// PCounts_like class: PCounts with log-likelihood
class PCounts_like : public PCounts
{
public:
  Loge log_likelihood;
  PCounts_like (const PScores& var_scores) : PCounts (var_scores), log_likelihood (0.0) { }
  void clear();
};

// template method definitions

template<class T>
void PVar_container<T>::show (ostream& o, const char* entries_descriptor) const
{
  int old_prec = o.precision(3);
  save_flags(o);
  right_align(o);
  
  o << "PVar " << entries_descriptor << ":\n";
  for (int g = 0; g < groups(); ++g)
    {
      for (int p = 0; p < (int) group[g].size(); ++p)
	{
	  PVar id (g, p);
	  o.width(6);
	  id.show (o, &group_suffix);
	  o << '=';
	  o.width(10);
	  o << group[g][p];
	  if (p < ((int) group[g].size()) - 1) o << ' ';
	}
      o << "\t" << group_name[g] << '\n';
    }
  restore_flags(o);
  o.precision (old_prec);
}

#endif
