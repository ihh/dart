#include <algorithm>

#include "util/sexpr.h"
#include "util/logfile.h"

// amount of context that's printed for missing parenthesis error messages
#define ERR_CONTEXT 40

SExpr::SExpr()
{ }

SExpr::SExpr (Ptr begin, Ptr end)
{
  init (begin, end);
}

SExpr::SExpr (const char* atom) : atom (atom)
{ }

void SExpr::init (Ptr begin, Ptr end)
{
  // set member vars
  atom.clear();
  child.clear();

  // set up stacks
  stack<SExprIter> ancestral_sexpr;
  stack<Ptr> ancestral_ptr;
  bool in_word = false, in_comment = false;
  Ptr word_begin;

  // parse text
  for (Ptr ptr = begin; ptr != end; ++ptr)
    {
      // get current SExpr
    GetCurrent:
      SExpr& current (ancestral_sexpr.size() ? *ancestral_sexpr.top() : *this);

      // parse next char
      if (in_comment)
	{
	  if (*ptr == '\n' || *ptr == '\r')
	    in_comment = false;
	}

      else if (in_word)
	switch (*ptr)
	  {
	  case ' ':
	  case '\t':
	  case '\n':
	  case '\r':
	  case '\0':
	  case '(':
	  case ')':
	  case ';':
	    current.atom = sstring (word_begin, ptr);
	    in_word = false;
	    pop_or_stop (ancestral_sexpr, ancestral_ptr, ptr, begin);
	    goto GetCurrent;

	  default:
	    break;
	  }

      else
	switch (*ptr)
	  {
	  case ' ':
	  case '\t':
	  case '\n':
	  case '\r':
	  case '\0':
	    break;

	  case '(':
	    current.child.push_back (SExpr());
	    ancestral_sexpr.push (--current.child.end());
	    ancestral_ptr.push (ptr);
	    break;

	  case ')':
	    pop_or_stop (ancestral_sexpr, ancestral_ptr, ptr, begin);
	    break;

	  case ';':
	    in_comment = true;
	    break;

	  default:
	    in_word = true;
	    word_begin = ptr;
	    current.child.push_back (SExpr());
	    ancestral_sexpr.push (--current.child.end());
	    ancestral_ptr.push (ptr);
	    break;
	  }
    }

  if (in_word)
    {
      SExpr& current (*ancestral_sexpr.top());
      current.atom = sstring (word_begin, end);
      pop_or_stop (ancestral_sexpr, ancestral_ptr, end, begin);  // tests if current is valid
    }

  if (ancestral_sexpr.size())
    THROWEXPR ("SExpr: Missing ')' at end of S-expression; corresponding left-bracket context '" << get_context (ancestral_ptr.top(), begin) << "'");
}

void SExpr::pop_or_stop (stack<SExprIter>& sexpr_stack, stack<Ptr>& ptr_stack, Ptr& ptr, Ptr& begin)
{
  if (!sexpr_stack.size())
    THROWEXPR ("SExpr: Missing '(', context '" << get_context (ptr, begin) << "'");
  sexpr_stack.pop();
  ptr_stack.pop();
}

sstring SExpr::get_context (Ptr& ptr, Ptr& begin)
{
  sstring context;
  int context_left;
  Ptr err_ptr;
  bool hide_space = true;
  for (err_ptr = ptr, context_left = ERR_CONTEXT; err_ptr >= begin && context_left > 0; --err_ptr)
    if (isspace (*err_ptr))
      {
	if (!hide_space)
	  context << ' ';
	hide_space = true;
      }
    else
      {
	context << *err_ptr;
	--context_left;
	hide_space = (*err_ptr == '(' || *err_ptr == ')');
      }
  if (err_ptr > begin)
    context << "...";
  reverse (context.begin(), context.end());
  return context;
}

void SExpr::swap (SExpr& sexpr)
{
  atom.swap (sexpr.atom);
  child.swap (sexpr.child);
}

ostream& operator<< (ostream& out, const SExpr& sexpr)
{
  if (sexpr.is_atom())
    out << sexpr.atom;
  else
    {
      bool word_boundary = false;
      for_const_contents (list<SExpr>, sexpr.child, c)
	{
	  if (word_boundary)
	    out << ' ';
	  if (c->is_atom())
	    out << c->atom;
	  else
	    out << '(' << *c << ')';
	  word_boundary = true;
	}
    }
  return out;
}

sstring& SExpr::tag()
{
  if (is_atom())
    THROWEXPR ("In SExpr (" << *this << "):\nAttempt to find name of atom");
  if (is_empty_list())
    THROWEXPR ("In SExpr (" << *this << "):\nAttempt to find name of empty list");
  if (!(*this)[0].is_atom())
    THROWEXPR ("In SExpr (" << *this << "):\nAttempt to call tag() method, which requires that first child is an atom");
  return (*this)[0].atom;
}

SExpr& SExpr::value()
{
  if (!has_value())
    THROWEXPR ("In SExpr (" << *this << "):\nMissing or multiple values");
  return (*this)[1];
}

vector<SExpr*> SExpr::values()
{
  if (!has_tag())
    THROWEXPR ("In SExpr: Missing tag");

  vector<SExpr*> vals;
  for (SExpr::SExprIter child_iter = ++child.begin();
       child_iter != child.end();
       ++child_iter)
    vals.push_back (&*child_iter);

  return vals;
}

SExpr* SExpr::find (const sstring& child_tag, int offset)
{
  for_contents (list<SExpr>, child, c)
    if (--offset < 0)
      if (c->has_tag())
	if (c->tag() == child_tag)
	return &*c;
  return 0;
}

SExpr* SExpr::find (const char* child_tag, int offset)
{
  return find (sstring (child_tag), offset);
}

SExpr& SExpr::find_or_die (const sstring& child_tag, int offset)
{
  SExpr* sexpr = find (child_tag, offset);
  if (!sexpr)
    THROWEXPR ("In (" << *this << "):\nCouldn't find tag (" << child_tag << ")");
  return *sexpr;
}

SExpr& SExpr::find_or_die (const char* child_tag, int offset)
{
  SExpr* sexpr = find (child_tag, offset);
  if (!sexpr)
    THROWEXPR ("In (" << *this << "):\nCouldn't find tag (" << child_tag << ")");
  return *sexpr;
}

vector<SExpr*> SExpr::find_all (const sstring& child_tag, int offset)
{
  vector<SExpr*> matches;
  for_contents (list<SExpr>, child, c)
    if (--offset < 0)
      if (c->has_tag())
	if (c->tag() == child_tag)
	  matches.push_back (&((SExpr&) *c));
  return matches;
}

vector<SExpr*> SExpr::find_all (const char* child_tag, int offset)
{
  return find_all (sstring (child_tag), offset);
}

const sstring& SExpr::get_atom() const
{
  if (!is_atom())
    THROWEXPR ("In SExpr (" << *this << "):\nSExpr is not an atom");
  return atom;
}

vector<sstring> SExpr::atoms_to_strings (int offset)
{
  if (is_atom())
    THROWEXPR ("In SExpr (" << *this << "):\nSExpr is not a list");
  vector<sstring> strings;
  for_contents (list<SExpr>, child, c)
    if (--offset < 0)
      {
	if (!c->is_atom())
	  THROWEXPR ("In SExpr (" << *this << "):\nSExpr (" << *c << ") is not an atom");
	strings.push_back (c->atom);
      }
  return strings;
}

SExpr& SExpr::operator[] (int n_child)
{
  const int for_err = n_child;
  SExprIter i;
  for (i = child.begin(); i != child.end() && n_child > 0; --n_child, ++i)
    { }
  if (i == child.end())
    THROWEXPR ("In SExpr (" << *this << "):\nChild [" << for_err << "] not found");
  return *i;
}

SExpr& SExpr::operator() (const sstring& child_tag, int offset)
{
  SExpr* c = find (child_tag, offset);
  if (!c)
    THROWEXPR ("In SExpr (" << *this << "):\nChild with tag (" << child_tag << ") not found");
  if (c->child.size() < 2)
    THROWEXPR ("In SExpr (" << *this << "):\nChild with tag (" << child_tag << ") has no value");
  if (c->child.size() > 2)
    CLOGERR << "WARNING -- in SExpr (" << *this << "):\nWARNING -- child with tag (" << child_tag << ") has multiple values; I'm just looking at the first value";
  return c->value();
}

SExpr& SExpr::operator() (const char* child_tag, int offset)
{

  return (*this) (sstring (child_tag), offset);
}

bool SExpr::operator== (const SExpr& sexpr) const
{
  if (is_atom())
    return sexpr.is_atom() && atom == sexpr.atom;
  return sexpr.is_list() && child == sexpr.child;
}

bool SExpr::operator!= (const SExpr& sexpr) const
{
  return !(*this == sexpr);
}

SExpr_file::SExpr_file (const vector<sstring>& filename)
{
  if (filename.size())
    for_const_contents (vector<sstring>, filename, fn)
      read_text_from_file (fn->c_str());
  else
    read_text_from_stdin();
  parse_text();
}

SExpr_file::SExpr_file (const char* filename)
{
  if (filename)
    read_text_from_file (filename);
  else
    read_text_from_stdin();
  parse_text();
}

void SExpr_file::read_text_from_file (const char* filename)
{
  ifstream in (filename);
  if (!in)
    THROWEXPR ("SExpr_file: file '" << filename << "' not found");
  read_text_from_stream (in);
}

void SExpr_file::read_text_from_stream (istream& in)
{
  sstring s;
  while (in && !in.eof())
    {	
      s.getline (in);
      text << s;
    }
}

void SExpr_file::read_text_from_stdin()
{
  CLOGERR << "[waiting for S-expression on standard input]\n";
  read_text_from_stream (cin);
}

void SExpr_file::parse_text()
{
  sexpr.init (text.begin(), text.end());
}

Regexp SExpr_validator::brackets_regexp ("\\([ \t]*([^ \t\\(\\)]+)[ \t]*\\)");
Regexp SExpr_validator::list_regexp ("([^ \t\\(\\)]+)\\*");
Regexp SExpr_validator::tagval_regexp ("\\([ \t]*([^ \t\\(\\)]+)[ \t]+([^ \t\\(\\)]+)[ \t]*\\)");
Regexp SExpr_validator::nonwhite_regexp ("[ \t]*([^ \t]+)[ \t]*");
Regexp SExpr_validator::first_nonterm_regexp ("([^ \t\\-]+)");

bool SExpr_validator::parse (SExpr_iterator begin, SExpr_iterator end, bool issue_warnings)
{
  if (first_nonterm_regexp.Match (grammar.c_str()))
    return parse (first_nonterm_regexp[1], begin, end, issue_warnings);
  CLOGERR << "Warning: couldn't find top-level nonterminal in grammar:\n" << grammar << "\n";
  return false;
}

bool SExpr_validator::parse (sstring nonterm, SExpr_iterator begin, SExpr_iterator end, bool issue_warnings)
{
  SExpr_iterator begin_plus_one = begin;
  ++begin_plus_one;
  if (nonterm == "atom") {
    const bool is_atom = end == begin_plus_one && begin->is_atom();
    if (!is_atom && issue_warnings)
      warn (nonterm, begin, end);
    return is_atom;
  } else if (nonterm == "wild") {
    return true;
  } else if (nonterm == "end") {
    const bool is_empty = begin == end;
    if (!is_empty && issue_warnings)
      warn (nonterm, begin, end);
    return is_empty;
  } else {
    sstring rule;
    rule << nonterm << "[ \t]*->([^;]+);";
    Regexp rule_regexp (rule.c_str());
    if (rule_regexp.Match(grammar.c_str())) {
      sstring rhs_string = rule_regexp[1];
      vector<sstring> rhs_options = rhs_string.split("|");
      for_const_contents (vector<sstring>, rhs_options, rhs_opt) {
	if (match_rhs (*rhs_opt, begin, end))
	  return true;
      }
      if (issue_warnings)
	warn (nonterm, begin, end);
      return false;
    } else {  // no grammar rule matched
      CLOGERR << "Nonterminal " << nonterm << " not found in grammar\n";
      return false;
    }
  }
  return false;
}

bool SExpr_validator::match_rhs (sstring rhs, SExpr_iterator begin, SExpr_iterator end)
{
  SExpr_iterator begin_plus_one = begin;
  ++begin_plus_one;
  const char* s = rhs.c_str();
  if (brackets_regexp.Match(s)) {
    if (end == begin_plus_one && begin->is_list())
      return parse (brackets_regexp[1], *begin, false);
    return false;
  } else if (list_regexp.Match(s)) {
    sstring list_nonterm = list_regexp[1];
    for (SExpr_iterator iter = begin; iter != end; ++iter)
      if (!parse (list_nonterm, *iter, false))
	return false;
    return true;
  } else if (tagval_regexp.Match(s)) {
    sstring keyword = tagval_regexp[1];
    sstring tagval_nonterm = tagval_regexp[2];
    if (end != begin_plus_one || !begin->is_list())
      return false;
    if (begin->child.size() != 2 || !begin->child.front().is_atom() || begin->child.front().atom != keyword)
      return false;
    return parse (tagval_nonterm, begin->child.back(), false);
  } else if (nonwhite_regexp.Match(s)) {
    // assume RHS represents a nonterminal
    return parse (nonwhite_regexp[1], begin, end, false);
  }
  // RHS unrecognizable
  CLOGERR << "Warning: can't understand RHS option in grammar rule: " << rhs << "\n";
  return false;
}

void SExpr_validator::warn (sstring nonterm, SExpr_iterator begin, SExpr_iterator end)
{
  if (!warned)
    {
      CLOGERR << "Warning: could not match the following S-expression to nonterminal " << nonterm << "\n";
      for (SExpr_iterator iter = begin; iter != end; ++iter)
	CLOGERR << " " << *iter;
      CLOGERR << "\n";
      warned = true;
    }
}

