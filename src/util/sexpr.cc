#include <algorithm>

#include "util/sexpr.h"
#include "util/logfile.h"

// amount of context that's printed for missing parenthesis error messages
#define ERR_CONTEXT 40

SExpr::SExpr()
{ }

SExpr::SExpr (Ptr begin, Ptr end, bool allow_quotes)
{
  init (begin, end, allow_quotes);
}

SExpr::SExpr (const char* atom) : atom (atom)
{ }

void SExpr::init (Ptr begin, Ptr end, bool allow_quotes)
{
  // set member vars
  atom.clear();
  child.clear();

  // set up stacks
  stack<SExprIter> ancestral_sexpr;
  stack<Ptr> ancestral_ptr;
  enum State { InAtom, InQuotes, EscapedInQuotes, InComment, Between } state = Between;
  Ptr word_begin;

  // parse text
  for (Ptr ptr = begin; ptr != end; ++ptr)
    {
      // get current SExpr
      SExpr& current (ancestral_sexpr.size() ? *ancestral_sexpr.top() : *this);

      // parse next char
      switch (state)
	{
	case InComment:
	  //	  CLOGERR << "InComment " << *ptr << '\n';
	  if (*ptr == '\n' || *ptr == '\r')
	    state = Between;
	  break;

	case InAtom:
	  //	  CLOGERR << "InAtom " << *ptr << '\n';
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
	      pop_or_stop (ancestral_sexpr, ancestral_ptr, ptr, begin);
	      state = Between;
	      --ptr;  // back up and consider this character again
	      break;
	      
	    default:
	      break;
	    }
	  break;

	case InQuotes:
	  //	  CLOGERR << "InQuotes " << *ptr << '\n';
	  if (!allow_quotes)
	    THROWEXPR ("SExpr: This program uses a limited S-expression format in which quoted atoms are not allowed.\nContext '" << get_context (ptr, begin) << "'");
	  switch (*ptr)
	    {
	    case '\\':
	      state = EscapedInQuotes;
	      break;

	    case '"':
	      state = Between;
	      pop_or_stop (ancestral_sexpr, ancestral_ptr, ptr, begin);
	      break;
	      
	    default:
	      current.atom.push_back (*ptr);
	      break;
	    }
	  break;

	case EscapedInQuotes:
	  //	  CLOGERR << "EscapedInQuotes " << *ptr << '\n';
	  current.atom.push_back (*ptr);
	  state = InQuotes;
	  break;

	case Between:
	  //	  CLOGERR << "Between " << *ptr << '\n';
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
	      state = InComment;
	      break;

	    case '"':
	      state = InQuotes;
	      current.child.push_back (SExpr());
	      ancestral_sexpr.push (--current.child.end());
	      ancestral_ptr.push (ptr);
	      break;

	    default:
	      state = InAtom;
	      word_begin = ptr;
	      current.child.push_back (SExpr());
	      ancestral_sexpr.push (--current.child.end());
	      ancestral_ptr.push (ptr);
	      break;
	    }
	  break;

	default:
	  THROWEXPR("Unreachable");
	  break;
	}
    }

  switch (state)
    {
    case InAtom:
      {
	SExpr& current (*ancestral_sexpr.top());
	current.atom = sstring (word_begin, end);
	pop_or_stop (ancestral_sexpr, ancestral_ptr, end, begin);  // tests if current is valid
      }
      break;

    case InQuotes:
    case EscapedInQuotes:
      THROWEXPR("SExpr: Missing end-quote '\"' at end of S-expression");
      break;

    default:
      break;
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

sstring SExpr::to_string() const
{
  sstring s;
  s << *this;
  return s;
}

sstring SExpr::to_parenthesized_string() const
{
  sstring s;
  if (is_list())
    s << '(' << *this << ')';
  else
    s << *this;
  return s;
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

Regexp special_char_regexp("[ \t\r\n;\"\\(\\)]");
ostream& operator<< (ostream& out, const SExpr_atom& atom)
{
  if (special_char_regexp.Match(atom.c_str()))
    {
      out << '"';
      for_const_contents (sstring, atom, chr)
	{
	  if (*chr == '"' || *chr == '\\')
	    out << '\\';
	  out << *chr;
	}
      out << '"';
    }
  else
    out << (sstring&) atom;
  return out;
}

vector<SExpr_atom> SExpr_atom::from_vector (const vector<sstring>& s)
{
  vector<SExpr_atom> v;
  v.reserve(s.size());
  for_const_contents (vector<sstring>, s, str)
    v.push_back (SExpr_atom (*str));
  return v;
}

SExpr_atom& SExpr::tag()
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
  sstring child_tag_str (child_tag);
  return find_or_die (child_tag_str, offset);
}

SExpr* SExpr::find_recursive (const sstring& tag, int max_depth)
{
  set<SExpr*> current, next;
  next.insert (this);
  while (next.size() && max_depth-- != 0)
    {
      current.swap(next);
      next.clear();
      for_const_contents (set<SExpr*>, current, sxpr)
	{
	  for_contents (list<SExpr>, (**sxpr).child, c)
	    {
	      if (c->has_tag())
		if (c->tag() == tag)
		  return &*c;
	      next.insert (&*c);
	    }
	}
    }
  return (SExpr*) 0;
}

SExpr* SExpr::find_recursive (const char* tag, int max_depth)
{
  sstring tag_str (tag);
  return find_recursive (tag_str, max_depth);
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

const SExpr_atom& SExpr::get_atom() const
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
    CLOGERR << "WARNING -- in SExpr (" << *this << "):\nWARNING -- child with tag (" << child_tag << ") has multiple values; I'm just looking at the first value\n";
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

SExpr_file::SExpr_file (const vector<sstring>& filename, bool allow_quotes)
{
  if (filename.size())
    for_const_contents (vector<sstring>, filename, fn)
      read_text_from_file (fn->c_str());
  else
    read_text_from_stdin();
  parse_text(allow_quotes);
}

SExpr_file::SExpr_file (const char* filename, bool allow_quotes)
{
  if (filename)
    read_text_from_file (filename);
  else
    read_text_from_stdin();
  parse_text(allow_quotes);
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

void SExpr_file::parse_text(bool allow_quotes)
{
  sexpr.init (text.begin(), text.end(), allow_quotes);
}

Regexp SExpr_validator::brackets_regexp ("^[ \t]*\\([ \t]*([^ \t\\(\\)]+)[ \t]*\\)[ \t]*$");
Regexp SExpr_validator::quote_regexp ("^[ \t]*'([^ \t\\(\\)]+)[ \t]*$");
Regexp SExpr_validator::list_regexp ("^[ \t]*([^ \t]+)(\\*|\\+)[ \t]*$");
Regexp SExpr_validator::tagval_regexp ("^[ \t]*\\([ \t]*([^ \t]+)[ \t]+([^ \t]+)[ \t]*\\)[ \t]*$");
Regexp SExpr_validator::leftemit_regexp ("^[ \t]*([^ \t]+)[ \t]+([^ \t]+)[ \t]*$");
Regexp SExpr_validator::nonwhite_regexp ("^[ \t]*([^ \t]+)[ \t]*$");
Regexp SExpr_validator::rule_regexp ("^[ \t]*([^ \t]+)[ \t]*->[ \t]*(.*)$");

SExpr_validator::SExpr_validator (const char* grammar_str)
  : warnings(0)
{
  const sstring gram (grammar_str);
  const vector<sstring> rule = gram.split(";");
  for_const_contents (vector<sstring>, rule, r) {
    if (rule_regexp.Match(r->c_str())) {
      const sstring lhs = rule_regexp[1];
      const sstring rhs = rule_regexp[2];
      const vector<sstring> rhs_opts = rhs.split("|");
      grammar[lhs] = rhs_opts;
      if (start_nonterm.size() == 0)
	start_nonterm = lhs;
      if (CTAGGING(2,SEXPR_SYNTAX))
	CL << "Rule #" << grammar.size() << ": " << lhs << " -> " << sstring::join(rhs_opts," | ") << ";\n";
    }
  }
}

bool SExpr_validator::parse (const SExpr& sexpr, bool issue_warnings)
{
  warnings = duplicate_warnings = 0;
  printed_warnings.clear();
  bool ok = parse (start_nonterm, sexpr, issue_warnings);
  if (duplicate_warnings > 0)
    CLOGERR << duplicate_warnings << " duplicate syntax warning" << (duplicate_warnings > 1 ? "s" : "") << '\n';
  return ok;
}

bool SExpr_validator::parse (sstring nonterm, SExpr_iterator begin, SExpr_iterator end, bool issue_warnings)
{
  if (CTAGGING(-1,SEXPR_VALIDATOR))
    {
      SExpr dump_sexpr;
      dump_sexpr.child.insert (dump_sexpr.child.begin(), begin, end);
      CL << "Attempting to match " << nonterm << " to " << dump_sexpr << "\n";
    }

  bool list_length_equals_one = false;
  if (begin != end) {
    SExpr_iterator begin_plus_one = begin;
    ++begin_plus_one;
    list_length_equals_one = end == begin_plus_one;
  }

  // implicitly recognized: Atom, Wild, End
  if (nonterm == "Atom") {
    return list_length_equals_one && begin->is_atom();
  } else if (nonterm == "Wild") {
    return true;
  } else if (nonterm == "End") {
    return begin == end;
  }

  // implicitly recognized: lists
  const char* s = nonterm.c_str();
  if (list_regexp.Match(s)) {  // X*, X+
    sstring list_nonterm = list_regexp[1];
    sstring list_star_or_plus = list_regexp[2];
    bool ok = true;
    if (begin == end && list_star_or_plus == "+")
      ok = false;
    int old_warnings = warnings, new_warnings = 0;
    for (SExpr_iterator iter = begin; iter != end; ++iter) {
      warnings = old_warnings;
      if (!parse (list_nonterm, *iter, issue_warnings)) {
	ok = false;
	if (!issue_warnings)  // if not giving warnings, don't bother trying to process the other list items: we have a fail, so bail
	  break;
	new_warnings += warnings - old_warnings;
      }
    }
    warnings = old_warnings + new_warnings;
    return ok;
  }

  // implicitly recognized: brackets
  if (brackets_regexp.Match(s)) {  // (X)
    if (list_length_equals_one && begin->is_list())
      return parse (brackets_regexp[1], begin->child.begin(), begin->child.end(), issue_warnings);
    if (issue_warnings)
      warn (nonterm, begin, end);
    return false;
  }

  // implicitly recognized: quoted atoms
  if (quote_regexp.Match(s)) {  // 'atom
    sstring atom = quote_regexp[1];
    return list_length_equals_one && begin->is_atom() && begin->atom == atom;
  }

  // try to match the nonterminal using the grammar
  Grammar::const_iterator grammar_iter = grammar.find(nonterm);
  if (grammar_iter != grammar.end()) {
    const vector<sstring>& rhs_options = grammar_iter->second;

    // first pass: shallow, without warnings (this pass only performed if issue_warnings==true)
    int shallow_matches = 0;
    if (issue_warnings)
      for_const_contents (vector<sstring>, rhs_options, rhs_opt) {
	if (match_rhs (*rhs_opt, begin, end, issue_warnings, true))
	  ++shallow_matches;
      }

    // second pass: deep, with warnings
    // if shallow_matches==1, then there is a unique shallow match, and so we can peg the warning lower down the tree
    for_const_contents (vector<sstring>, rhs_options, rhs_opt) {
      if (match_rhs (*rhs_opt, begin, end, issue_warnings && shallow_matches == 1, false))
	return true;
    }

    // fail
    if (issue_warnings)
      warn (nonterm, begin, end);
    return false;
  }

  // no grammar rule matched
  CLOGERR << "Nonterminal " << nonterm << " not found in grammar\n";
  return false;
}

bool SExpr_validator::match_rhs (sstring rhs, SExpr_iterator begin, SExpr_iterator end, bool issue_warnings, bool shallow)
{
  bool list_length_equals_one = false;
  SExpr_iterator begin_plus_one = begin;
  if (begin != end) {
    ++begin_plus_one;
    list_length_equals_one = end == begin_plus_one;
  }

  // resolve deterministic transition paths X -> ... -> Y
  Grammar::const_iterator gram_iter;
  while ((gram_iter = grammar.find(rhs)) != grammar.end() && gram_iter->second.size() == 1)
    rhs = gram_iter->second.front();

  const char* s = rhs.c_str();
  if (tagval_regexp.Match(s)) {  // (tag X)
    sstring tag = tagval_regexp[1];
    sstring val = tagval_regexp[2];
    if (!list_length_equals_one || !begin->is_list())
      return false;
    if (begin->child.size() == 0)
      return false;

    // shallow-match quoted tags if we can
    if (quote_regexp.Match(tag.c_str()) && (!begin->child.front().is_atom() || begin->child.front().atom != quote_regexp[1]))
      return false;

    SExpr_iterator child_begin_plus_one = begin->child.begin();
    ++child_begin_plus_one;

    return shallow
      ? true
      : (parse (tag, begin->child.front(), issue_warnings)
	 && parse (val, child_begin_plus_one, begin->child.end(), issue_warnings));

  } else if (leftemit_regexp.Match(s)) {  // tag X
    sstring tag = leftemit_regexp[1];
    sstring val = leftemit_regexp[2];
    if (end == begin)
      return false;

    // shallow-match quoted tags if we can
    if (quote_regexp.Match(tag.c_str()) && (!begin->is_atom() || begin->atom != quote_regexp[1]))
      return false;

    return shallow
      ? true
      : (parse (tag, *begin, issue_warnings)
	 && parse (val, begin_plus_one, end, issue_warnings));

  } else if (nonwhite_regexp.Match(s)) {  // X
    sstring nonterm = nonwhite_regexp[1];
    return shallow ? true : parse (nonterm, begin, end, issue_warnings);
  }

  // RHS unrecognizable
  CLOGERR << "Warning: can't understand RHS option in syntax validation grammar rule: " << rhs << "\n";
  return false;
}

void SExpr_validator::warn (sstring nonterm, SExpr_iterator begin, SExpr_iterator end)
{
  if (warnings == 0)
    {
      SExpr dump_sexpr;
      dump_sexpr.child.insert (dump_sexpr.child.begin(), begin, end);
      sstring warning;
      warning << "Syntax warning: the following does not look like a " << nonterm << ": " << dump_sexpr;
      if (printed_warnings.find (warning) == printed_warnings.end()) {
	CLOGERR << warning << "\n";
	printed_warnings.insert(warning);
      } else
	++duplicate_warnings;
    }
  ++warnings;
}
