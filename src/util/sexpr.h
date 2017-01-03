#ifndef SEXPR_INCLUDED
#define SEXPR_INCLUDED

#include <list>
#include <stack>
#include <map>
#include <vector>
#include <set>

#include "util/sstring.h"
#include "util/Regexp.h"

struct SExpr_atom : sstring
{
  // data
  bool is_quoted;  // if true, then the atom is in double quotes, i.e. a Scheme string
  // constructors
  SExpr_atom() : sstring(), is_quoted(false) { }
  SExpr_atom(const sstring& s) : sstring(s), is_quoted(false) { }
  SExpr_atom(const char* s) : sstring(s), is_quoted(false) { }
  SExpr_atom(const SExpr_atom& s) : sstring(s), is_quoted(s.is_quoted) { }
  SExpr_atom& operator= (const SExpr_atom& s) { ((sstring&)*this) = s; return *this; }  // prevent "copy assignment operator implicitly deleted" errors
  SExpr_atom& operator= (const sstring& s) { ((sstring&)*this) = s; return *this; }
  SExpr_atom& operator= (const char* s) { ((sstring&)*this) = s; return *this; }
  // vector of strings
  static vector<SExpr_atom> from_vector (const vector<sstring>& s);
  // quoted output
  friend ostream& operator<< (ostream& out, const SExpr_atom& atom);
};

struct SExpr
{
  // typedefs
  typedef basic_string<char>::const_iterator Ptr;
  typedef list<SExpr>::iterator SExprIter;

  // data
  // note: it should never be the case that both atom.size() and child.size() are nonzero
  SExpr_atom atom;
  list<SExpr> child;

  // constructors
  // (default copy constructor automatically does a deep copy)
  SExpr();  // default constructor, creates an empty SExpr: is_empty_list() == true
  SExpr (Ptr begin, Ptr end, bool allow_quotes = true);  // parses a text block into an SExpr
  SExpr (const char* atom);  // creates an atomic SExpr

  // initialiser
  void init (Ptr begin, Ptr end, bool allow_quotes = true);

  // swap method
  void swap (SExpr& sexpr);

  // accessors
  inline bool is_list() const;
  inline bool is_empty_list() const;
  inline bool is_atom() const;

  inline bool has_tag();  // true if there are at least two children and the first child is an atom
  inline bool has_value();  // true if there are exactly two children, or is_scheme_pair()
  inline bool is_scheme_pair();  // true if there are exactly three children and the middle one is the atom "."

  // tag() returns the atom of the first child, or an error if that doesn't exist
  SExpr_atom& tag();

  // value() returns the second child (if there are exactly 2 children), third child (if is_scheme_pair), or an error
  SExpr& value();

  // values() returns all children but the first, or an error if there are no children
  vector<SExpr*> values();

  // find() returns iterator for first child past offset with a given tag (or tags), or null if not found
  SExpr* find (const sstring& child_tag, int offset = 0);
  SExpr* find (const char* child_tag, int offset = 0);

  // find_any() allows multiple tags
  SExpr* find_any (const vector<sstring>& child_tags, int offset = 0);
  SExpr* find_any (const char* child_tag_string, int offset = 0);  // child_tag_string is space-separated

  // find_all() returns a vector of iterators for all children past offset with a given name
  vector<SExpr*> find_all (const sstring& child_tag, int offset = 0);
  vector<SExpr*> find_all (const char* child_tag, int offset = 0);

  // find_or_die() wraps find() with a null pointer test
  SExpr& find_or_die (const sstring& child_tag, int offset = 0);
  SExpr& find_or_die (const char* child_tag, int offset = 0);

  // find_recursive() does a breadth-first search looking for child_tag
  // with max_depth==1, it is equivalent to find()
  SExpr* find_recursive (const sstring& tag, int max_depth = -1);
  SExpr* find_recursive (const char* tag, int max_depth = -1);

  // get_atom() throws an error if this is not an atom
  const SExpr_atom& get_atom() const;

  // atoms_to_strings() converts childrens' atoms into a vector of strings
  vector<sstring> atoms_to_strings (int offset = 0);

  // operator[] returns indexed child, or throws an error
  SExpr& operator[] (int n_child);

  // operator() returns value of first child past offset with given tag, or throws an error
  // it's equivalent to find_or_die(child_tag,offset).value()
  SExpr& operator() (const sstring& child_tag, int offset = 0);
  SExpr& operator() (const char* child_tag, int offset = 0);

  // operator==() is a recursive test for string equality
  bool operator== (const SExpr& sexpr) const;
  bool operator!= (const SExpr& sexpr) const;

  // build helpers
  static void pop_or_stop (stack<SExprIter>& sexpr_stack, stack<Ptr>& ptr_stack, Ptr& ptr, Ptr& begin);
  static sstring get_context (Ptr& ptr, Ptr& begin);

  // output methods
  friend ostream& operator<< (ostream& out, const SExpr& sexpr);  // does NOT add enclosing parentheses to lists
  sstring to_string() const;  // simple wrapper for operator<<
  sstring to_parenthesized_string() const;  // wrapper for to_string() that DOES add enclosing parentheses to lists
};

// SExpr input from file, list of files, or standard input
struct SExpr_file
{
  // contents of file
  sstring text;
  SExpr sexpr;
  // constructors
  SExpr_file (const char* filename, bool allow_quotes = true);  // null pointer will use stdin
  SExpr_file (const vector<sstring>& filename, bool allow_quotes = true);  // empty list will use stdin
  // builders -- treat as private
  void read_text_from_file (const char* filename);
  void read_text_from_stream (istream& in);
  void read_text_from_stdin();
  void parse_text(bool allow_quotes = true);
};

// SExpr syntax validation grammar
struct SExpr_validator
{
  // typedefs
  typedef list<SExpr>::const_iterator SExpr_iterator;
  typedef map<sstring,vector<sstring> > Grammar;
  // statics
  static Regexp rule_regexp, brackets_regexp, list_regexp, quote_regexp, tagval_regexp, leftemit_regexp, nonwhite_regexp;
  // members
  Grammar grammar;  // maps lhs to rhs
  sstring start_nonterm;
  int warnings, duplicate_warnings;
  set<sstring> printed_warnings;
  // constructor
  SExpr_validator (const char* grammar_str);
  // methods
  bool parse (const SExpr& sexpr, bool issue_warnings = true);  // top-level parse method
  bool parse (sstring nonterm, SExpr_iterator begin, SExpr_iterator end, bool issue_warnings);
  bool parse (sstring nonterm, const SExpr& sexpr, bool issue_warnings)
  {
    const list<SExpr> sl (1, sexpr);
    return parse (nonterm, sl.begin(), sl.end(), issue_warnings);
  }
  bool match_rhs (sstring rhs, SExpr_iterator begin, SExpr_iterator end, bool issue_warnings, bool shallow);  // does recursive descent if shallow==false
  void warn (sstring nonterm, SExpr_iterator begin, SExpr_iterator end);
};


// Inline method defs

bool SExpr::is_list() const
{
  return atom.size() == 0;
}

bool SExpr::is_atom() const
{
  return !is_list();
}

bool SExpr::is_empty_list() const
{
  return is_list() && child.begin() == child.end();
}

bool SExpr::has_tag()
{
  return is_list() && !is_empty_list() && (*this)[0].is_atom();
}

bool SExpr::has_value()
{
  return child.size() == 2 || is_scheme_pair();
}

bool SExpr::is_scheme_pair()
{
  return child.size() == 3 && (*this)[1].is_atom() && (*this)[1].get_atom() == ".";
}

#endif /* SEXPR_INCLUDED */
