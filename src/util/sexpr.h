#ifndef SEXPR_INCLUDED
#define SEXPR_INCLUDED

#include <list>
#include <stack>
#include <map>
#include <vector>
#include "util/sstring.h"
#include "util/Regexp.h"

struct SExpr
{
  // typedefs
  typedef basic_string<char>::const_iterator Ptr;
  typedef list<SExpr>::iterator SExprIter;

  // data
  // note: it should never be the case that both atom.size() and child.size() are nonzero
  sstring atom;
  list<SExpr> child;

  // constructors
  // (default copy constructor automatically does a deep copy)
  SExpr();  // creates an empty SExpr
  SExpr (Ptr begin, Ptr end);  // parses a text block into an SExpr
  SExpr (const char* atom);  // creates an atomic SExpr

  // initialiser
  void init (Ptr begin, Ptr end);

  // swap method
  void swap (SExpr& sexpr);

  // accessors
  inline bool is_list() const;
  inline bool is_empty_list() const;
  inline bool is_atom() const;

  // tag() returns the atom of the first child, or an error if that doesn't exist
  inline bool has_tag();
  sstring& tag();

  // value() returns the second child, or an error if there aren't exactly two children
  inline bool has_value();
  SExpr& value();

  // values() returns all children but the first, or an error if there are no children
  vector<SExpr*> values();

  // find() returns iterator for first child past offset with a given tag, or null if not found
  SExpr* find (const sstring& child_tag, int offset = 0);
  SExpr* find (const char* child_tag, int offset = 0);

  // find_all() returns a vector of iterators for all children past offset with a given name
  vector<SExpr*> find_all (const sstring& child_tag, int offset = 0);
  vector<SExpr*> find_all (const char* child_tag, int offset = 0);

  // find_or_die() wraps find() with a null pointer test
  SExpr& find_or_die (const sstring& child_tag, int offset = 0);
  SExpr& find_or_die (const char* child_tag, int offset = 0);

  // get_atom() throws an error if this is not an atom
  const sstring& get_atom() const;

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

  // output method
  friend ostream& operator<< (ostream& out, const SExpr& sexpr);
};

// SExpr input from file, list of files, or standard input
struct SExpr_file
{
  // contents of file
  sstring text;
  SExpr sexpr;
  // constructors
  SExpr_file (const char* filename);  // null pointer will use stdin
  SExpr_file (const vector<sstring>& filename);  // empty list will use stdin
  // builders -- treat as private
  void read_text_from_file (const char* filename);
  void read_text_from_stream (istream& in);
  void read_text_from_stdin();
  void parse_text();
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
  int warnings;
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
  return child.size() == 2;
}

#endif /* SEXPR_INCLUDED */
