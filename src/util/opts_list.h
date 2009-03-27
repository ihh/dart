#ifndef OPTS_LIST_INCLUDED
#define OPTS_LIST_INCLUDED

#include <iostream>
#include <math.h>
#include <set>
#include <map>
#include <list>
#include "util/sstring.h"
#include "util/dexception.h"

// parser class for command line options & arguments of an executable
struct Opts_list {

  // syntax error exception class... this is really archaic stuff
  struct Syntax_exception : Dart_exception
  {
    const Opts_list&  opts;  // pointer to parent class
    sstring           temp;  // temporary variable used by Dart_exception method
    sstring           info;  // used to hold the error message

    Syntax_exception (const Opts_list& opts, const char* msg) :
      Dart_exception(), opts(opts), info(msg) {}
    Syntax_exception (const Opts_list& opts, const char* msg1, const char* msg2) :
      Dart_exception(), opts(opts), info(msg1) { info.append(msg2); }

    // overrides Dart_exception method
    const char* details() const
      {
	(sstring&) temp = "";
	((sstring&) temp).append("While parsing command line: ").append(info).append("\n");
	return temp.c_str();
      }
  };

  // typedef for "option handlers" (callback member functions)
  typedef bool (*Option_handler) (Opts_list*);  // returns TRUE if parsed OK, FALSE if error

  // type of option taking values "yes", "no" and "auto" (used in one piece of code dating from 1999, AFAIK...)
  enum Expect_flag { YES, NO, AUTO };
  void set_expect_flag (Expect_flag& flag, const char* value_str);

  // member variables

  // "current" argc and argv - these are used as pointers and are changed by the option parse code
  int    argc;
  char** argv;

  // aliased arguments
  list<sstring> alias_args;
  list<sstring>::iterator next_alias_arg;

  // method to check if there are more args remaining
  bool more_args() const;

  // initial values of argc and argv
  int    init_argc;
  char** init_argv;

  int expect_args;   // expected number of args: -1 for "any", otherwise will complain if not equal to this

  // all options
  set<sstring> all_opts;
  void new_opt (const sstring& s);
  static sstring neg_opt (const sstring& s);  // puts a "no" in front of an option

  // options by type
  map <sstring, bool*>          bool_opts;
  map <sstring, bool*>          bool_no_opts;
  map <sstring, int*>           int_opts;
  map <sstring, double*>        double_opts;
  map <sstring, sstring*>       string_opts;
  map <sstring, Option_handler> callback_opts;
  map <sstring, sstring>        alias_opts;

  // post-parse callback hooks
  vector <Option_handler>       post_parse_callback;

  // various bits of help text
  sstring program_name;
  sstring short_description;
  sstring version_info;
  sstring syntax;
  sstring short_help_text;
  sstring options_help_text;

  // arguments (extracted from the command line & stuck into this vector by the option-parsing code)
  vector<sstring> args;

  // member functions

  // constructor
  Opts_list (int argc, char** argv);

  // virtual destructor
  virtual ~Opts_list();

  // builder methods to add options, with comments in help text
  void add (const char* opt, bool& var, const char* desc = 0, bool show_default = 1, const char* negopt = 0, const char* negdesc = 0);
  void add (const char* opt, int& var,     const char* desc = 0, bool show_default = true);
  void add (const char* opt, double& var,  const char* desc = 0, bool show_default = true);
  void add (const char* opt, sstring& var, const char* desc = 0, bool show_default = true);
  void add (const char* opt, const char* alias, const char* desc = 0, bool show_alias = true);
  void add (const char* opt, Option_handler callback, const char* desc = 0);

  // builder methods to add comments to help text
  void print (const char* text);
  void newline();
  void print_title (const char* text);

  // parser methods
  bool parse();  // returns TRUE if parsed OK
  void parse_or_die();

  // virtual parse method
  enum Parse_status { UNPARSED = 0, PARSE_OK = 1, PARSE_NOT_OK = 2 };
  virtual Parse_status parse_opt (const sstring& opt);

  // parser helper methods
  double next_double();
  int    next_int();
  char*  next_string();

  // help text accessors
  virtual sstring short_help() const;  // prints program_name/short_description/syntax, short_help_text
  virtual sstring help() const;  // prints program_name/short_description/syntax, options_help_text

  // option handler functions to display various bits of help text
  static bool display_help (Opts_list* ol);
  static bool display_version (Opts_list* ol);
};

// build macros
#define SET_VERSION_INFO(OPTS) (OPTS).version_info.clear(); (OPTS).version_info << "DART version " << DART_VERSION_STRING << " (" << DART_DEBUG_MODE_STRING << ") compiled " << __DATE__ << " " << __TIME__ << "\n";

#define INIT_CONSTRUCTED_OPTS_LIST(OPTS,ARGS,SYNTAX,SHORTDESC) OPTS.short_description = (SHORTDESC); OPTS.syntax = (SYNTAX); OPTS.expect_args = (ARGS); SET_VERSION_INFO(OPTS); Rnd::add_opts(OPTS); OPTS.newline(); Log_stream::add_opts(OPTS);

#define INIT_TYPED_OPTS_LIST(OPTS_TYPE,OPTS,ARGC,ARGV,ARGS,SYNTAX,SHORTDESC) OPTS_TYPE OPTS(ARGC,ARGV); INIT_CONSTRUCTED_OPTS_LIST(OPTS,ARGS,SYNTAX,SHORTDESC)

#define INIT_OPTS_LIST(OPTS,ARGC,ARGV,ARGS,SYNTAX,SHORTDESC) INIT_TYPED_OPTS_LIST(Opts_list,OPTS,ARGC,ARGV,ARGS,SYNTAX,SHORTDESC)

#endif
