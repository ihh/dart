#ifndef LOGFILE_INCLUDED
#define LOGFILE_INCLUDED

#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <fstream>
#include <list>
#include <map>
#include <stack>

#include "util/sstring.h"
#include "util/macros.h"
#include "util/opts_list.h"
#include "util/Regexp.h"

#define CLOGXMLTAGNAME    "log"
#define DEBUG_LOG_LEVEL   0
#define INITIAL_LOG_LEVEL 5
#define QUIET_LOG_LEVEL   9
#define ERROR_LOG_LEVEL   10
#define ERROR_TAG         "ERROR"

#define CLOGSTREAM             clog_stream
#define CLOGSTREAMBUF          clog_streambuf
#define CLOG_FILE_LINE         file_tag(__FILE__,__LINE__,"")
#define CLOG_FILE_LINE_TAGS(T) file_tag(__FILE__,__LINE__,"" #T "")

#define CLOGREQUEST(n,STDOUT_FLAG)      CLOGSTREAM.request_logging(n,__FILE__ " ",__LINE__,STDOUT_FLAG)
#define CTAGREQUEST(n,TAGS,STDOUT_FLAG) CLOGSTREAM.request_logging(n,__FILE__ " " #TAGS " ",__LINE__,STDOUT_FLAG)

#define CLOG(n)          CLOGREQUEST(n,0).CLOG_FILE_LINE
#define CTAG(n,TAGS)     CTAGREQUEST(n,TAGS,0).CLOG_FILE_LINE_TAGS(TAGS)
#define CLOGERR          CLOGREQUEST(ERROR_LOG_LEVEL,0).CLOG_FILE_LINE_TAGS(ERROR_TAG)
#define CLOGOUT          CLOGREQUEST(ERROR_LOG_LEVEL,1).CLOG_FILE_LINE
#define CL               CLOGSTREAM
#define CLOGGING(n)      (CLOGREQUEST(n,0).logging() ? CLOGSTREAM.CLOG_FILE_LINE.logging() : 0)
#define CTAGGING(n,TAGS) (CTAGREQUEST(n,TAGS,0).logging() ? CLOGSTREAM.CLOG_FILE_LINE_TAGS(TAGS).logging() : 0)

// max buffer sizes
#define TBUF_MAX_SZ 60
#define FILETAG_MAX_SZ 400

// helpful exception stuff

#define THROW              CLOG(7) << "Exception thrown\n", throw
#define THROWSTR           "", THROW Standard_exception (ERRSTR.c_str())
#define THROWEXPR(EXPR)    { ERRSTR.clear(); ERRSTR << EXPR; THROW Standard_exception (ERRSTR.c_str()); }
#define THROWIFTRUE(EXPR)  { if (EXPR) { ERRSTR.clear(); ERRSTR << #EXPR; THROW Standard_exception (ERRSTR.c_str()); } }
#define THROWASSERT(EXPR)  { if (!(EXPR)) { ERRSTR.clear(); ERRSTR << "Failed assertion '" << #EXPR << "'"; THROW Standard_exception (ERRSTR.c_str()); } }

// Logfile system variables
class Log_streambuf_vars
{
protected:
  bool     logfile_open;
  sstring  last_filename;
  bool     bare_newline;
  sstring  file_tag;
  sstring  last_file_tag;

public:
  bool     verbose_on_stderr;
  bool     log_to_stderr;
  bool     log_to_stdout;
  bool     log_to_file;
  bool     log_in_XML;
};

// streambuf class
class Log_streambuf : public streambuf, public Log_streambuf_vars
{
 protected:
  int overflow(int c);

 private:

  ofstream logfile_stream;
  char*    datebuf;  // date string buffer
  char*    timebuf;  // time string buffer
  sstring  newtag_buf;

  stack<Log_streambuf_vars> saved_logfile_state;

  inline void prepare_date_and_time()
  {
    const time_t t = time(NULL);
    strftime (datebuf, TBUF_MAX_SZ, "%Y/%d/%m", localtime(&t));
    strftime (timebuf, TBUF_MAX_SZ, "%H:%M:%S", localtime(&t));
  }
  
 public:

  void save_logfile_state();
  void restore_logfile_state();

  inline void new_file_tag (const char* s)
    {
      // print a newline
      if (!bare_newline)
	{
	  if (log_to_file && logfile_open) logfile_stream << '\n';
	  if (log_to_stderr) cerr << '\n';
	  if (log_to_stdout) cout << '\n';
	}
      // update vars
      file_tag = s;
      bare_newline = TRUE;
    }

  void print_long_XML_tag();
  void print_short_XML_tag();
  void print_tabbed_tag();
  void print_newtag_buf();

  inline void close_logfile()
    {
      if (logfile_open)
	logfile_stream.close();
      logfile_open = FALSE;
    }

  inline void open_logfile(const sstring& filename)
    {
      if (filename != last_filename)
	{
	  close_logfile();
	  logfile_stream.open (filename.c_str(), ios_base::app);
	  logfile_open = 1;
	  last_filename = filename;
	}
    }

  Log_streambuf();
  ~Log_streambuf();
};

// vars for Log_stream
struct Log_stream_vars
{
  bool log_to_stderr_mask;
  bool log_to_file_mask;
  int  lowest_log_level;
  int  current_log_level;
};

// ostream class
class Log_stream : public ostream, public Log_stream_vars
{
 private:
  Log_streambuf& log_streambuf;

  char* file_tag_buf;

 public:

  stack<Log_stream_vars> saved_logfile_state;

  struct Log_directive
  {
    sstring dirstr_match;
    Regexp  filename_pattern;
    int     start_line;
    int     end_line;
    int     lowest_level;
    int     highest_level;
    int     stderr_or_logfile;   // 1 means directive applies to stderr, 2 means it applies to logfile, 3 means both
    bool    on;                  // true if this directive turns logging on & not off
    bool    process (int level, const char* src_file, int src_line, int& stderr_ret, int& file_ret);
    Log_directive (const sstring& directive_string);

    static bool syntax_help (Opts_list* ol);
    static bool tags_help (Opts_list* ol);
    static bool tags_long_help (Opts_list* ol);
  };
  
  list<Log_directive> directives;

  void add_new_directive (const char* dir);       // can be called from the debugger

  void save_logfile_state();
  void restore_logfile_state();
  
  void log_to_stderr_only() { log_to_stderr_mask = 1; log_to_file_mask = 0; }
  void log_to_file_only() { log_to_stderr_mask = 0; log_to_file_mask = 1; }
  void log_everywhere() { log_to_stderr_mask = log_to_file_mask = 1; }

  inline Log_stream& request_logging (int level, const char* src_file, int src_line, bool log_to_stdout = 0)
    {
      if (directives.empty())
	{
	  log_streambuf.log_to_stderr = log_to_stderr_mask && level >= lowest_log_level && !log_to_stdout;
	  log_streambuf.log_to_file   = log_to_file_mask   && level >= lowest_log_level;
	  log_streambuf.log_to_stdout = log_to_stdout;
	}
      else
	{
	  src_file = strip_path (src_file);
	  int stderr_flag = -1;
	  int file_flag = -1;
	  for_iterator (list<Log_directive>::reverse_iterator, dir, directives.rbegin(), directives.rend())     // loop backwards through directives list
	    {
	      dir->process (level, src_file, src_line, stderr_flag, file_flag);
	      if (stderr_flag != -1 && file_flag != -1)
		{
		  if (dir->dirstr_match.size())
		    {
		      for_const_contents (sstring, dir->dirstr_match, c)
			if (*c == ' ' || *c == '\t')  // chop off whitespace
			  break;
		    }
		  break;
		}
	    }
	  log_streambuf.log_to_stderr = log_to_stderr_mask && (stderr_flag == -1 ? level >= lowest_log_level : stderr_flag) && !log_to_stdout;
	  log_streambuf.log_to_file   = log_to_file_mask   && (file_flag == -1   ? level >= lowest_log_level : file_flag);
	  log_streambuf.log_to_stdout = log_to_stdout;
	}
      current_log_level = level;
      return *this;
    }
  
  inline bool logging() { return log_streambuf.log_to_stderr || log_streambuf.log_to_file || log_streambuf.log_to_stdout; }

  static inline const char* strip_path (const char* file_with_path)
    {
      const char* s;
      while ((s = strchr (file_with_path, DIR_SEP_CHAR)) != NULL) file_with_path = s + 1;
      return file_with_path;
    }

  inline Log_stream& file_tag (const char* src_file, int src_line, const char* tags)
    {
      if (log_streambuf.log_in_XML)
	{
	  if (tags ? (strlen(tags) > 0) : false)
	    sprintf (file_tag_buf, "file=\"%s\" line=%d level=%d tags=\"%s\"", strip_path (src_file), src_line, current_log_level, tags);
	  else
	    sprintf (file_tag_buf, "file=\"%s\" line=%d level=%d", strip_path (src_file), src_line, current_log_level);
	}
      else
	{
	  if (tags ? (strlen(tags) > 0) : false)
	    sprintf (file_tag_buf, "<%s#%d> <%d %s>", strip_path (src_file), src_line, current_log_level, tags);
	  else
	    sprintf (file_tag_buf, "<%s#%d> <%d>", strip_path (src_file), src_line, current_log_level);
	}
      log_streambuf.new_file_tag (file_tag_buf);
      return *this;
    }

  // constructor
  Log_stream (Log_streambuf& log_streambuf, int log_level);

  // destructor
  ~Log_stream();

  static bool clog_to_file (Opts_list* ol);
  static bool clog_to_stderr (Opts_list* ol);
  static bool clog_everywhere (Opts_list* ol);
  static bool clog_directive (Opts_list* ol);
  static bool clog_directive (const sstring& directive_string);
  static bool clog_negated_directive (Opts_list* ol);

  static void add_opts (Opts_list& ol);
};

extern Log_streambuf clog_streambuf;
extern Log_stream    clog_stream;

// the following string does nothing special except serve as a buffer for exception-throwing via the THROWSTR macro:
extern sstring ERRSTR;

#endif
