#include "util/vector_output.h"
#include "util/logfile.h"
#include "util/map_keys.h"
#include "util/logtags.h"

Log_streambuf clog_streambuf;
Log_stream    clog_stream (clog_streambuf, 9);
sstring       ERRSTR;   /* used as a temporary variable by THROWEXPR() macro */

Log_streambuf::Log_streambuf()
  : streambuf()
{
  bare_newline = true;
  verbose_on_stderr = false;
  log_to_stderr = true;
  log_to_file = false;
  
  datebuf = new char[TBUF_MAX_SZ + 1];
  timebuf = new char[TBUF_MAX_SZ + 1];
}

Log_streambuf::~Log_streambuf()
{
  delete[] datebuf;
  delete[] timebuf;
}

int Log_streambuf::overflow (int c)
{
  // print tags
  if (bare_newline)
    {
      if (log_in_XML)
	if (last_file_tag != file_tag)
	  print_long_XML_tag();
	else
	  print_short_XML_tag();
      else
	print_tabbed_tag();
      last_file_tag = file_tag;
    }
  // log the character
  if (log_to_stdout) cout.put(c);
  if (log_to_stderr) cerr.put(c);
  if (log_to_file && logfile_open) logfile_stream.put(c);
  // set (or unset) bare_newline flag
  if (c == '\n')
    {
      bare_newline = TRUE;
      if (log_to_stdout) cout.flush();
      if (log_to_stderr) cerr.flush();
      if (log_to_file && logfile_open) logfile_stream.flush();
    }
  else
    bare_newline = FALSE;
  // return
  return 0;
}

void Log_streambuf::save_logfile_state()
{
  saved_logfile_state.push (Log_streambuf_vars (*this));
}

void Log_streambuf::restore_logfile_state()
{
  ((Log_streambuf_vars&) *this) = saved_logfile_state.top();
  saved_logfile_state.pop();
}

void Log_streambuf::print_long_XML_tag()
{
  prepare_date_and_time();
  newtag_buf.clear();
  newtag_buf << '<' << CLOGXMLTAGNAME;
  newtag_buf << " date=\"" << datebuf << "\" time=\"" << timebuf << "\" " << file_tag << "/>\n";
  print_newtag_buf();
}

void Log_streambuf::print_short_XML_tag()
{
  const sstring old_datebuf (datebuf);
  const sstring old_timebuf (timebuf);
  prepare_date_and_time();
  const bool date_changed = strcmp (datebuf, old_datebuf.c_str()) != 0;
  const bool time_changed = strcmp (timebuf, old_timebuf.c_str()) != 0;
  if (date_changed || time_changed)
    {
      newtag_buf.clear();
      newtag_buf << '<' << CLOGXMLTAGNAME;
      if (date_changed) newtag_buf << " date=\"" << datebuf << "\"";
      if (time_changed) newtag_buf << " time=\"" << timebuf << "\"";
      newtag_buf << "/>\n";
      print_newtag_buf();
    }
}

void Log_streambuf::print_tabbed_tag()
{
  prepare_date_and_time();
  newtag_buf.clear();
  newtag_buf << '[' << datebuf << ',' << timebuf << "] " << file_tag << "\t";
  print_newtag_buf();
}

void Log_streambuf::print_newtag_buf()
{
  if (log_to_file && logfile_open) logfile_stream << newtag_buf;
  if (log_to_stderr && verbose_on_stderr) cerr << newtag_buf;
}

Regexp stream_qualifier_re ("^(.*);([^;]*)$");
Regexp bad_re ("[!#,:;]");
Regexp level_range_re ("^(.*),([^,]*)$");
Regexp range_re ("^(-?[0123456789]+):(-?[0123456789]+)$");
Regexp above_re ("^(-?[0123456789]+)$");
Regexp below_re ("^:(-?[0123456789]+)$");
Regexp line_range_re ("^(.*)#([0123456789]+)-([0123456789]+)$");
Regexp single_line_re ("^(.*)#([0123456789]+)$");
Regexp pling_re ("^!(.*)$");
Regexp filename_re ("^[ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_]+\\.[ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]+$");
Regexp tag_re ("^[ABCDEFGHIJKLMNOPQRSTUVWXYZ][ABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789]+$");

Log_stream::Log_directive::Log_directive (const sstring& directive_string)
{
  sstring d = directive_string;

  // chop off the ';f' or the ';e'
  //
  if (stream_qualifier_re.Match (d.c_str()))
    {
      d = stream_qualifier_re[1];
      const sstring qualifier = stream_qualifier_re[2];
      if (qualifier.size() != 1) THROWEXPR ("Bad qualifier char: '" << qualifier << "' in '" << d << "'");
      const char q = qualifier[0];
      if (q != 'f' && q != 'e') { syntax_help(0); THROW Standard_exception ("Syntax error in logging directive"); }
      stderr_or_logfile = q == 'e' ? 1 : 2;
    }
  else
    stderr_or_logfile = 3;

 // chop off the ',x:y' log-level range
  //
  if (level_range_re.Match (d.c_str()))
    {
      d = level_range_re[1];
      const sstring level_range = level_range_re[2];
      if (range_re.Match (level_range.c_str())) { sstring l = range_re[1]; sstring h = range_re[2]; lowest_level = atoi(l.c_str()); highest_level = atoi(h.c_str()); }
      else if (above_re.Match (level_range.c_str())) { sstring l = above_re[1]; lowest_level = atoi(l.c_str()); highest_level = +123456789; }
      else if (below_re.Match (level_range.c_str())) { sstring h = below_re[1]; lowest_level = -123456789; highest_level = atoi(h.c_str()); }
      else { syntax_help(0); THROW Standard_exception ("Syntax error in logging directive"); }
    }
  else
    {
      lowest_level  = -123456789;
      highest_level = +123456789;
    }
  if (lowest_level > highest_level) { syntax_help(0); THROW Standard_exception ("Syntax error in logging directive"); }

  // chop off the '#x-y' line-number range
  //
  if (line_range_re.Match (d.c_str())) { d = line_range_re[1]; sstring s = line_range_re[2]; start_line = atoi(s.c_str()); sstring e = line_range_re[3]; end_line = atoi(e.c_str()); }
  else if (single_line_re.Match (d.c_str())) { d = single_line_re[1]; sstring s = single_line_re[2]; start_line = end_line = atoi(s.c_str()); }
  else { start_line = 0; end_line = +123456789; }
  if (start_line > end_line) { syntax_help(0); THROW Standard_exception ("Syntax error in logging directive"); }

  // strip the '!' off the front
  //
  if (pling_re.Match (d.c_str())) { d = pling_re[1]; on = 0; }
  else on = 1;

  // if any special characters remain, then something's wrong
  //
  if (bad_re.Match (d.c_str())) { syntax_help(0); THROW Standard_exception ("Syntax error in logging directive"); }

  // convert '^' and '$' to whitespace
  for_contents (sstring, d, c) if (*c == '^' || *c == '$') *c = ' ';

  // if the remaining pattern fits
  //      /^[A-Za-z_]+\.[A-Za-z]+$/ (has a dot suffix => is a filename)
  //   or /^[A-Z][A-Z_0-9]+$/       (upper case variable name => is a tag),
  // then add a space to the start and the end, to force a complete match
  // (or in the case of filenames, add a '^' to the start and a space to the end)
  //
  // also add parentheses at this stage
  //
  const char* dc = d.c_str();
  if (filename_re.Match(dc)) { d.insert ((sstring::size_type) 0, "^()"); d.append (" "); }
  else if (tag_re.Match(dc)) { d.insert ((sstring::size_type) 0, " ("); d.append (") "); }
  else                       { d.insert ((sstring::size_type) 0, "()"); }

  filename_pattern = Regexp (d.c_str());
}

bool Log_stream::Log_directive::syntax_help (Opts_list* ol)
{
  cerr << "\nLogging directives have a flexible syntax. The simplest kind of directive is a tag in upper-case,\nsuch as '-log INIT', or a priority level, e.g. '-log 6'; but you can also turn on (or off) messages\nby source file, by line, and by combinations of these.\n\n";
  cerr << "Examples include:\n\n";
  cerr << "-log JAM               turns on logging for all messages with tag 'JAM' (see -logtags or -logtaginfo for list)\n";
  cerr << "-log 3                 turns on logging for messages of priority level 3 or higher\n";
  cerr << "-log cabin,3           turns on logging for messages of priority level 3 or higher in source files matching '*cabin*' (eg cabin.h, mycabin.cc)\n";
  cerr << "-log cabin.h,3         turns on logging for messages of priority level 3 or higher in source file 'cabin.h' only\n";
  cerr << "-log 'lady,3:6'        turns on logging for messages from priority levels 3 to 6 in source files '*lady*'\n";
  cerr << "                        (NB excluding high log levels may have the side-effect of blocking lower ones due to nested 'if' clauses in C++ code)";
  cerr << "\n-log 'lady#1-100,3:6'  turns on logging for messages from levels 3 to 6 for first 100 lines of files '*lady*'";
  cerr << "\n-log 'lady,3:6;f'      turns on logging for messages from levels 3 to 6 in files '*lady*' for logfile only\n";
  cerr << "                        (suffix ';e' applies to standard error only)\n";
  cerr << "-log '!^fire'          turns off logging for all messages in source files starting with 'fire' (eg firefight.f77, fired.exe)\n";
  cerr << "-nolog fire            same as -log '!fire'\n\n";
  exit(0);
  return 0;
}

bool Log_stream::Log_directive::tags_help (Opts_list* ol)
{
  cout << all_log_tags;
  exit(0);
  return 1;
}

bool Log_stream::Log_directive::tags_long_help (Opts_list* ol)
{
  cout << all_log_tag_info;
  exit(0);
  return 1;
}

bool Log_stream::Log_directive::process (int level, const char* src_file, int src_line, int& stderr_ret, int& file_ret)
{
  if (filename_pattern.Match (src_file) && level >= lowest_level && level <= highest_level && src_line >= start_line && src_line <= end_line)
    {
      if ((stderr_or_logfile & 1) && stderr_ret == -1) stderr_ret = on ? 1 : 0;
      if ((stderr_or_logfile & 2) && file_ret == -1) file_ret = on ? 1 : 0;
      dirstr_match = filename_pattern[1];
      return stderr_or_logfile ? 1 : 0;
    }
  return 0;
}

void Log_stream::add_new_directive (const char* dir)
{
  sstring d (dir);
  directives.push_back (Log_directive (d));
}

void Log_stream::save_logfile_state()
{
  saved_logfile_state.push (Log_stream_vars (*this));
  log_streambuf.save_logfile_state();
}

void Log_stream::restore_logfile_state()
{
  ((Log_stream_vars&) *this) = saved_logfile_state.top();
  saved_logfile_state.pop();
  log_streambuf.restore_logfile_state();
}

Log_stream::Log_stream (Log_streambuf& log_streambuf, int log_level) :
  ostream (&log_streambuf),
  log_streambuf (log_streambuf)
{
  log_to_stderr_mask = true;
  log_to_file_mask = false;
  lowest_log_level = log_level;
  current_log_level = INITIAL_LOG_LEVEL;

  file_tag_buf = new char[FILETAG_MAX_SZ + 1];
  this->request_logging(0,__FILE__,__LINE__).CLOG_FILE_LINE << "Log opened\n";
}

Log_stream::~Log_stream()
{
  this->request_logging(0,__FILE__,__LINE__).CLOG_FILE_LINE << "Log closed\n";
  clog_streambuf.close_logfile();
  delete[] file_tag_buf;
}


bool Log_stream::clog_to_file (Opts_list* ol)
{
  sstring filename = ol->next_string();
  clog_streambuf.open_logfile (filename);
  clog_stream.log_to_file_only();
  CLOG(0) << "Opened logfile " << filename << "\n";
  return 1;
}

bool Log_stream::clog_to_stderr (Opts_list* ol)
{
  clog_stream.log_to_stderr_only();
  return 1;
}

bool Log_stream::clog_everywhere (Opts_list* ol)
{
  sstring filename = ol->next_string();
  clog_streambuf.open_logfile (filename);
  clog_stream.log_everywhere();
  CLOG(0) << "Opened logfile " << filename << "\0";
  return 1;
}

bool Log_stream::clog_directive (Opts_list* ol)
{
  sstring directive_string (ol->next_string());
  Regexp numeric_re ("^-?[0123456789]+$");
  if (numeric_re.Match (directive_string.c_str())) clog_stream.lowest_log_level = atoi (directive_string.c_str());
  else clog_stream.directives.push_back (Log_directive (directive_string));
  return 1;
}

bool Log_stream::clog_negated_directive (Opts_list* ol)
{
  sstring directive_string (ol->next_string());
  clog_stream.directives.push_back (Log_directive (directive_string));
  clog_stream.directives.back().on = !clog_stream.directives.back().on;
  return 1;
}

void Log_stream::add_opts (Opts_list& ol)
{
  ol.print_title ("Logging options");
  ol.add ("log", &Log_stream::clog_directive, " <string>\tturn on diagnostic logging (-loghelp shows syntax)");
  ol.add ("logfile", &Log_stream::clog_to_file, " <file>\tlog to file");
  ol.add ("logcopy", &Log_stream::clog_everywhere, " <file>\tlog to file and standard error");
  ol.add ("logtime", clog_streambuf.verbose_on_stderr = 0, "timestamp standard error (logfile stamped automatically)");
  ol.add ("logxml",  clog_streambuf.log_in_XML = TRUE, "add XML timestamps");
  ol.add ("logerr",  &Log_stream::clog_to_stderr, "\tlog on standard error (default)");
  ol.add ("nolog",   &Log_stream::clog_negated_directive);
  ol.add ("loghelp", &Log_stream::Log_directive::syntax_help);
  ol.add ("logtags", &Log_stream::Log_directive::tags_help);
  ol.add ("logtaginfo", &Log_stream::Log_directive::tags_long_help);
  ol.newline();
}
