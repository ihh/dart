#include <algorithm>
#include "util/sstring.h"
#include "util/macros.h"
#include "util/Regexp.h"
#include "util/logfile.h"
#include "util/math_fn.h"

// minimum buffer size for sstring::getline()
#define DART_GETLINE_BUF_SIZE 100

int string_with_streambuf::DART_string_streambuf::overflow(int c) { owner.push_back ((char) c); return c; }

void sstring::to_lower()
{
  for_contents (sstring, *this, c)
    if (*c >= 'A' && *c <= 'Z') *c += 'a' - 'A';
}

void sstring::to_upper()
{
  for_contents (sstring, *this, c)
    if (*c >= 'a' && *c <= 'z') *c += 'A' - 'a';
}

vector<sstring> sstring::split (const char* split_chars, bool skip_empty_fields, int max_fields) const
{
  sstring expr;
  expr << "^([^" << split_chars << "]*)([" << split_chars << "])(.*)";
  Regexp parser (expr.c_str());
  
  vector<sstring> result;
  sstring tmp = *this;
  while (--max_fields != 0 && parser.Match (tmp.c_str()))
    {
      if (parser[1].size() || !skip_empty_fields) result.push_back (parser[1]);
      tmp = parser[3];
    }
  if (tmp.size() || !skip_empty_fields) result.push_back(tmp);
  
  return result;
}

sstring sstring::join (const vector<sstring>& v, const char* sep)
{
  sstring result;
  for (int i = 0; i+1 < (int) v.size(); i++) result.append(v[i]).append(sep);
  if (v.size()) result.append(v.back());
  return result;
}

sstring& sstring::getline (istream& is, size_type max_size)
{
  clear();
  size_type buf_sz = max ((size_t) capacity(), (size_t) DART_GETLINE_BUF_SIZE);
  const char term = '\n';
  while (1)
    {
      // read from the stream, until we reach the delimiter character, or fill the buffer.
      char buf [buf_sz];
      is.get (buf, buf_sz, term);
      append (buf);
      if (is.eof())
	break;
      is.clear();
      if (is.fail()) THROWEXPR ("Couldn't clear fail bit");
      // read next character. If it's the delimiter, stop; otherwise, enlarge the buffer.
      char c = '\0';
      is.get (c);
      push_back (c);
      if (is.fail()) THROW Format_exception (is, "Read failure\n");
      if (is.eof())
	break;
      if (c == term)
	{
	  is.clear();
	  break;
	}
      if (buf_sz >= max_size) THROW Format_exception (is, "Line too long\n");
      buf_sz = min (buf_sz * 2, max_size);
    }
  return *this;
}

Regexp int_regexp ("^\\-?[0-9]+$");
int sstring::to_int_strict (const char* err_prefix) const
{
  if (!int_regexp.Match (c_str()))
    THROWEXPR (err_prefix << ": " << *this);
  return to_int();
}

Regexp nonneg_int_regexp ("^[0-9]+$");
int sstring::to_nonneg_int_strict (const char* err_prefix) const
{
  if (!nonneg_int_regexp.Match (c_str()))
    THROWEXPR (err_prefix << ": " << *this);
  return to_int();
}

Regexp double_regexp ("^\\-?([0-9]*\\.)?[0-9]+([eE][\\+\\-]?[0-9]+)?$");
double sstring::to_double_strict (const char* err_prefix) const
{
  if (!double_regexp.Match (c_str()))
    THROWEXPR (err_prefix << ": " << *this);
  return to_double();
}

Regexp nonneg_double_regexp ("^([0-9]*\\.)?[0-9]+([eE][\\+\\-]?[0-9]+)?$");
double sstring::to_nonneg_double_strict (const char* err_prefix) const
{
  if (!nonneg_double_regexp.Match (c_str()))
    THROWEXPR (err_prefix << ": " << *this);
  return to_double();
}

sstring sstring::substr (int start, int len) const
{
  const int real_start = minmax (start, 0, (int) size());
  const int real_end = minmax (start + len, real_start, (int) size());
  sstring sub (begin() + real_start, begin() + real_end);
  return sub;
}

char_string::char_string (char c) : sstring()
{
  add_char(c);
}


char_string::char_string (const char* s) : sstring()
{
  reserve (strlen (s));
  for (; *s != '\0'; ++s) add_char(*s);
}

void char_string::add_char (char c)
{
  if (isprint(c))
    *this << c;
  else
    {
      char cstr[12];
      sprintf (cstr, "\\%.3o", c);
      *this << cstr;
    }
}

void sstring::read_entire_file (const char* filename)
{
  ifstream sstring_stream (filename);
  if (!sstring_stream) THROWEXPR ("File '" << filename << "' not found");
  std::copy (istream_iterator<char> (sstring_stream), istream_iterator<char>(), back_inserter (*this));
}
