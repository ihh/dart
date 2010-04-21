// Commented out as it doesn't work on OSX 10.4 - IH 4/20/2010
// #include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>

#include "util/logfile.h"

#define MAX_BACKTRACE_DEPTH 10
Dart_exception::Dart_exception()
{
  stack_trace.clear();
  // Commented out as it doesn't work on OSX 10.4 - IH 4/20/2010
  /*
  void *array[MAX_BACKTRACE_DEPTH];
  size_t size = backtrace (array, MAX_BACKTRACE_DEPTH);
  char** strings = backtrace_symbols (array, size);

  for (size_t i = 0; i < size; i++)
    stack_trace << strings[i] << '\n';

  free (strings);
  */
}

Dart_exception::~Dart_exception() { }

const char* Dart_exception::what() const
{
  sstring& what = (sstring&) _what;  // cast away const
  what.clear();
  what << stack_trace;
  what << details();
  return what.c_str();
}

const char* Dart_exception::details() const
{
  return "unknown exception\n";
}

Standard_exception::Standard_exception (const char* m) : Dart_exception(), msg(m)
{
  msg.append("\n");
}

Standard_exception::Standard_exception (const sstring& m) : Dart_exception(), msg(m)
{
  msg.append("\n");
}

const char* Standard_exception::details() const
{
  return msg.c_str();
}

Format_exception::Format_exception (istream& in) : Dart_exception() { info = ""; setup_info (in); }
Format_exception::Format_exception (istream& in, const char* prefix) : Dart_exception() { info << prefix << '\n'; setup_info (in); }
Format_exception::Format_exception (istream& in, const sstring& prefix) : Dart_exception() { info << prefix << '\n'; setup_info (in); }

void Format_exception::setup_info (istream& in)
{
  info << "Bad input format, somewhere before ";
  if (in.eof()) info.append("EOF\n");
  else
    {
      info.append("'");
      char c;
      for (int i = 0; i < 20 && !in.eof(); i++)
	if ((c = in.get()) != '\n') { if (!in.eof()) info.push_back(c); } else info.append("\\n");
      info.append("'\n");
    }
}

const char* Format_exception::details() const
{
  return info.c_str();
}

const char* String_exception::details() const
{
  return msg.c_str();
}
