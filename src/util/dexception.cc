#include <stdio.h>
#include "util/logfile.h"

// uncomment the following #define to allow debugging display of the function stack
// (finally decided against leaving this as the default behaviour: it garbled the error message too many times - ihh, 8/15/03)
// #define __SHOW_FUNCTION_STACK__

Dart_exception::Dart_exception()
{
  stack_trace.clear();  // dummy, redundant line of code just so we have some place to set a breakpoint
#ifdef __SHOW_FUNCTION_STACK__
#ifdef __GNUC__
  // this code is VERY messy and doesn't work right, but a stack trace is a nice thing to have, even a half-working one
  void* ret_addr[20];
  int d = 0;
  void* fa;
  fa = __builtin_frame_address(0);
  if (fa) { ret_addr[d++] = __builtin_return_address(0); fa = __builtin_frame_address(1); }
  if (fa) { ret_addr[d++] = __builtin_return_address(1); fa = __builtin_frame_address(2); }
  if (fa) { ret_addr[d++] = __builtin_return_address(2); fa = __builtin_frame_address(3); }
  if (fa) { ret_addr[d++] = __builtin_return_address(3); fa = __builtin_frame_address(4); }
  if (fa) { ret_addr[d++] = __builtin_return_address(4); fa = __builtin_frame_address(5); }
  if (fa) { ret_addr[d++] = __builtin_return_address(5); fa = __builtin_frame_address(6); }
  if (fa) { ret_addr[d++] = __builtin_return_address(6); fa = __builtin_frame_address(7); }
  if (fa) { ret_addr[d++] = __builtin_return_address(7); fa = __builtin_frame_address(8); }
  if (fa) { ret_addr[d++] = __builtin_return_address(8); fa = __builtin_frame_address(9); }
  if (fa) { ret_addr[d++] = __builtin_return_address(9); fa = __builtin_frame_address(10); }
  if (fa) { ret_addr[d++] = __builtin_return_address(10); fa = __builtin_frame_address(11); }
  if (fa) { ret_addr[d++] = __builtin_return_address(11); fa = __builtin_frame_address(12); }
  if (fa) { ret_addr[d++] = __builtin_return_address(12); fa = __builtin_frame_address(13); }
  if (fa) { ret_addr[d++] = __builtin_return_address(13); fa = __builtin_frame_address(14); }
  if (fa) { ret_addr[d++] = __builtin_return_address(14); fa = __builtin_frame_address(15); }
  if (fa) { ret_addr[d++] = __builtin_return_address(15); fa = __builtin_frame_address(16); }
  if (fa) { ret_addr[d++] = __builtin_return_address(16); fa = __builtin_frame_address(17); }
  if (fa) { ret_addr[d++] = __builtin_return_address(17); fa = __builtin_frame_address(18); }
  if (fa) { ret_addr[d++] = __builtin_return_address(18); fa = __builtin_frame_address(19); }
  if (fa) ret_addr[d++] = __builtin_return_address(19);
  // the following line was added to work around mysterious, intermittent corruption of stack_trace...
  // this now seems to have disappeared, but it did that once before and came back... perhaps the code is haunted?
  stack_trace << "*** exception *** exception *** exception *** exception *** exception *** exception *** exception *** exception ***\n";
  stack_trace << "Function stack (depth=" << d;
  if (d == 20) stack_trace << '+';
  stack_trace << "):";
  for (int level = 0; level < d; ++level)
    {
      char hex[100];   // i assume this will never run on a >400-bit machine
      sprintf (hex, "%.8x", ret_addr[level]);
      stack_trace << " *0x" << hex;
    }
  if (d == 20) stack_trace << " ...";
  stack_trace << '\n';
  stack_trace << "*** exception *** exception *** exception *** exception *** exception *** exception *** exception *** exception ***\n";
  stack_trace << '\n';
#endif  /* __GNUC__ */  
#endif  /* __SHOW_FUNCTION_STACK__ */
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
