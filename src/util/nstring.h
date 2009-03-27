#ifndef NSTRING_INCLUDED
#define NSTRING_INCLUDED

// numeric-to-sstring adaptors
// these classes are very much obsoleted by the sstring class which allows ostream-like '<<' operators.
// (although these don't seem to work sometimes, so maybe nstring's are still useful...)

#include <stdlib.h>
#include "util/sstring.h"

class int_string : public sstring
{
 public:
  int_string (int n) : sstring()
    {
      char buf[100];
      sprintf (buf, "%d", n);
      append (buf);
    }
  int_string (const char* s, int n) : sstring()
    {
      char* buf = new char[strlen(s) + 100];
      sprintf (buf, s, n);
      append (buf);
      delete[] buf;
    }
};

class double_string : public sstring
{
 public:
  double_string (double n) : sstring()
    {
      char buf[100];
      sprintf (buf, "%g", n);
      append (buf);
    }
  double_string (const char* s, double n) : sstring()
    {
      char* buf = new char[strlen(s) + 100];
      sprintf (buf, s, n);
      append (buf);
      delete[] buf;
    }
};

#endif
