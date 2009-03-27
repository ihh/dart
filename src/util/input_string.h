#ifndef INPUT_STRING_INCLUDED
#define INPUT_STRING_INCLUDED

#include <iostream>
#include "util/sstring.h"

class Input_string : public sstring
{
 public:
  sstring terminator_chars;

  Input_string (const char* terminator_chars = "\n") :  sstring(), terminator_chars(terminator_chars) {}
  Input_string (const sstring& terminator_chars) :  sstring(), terminator_chars(terminator_chars) {}

  Input_string& read (istream& in)
    {
      char c;
      while (!in.eof())
	{
	  in.get(c);
	  if (terminator_chars.find(c) != terminator_chars.npos) break;
	  push_back(c);
	}
      return *this;
    }

  friend istream& operator>> (istream& in, Input_string& s) { s.read (in); return in; }
};

#endif

