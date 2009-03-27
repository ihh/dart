#ifndef NULL_STREAM_INCLUDED
#define NULL_STREAM_INCLUDED

#include <iostream>

using namespace std;

class Null_streambuf : public streambuf
{
protected:
  int overflow (int c);
};

extern Null_streambuf singleton_null_streambuf;

class Null_ostream : public ostream
{
public:
  Null_ostream()
    : ostream (&singleton_null_streambuf)
  { }
};

#endif /* NULL_STREAM_INCLUDED */
