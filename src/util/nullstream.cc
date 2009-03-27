#include "util/nullstream.h"

int Null_streambuf::overflow (int c)
{
  return c;
}

Null_streambuf singleton_null_streambuf;
