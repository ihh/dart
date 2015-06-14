#ifndef STREAM_SAVER_INCLUDED
#define STREAM_SAVER_INCLUDED

#include <iostream>
#include <stack>

using namespace std;

struct Stream_saver
{
  typedef iostream::fmtflags FlagType;

  static stack<FlagType> old_stream_flags;

  static void save_flags (ostream& o);
  static void left_align (ostream& o);
  static void right_align (ostream& o);
  static void restore_flags (ostream& o);
};

#endif /* STREAM_SAVER_INCLUDED */
