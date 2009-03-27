#include "util/strsaver.h"
#include "util/logfile.h"

stack<Stream_saver::FlagType> Stream_saver::old_stream_flags;

void Stream_saver::save_flags (ostream& o)
{
  FlagType old_flags = o.flags();
  old_stream_flags.push (old_flags);
}

void Stream_saver::left_align (ostream& o)
{
  FlagType f = o.flags();
  f = (f & ~ios::right) | ios::left;
  o.flags (f);
}

void Stream_saver::right_align (ostream& o)
{
  FlagType f = o.flags();
  f = (f & ~ios::left) | ios::right;
  o.flags (f);
}

void Stream_saver::restore_flags (ostream& o)
{
  if (old_stream_flags.size() == 0)
    THROWEXPR ("Tried to restore stream flags from an empty stack");
  o.flags (old_stream_flags.top());
  old_stream_flags.pop();
}
