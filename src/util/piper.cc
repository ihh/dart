#include <errno.h>
#include "util/piper.h"
#include "util/logfile.h"

#if (defined __GNUC__ && (__GNUC__ < 3 || (__GNUC__ == 3 && defined __GNUC_MINOR__ && __GNUC_MINOR__ < 4)))
#undef FILEBUF_TAKES_3_ARGS
#else
#define FILEBUF_TAKES_3_ARGS
#endif

istream_cfile::istream_cfile (int fd) :
#ifdef FILEBUF_TAKES_3_ARGS
  filebuf (fd, ios_base::in, PIPE_BUF_SIZE),
#else
  filebuf (fd, ios_base::in, false, PIPE_BUF_SIZE),
#endif
  in (&filebuf)
{ }

ostream_cfile::ostream_cfile (int fd) :
#ifdef FILEBUF_TAKES_3_ARGS
  filebuf (fd, ios_base::out, PIPE_BUF_SIZE),
#else
  filebuf (fd, ios_base::out, false, PIPE_BUF_SIZE),
#endif
  out (&filebuf)
{ }

void ostream_cfile::close()
{
  filebuf.close();
}

Piper::Piper (int max_fork)
{
  init (max_fork);
}

void Piper::init (int new_max_fork)
{
  max_fork = new_max_fork;
  fd = vector<vector<int> > (max_fork, vector<int> (2));
  pid = vector<pid_t> (max_fork);
}

void Piper::open_pipes()
{
  for (int j = 0; j < max(max_fork,1); ++j)
    if (pipe (&fd[j][0]) == -1)
      THROWEXPR ("Couldn't open pipe: error code " << errno);
}

void Piper::close_pipes() const
{
  for (int j = 0; j < max(max_fork,1); ++j)
    {
      close (fd[j][0]);
      close (fd[j][1]);
    }
}
