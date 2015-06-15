#ifndef PIPER_INCLUDED
#define PIPER_INCLUDED

#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <iostream>

using namespace std;

// #include's for C++ stream wrappers for C files
#if defined(__GNUC__) && !defined(__clang__)
 #include <ext/stdio_filebuf.h>
 typedef __gnu_cxx::stdio_filebuf<char> Stdio_filebuf;
 using namespace __gnu_cxx;
#else
 // clang doesn't have stdio_filebuf, so create a dummy object and throw runtime exceptions if this functionality is used (EWW: DISGUSTING HACK)
 struct Stdio_filebuf : basic_streambuf<char> {
   Stdio_filebuf (int __fd, std::ios_base::openmode __mode) { }
   Stdio_filebuf (int __fd, std::ios_base::openmode __mode, size_t __size) { }
   void close() { }
 };
#endif
#define PIPE_BUF_SIZE 100   /* number of bytes for pipe buffer */

// Piper class handles interprocess pipes between parent & forked children
struct Piper
{
  // data
  int max_fork;  // number of child processes to run
  vector<vector<int> > fd;  // file descriptors for pipes
  vector<pid_t> pid;  // pid's for child processes

  // constructor, initialiser
  Piper (int max_fork);
  void init (int new_max_fork);

  // pipe & fork methods
  void open_pipes();
  void close_pipes() const;
  int read_fd (int i) const { return fd[i][0]; }  // pipe file descriptor for read
  int write_fd (int i) const { return fd[i][1]; }  // pipe file descriptor for write
  bool fork_child (int n) { return !(pid[n] = fork()); }
  void wait_for_child (int n) const { waitpid (pid[n], NULL, 0); }
};

// C++ stream wrappers for C file descriptors (e.g. for pipes)

struct istream_cfile
{
  // data
  Stdio_filebuf filebuf;
  istream in;
  // constructor
  istream_cfile (int fd);
};

struct ostream_cfile
{
  // data
  Stdio_filebuf filebuf;
  ostream out;
  // constructor
  ostream_cfile (int fd);
  // close method
  void close();
};

#endif /* PIPER_INCLUDED */
