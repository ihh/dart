#include <stdio.h>
#include <errno.h>
#include "util/piper.h"
#include "handel/talike.h"

#define MAX_FORK_ATTEMPTS 5

void Likelihood_executable::add_help (Opts_list* ol)
{
  ol->add ("xlh -execlike-help", &Likelihood_executable::display_help, "\thelp on alignment-likelihood executables");
}

bool Likelihood_executable::display_help (Opts_list* ol)
{
  cout << "An alignment-likelihood executable should read a Stockholm-format alignment+tree\n";
  cout << "from standard input, print a bit-score on standard output, and exit.\n";
  cout << "The program should thus evaluate the function log(P(A,T))/log(2),\n";
  cout << "where P(A,T) is the likelihood of the (Alignment,Tree) pair.\n";
  exit(0);
  return false;
}

sstring Likelihood_executable::shell_executable ("/bin/sh");  // TODO: use autoconf to set this in config.h
sstring Likelihood_executable::shell_command_switch ("-c");  // TODO: set this in config.h as well

Likelihood_executable::Likelihood_executable (const Alphabet& alph,
					      const char* command)
  : shell_command (command), alphabet (alph)
{ }

#define PIPE_READ_END  0
#define PIPE_WRITE_END 1
Loge Likelihood_executable::loglike (const Tree_alignment& tree_align)
{
  Log2 bit_score = 0.;

  // Create two pipes.
  // Alignment goes down one pipe (from parent to child),
  // score comes back up t'other (from child to parent).
  int fd_align[2], fd_score[2];
  if (pipe(fd_align) != 0 || pipe(fd_score) != 0)
    THROWEXPR ("Pipe creation error");
  CTAG(5,MCMC_PIPE) << "Pipe file numbers: alignment " << fd_align[PIPE_READ_END] << " (in), " << fd_align[PIPE_WRITE_END] << " (out); score " << fd_score[PIPE_READ_END] << " (in), " << fd_score[PIPE_WRITE_END] << " (out)\n";

  pid_t pid = -1;
  int fork_attempt = 1;
  while ((pid = fork()) == -1)   /* try a few times before giving up */
    if (fork_attempt++ >= MAX_FORK_ATTEMPTS)
      THROWEXPR ("Failed " << MAX_FORK_ATTEMPTS << " fork attempts");

  // Fork.
  if (pid == 0)
    {
      // Child process.
      // First close the pipe ends we don't need.
      close(fd_score[PIPE_READ_END]);
      close(fd_align[PIPE_WRITE_END]);

      // Tie stdout to score pipe and stdin to alignment pipe.
      if (dup2 (fd_score[PIPE_WRITE_END], fileno (stdout)) < 0)
	  THROWEXPR ("Unable to redirect child STDOUT to write end of score pipe");

      if (dup2 (fd_align[PIPE_READ_END], fileno (stdin)) < 0)
        THROWEXPR ("Unable to redirect child STDIN to read end of alignment pipe");

      // run shell command
      CTAG(3,MCMC_EXEC) << "Running \"" << shell_executable << ' ' << shell_command_switch << ' ' << shell_command << "\"\n";
      execlp (shell_executable.c_str(), shell_executable.c_str(), shell_command_switch.c_str(), shell_command.c_str(), (char*) 0);

      // should not get here
      THROWEXPR ("execlp() failed with errno " << errno);
    }
  else
    {
      // Parent process.
      CTAG(5,MCMC_FORK) << "Forked on attempt #" << fork_attempt << "; child pid " << pid << "\n";

      // First close the pipe ends we don't need.
      close (fd_score[PIPE_WRITE_END]);
      close (fd_align[PIPE_READ_END]);

      // Create C++ stream wrappers for the pipes.
      ostream_cfile align_stream (fd_align[PIPE_WRITE_END]);
      istream_cfile score_stream (fd_score[PIPE_READ_END]);

      // Write the alignment into the alignment pipe.
      // The child will run it through a shell command.
      tree_align.write_Stockholm (align_stream.out, alphabet);
      align_stream.close();
      close (fd_align[PIPE_WRITE_END]);  // close write end of pipe... hope this doesn't happen twice... hmm.... need to think about this

      // Read the output of the shell command from the score pipe.
      score_stream.in >> bit_score;
      close (fd_score[PIPE_READ_END]);  // close read end of pipe... hope this doesn't happen twice

      // wait for child to die
      CTAG(5,MCMC_WAIT) << "Waiting for child pid " << pid << "\n";
      int status = -1;
      waitpid (pid, &status, 0);
      if (status < 0)
	THROWEXPR ("waitpid() failed with error " << errno);
    }

  // return the score
  return Bits2Nats (bit_score);
}
