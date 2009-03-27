#ifndef TREE_ALIGNMENT_LIKELIHOOD
#define TREE_ALIGNMENT_LIKELIHOOD

#include "util/score.h"
#include "tree/tree_alignment.h"

// Tree_alignment_likelihood abstract base class
struct Tree_alignment_likelihood
{
  // override loglike() with an alignment likelihood function
  virtual Loge loglike (const Tree_alignment& tree_align) = 0;
  // virtual destructor
  virtual ~Tree_alignment_likelihood() { }
};

// Executable Tree_alignment_likelihood class
// forks off a separate process and exec's a specified command
struct Likelihood_executable : Tree_alignment_likelihood
{
  // data
  static sstring shell_executable, shell_command_switch;  // "/bin/sh", "-c"
  sstring shell_command;    // likelihood command (takes alignment on standard input)
  const Alphabet& alphabet;
  // constructor
  Likelihood_executable (const Alphabet& alph, const char* command);
  // loglike method
  Loge loglike (const Tree_alignment& tree_align);
  // format help methods
  static void add_help (Opts_list* ol);
  static bool display_help (Opts_list* ol);
};


#endif /* TREE_ALIGNMENT_LIKELIHOOD */
