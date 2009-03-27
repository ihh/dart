#include "stemloc/slmain.h"

// Main program: creates & runs a Stemloc object from cmdline args
int main (int argc, char** argv)
{
  Stemloc stemloc (argc, argv);
  return stemloc.run();
}
