#ifndef RECORDER_INCLUDED
#define RECORDER_INCLUDED

#include "handel/movement.h"

// general class to save sequentially increasing numbered files in a directory
struct Numbered_file_opener
{
  // data
  sstring directory, prefix, suffix;  // filename info
  int n_files;  // count of files saved to this directory
  // method to test whether recorder is "on"
  bool active() const { return directory.size() > 0; }
  // constructor
  Numbered_file_opener() : n_files (0) { }
  // method to open a file
  ofstream* open_file();
};

// class for losers (ECM, 6/26/'06)
// records graph representations of transducers to GraphViz dot-format files in a given directory
struct Transducer_dotfile_recorder : Numbered_file_opener
{
  // constructor
  Transducer_dotfile_recorder();
  // dotfile save method
  void save (const Transducer<Score>& trans, const char* comment = 0);
};

// class for recording S-expression
struct Composition_recorder : Numbered_file_opener
{
  // constructor
  Composition_recorder();
  // composition save method
  void save (Handel_movement& move);
};

#endif /* RECORDER_INCLUDED */
