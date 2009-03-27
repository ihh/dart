#ifndef FOLDFILE_INCLUDED
#define FOLDFILE_INCLUDED

#include "seq/biosequence.h"

class Fold_file
{
private:
  Named_profile* np_ptr;
public:
  // data
  vector<sstring> fold_string;
  vector<sstring> fold_label;
  // constructor
  Fold_file (const char* filename, Sequence_database& seq_db);
  // accessors
  Named_profile& np();
};

#endif
