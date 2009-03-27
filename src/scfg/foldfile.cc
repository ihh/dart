#include <fstream>
#include "scfg/foldfile.h"

Fold_file::Fold_file (const char* filename, Sequence_database& seq_db) : np_ptr (0)
{
  ifstream infile (filename);
  if (!infile) THROWEXPR ("Couldn't open fold file '" << filename << "'");
  sstring instr;
  instr.getline (infile);
  vector<sstring> f = instr.split();
  if (f.size() != 2) THROWEXPR ("Bad first line in fold file");
  seq_db.push_back (Named_profile());
  Named_profile& the_np = seq_db.back();
  the_np.name = f[0];
  the_np.seq = f[1];
  np_ptr = &the_np;
  while (1)
    {
      instr.getline (infile);
      f = instr.split();
      if (f.size() == 0) break;
      if (f.size() != 2) THROWEXPR ("Bad line in fold file");
      fold_label.push_back (f[0]);
      fold_string.push_back (f[1]);
    }
}

Named_profile& Fold_file::np()
{
  return *np_ptr;
}
