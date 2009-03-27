#ifndef LOCAL_INCLUDED
#define LOCAL_INCLUDED

#include "seq/alignment.h"

struct Local_alignment : Alignment
{
  vector<int> start_coord;

  Local_alignment (int rows = 0) : Alignment (rows), start_coord (rows, (int) 0) { }
  Local_alignment (const Alignment_path& path) : Alignment (path), start_coord (path.rows(), (int) 0) { }
  Local_alignment (const Alignment& align, const vector<int>& row_set, bool collapse_empty_columns) : Alignment (align, row_set, collapse_empty_columns), start_coord (row_set.size(), (int) 0) { }

  void reset (int rows = 0) { Alignment::reset(rows); start_coord = vector<int> (rows, (int) 0); }
  void clear() { reset(0); }

  int  end_coord (int row) const;

  void read_MUL (istream& in, Sequence_database& db);                                                // discards co-ordinates
  void write_MUL (ostream& out, const Alphabet& alphabet, const sstring* column_labels = 0) const;   // adds co-ordinates
};

#endif

