#include <functional>
#include <algorithm>

#include "util/macros.h"
#include "util/istream_get_iterator.h"
#include "seq/alignment.h"
#include "util/Regexp.h"
#include "util/logfile.h"

Alignment_path::Alignment_path (int rows, int cols)
  : path_data (rows, Row (cols, (bool) 0))
{ }

bool Alignment_path::operator== (const Alignment_path& a) const
{
  return path_data == a.path_data;
}

bool Alignment_path::operator!= (const Alignment_path& a) const
{
  return !(*this == a);
}

void Alignment_path::make_flush (int n_rows)
{
  if (rows() > n_rows) THROWEXPR ("Too many rows");
  while (rows() < n_rows) path_data.push_back (Row());
  make_flush();
}

void Alignment_path::make_flush()
{
  Row::size_type cmax = 0;
  for (int r = 0; r < rows(); r++) cmax = max (cmax, row(r).size());
  for (int r = 0; r < rows(); r++) row(r).insert (row(r).end(), cmax - row(r).size(), 0);
}

void Alignment_path::insert_columns (int col, int n)
{
  for (int i = 0; i < rows(); i++) row(i).insert (row(i).begin() + col, n, (bool) 0);
}

void Alignment_path::insert_columns (int col, const vector<bool>& col_data, int n)
{
  for (int i = 0; i < rows(); i++) row(i).insert (row(i).begin() + col, n, col_data[i]);
}

void Alignment_path::erase_columns (int col, int n)
{
  for (int i = 0; i < rows(); i++) row(i).erase (row(i).begin() + col, row(i).begin() + col + n);
}

bool Alignment_path::column_is_empty (int col) const
{
  bool is_empty = 1;
  for (int row = 0; row < rows(); row++)
    if (path_data[row][col]) { is_empty = 0; break; }
  return is_empty;
}

void Alignment_path::erase_empty_columns()
{
  int ncol = 0;
  for (int col = 0; col < columns(); col++)
    if (!column_is_empty (col))
      {
	if (ncol != col)
	  for (int row = 0; row < rows(); row++)
	    path_data[row][ncol] = path_data[row][col];
	++ncol;
      }
  erase_columns (ncol, columns() - ncol);
}

void Alignment_path::insert_rows (int row, int n)
{
  path_data.insert (path_data.begin() + row, n, create_empty_row());
}

void Alignment_path::insert_rows (int row, const Row& row_data, int n)
{
  path_data.insert (path_data.begin() + row, n, row_data);
}

void Alignment_path::erase_rows (int row, int n)
{
  path_data.erase (path_data.begin() + row, path_data.begin() + row + n);
}

void Alignment_path::erase_rows (const Row_index_set& row_set)
{
  vector<bool> keep_row (rows(), (bool) 1);
  for_const_contents (Row_index_set, row_set, r) keep_row[*r] = 0;
  int j = 0;
  for (int i = 0; i < rows(); i++)
    if (keep_row[j])
      {
	if (i != j) row(j).swap (row(i));
	j++;
      }
  erase_rows (j, rows() - j);
}

bool Alignment_path::is_ungapped() const
{
  for_const_contents (vector<Row>, path_data, row)
    for_const_contents (Row, *row, c)
    if (!*c)
      return false;
  return true;
}

Alignment_path::Sequence_coords Alignment_path::seq_coords_begin() const
{
  return create_seq_coords();
}

Alignment_path::Sequence_coords Alignment_path::seq_coords_end() const
{
  Sequence_coords seq_coords = create_seq_coords();
  for (int col = 0; col < columns(); ++col) inc_seq_coords (seq_coords, col);
  return seq_coords;
}

void Alignment_path::zero_seq_coords (Sequence_coords& seq_coords) const
{
  if ((int) seq_coords.size() != rows()) THROWEXPR ("Sequence_coords size doesn't match number of rows in alignment");
  for (int i = 0; i < rows(); i++) seq_coords[i] = 0;
}

void Alignment_path::inc_seq_coords (Sequence_coords& seq_coords, int col) const
{
  if ((int) seq_coords.size() != rows()) THROWEXPR ("Sequence_coords size doesn't match number of rows in alignment");
  for (int i = 0; i < rows(); i++) seq_coords[i] += path_data[i][col];
}

void Alignment_path::dec_seq_coords (Sequence_coords& seq_coords, int col) const
{
  if ((int) seq_coords.size() != rows()) THROWEXPR ("Sequence_coords size doesn't match number of rows in alignment");
  for (int i = 0; i < rows(); i++) seq_coords[i] -= path_data[i][col];
}

void Alignment_path::explode_row (int row, const Alignment_path& exploded_path)
{
  if (exploded_path.rows() == 0)
    erase_rows (row);
  else
    {
      path_data.insert (path_data.begin() + row + 1, exploded_path.rows() - 1, create_empty_row());
      int expl_col = 0;
      for (int col = 0; col < columns(); col++)
	if (path_data[row][col])
	  {
	    for (int expl_row = 0; expl_row < exploded_path.rows(); expl_row++)
	      path_data[row+expl_row][col] = exploded_path(expl_row,expl_col);
	    ++expl_col;
	  }
      if (expl_col != exploded_path.columns())
	THROW Standard_exception ("Number of steps in exploding row doesn't match number of columns in exploded path");
    }
}

void Alignment_path::explode (const vector<Row_index_set>& row_sets, const vector<Alignment_path>& exploded_paths)
{
  vector<int> new_row_pos;
  for (int r = 0; r < (int) row_sets.size(); r++)
    {
      if ((int) row_sets[r].size() != exploded_paths[r].rows())
	THROW Standard_exception ("Parameter set mismatch during alignment explosion");
      explode_row (new_row_pos.size(), exploded_paths[r]);
      for_const_contents (Row_index_set, row_sets[r], s) new_row_pos.push_back (*s);
    }
  vector<Row> new_path_data (path_data.size(), Row());
  for (int r = 0; r < (int) new_row_pos.size(); r++)
    new_path_data[new_row_pos[r]].swap (path_data[r]);
  path_data.swap (new_path_data);
}

vector<Alignment_path> Alignment_path::implode (const vector<Row_index_set>& row_sets)
{
  make_flush();

  vector<bool> accounted_for (rows(), (bool) 0);
  vector<Alignment_path> v;
  v.reserve (row_sets.size());
  
  for (int i = 0; i < (int) row_sets.size(); i++)
    {
      v.push_back (Subalignment_path (*this, row_sets[i], 1));
      for_const_contents (Row_index_set, row_sets[i], row)
	{
	  if (accounted_for[*row]) THROW Standard_exception ("Attempt to implode row into more than one set");
	  accounted_for[*row] = 1;
	}
    }
  for_contents (vector<bool>, accounted_for, a)
    if (!*a) THROW Standard_exception ("Must account for all rows when doing an alignment implosion");

  Alignment_path new_path (row_sets.size(), columns());
  for (int i = 0; i < (int) row_sets.size(); i++)
    for (int col = 0; col < columns(); col++)
      {
	bool nonempty = 0;
	for_const_contents (Row_index_set, row_sets[i], r) if (path_data[*r][col]) { nonempty = 1; break; }
	new_path[i][col] = nonempty;
      }
  
  path_data.swap (new_path.path_data);
  return v;
}


Alignment_path::Decomposition Alignment_path::decompose (const vector<Alignment_path::Row_pair>& pairs) const
{
  Decomposition decomp;
  for_const_contents (vector<Row_pair>, pairs, pair)
    decomp[*pair] = Pairwise_path (*this, (*pair).first, (*pair).second, 1);
  return decomp;
}

void Alignment_path::compose_and_log (const Decomposition& decomp, bool squash_columns, bool transducer_gap_ordering)
{
  if (CTAGGING(1,COMPOSE))
    {
      CL << "Composing an alignment path from the following set of pairwise paths:\n";
      for_const_contents (Decomposition, decomp, d)
	{
	  int parent = (*d).first.first, child = (*d).first.second;
	  vector<sstring> row_name (2);
	  row_name[0] << parent;
	  row_name[1] << child;
	  (*d).second.show (CL, row_name);
	  CL << "\n";
	}
    }
  compose (decomp, squash_columns, transducer_gap_ordering);
  if (CTAGGING(0,COMPOSE))
    {
      CL << "Composed alignment:\n";
      vector<sstring> row_name (rows());
      for (int r = 0; r < rows(); ++r) row_name[r] << r;
      show (CL, row_name);
    }
}

void Alignment_path::compose (const Decomposition& decomp, bool squash_columns, bool transducer_gap_ordering)
{
  int new_rows = 0;
  for_const_contents (Decomposition, decomp, d)
    new_rows = max (new_rows, 1 + max ((*d).first.first, (*d).first.second));

  vector<Row_index_set> children (new_rows, Row_index_set());
  vector<bool>          has_child  (new_rows, (bool) 0);
  vector<bool>          has_parent (new_rows, (bool) 0);
  Row_pair_cursor       cursor;
  Alignment_path        buffer (new_rows);

  // identify "root" nodes (i.e. unparented ones)
  //
  for_const_contents (Decomposition, decomp, d)
    {
      int parent = (*d).first.first, child = (*d).first.second;
      children[parent].push_back (child);
      has_child[parent] = has_parent[child] = 1;
      cursor[(*d).first] = -1;
    }

  // if transducer_gap_ordering is true, then sort all children[] vectors in reverse
  for_contents (vector<Row_index_set>, children, c)
    {
      sort (c->begin(), c->end());
      if (transducer_gap_ordering)
	reverse (c->begin(), c->end());
    }

  // advance through each "root" sequence
  //
  for (int r = 0; r < buffer.rows(); r++)
    if (has_child[r] && !has_parent[r])
      while (buffer.add_child_step_to_buffer (Row_pair(-1,r), decomp, children, squash_columns, cursor)) {}

  // check that every cursor is at the end of its path
  //
  for_const_contents (Decomposition, decomp, d)
    {
      int s;
      if ((s = (*d).second.columns() - cursor[(*d).first]) != 0)
	{
	  verify_decomposition (decomp);  // see if we can pin down the error
	  sstring e;
	  e << "During alignment composition, pairwise cursor for rows " << (*d).first.first << " and " << (*d).first.second << " was " << s << " columns short of the end of the row";
	  THROW Standard_exception (e);
	}
    }

  buffer.make_flush();
  path_data.swap (buffer.path_data);
}

void Alignment_path::verify_decomposition (const Decomposition& decomp)
{
  map<int,int> row_len;
  for_const_contents (Decomposition, decomp, di)
    {
      const Row_pair& row_pair = di->first;
      const Alignment_path& path = di->second;
      if (path.rows() != 2) THROWEXPR ("Bad path for row pair (" << row_pair.first << "," << row_pair.second << ")\n");

      const int row0 = row_pair.first;
      const int row1 = row_pair.second;

      const int row0_len = path.count_steps_in_row (0);
      const int row1_len = path.count_steps_in_row (1);

      map<int,int>::iterator i0 = row_len.find (row0);
      map<int,int>::iterator i1 = row_len.find (row1);

      if (i0 != row_len.end())
	{
	  if (i0->second != row0_len)
	    THROWEXPR ("Inconsistent lengths for row " << row0 << " (" << i0->second << ", " << row0_len << ")");
	}
      else
	row_len[row0] = row0_len;

      if (i1 != row_len.end())
	{
	  if (i1->second != row1_len)
	    THROWEXPR ("Inconsistent lengths for row " << row1 << " (" << i1->second << ", " << row1_len << ")");
	}
      else
	row_len[row1] = row1_len;
    }
}

bool Alignment_path::advance_cursor_to_next_parent_step (const Row_pair& pair, const Decomposition& decomp, const vector<Row_index_set>& children, bool squash_columns, Row_pair_cursor& cursor, Nascent_column& col_data)
{
  const int             child   = pair.second;
  const Alignment_path& ppath   = (*decomp.find(pair)).second;
  int&                  pcursor = cursor[pair];

  while (1)
    {
      //      CLOG(1) << "p" << parent << ":" << pcursor << " ";

      if (++pcursor == ppath.columns())                                                                   // is final column
	{
	  advance_child_dependent_cursors (pair, decomp, children, squash_columns, cursor, col_data);
	  break;
	}
      
      bool parent_flag = ppath[0][pcursor];
      bool child_flag  = ppath[1][pcursor];

      if (child_flag && !parent_flag)                                                                     // is an insertion
	add_child_step_to_buffer (pair, decomp, children, squash_columns, cursor);

      if (child_flag && parent_flag)                                                                      // is a match
	advance_child_dependent_cursors (pair, decomp, children, squash_columns, cursor, col_data);

      if (parent_flag) { col_data[child] = child_flag; break; }
    }

  return pcursor < ppath.columns();
}

bool Alignment_path::advance_child_dependent_cursors (const Row_pair& pair, const Decomposition& decomp, const vector<Row_index_set>& children, bool squash_columns, Row_pair_cursor& cursor, Nascent_column& col_data)
{
  const int child  = pair.second;
  
  bool found_dependents        = 0;
  bool dependents_still_active = 0;

  //  CLOG(1) << "c" << child << "(";
  
  for_const_contents (Row_index_set, children[child], c)
    {
      Row_pair child_pair (child, *c);
      if (decomp.find (child_pair) != decomp.end())
	{
	  found_dependents = 1;
	  if (advance_cursor_to_next_parent_step (child_pair, decomp, children, squash_columns, cursor, col_data))
	    dependents_still_active = 1;
	}
    }

  //  CLOG(1) << ")";
  
  return !found_dependents || dependents_still_active;}

bool Alignment_path::add_child_step_to_buffer (const Row_pair& pair, const Decomposition& decomp, const vector<Row_index_set>& children, bool squash_columns, Row_pair_cursor& cursor)
{
  const int parent = pair.first;
  const int child  = pair.second;

  Nascent_column col_data;
  col_data[child] = 1;
  if (parent != -1) col_data[parent] = 0;

  //  CLOG(1) << "b" << child << ":" << cursor[pair] << " ";
  
  bool ok = advance_child_dependent_cursors (pair, decomp, children, squash_columns, cursor, col_data);
  if (ok)
    {
      Row::size_type col = 0;
      for (int row = 0; row < rows(); row++)
	if (!squash_columns || col_data.find(row) != col_data.end())
	  col = max (col, path_data[row].size());
      
      for (int row = 0; row < rows(); row++)
	if (col_data.find(row) != col_data.end())
	  {
	    path_data[row].insert (path_data[row].end(), col + 1 - path_data[row].size(), 0);
	    path_data[row].back() = col_data[row];
	  }
	else if (!squash_columns)
	  {
	    path_data[row].insert (path_data[row].end(), col + 1 - path_data[row].size(), 0);
	    path_data[row].back() = 0;
	  }
    }
  return ok;
}


Alignment_path::Row Alignment_path::flag_nonempty_columns() const
{
  Row v = create_empty_row();
  for (int col = 0; col < columns(); col++) v[col] = !column_is_empty (col);
  return v;
}

Alignment_path Alignment_path::realign_pair (const Row_pair& pair, const Alignment_path& new_align, const vector<Row_pair>& dependencies)
{
  if (new_align.rows() != 2)
    THROW Standard_exception ("While re-aligning pair: new alignment must be a pairwise path");

  if (CLOGGING(2))
    {
      vector<sstring> rn (2, sstring("row_"));
      rn[0] << pair.first;
      rn[1] << pair.second;
      CL << "Changing pairwise alignment between rows " << pair.first << " and " << pair.second << " to:\n";
      new_align.show (CL, rn);
    }

  if (new_align.count_steps_in_row(0) != count_steps_in_row(pair.first))
    THROW Standard_exception ("While re-aligning pair: new pair alignment did not match original sequence lengths");

  if (new_align.count_steps_in_row(1) != count_steps_in_row(pair.second))
    THROW Standard_exception ("While re-aligning pair: new pair alignment did not match original sequence lengths");
  
  int reference_row = rows();
  append_row (create_full_row());

  vector<Row_pair> dep = dependencies;
  dep.push_back (Row_pair (pair.first, reference_row));

  Decomposition decomp = decompose (dep);
  if (decomp.find(pair) == decomp.end())
    {
      Row_pair tmp_pair = pair;
      swap (tmp_pair.first, tmp_pair.second);
      Alignment_path& tmp_align = (decomp[tmp_pair] = new_align);
      tmp_align[0].swap (tmp_align[1]);
    }
  else
    decomp[pair] = new_align;
  compose_and_log (decomp, 1);
  
  Pairwise_path align_delta;
  align_delta[0] = row(reference_row);
  erase_rows (reference_row);
  align_delta[1] = flag_nonempty_columns();
  erase_empty_columns();

  align_delta.erase_empty_columns();
  return align_delta;
}

Alignment_path Alignment_path::realign_row (int r, const Alignment_path& align_to_new_row_map, const vector<Row_pair>& dependencies)
{
  if (align_to_new_row_map.rows() != 2)
    THROW Standard_exception ("While re-aligning row: deltas must be pairwise paths");

  if (CLOGGING(2))
    {
      vector<sstring> rn;
      rn.push_back (sstring ("Old_align"));
      rn.push_back (sstring ("Old_row"));
      rn.push_back (sstring ("New_row"));
      CL << "Changing row " << r << "\n";
      Alignment_path::Decomposition decomp;
      Pairwise_path align_to_old_row_map;
      align_to_old_row_map[0] = create_full_row();
      align_to_old_row_map[1] = row(r);
      decomp[Alignment_path::Row_pair (0, 1)] = align_to_old_row_map;
      decomp[Alignment_path::Row_pair (0, 2)] = align_to_new_row_map;
      Alignment_path display_path;
      display_path.compose (decomp, 1);
      display_path.show (CL, rn);
    }

  if (align_to_new_row_map.count_steps_in_row(0) != columns())
    THROW Standard_exception ("While re-aligning row: row delta did not match original alignment width");

  row(r) = create_full_row();

  int new_row = rows();
  Decomposition decomp = decompose (dependencies);
  decomp[Row_pair (r, new_row)] = align_to_new_row_map;
  compose_and_log (decomp, 1);
  row(r).swap (row(new_row));

  Pairwise_path align_delta;
  align_delta[0].swap (row(new_row));
  erase_rows (new_row);
  align_delta[1] = flag_nonempty_columns();
  erase_empty_columns();

  align_delta.erase_empty_columns();
  return align_delta;
}

Alignment_path Alignment_path::realign_subpath (const Row_index_set& subpath_rows, const Alignment_path& new_subpath, const Alignment_path& old_to_new_subpath_map)
{
  Subalignment_path old_subpath (*this, subpath_rows, 1);

  if (new_subpath.rows() != (int) subpath_rows.size())
    THROW Standard_exception ("While re-aligning subpath: new subpath different size from old subpath");

  if (old_to_new_subpath_map.rows() != 2)
    THROW Standard_exception ("While re-aligning subpath: deltas must be pairwise paths");

  if (old_to_new_subpath_map.count_steps_in_row(0) != old_subpath.columns())
    THROW Standard_exception ("While re-aligning subpath: subpath delta did not match old subpath");

  if (old_to_new_subpath_map.count_steps_in_row(1) != new_subpath.columns())
    THROW Standard_exception ("While re-aligning subpath: subpath delta did not match new subpath");

  if (CLOGGING(2))
    {
      Alignment_path display_path = old_to_new_subpath_map;
      display_path.explode_row (1, new_subpath);
      display_path.explode_row (0, old_subpath);
      vector<sstring> rn;
      for_const_contents (Row_index_set, subpath_rows, r)
	{ sstring tmp ("Old_row_"); tmp << *r; rn.push_back (tmp); }
      for_const_contents (Row_index_set, subpath_rows, r)
	{ sstring tmp ("New_row_"); tmp << *r; rn.push_back (tmp); }
      CL << "Changing subpath:\n";
      display_path.show (CL, rn);
    }

  // create three row sets:
  // set 0 contains all the rows in the subpath that's being changed
  // set 1 contains all the rows that aren't in the subpath that's being changed
  // set 2 contains a single reference row with a 1 in every column
  //
  vector<Row_index_set> row_sets (3, Row_index_set());
  row_sets[0] = subpath_rows;

  vector<bool> in_set_zero (rows(), (bool) 0);
  for_const_contents (Row_index_set, subpath_rows, r) in_set_zero[*r] = 1;
  for (int r = 0; r < rows(); r++) if (!in_set_zero[r]) row_sets[1].push_back (r);

  int reference_row = rows();
  append_row (create_full_row());
  row_sets[2].push_back (reference_row);

  vector<Alignment_path> old_paths = implode (row_sets);

  // tie row sets 0 & 1 to the reference row
  //
  vector<Row_pair> dep;
  dep.push_back (Row_pair (2, 0));
  dep.push_back (Row_pair (2, 1));

  // replace row set 0 with the new subpath
  //
  Decomposition decomp = decompose (dep);
  decomp[Row_pair (0, 3)] = old_to_new_subpath_map;
  compose_and_log (decomp, 0);

  row(0).swap (row(3));
  old_paths[0] = new_subpath;
  erase_rows (3);
  explode (row_sets, old_paths);
  
  // calculate the old-to-new alignment map, then erase the reference row
  //
  Pairwise_path align_delta;
  align_delta[0] = row(reference_row);
  erase_rows (reference_row);
  align_delta[1] = flag_nonempty_columns();
  erase_empty_columns();

  align_delta.erase_empty_columns();
  return align_delta;
}


void Alignment_path::show (ostream& o, const vector<sstring>& row_names) const
{
  if ((int) row_names.size() != rows())
    THROW Standard_exception ("While trying to show an Alignment_path: row_names vector size mismatch");

  sstring::size_type max_name_len = 0;
  for (int r = 0; r < rows(); r++) max_name_len = max (max_name_len, row_names[r].size());

  for (int r = 0; r < rows(); r++)
    {
      o << row_names[r];
      for (int pad = (int) row_names[r].size(); pad < (int) max_name_len+1; pad++) o << ' ';
      for (int col = 0; col < (int) row(r).size(); col++) o << (path_data[r][col] ? '*' : '-');
      o << "\n";
    }
}

Subalignment_path::Subalignment_path (const Alignment_path& path, const Row_index_set& rows, bool collapse_empty_columns) : Alignment_path (0)
{
  path_data.reserve (rows.size());
  for_const_contents (Row_index_set, rows, r) path_data.push_back (path.row(*r));
  if (collapse_empty_columns) erase_empty_columns();
}

Subalignment_path::Subalignment_path (const Alignment_path& path, const Row_index_set& rows, Pairwise_path& align_subalign_map) : Alignment_path (rows.size(), path.columns())
{
  align_subalign_map.reset_rows();
  int mycol = 0;
  for (int col = 0; col < path.columns(); col++)
    {
      bool nonempty = 0;
      for (int r = 0; r < (int) rows.size(); r++)
	{
	  path_data[r][mycol] = path(rows[r],col);
	  nonempty = nonempty || path_data[r][mycol];
	}
      if (nonempty) { ++mycol; align_subalign_map.append_column (1,1); }
      else align_subalign_map.append_column (1,0);
    }
  erase_columns (mycol, columns() - mycol);
}

Pairwise_path::Pairwise_path() : Alignment_path (2)
{ }

Pairwise_path::Pairwise_path (const Alignment_path& path) : Alignment_path (path)
{
  if (rows() != 2) THROW Standard_exception ("Attempt to construct a Pairwise_path with rows != 2");
}

Pairwise_path::Pairwise_path (const Alignment_path& path, int row0, int row1, bool collapse_empty_columns) : Alignment_path (0)
{
  path_data.reserve (2);
  path_data.push_back (path.row(row0));
  path_data.push_back (path.row(row1));
  if (collapse_empty_columns) erase_empty_columns();
}

int Pairwise_path::match_columns() const
{
  int c = 0;
  for (int i = 0; i < columns(); ++i)
    if (parent()[i] && child()[i])
      ++c;
  return c;
}

int Pairwise_path::match_overlap (const Pairwise_path& p) const
{
  int overlap = 0;
  Sequence_coords coords_i = create_seq_coords();
  Sequence_coords coords_j = p.create_seq_coords();
  int col_i = 0;
  for (int col_j = 0; col_j < p.columns(); p.inc_seq_coords (coords_j, col_j++))
    {
      if (p.parent()[col_j] && p.child()[col_j])
	{
	  while (col_i < columns())
	    if (coords_i[0] < coords_j[0] || coords_i[1] < coords_j[1])
	      inc_seq_coords (coords_i, col_i++);
	    else
	      break;
	  if (coords_i == coords_j)
	    ++overlap;
	}
    }
  return overlap;
}

double Pairwise_path::mean_match_displacement() const
{
  int total_d = 0;
  int matches = 0;
  Sequence_coords coords = create_seq_coords();
  for (int col = 0; col < columns(); inc_seq_coords (coords, col++))
    if (parent()[col] && child()[col])
      {
	total_d += coords[0] - coords[1];
	++matches;
      }
  return matches == 0 ? 0 : ((double) total_d) / (double) matches;
}

void Pairwise_path::sort_indels()
{
  // pre log
  if (CTAGGING(2,SORTED_PAIRWISE_PATH))
    {
      CL << "Sorting the following pairwise path:\n";
      vector<sstring> row_name (2);
      row_name[0] << "Anc";
      row_name[1] << "Des";
      show (CL, row_name);
    }      

  // move I's before D's
  int last_nondel_col = -1;
  int n_shunts = 0;
  for (int col = 0; col < columns(); ++col)
    {
      const bool is_insert = !parent()[col] && child()[col];
      const bool is_delete = parent()[col] && !child()[col];
      if (is_insert && last_nondel_col >= 0 && last_nondel_col < col - 1)  // try to shunt this column left?
	{
	  for (int r = 0; r < rows(); ++r)
	    {
	      // have to do the swap explicitly, since STL's swap<> template apparently can't handle std::_Bit_reference as a type... sigh (IH, 5/20/09)
	      const bool tmp = row(r)[last_nondel_col+1];
	      row(r)[last_nondel_col+1] = row(r)[col];
	      row(r)[col] = tmp;
	    }
	  ++last_nondel_col;
	  ++n_shunts;
	}
      else if (!is_delete)
	last_nondel_col = col;
    }

  // post log
  if (n_shunts && CTAGGING(4,SORTED_PAIRWISE_PATH))
    CL << "Shunted " << n_shunts << " insert columns leftwards past delete columns\n";

  if (CTAGGING(2,SORTED_PAIRWISE_PATH))
    {
      CL << "Sorted pairwise path (" << n_shunts << " insert columns moved):\n";
      vector<sstring> row_name (2);
      row_name[0] << "Anc";
      row_name[1] << "Des";
      show (CL, row_name);
    }
}

Pairwise_path& Pairwise_path::operator*= (const Pairwise_path& p2)
{
  Decomposition decomp;
  decomp[Row_pair (0, 1)] = *this;
  decomp[Row_pair (1, 2)] = p2;
  Alignment_path align;
  align.compose (decomp, 0);    // NB don't squash columns ...
  row(0) = align.row(0);
  row(1) = align.row(2);
  erase_empty_columns();
  return *this;
}

Pairwise_path Pairwise_path::operator*  (const Pairwise_path& p2) const
{
  Pairwise_path p (*this);
  return p *= p2;
}

void Alignment_path::Output_mask::initialise_global (const Alignment_path& path)
{
  last_match = first_match = Sequence_coords (path.rows(), 0);
  for (int row = 0; row < path.rows(); ++row)
    {
      const int s = path.count_steps_in_row (row);
      if (s >= 0)
	last_match[row] = s - 1;
      else
	first_match[row] = last_match[row] = -1;
    }
  printable_columns.clear();
  for (int col = 0; col < path.columns(); ++col)
    printable_columns.push_back (col);
}

void Alignment_path::Output_mask::initialise_local (const Alignment_path& path)
{
  // initialise vectors
  last_match = first_match = Sequence_coords (path.rows(), -1);
  printable_columns.clear();

  // first pass: fill first_match[] and last_match[]
  Alignment_path::Sequence_coords seq_coords (path.create_seq_coords());
  Alignment_path::Row_index_set match_rows;
  match_rows.reserve (path.rows());
  for (int col = 0; col < path.columns(); ++col)
    {
      // first count matches in this column
      match_rows.clear();
      for (int row = 0; row < path.rows(); ++row)
	if (path (row, col))
	  match_rows.push_back (row);
      // now update rows involved in matches
      if (match_rows.size() > 1)
	for_const_contents (Alignment_path::Row_index_set, match_rows, row)
	  if (path (*row, col))
	    {
	      last_match[*row] = seq_coords[*row];
	      if (first_match[*row] < 0)
		first_match[*row] = seq_coords[*row];
	    }
      // increment sequence coords
      path.inc_seq_coords (seq_coords, col);
    }

  // second pass: set printable_columns[]
  update_printable_columns (path);
}

void Alignment_path::Output_mask::update_printable_columns (const Alignment_path& path)
{
  printable_columns.clear();
  Alignment_path::Sequence_coords seq_coords (path.create_seq_coords());
  for (int col = 0; col < path.columns(); ++col)
    {
      for (int row = 0; row < path.rows(); ++row)
	if (path (row, col))
	  {
	    const int f = first_match[row];
	    const int l = last_match[row];
	    const int c = seq_coords[row];
	    if (f >= 0 && c >= f && c <= l)
	      {
		printable_columns.push_back (col);
		break;
	      }
	  }
      path.inc_seq_coords (seq_coords, col);
    }
}

void Alignment::set_gap_chars (const sstring& gap_chars)
{
  if (gap_chars.size() == 0)
    THROWEXPR ("Must have at least one gap character");
  primary_gap_char = ((sstring&)gap_chars).front();  // cast away const
  char_is_gap = vector<int> (256, (int) 0);   // caution: hardwired 256!
  for_const_contents (sstring, gap_chars, gc)
    char_is_gap[*gc] = 1;
}

sstring Alignment::get_gap_chars()
{
  sstring gc (1, (char) primary_gap_char);
  for (unsigned char c = 0; c < (unsigned char) char_is_gap.size(); ++c)
    if (is_gap_char(c) && c != primary_gap_char)
      gc.push_back(c);
  return gc;
}

Alignment_initializer::Alignment_initializer()
{
  Alignment::set_gap_chars (sstring (DEFAULT_GAP_CHARS));
}

char Alignment::primary_gap_char = '\0';
vector<int> Alignment::char_is_gap (256, (int) 0);  // caution: hardwired 256!
Alignment_initializer singleton_alignment_initializer;

void Alignment::read_MUL (istream& in, Sequence_database& db)
{
  clear();

  vector<sstring> gapped_seq;
  Regexp re_ignore = "^[ \t]*#.*$";
  Regexp re = "^[ \t]*([^ \t]+)[ \t]+([^ \t]*)[ \t]*$";
  while (in && !in.eof())
    {
      sstring input;
      input.getline(in).chomp();
      if (!re_ignore.Match(input.c_str()))
	{
	  if (re.Match(input.c_str()) && re.SubStrings() == 2)
	    {
	      Named_profile tmp_np;
	      db.push_back(tmp_np);
	      Named_profile& np_cur = db.back();
	    
	      np_cur.name = re[1];
	      gapped_seq.push_back (re[2]);
	    
	      Biosequence& g = gapped_seq.back();
	      remove_copy_if (g.begin(), g.end(), back_inserter(np_cur.seq), Alignment::is_gap_char);
	    
	      vector<bool> seq_step_flags;
	      transform (g.begin(), g.end(), back_inserter(seq_step_flags), not1(ptr_fun(Alignment::is_gap_char)));
	    
	      prof.push_back (&np_cur.prof_sc);
	      row_name.push_back (np_cur.name);
	      path.append_row (seq_step_flags);
	    }
	  else if (input.size())
	    CLOGERR << "Warning: couldn't parse the following alignment input line: \"" << input << "\"\n";
	}
    }
  path.make_flush();
}

void Alignment::write_MUL (ostream& out, const Alphabet& alphabet, const sstring* column_labels, bool print_empty_rows) const
{
  if (path.rows() == 0) return;

  save_flags (out);
  out.setf (ios::left, ios::adjustfield);

  sstring::size_type max_name_width = 0;
  for (int r = 0; r < rows(); r++)
    max_name_width = max (max_name_width, row_name[r].size());

  if (column_labels != 0)
    {
      sstring rf ("#=RF");
      max_name_width = max (max_name_width, rf.size());
      out.width (max_name_width + 1);
      out << rf;
      out << *column_labels << "\n";
    }
  
  const Symbol_score_map wildcard_score_map = alphabet.flat_score_map (0);
  const char wildcard_char = alphabet.score2char (wildcard_score_map);
  for (int row = 0; row < path.rows(); row++)
    if (print_empty_rows || path.count_steps_in_row(row))  // omit empty rows, unless so requested
      {
	out.width (max_name_width + 1);
	out << row_name[row];
	int pos = 0;
	for (int col = 0; col < path.columns(); col++)
	  {
	    out << (path(row,col) ? (prof[row] ? alphabet.score2char((*prof[row])[pos++]) : wildcard_char) : gap_char());
	  }
	out << '\n';
      }
  restore_flags (out);
}

void Alignment::reset (int rows, int cols)
{
  path = Alignment_path (rows, cols);
  prof = vector<const Score_profile*> (rows, (const Score_profile*) 0);
  row_name = vector<sstring> (rows, sstring());
}

Alignment::Alignment (int rows)
  : Named_rows (rows), path (rows), prof (rows, (const Score_profile*) 0) { }

Alignment::Alignment (const Alignment_path& p)
  : Named_rows (p.rows()), path (p), prof (p.rows(), (const Score_profile*) 0) { }

Alignment::Alignment (const Pairwise_path& path, const Named_profile& xseq, const Named_profile& yseq)
  : Named_rows(2), path(path), prof(2)
{
  prof[0] = &xseq.prof_sc;
  prof[1] = &yseq.prof_sc;
  row_name[0] = xseq.name;
  row_name[1] = yseq.name;
}

Alignment::Alignment (const Alignment& align, const Row_index_set& row_set, bool collapse_empty_columns)
  : Named_rows (row_set.size()), path (row_set.size()), prof (row_set.size(), (const Score_profile*) 0)
{
  path = Subalignment_path (align.path, row_set, collapse_empty_columns);
  for (int i = 0; i < (int) row_set.size(); i++)
    {
      prof[i] = align.prof[row_set[i]];
      row_name[i] = align.row_name[row_set[i]];
    }
}

Alignment::Alignment (const Sequence_database_index& seqdb_index)
  : Named_rows (seqdb_index.size()),
    path (seqdb_index.size()),
    prof (seqdb_index.size())
{
  for (int i = 0; i < rows(); i++)
    {
      prof[i] = &seqdb_index.profile[i]->prof_sc;
      row_name[i] = seqdb_index.name[i];
    }
}

void Alignment::assert_rows_fit_sequences() const
{
  for (int r = 0; r < rows(); ++r)
    if (prof[r] != 0)
      if (path.count_steps_in_row(r) != (int) prof[r]->size())
	THROW Standard_exception ("Alignment is not global");
}

void Alignment::discard_wild_sequences (const Alphabet& alphabet)
{
  for (int r = 0; r < rows(); ++r)
    if (prof[r] != 0)
      if (alphabet.all_wild (*prof[r])) prof[r] = 0;
}

Named_rows::Named_rows (int rows)
  : row_name (rows, sstring())
{ }

Named_rows::Row_index_set Named_rows::find_rows_by_name (const char* name) const
{
  Row_index_set result;
  const sstring namestr (name);
  for (int r = 0; r < rows(); ++r)
    if (((Named_rows&)*this).row_name[r] == namestr)
      result.push_back(r);
  return result;
}

int Named_rows::find_row_by_name (const char* name) const
{
  Row_index_set rows = find_rows_by_name (name);
  if (rows.size() < 1) { sstring e; e << "Couldn't find any rows named '" << name << "'"; THROW Standard_exception(e); }
  if (rows.size() > 1) { sstring e; e << "Multiple rows named '" << name << "'"; THROW Standard_exception(e); }
  return rows[0];
}

bool Named_rows::names_unique() const
{
  set<sstring> names;
  for (int r = 0; r < rows(); ++r)
    {
      const sstring& name = row_name[r];
      if (names.find (name) != names.end())
	return false;
      names.insert (name);
    }
  return true;
}

void Named_rows::assert_names_unique() const
{
  if (!names_unique())
    THROWEXPR ("Named_rows (base class of Alignment) failed assert_names_unique()");
}

Score Alignment::null_emit_score (const vector<Score>& null_model) const
{
  Score sc = 0;
  for_const_contents (vector<const Score_profile*>, prof, row)
    if (*row)
      ScorePMulAcc (sc, (*row)->null_score (null_model));
  return sc;
}

const Symbol_score_map& Alignment::get_ssm (int row, int col) const
{
  if (is_gap(row,col)) THROW Standard_exception ("Attempt to retrieve gap character"); 
  return (*prof[row]) [path.get_seq_pos (row, col)];
}

int Alignment::residue_pairs() const
{
  int pairs = 0;
  for (int col = 0; col < path.columns(); ++col)
    {
      int rows = 0;
      for (int row = 0; row < path.rows(); ++row)
	if (path (row, col))
	  ++rows;
      pairs += rows * (rows - 1) / 2;
    }
  return pairs;
}

int Alignment::count_pairwise_insertions() const
{
	int inserts = 0;
	Pairwise_gap_set my_gap_set = pairwise_insertions();
	Pairwise_gap_set::iterator it;	
	for(it = my_gap_set.begin(); it != my_gap_set.end(); it++)
	{
		inserts += it->second.size();
	}
	return inserts;
}

int Alignment::count_pairwise_deletions() const
{
	int deletes = 0;
	Pairwise_gap_set my_gap_set = pairwise_deletions();
	Pairwise_gap_set::iterator it;	
	for(it = my_gap_set.begin(); it != my_gap_set.end(); it++)
	{
		deletes += it->second.size();
	}
	return deletes;
}

int Alignment::residue_pair_overlap (const Alignment& align) const
{
  Residue_pair_set my_pair_set = residue_pair_set();
  Residue_pair_set your_pair_set = align.residue_pair_set();
  return residue_pair_overlap (my_pair_set, your_pair_set);
}

int Alignment::residue_pair_overlap (const Residue_pair_set& my_pair_set, const Residue_pair_set& your_pair_set)
{
  int overlap = 0;
  for_const_contents (Residue_pair_set, my_pair_set, my_name_pset)
    {
      Residue_pair_set::const_iterator your_name_pset = your_pair_set.find (my_name_pset->first);
      if (your_name_pset != your_pair_set.end())
	for_const_contents (Pairform_alignment, my_name_pset->second, my_pair)
	  if (your_name_pset->second.find (*my_pair) != your_name_pset->second.end())
	    ++overlap;
    }
  return overlap;
}

int Alignment::pairwise_deletion_overlap (const Alignment& align) const
{
	Pairwise_gap_set my_gap_set = pairwise_deletions();
	Pairwise_gap_set your_gap_set = align.pairwise_deletions();
	return pairwise_gap_overlap(my_gap_set, your_gap_set);
}

int Alignment::pairwise_insertion_overlap (const Alignment& align) const
{
	Pairwise_gap_set my_gap_set = pairwise_insertions();
	Pairwise_gap_set your_gap_set = align.pairwise_insertions();
	return pairwise_gap_overlap(my_gap_set, your_gap_set);
}

int Alignment::pairwise_gap_overlap (const Pairwise_gap_set& my_gap_set, const Pairwise_gap_set& your_gap_set)
{
	int overlap = 0;
	for_const_contents (Pairwise_gap_set, my_gap_set, my_name_pset)
	{
		Pairwise_gap_set::const_iterator your_name_pset = your_gap_set.find (my_name_pset->first);
		if (your_name_pset != your_gap_set.end())
			for_const_contents (Gap_coord_set, my_name_pset->second, my_coord)
				if(your_name_pset->second.find (*my_coord) != your_name_pset->second.end())
					++overlap;
	}
	return overlap;
} 

float Alignment::ama_similarity (const Alignment &align) const
{
	float score_builder = 0;
	int length_sum = 0;
	for (int x = 0; x < rows(); ++x)
	{
		length_sum += (int) path.count_steps_in_row(x);
	}
	length_sum = length_sum * (rows() - 1);
	int HH_HI_HD = 2 * residue_pair_overlap(align) + pairwise_insertion_overlap(align) + pairwise_deletion_overlap(align);
	score_builder = length_sum - HH_HI_HD; 
	score_builder = 1 - (score_builder / length_sum); 
	return score_builder;
}	

int Alignment::column_overlap (const Alignment& align) const
{
 Col_set my_col_set = col_set();
 Col_set your_col_set = align.col_set();
 return column_overlap (my_col_set, your_col_set); 
}

int Alignment::column_overlap (const Col_set& my_col_set, const Col_set& your_col_set)
{
   int overlap = 0;
   for_const_contents(Col_set, my_col_set, my_coord_vector)
      if(your_col_set.find(*my_coord_vector) != your_col_set.end())
         ++overlap;
   return overlap;
} 	

Alignment::Residue_pair_set Alignment::residue_pair_set() const
{
  Residue_pair_set pair_set;
  Row_index_set rows;
  Alignment_path::Sequence_coords seq_coords = path.create_seq_coords();
  for (int col = 0; col < path.columns(); path.inc_seq_coords (seq_coords, col), ++col)
    {
      rows.clear();
      for (int row = 0; row < path.rows(); ++row)
	if (path (row, col))
	  rows.push_back (row);
      for (int i = 0; i < (int) rows.size(); ++i)
	for (int j = i + 1; j < (int) rows.size(); ++j)
	  {
	    sstring name1 = row_name[rows[i]];
	    sstring name2 = row_name[rows[j]];
	    int coord1 = seq_coords[rows[i]];
	    int coord2 = seq_coords[rows[j]];
	    if (name2 < name1) { swap (name1, name2); swap (coord1, coord2); }
	    pair_set [Row_name_pair (name1, name2)].insert (Pairform_match (coord1, coord2));
	  }
    }
  return pair_set;
}

Alignment::Col_set Alignment::col_set() const
{
  Col_set column_set;
  Alignment_path::Sequence_coords seq_coords = path.create_seq_coords();
  Col_coord_vector coords;
  Row_index_set rows = sorted_rows();
  vector<int> last_coord;
  for (int col = 0; col < path.columns(); path.inc_seq_coords (seq_coords, col), ++col)
  {
     coords.clear();
     for (unsigned int i = 0; i < rows.size(); ++i)
     {
	if(not_gap(rows[i] , col))	
           coords.push_back (seq_coords[rows[i]]);
        else
           coords.push_back (-1);
     }
     column_set.insert(coords);
  }
  return column_set;
}

const Alignment::Row_index_set Alignment::sorted_rows() const
{
   Row_index_set rows;
   for(int row = 0; row < path.rows(); ++row)   // insertion sort implemented to fix bug where different alignment orderings produce different 
   {                                               // AMA scores... probably a more efficient way to do this? -LEB 05/19/2008
      if(rows.empty()){
         rows.push_back(row);
      }
      else{
         Row_index_set::iterator it;
         bool inserted = false;
         sstring new_row_name = row_name[row];
         for( it = rows.begin(); it != rows.end(); it++)
         {
            sstring old_row_name = row_name[*it];
            if (new_row_name < old_row_name){
               rows.insert(it, row);
               inserted = true;
               break;
            }
         }
         if(!inserted){
            rows.insert(rows.end(), row);
         }
      }
   }
   return rows;
}


Alignment::Pairwise_gap_set Alignment::pairwise_deletions() const
{
	Row_index_set rows = sorted_rows();
	Pairwise_gap_set pair_set;
	Alignment_path::Sequence_coords seq_coords = path.create_seq_coords();
	for  (int col = 0; col < path.columns(); path.inc_seq_coords(seq_coords, col), ++col)
	{	
		for(int i = 0; i < (int) rows.size(); i++)
			for(int j = i + 1; j < (int) rows.size(); j++)
				{
					if((path(rows[i], col)) && (!path (rows[j], col)))
					{
						sstring name1 = row_name[rows[i]];
						sstring name2 = row_name[rows[j]];
						int coord1 = seq_coords[rows[i]];
						pair_set [Row_name_pair (name1, name2)].insert (coord1);
					}
				}
			
	}
	return pair_set;
}

Alignment::Pairwise_gap_set Alignment::pairwise_insertions() const
{
	Row_index_set rows = sorted_rows();
	Pairwise_gap_set pair_set;
	Alignment_path::Sequence_coords seq_coords = path.create_seq_coords();
	for  (int col = 0; col < path.columns(); path.inc_seq_coords(seq_coords, col), ++col)
	{
		for(int i = 0; i < (int) rows.size(); i++)
				for(int j = i + 1; j < (int) rows.size(); j++)
				{
					if((!path (rows[i], col)) && (path (rows[j], col)))
					{
						sstring name1 = row_name[rows[i]];
						sstring name2 = row_name[rows[j]];
						int coord2 = seq_coords[rows[j]]; 
						pair_set [Row_name_pair (name1, name2)].insert (coord2);
					}
				}
			
	}
	return pair_set;
}

void Alignment::swap_rows (int row1, int row2)
{
  swap (row_name[row1], row_name[row2]);
  swap (path[row1], path[row2]);
  swap (prof[row1], prof[row2]);
}	

void Aligned_score_profile::init (const Alignment& align, const vector<Named_profile*>& np, const Alphabet& alph)
{
  if (CTAGGING(3,ALLOC))
    CL << "Allocating " << sizeof(Symbol_score_map*)*align.rows()*align.columns() << " bytes for alignment ("
       << align.rows() << " rows * " << align.columns() << " columns)\n";
  resize (align.rows(), align.columns(), (Symbol_score_map*) 0);
  if (rows() != (int) np.size())
    THROWEXPR ("Aligned_score_profile: supplied vector of Named_profile*'s has wrong size");
  Alignment_path::Sequence_coords seq_coords = align.path.create_seq_coords();
  for (int col = 0; col < columns(); ++col) {
    for (int row = 0; row < rows(); ++row) {
      if (align.prof[row] != 0 && align.not_gap (row, col))
	(*this) (row, col) = &alph.char2score (np[row]->seq[seq_coords[row]]);
    }
    align.path.inc_seq_coords (seq_coords, col);
  }
}

Aligned_score_profile Aligned_score_profile::window (int start, int len) const
{
  Aligned_score_profile asp;
  asp.resize (rows(), len, (Symbol_score_map*) 0);
  for (int col = 0; col < len; ++col)
    for (int row = 0; row < rows(); ++row)
      asp (row, col) = (*this) (row, col + start);
  return asp;
}

void Aligned_score_profile::show (ostream& out) const
{
  Symbol_score_iterator ssi;
  const int cols = ysize(), rows = xsize();
  out << "Aligned_score_profile:\n";
  for (int col = 0; col < cols; col++)
  {
	for (int row = 0; row < rows; row++)
	{
	  out << "Column " << col << " Row " << row << " ";
	  const Symbol_score_map* ssm = (*this) (row,col);
	  if (ssm != 0)
	  {
		for (ssi = ssm->begin(); ssi != ssm->end(); ssi++)
		  out << ssi->first << ":" << ssi->second << " ";
		out << "\n";
	  }
	  else
		out << "Null\n";
	}
  }
}
