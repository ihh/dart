#include <stdlib.h>
#include "seq/local.h"
#include "util/Regexp.h"

void Local_alignment::read_MUL (istream& in, Sequence_database& db)
{
  clear();
  db.clear();
  start_coord.clear();

  vector<sstring> gapped_seq;
  Regexp re_ignore = "^[ \t]*#.*$";
  Regexp re = "^[ \t]*([^ \t/]+)[/0123456789-]*[ \t]+([^ \t]*)[ \t]*$";
  while (!in.eof())
    {
      sstring input;
      input.getline(in).chomp();
      if (!re_ignore.Match(input.c_str()) && re.Match(input.c_str()) && re.SubStrings() == 2)
	{
	  db.push_back(Named_profile());
	  Named_profile& np_cur = db.back();

	  np_cur.name = re[1];
	  gapped_seq.push_back (re[2]);

	  Biosequence& g = gapped_seq.back();
	  remove_copy_if (g.begin(), g.end(), back_inserter(np_cur.seq), Local_alignment::is_gap_char);

	  vector<bool> seq_step_flags;
	  transform (g.begin(), g.end(), back_inserter(seq_step_flags), not1(ptr_fun(Local_alignment::is_gap_char)));

	  prof.push_back (&np_cur.prof_sc);
	  row_name.push_back (np_cur.name);
	  path.append_row (seq_step_flags);
	  start_coord.push_back (0);
	}
      else if (input.size())
	CLOGERR << "Warning: couldn't parse the following alignment input line: \"" << input << "\"\n";
    }
  path.make_flush();
}

void Local_alignment::write_MUL (ostream& out, const Alphabet& alphabet, const sstring* column_labels) const
{
  if (path.rows() == 0) return;

  save_flags (out);
  left_align (out);

  vector<sstring> row_name_with_coords (rows(), sstring());
  sstring::size_type max_name_width = 0;
  for (int r = 0; r < rows(); r++)
    {
      row_name_with_coords[r] << row_name[r] << "/" << start_coord[r] << "-" << end_coord(r);
      max_name_width = max (max_name_width, row_name_with_coords[r].size());
    }
  
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
  for (int seq = 0; seq < path.rows(); seq++)
    {
      out.width (max_name_width + 1);
      out << row_name_with_coords[seq];
      
      int pos = 0;
      for (int col = 0; col < path.columns(); col++)
	out << (path(seq,col) ? (prof[seq] ? alphabet.score2char((*prof[seq])[start_coord[seq] + pos++]) : wildcard_char) : '-');
      out << '\n';
    }
  
  restore_flags (out);
}

int Local_alignment::end_coord (int row) const
{
  return ((Local_alignment&) *this) . start_coord[row] + path.count_steps_in_row (row) - 1;
}
