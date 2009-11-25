// pk 3/05 accept long input lines

#include "seq/stockholm.h"

Stockholm::Stockholm (int rows, int cols) : Alignment()
{
  reset (rows, cols);
}

Stockholm* Stockholm::deep_copy (const Stockholm& stock, Sequence_database& seqdb)
{
  Stockholm* new_stock = new Stockholm (stock);
  for (unsigned int r = 0; r < stock.np.size(); ++r)
    if (stock.np[r]) {
      Named_profile np (*stock.np[r]);
      seqdb.push_back (np);
      new_stock->np[r] = &seqdb.back();
    }
  return new_stock;
}

void Stockholm::reset (int rows, int cols)
{
  Alignment::reset (rows, cols);
  row_index.clear();
  gf_annot.clear();
  gf_index.clear();
  gc_annot.clear();
  gs_annot.clear();
  gr_annot.clear();
  np = vector<Named_profile*> (rows);
}

void Stockholm::update_row_index()
{
  row_index.clear();
  for (int i = 0; i < (int) row_name.size(); ++i)
    row_index[row_name[i]] = i;
}

void Stockholm::read_Stockholm (istream& in, Sequence_database& db, int max_line)	// increase max_line to accept long input lines
{
  clear();

  vector<sstring> gapped_seq;
  Regexp re_stock = "^[ \t]*#[ \t]*" Stockholm_format_identifier "[ \t]*" Stockholm_version_identifier "[ \t]*$";  // format & version identifiers
  Regexp re_sep = "^[ \t]*" Stockholm_alignment_separator "[ \t]*$";  // alignment separator lines, "//"
  Regexp re_gf = "^[ \t]*" Stockholm_file_annotation "[ \t]+([^ \t]+)[ \t]+(.*)$";  // #=GF lines
  Regexp re_gc = "^[ \t]*" Stockholm_column_annotation "[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]*$";  // #=GC lines
  Regexp re_gs = "^[ \t]*" Stockholm_sequence_annotation "[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]+(.*)$";  // #=GS lines
  Regexp re_gr = "^[ \t]*" Stockholm_sequence_column_annotation "[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]*$";  // #=GR lines
  Regexp re_row = "^[ \t]*([^ \t#]+)[ \t]+([^ \t]*)[ \t]*$";  // alignment row data
  Regexp re_nonwhite = "[^ \t]";  // whitespace

  bool found_format_version_id = FALSE;
  while (in && !in.eof())
    {
      sstring input;
	  if (max_line > 0)							// pk if not using default input line length
		input.getline(in, max_line).chomp();
	  else
		input.getline(in).chomp();
      if (re_sep.Match(input.c_str()))
	break;

      if (re_gf.Match(input.c_str()) && re_gf.SubStrings() == 2)  // #=GF
	{
	  const sstring& gf_tag = re_gf[1];
	  const sstring& gf_data = re_gf[2];
	  gf_annot.push_back (Tag_value (gf_tag, gf_data));
	  gf_index[gf_tag].insert (gf_annot.size() - 1);
	}

      else if (re_gc.Match(input.c_str()) && re_gc.SubStrings() == 2)  // #=GC
	gc_annot[re_gc[1]].append (re_gc[2]);

      else if (re_gs.Match(input.c_str()) && re_gs.SubStrings() == 3)  // #=GS
	gs_annot[re_gs[1]][re_gs[2]] = re_gs[3];

      else if (re_gr.Match(input.c_str()) && re_gr.SubStrings() == 3)  // #=GR
	gr_annot[re_gr[1]][re_gr[2]].append (re_gr[3]);

      else if (re_row.Match(input.c_str()) && re_row.SubStrings() == 2)  // row data
	{
	  const sstring& name = re_row[1];
	  const sstring& row_data = re_row[2];

	  const bool is_new_row = row_index.find(name) == row_index.end();

	  Named_profile* np_cur;
	  Biosequence* g;
	  vector<bool>* p;

	  if (is_new_row)
	    {
	      Named_profile tmp_np;
	      db.push_back (tmp_np);
	      np_cur = &db.back();
	      np_cur->name = name;

	      np.push_back (np_cur);
	      prof.push_back (&np_cur->prof_sc);
	      row_name.push_back (np_cur->name);

	      row_index[name] = rows() - 1;  // update row_index manually

	      if (gs_annot.find (name) == gs_annot.end())
		gs_annot[name] = Annotation();

	      if (gr_annot.find (name) == gr_annot.end())
		gr_annot[name] = Annotation();

	      gapped_seq.push_back (Biosequence());
	      g = &gapped_seq.back();

	      path.append_row (vector<bool>());
	      p = &path[path.rows()-1];
	    }
	  else
	    {
	      const int r = row_index[name];
	      np_cur = np[r];
	      g = &gapped_seq[r];
	      p = &path[r];
	    }

	  // horrible functional-style methods that i thought were cool when I first learned STL and now can't be arsed to excise
	  g->append (row_data);
	  remove_copy_if (row_data.begin(), row_data.end(), back_inserter(np_cur->seq), Alignment::is_gap_char);
	  transform (row_data.begin(), row_data.end(), back_inserter(*p), not1(ptr_fun(Alignment::is_gap_char)));
	}
      
      else if (re_stock.Match (input.c_str()))  // ignore format-version identifiers
	found_format_version_id = TRUE;

      else if (input.size() && re_nonwhite.Match (input.c_str()))
	CLOGERR << "Warning: couldn't parse the following alignment input line: \"" << input << "\"\n";
    }

  path.make_flush();

  CTAG(3,STOCKHOLM) << "Read Stockholm alignment: " << rows() << " rows, " << columns() << " columns\n";
}

void Stockholm::write_Stockholm_NSE (ostream& out, char annotation_wildcard_char) const
{
  // create local output mask
  Output_mask out_mask;
  out_mask.initialise_local (path);
  // extend the output mask using annotations (hacky, intended to fix bug where unaligned fold chars are sometimes chopped off; 2/27/2008)
  Alignment_path::Sequence_coords seq_coords (path.create_seq_coords());
  for (int col = 0; col < path.columns(); ++col)
    {
      for (int row = 0; row < path.rows(); ++row)
	if (path (row, col))
	  {
	    bool found_annot = false;
	    for_const_contents (Annotation, ((Stockholm&)*this).gr_annot[row_name[row]], tag_val)
	      if (tag_val->second[col] != annotation_wildcard_char)
		{
		  found_annot = true;
		  break;
		}
	    if (found_annot)
	      {
		const int seq_pos = seq_coords[row];
		if (seq_pos > out_mask.last_match[row])
		  out_mask.last_match[row] = seq_pos;
		if (seq_pos < out_mask.first_match[row])
		  out_mask.first_match[row] = seq_pos;
	      }
	  }
      // increment sequence coords
      path.inc_seq_coords (seq_coords, col);
    }
  out_mask.update_printable_columns (path);
  // print
  write_Stockholm_header (out);
  write_Stockholm_body (out, out_mask, TRUE);
  write_Stockholm_separator (out);
}

void Stockholm::write_Stockholm (ostream& out) const
{
  // create global output mask
  Output_mask out_mask;
  out_mask.initialise_global (path);
  // print
  write_Stockholm_header (out);
  write_Stockholm_body (out, out_mask, FALSE);
  write_Stockholm_separator (out);
}

void Stockholm::write_Stockholm_body (ostream& out, const Output_mask& out_mask, bool use_NSE_names) const
{
  Stockholm& stock = (Stockholm&) *this;  // cast away const to allow operator[] access to maps... hacky
  save_flags (out);
  left_align (out);

  // create the NSE names
  vector<sstring> nse_row_name = row_name;
  if (use_NSE_names)
    for (int r = 0; r < rows(); r++)
      nse_row_name[r] << '/' << (int) (out_mask.first_match[r] + 1) << '-' << (int) (out_mask.last_match[r] + 1);
  
  // get the field width for the tag columns
  // max_name_width takes into account sequence names, GF/GS/GC/GR annotations,
  // and GS annotations for internal nodes...
  // scroll down to see where it gets its final value
  int max_name_width = 0;
  for (int r = 0; r < rows(); r++)
    {
      const sstring& name = row_name[r];
      const sstring& nse_name = nse_row_name[r];
      const int nse_name_len = nse_name.size();
      max_name_width = max (max_name_width, nse_name_len);
      for_const_contents (Annotation, stock.gs_annot[name], na)
	max_name_width = max (max_name_width, (int) na->first.size() + nse_name_len + 6);  // 6 chars for "#=GS name "
      for_const_contents (Annotation, stock.gr_annot[name], na)
	max_name_width = max (max_name_width, (int) na->first.size() + nse_name_len + 6);  // 6 chars for "#=GR name "
    }
  for_const_contents (vector<Tag_value>, stock.gf_annot, na)
    max_name_width = max (max_name_width, (int) na->first.size() + 5);  // 5 chars for "#=GF "
  for_const_contents (Annotation, stock.gc_annot, na)
    max_name_width = max (max_name_width, (int) na->first.size() + 5);  // 5 chars for "#=GC "

  // go ahead and set up some GS stuff
  // (we need to go through the internal nodes anyways to set max_name_width)
  const set<sstring> row_names (row_name.begin(), row_name.end());
  vector<sstring> gs_name (row_name), nse_gs_name (nse_row_name);
  for_const_contents (Row_annotation, stock.gs_annot, ra)
    if (row_names.find (ra->first) == row_names.end())
      {
	gs_name.push_back (ra->first);
	nse_gs_name.push_back (ra->first);
	// resize max_name_width if necessary
	max_name_width = max (max_name_width, (int) ra->first.size() + (int) ra->first.size() + 6);  // 6 chars for "#=GS name "
      }
  // now max_name_width has its final value

  // print it out...
  // #=GF lines
  for_const_contents (vector<Tag_value>, stock.gf_annot, na)
    {
      sstring gf;
      gf << Stockholm_file_annotation << ' ' << na->first;
      out.width (max_name_width + 1);
      out << gf << na->second << "\n";
    }

  // #=GS lines
  for (int r = 0; r < (int) gs_name.size(); ++r)
    {
      const sstring& name = gs_name[r];
      const sstring& nse_name = nse_gs_name[r];
      for_const_contents (Annotation, stock.gs_annot[name], na)
	{
	  sstring gs;
	  gs << Stockholm_sequence_annotation << ' ' << nse_name << ' ' << na->first;
	  const vector<sstring> gs_vals = na->second.split ("\n");
	  for_const_contents (vector<sstring>, gs_vals, gs_val)
	    {
	      out.width (max_name_width + 1);
	      out << gs << *gs_val << '\n';
	    }
	}
    }

  // main alignment
  for (int row = 0; row < path.rows(); row++)
    {
      const sstring& name = row_name[row];
      const sstring& nse_name = nse_row_name[row];
      // row data
      if ((int) path[row].size() != columns())
	CLOGERR << "Warning: skipping row '" << nse_name << "' as it has the wrong number of columns ("
		<< path[row].size() << " instead of " << columns() << ")\n";
      else
	{
	  out.width (max_name_width + 1);
	  out << nse_name;
	  if (path.count_steps_in_row (row))
	    {
	      int pos = out_mask.first_match[row];
	      for_const_contents (Column_index_set, out_mask.printable_columns, col)
		{
		  const char c =
		    path (row, *col)
		    ? (np[row]
		       ? np[row]->seq[pos++]  // pos gets incremented here
		       : Default_alignment_wildcard_char)
		    : gap_char();
		  out << c;
		}
	    }
	  else
	    for_const_contents (Column_index_set, out_mask.printable_columns, col)
	      out << gap_char();
	  out << '\n';
	}

      // #=GR lines
      for_const_contents (Annotation, stock.gr_annot[name], na)
	if ((int) na->second.size() != columns())
	  CLOGERR << "Warning: skipping '" << Stockholm_sequence_column_annotation << ' ' << name << ' ' << na->first << "' as it has the wrong number of columns ("
		  << na->second.size() << " instead of " << columns() << ")\n";
	else
	  {
	    sstring gr;
	    gr << Stockholm_sequence_column_annotation << ' ' << nse_name << ' ' << na->first;
	    out.width (max_name_width + 1);
	    out << gr;
	    for_const_contents (Column_index_set, out_mask.printable_columns, col)
	      out << na->second[*col];
	    out << "\n";
	    }
    }

  // #=GC lines
  for_const_contents (Annotation, stock.gc_annot, na)
    if ((int) na->second.size() != columns())
      CLOGERR << "Warning: skipping '" << Stockholm_sequence_column_annotation << ' ' << na->first << "' as it has the wrong number of columns ("
	      << na->second.size() << " instead of " << columns() << ")\n";
    else
      {
	sstring gc;
	gc << Stockholm_column_annotation << ' ' << na->first;
	out.width (max_name_width + 1);
	out << gc;
	for_const_contents (Column_index_set, out_mask.printable_columns, col)
	  out << na->second[*col];
	out << "\n";
      }

  restore_flags (out);
}

sstring Stockholm::get_row_as_string (int row)
{
  sstring s;
  if (path.count_steps_in_row (row))
    {
      int pos = 0;
      for (int col = 0; col < columns(); ++col)
	{
	  const char c =
	    path (row, col)
	    ? (np[row]
	       ? np[row]->seq[pos++]  // pos gets incremented here
	       : Default_alignment_wildcard_char)
	    : gap_char();
	  s << c;
	}
    }
  else
    for (int col = 0; col < columns(); ++col)
      s << gap_char();
  return s;
}

void Stockholm::write_Stockholm_header (ostream& out)
{
    out << Stockholm_header;
}

void Stockholm::write_Stockholm_separator (ostream& out)
{
  out << Stockholm_alignment_separator << '\n';
}

void Stockholm::set_np (int row, const Named_profile& row_np)
{
  np[row] = (Named_profile*) &row_np;
  prof[row] = &row_np.prof_sc;
  row_name[row] = row_np.name;
  update_row_index();
}

void Stockholm::discard_wild_sequences (const Alphabet& alphabet)
{
  Alignment::discard_wild_sequences (alphabet);
  for (int r = 0; r < rows(); ++r)
    if (!prof[r])
      np[r] = 0;
}

int Stockholm::add_row()
{
  return add_row (path.create_empty_row(), sstring(), (Named_profile*) 0);
}

int Stockholm::add_row (const vector<bool>& new_row, const sstring& new_row_name, Named_profile* new_np)
{
  np.push_back (new_np);
  prof.push_back (new_np ? &new_np->prof_sc : (Score_profile*) 0);
  row_name.push_back (new_row_name);
  path.append_row (new_row);
  const int new_row_index = rows() - 1;
  row_index[new_row_name] = new_row_index;
  return new_row_index;
}

sstring Stockholm::get_gf_annot (const sstring& tag) const
{
  sstring s;
  const Tag_index_map::const_iterator gf_iter = gf_index.find (tag);
  if (gf_iter != gf_index.end())
    for_const_contents (set<int>, gf_iter->second, gfi)
      s.append (gf_annot[*gfi].second);
  return s;
}

sstring Stockholm::get_gc_annot (const sstring& tag) const
{
  sstring s;
  const Annotation::const_iterator gc_iter = gc_annot.find (tag);
  if (gc_iter != gc_annot.end())
    s = gc_iter->second;
  return s;
}

sstring Stockholm::get_gr_annot (const sstring& seqname, const sstring& tag) const
{
  sstring s;
  const Row_annotation::const_iterator gr_iter = gr_annot.find (seqname);
  if (gr_iter != gr_annot.end())
    {
      const Annotation::const_iterator row_iter = gr_iter->second.find (tag);
      if (row_iter != gr_iter->second.end())
	s = row_iter->second;
    }
  return s;
}

sstring Stockholm::get_gs_annot (const sstring& seqname, const sstring& tag) const
{
  sstring s;
  const Row_annotation::const_iterator gs_iter = gs_annot.find (seqname);
  if (gs_iter != gs_annot.end())
    {
      const Annotation::const_iterator row_iter = gs_iter->second.find (tag);
      if (row_iter != gs_iter->second.end())
	s = row_iter->second;
    }
  return s;
}

void Stockholm::clear_gf_annot()
{
  gf_annot.clear();
  gf_index.clear();
}

void Stockholm::clear_gf_annot (const sstring& tag)
{
  vector<Tag_value> new_gf_annot;
  Tag_index_map new_gf_index;
  for (int old_idx = 0; old_idx < (int) gf_annot.size(); ++old_idx)
    if (gf_annot[old_idx].first != tag)
      {
	new_gf_index[gf_annot[old_idx].first].insert (new_gf_annot.size());
	new_gf_annot.push_back (gf_annot[old_idx]);
      }
  swap (gf_annot, new_gf_annot);
  swap (gf_index, new_gf_index);
}

void Stockholm::add_gf_annot (const sstring& tag, const sstring& val)
{
  gf_annot.push_back (Tag_value (tag, val));
  gf_index[tag].insert ((int) gf_annot.size() - 1);
}

void Stockholm::set_gc_annot (const sstring& tag, const sstring& val)
{
  gc_annot[tag] = val;
}

void Stockholm::set_gs_annot (const sstring& seqname, const sstring& tag, const sstring& val)
{
  gs_annot[seqname][tag] = val;
}

void Stockholm::set_gr_annot (const sstring& seqname, const sstring& tag, const sstring& val)
{
  gr_annot[seqname][tag] = val;
}

void Stockholm::erase_empty_columns()
{
  assert_flush();

  int ncol = 0;
  for (int col = 0; col < columns(); col++)
    if (!path.column_is_empty (col))
      {
	if (ncol != col)
	  {
	    for (int row = 0; row < rows(); row++)
	      path.row(row)[ncol] = path.row(row)[col];
	    for_contents (Annotation, gc_annot, gc)
	      gc->second[ncol] = gc->second[col];
	    for_contents (Row_annotation, gr_annot, gr)
	      for_contents (Annotation, gr->second, gr_row)
	      gr_row->second[ncol] = gr_row->second[col];
	  }
	++ncol;
      }

  const int cols_to_erase = columns() - ncol;
  path.erase_columns (ncol, cols_to_erase);
  for_contents (Annotation, gc_annot, gc)
    gc->second.erase (ncol, cols_to_erase);
  for_contents (Row_annotation, gr_annot, gr)
    for_contents (Annotation, gr->second, gr_row)
    gr_row->second.erase (ncol, cols_to_erase);
}

void Stockholm::assert_flush() const
{
  const int cols = columns();
  for (int row = 0; row < rows(); row++)
    if ((int) path.row(row).size() != cols)
      THROWEXPR ("Row #" << row << " has " << (int) path.row(row).size() << " columns; expected " << cols);
  for_const_contents (Annotation, gc_annot, gc)
    if ((int) gc->second.size() != cols)
      THROWEXPR ("Annotation '" << gc->first << "' has " << (int) gc->second.size() << " columns; expected " << cols);
  for_const_contents (Row_annotation, gr_annot, gr)
    for_const_contents (Annotation, gr->second, gr_row)
    if ((int) gr_row->second.size() != cols)
      THROWEXPR ("Annotation '" << gr_row->first << "' for sequence '" << gr->first
		 << "' has " << (int) gr_row->second.size() << " columns; expected " << cols);
}

double Stockholm::get_alignment_weight() const
{
  double w = 1.;
  const sstring wt = get_gf_annot (Stockholm_Dart_alignment_weight_tag);
  if (wt.size()) w = wt.to_double();
  return w;
}

sstring Stockholm::get_fold_string (const sstring& seqname) const
{
  sstring fold_string;
  const sstring ss_annot = get_gr_annot (seqname, sstring (Stockholm_secondary_structure_tag));
  if (ss_annot.size())
    {
      if ((int) ss_annot.size() != columns())
	THROWEXPR ("Secondary structure for '" << seqname << "' has wrong number of columns");
      const Phonebook::const_iterator row_iter = row_index.find (seqname);
      if (row_iter == row_index.end())
	{
	  THROWEXPR ("Alignment row index doesn't contain name '" << seqname << "'\n");
	}
      else
	{
	  fold_string.reserve (columns());
	  const int row = row_iter->second;
	  for (int col = 0; col < columns(); ++col)
	    if (path(row,col))
	      fold_string.push_back (ss_annot[col]);
	}
    }
  else
    CTAG(2,STOCKHOLM) << "No secondary structure annotation found for sequence '" << seqname << "'\n";
  return fold_string;
}

void Stockholm::set_fold_string (const sstring& seqname, const sstring& fold)
{
  if (row_index.find (seqname) == row_index.end())
    THROWEXPR ("While trying to set fold string: sequence '" << seqname << "' not in alignment");
  const int r = row_index[seqname];
  if (path.count_steps_in_row(r) != (int) fold.size())
    THROWEXPR ("Fold string '" << fold << "' is wrong length for alignment row '" << seqname << "'");
  sstring gapped_fold (columns());
  int x = 0;
  for (int c = 0; c < columns(); ++c)
    gapped_fold[c] = path(r,c) ? fold[x++] : Alignment::gap_char();
  gr_annot[seqname][sstring (Stockholm_secondary_structure_tag)] = gapped_fold;
}

void Stockholm::propagate_consensus_fold (bool override_row_folds)
{
  // get secondary structure
  const sstring ss_cons_tag = StockholmConsensus (Stockholm_secondary_structure_tag);
  const sstring ss_tag = Stockholm_secondary_structure_tag;
  Annotation::const_iterator ss_iter = gc_annot.find (ss_cons_tag);
  if (ss_iter == gc_annot.end())
    {
      CTAG(2,STOCKHOLM) << "Alignment has no '" << ss_cons_tag << "' tag, so can't propagate consensus structure\n";
      return;
    }
  const sstring& cons_ss = ss_iter->second;
  // map secondary structure to a set of column pairs
  typedef pair<int,int> Pair_index;
  set<Pair_index> pairs;
  set<int> singles;
  stack<int> left_pair_index;
  for (int pos = 0; pos < (int) cons_ss.size(); ++pos)
    {
      const char c = cons_ss[pos];
      if (Fold_char_enum::is_lchar (c))
	left_pair_index.push (pos);
      else if (Fold_char_enum::is_rchar (c))
	{
	  if (left_pair_index.empty())
	    THROWEXPR ("Too many >'s!");
	  pairs.insert (Pair_index (left_pair_index.top(), pos));
	  left_pair_index.pop();
	}
      else
	singles.insert (pos);
    }
  if (!left_pair_index.empty())
    THROWEXPR ("Too many <'s!");
  // loop through rows
  for (int row = 0; row < rows(); ++row)
    {
      const sstring& rname = row_name[row];
      Annotation& annot = gr_annot[rname];
      if (annot.find (ss_tag) != annot.end())
	{
	  if (override_row_folds)
	    CTAG(3,STOCKHOLM) << "Overwriting '" << ss_tag << "' tag for sequence '" << rname[row] << "'\n";
	  else
	    {
	      CTAG(3,STOCKHOLM) << "Sequence '" << rname[row] << "' already has a '" << ss_tag << "' tag; skipping\n";
	      continue;
	    }
	}
      sstring ss (cons_ss.size(), gap_char());
      for_const_contents (set<int>, singles, s)
	if (path(row,*s))
	  ss[*s] = cons_ss[*s];
      for_const_contents (set<Pair_index>, pairs, p)
	if (path(row,p->first) && path(row,p->second))
	  {
	    ss[p->first] = cons_ss[p->first];
	    ss[p->second] = cons_ss[p->second];
	  }
      annot[ss_tag] = ss;
    }
}

void Stockholm::make_consensus_fold()
{
  // get all fold strings
  vector<sstring> fold_seqname;
  vector<sstring> fold_string;
  for_const_contents (vector<sstring>, row_name, seqname)
    {
      const sstring ss_annot = get_gr_annot (*seqname, sstring (Stockholm_secondary_structure_tag));
      if (ss_annot.size())
	{
	  if ((int) ss_annot.size() != columns())
	    THROWEXPR ("Fold string for '" << *seqname << "' has " << ss_annot.size() << " columns, whereas alignment has " << columns());
	  fold_seqname.push_back (*seqname);
	  fold_string.push_back (ss_annot);
	}
    }

  // make consensus fold
  sstring cons (columns(), NOFOLD_CHAR);
  if (fold_string.size())  // if no fold strings, consensus fold is null
    {
      // liberally assign basepairs, checking for conflicts
      vector<int> paired_col (columns(), -1);
      vector<int> paired_source_row (columns(), -1);  // keep track of which rows contributed basepairs, so we can identify conflicts
      for (int i = 0; i < (int) fold_string.size(); ++i)
	{
	  stack<int> lcol_stack;
	  for (int col = 0; col < columns(); ++col)
	    {
	      // parse fold string character
	      const char c = fold_string[i][col];
	      if (is_lchar (c))  // fold string '<'
		lcol_stack.push (col);
	      else if (is_rchar (c))  // fold string '>'
		{
		  // check fold string validity
		  if (lcol_stack.empty())
		    THROWEXPR ("Too many >'s at column " << col << " of fold string for sequence '" << fold_seqname[i] << "'");
		  // pop paired column from stack
		  const int lcol = lcol_stack.top();
		  lcol_stack.pop();
		  // check for mismatches
		  if (paired_col[lcol] != -1 && paired_col[lcol] != col)
		    THROWEXPR ("Paired columns (" << lcol << "," << col << ") of sequence '" << fold_seqname[i] << " don't match paired columns (" << lcol << "," << paired_col[lcol] << ") of sequence '" << fold_seqname[paired_source_row[lcol]] << "'");
		  if (paired_col[col] != -1 && paired_col[col] != lcol)
		    THROWEXPR ("Paired columns (" << lcol << "," << col << ") of sequence '" << fold_seqname[i] << " don't match paired columns (" << paired_col[col] << "," << col << ") of sequence '" << fold_seqname[paired_source_row[col]] << "'");
		  // write paired columns and contributing row
		  paired_col[lcol] = col;
		  paired_col[col] = lcol;
		  paired_source_row[col] = paired_source_row[lcol] = i;
		}
	    }
	  // check fold string validity
	  if (!lcol_stack.empty())
	    THROWEXPR ("Too many <'s in fold string for sequence '" << fold_seqname[i] << "'");
	}

      // write consensus fold for all paired columns with no unpaired gaps
      for (int lcol = 0; lcol < columns(); ++lcol)
	{
	  const int rcol = paired_col[lcol];
	  if (rcol > lcol)
	    {
	      bool no_unpaired_gaps = TRUE;
	      for (int row = 0; row < rows(); ++row)
		if (is_gap (row, lcol) != is_gap (row, rcol))
		  {
		    no_unpaired_gaps = FALSE;
		    break;
		  }
	      if (no_unpaired_gaps)
		{
		  // check that more than 50% of sequences have basepairs in these two columns
		  int n_paired = 0;
		  for_const_contents (vector<sstring>, fold_string, fs)
		    if ((*fs)[lcol] == FOLD_LCHAR)
		      ++n_paired;
		  if (n_paired > (int) fold_string.size() / 2)
		    {
		      cons[lcol] = FOLD_LCHAR;
		      cons[rcol] = FOLD_RCHAR;
		    }
		}
	    }
	}
    }

  // store consensus fold
  gc_annot[StockholmConsensus (Stockholm_secondary_structure_tag)] = cons;
}

void Stockholm::add_gff (const GFF& gff)
{
  const sstring gff_tag (Stockholm_GFF_tag);
  sstring gff_val;
  gff_val << gff;
  gff_val.chomp();
  add_gf_annot (gff_tag, gff_val);
}

void Stockholm::add_gff (const GFF_list& gff_list)
{
  for_const_contents (GFF_list, gff_list, gff)
    add_gff (*gff);
}

sstring Stockholm::get_ID() const
{
  return get_gf_annot (sstring (Stockholm_identifier_tag));
}

sstring Stockholm::get_AC() const
{
  return get_gf_annot (sstring (Stockholm_accession_tag));
}

sstring Stockholm::get_name() const
{
  sstring name;
  name = get_ID();
  if (!name.size())
    name = get_AC();
  return name;
}

Stockade::Stockade (int rows, int cols)
  : align (rows, cols)
{
  Named_profile tmp_np;
  np = vector<Named_profile> (rows, tmp_np);
  for (int i = 0; i < (int) np.size(); ++i)
    align.np[i] = &np[i];
}

Stockade::Stockade (const Stockade& s)
  : align(), np()
{
  *this = s;
}

Stockade::Stockade (const Stockholm& s)
  : align (s)
{
  Named_profile tmp_np;
  np = vector<Named_profile> (s.rows(), tmp_np);
  for (int i = 0; i < s.rows(); ++i)
    {
      np[i] = Named_profile (*s.np[i]);
      align.np[i] = &np[i];
    }
}

Stockade& Stockade::operator= (const Stockade& s)
{
  align = s.align;
  np = s.np;
  for (int i = 0; i < (int) np.size(); ++i)
    align.np[i] = &np[i];
  return *this;
}

void Stockade::add_np (const Named_profile& prof)
{
  // add a Named_profile to the Stockade
  // then update all pointers in alignment, in case addresses changed when Named_profile was added
  np.push_back (prof);
  align.np.push_back (0);  // this gets updated below
  for (int i = 0; i < (int) np.size(); ++i)
    align.np[i] = &np[i];
  // update row_name, path & row_index
  align.row_name.push_back (prof.name);
  align.path.insert_rows (align.path.rows());
  align.update_row_index();
}

void Stockholm_database::read (istream& in, Sequence_database& seq_db, int max_line)		// pk accept long input lines
{
  int n_align = 0;
  while (in && !in.eof())
    {
      Stockholm stock;
      stock.read_Stockholm (in, seq_db, max_line);		// pk
      if (stock.rows())
	{
	  align.push_back (stock);
	  ++n_align;
	}
      else
	break;
    }
  update_index();
  CTAG(5,STOCKHOLM) << "Read " << n_align << " alignments\n";
}

void Stockholm_database::write (ostream& out) const
{
  for_const_contents (list<Stockholm>, align, a)
    a->write_Stockholm (out);
}

void Stockholm_database::read_Stockholm_or_FASTA (istream& in, Sequence_database& seq_db)
{
  if (Named_profile::detect_FASTA (in))
    {
      // read new sequences
      FASTA_sequence_database new_seqs;
      new_seqs.read_FASTA (in);
      Sequence_database::iterator first_new_seq = seq_db.append_sequence_database (new_seqs);
      // create a new alignment for every new sequence
      for_iterator (Sequence_database::iterator, seq, first_new_seq, seq_db.end())
	{
	  Stockholm stock (1, seq->size());
	  stock.row_name[0] = seq->name;
	  stock.prof[0] = &seq->prof_sc;
	  stock.np[0] = (Named_profile*) &*seq;  // cast away const
	  for (int col = 0; col < seq->size(); ++col)
	    stock.path[0][col] = 1;
	  align.push_back (stock);
	}
      // update index
      update_index();
    }
  else
    read (in, seq_db);  // calls update_index()
}

void Stockholm_database::update_index()
{
  align_row.clear();
  align_index.clear();
  Phonebook name_count;
  for_contents (list<Stockholm>, align, a)
    {
      const int nAlign = align_index.size();
      align_index.push_back (&*a);  // update align_index
      for (int nRow = 0; nRow < a->rows(); ++nRow)
	{
	  const sstring& row_name = a->row_name[nRow];
	  if (row_name.size())
	    {
	      if (align_row.find (row_name) == align_row.end())
		align_row[row_name] = Align_row (nAlign, nRow);  // update align_row
	      ++name_count[row_name];
	    }
	}
    }
  for_const_contents (Phonebook, name_count, nc)
    if (nc->second > 1)
      CTAG(8,STOCKHOLM) << "Warning: duplicate sequence name " << nc->first << " (occurs " << nc->second << " times)\n";
}

void Stockholm_database::propagate_consensus_folds (bool override_row_folds)
{
  for_contents (list<Stockholm>, align, stock)
    stock->propagate_consensus_fold (override_row_folds);
}

void Stockholm_database::add (const Stockholm& stock)
{
  align.push_back (stock);
  update_index();
}

void Stockholm_database::add (const Stockholm_database& stock_db)
{
  for_const_contents (list<Stockholm>, stock_db.align, stock)
    align.push_back (*stock);
  update_index();
}
