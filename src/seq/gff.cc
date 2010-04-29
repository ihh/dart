#include <algorithm>
#include "seq/gff.h"
#include "util/logfile.h"
#include "util/vector_output.h"

#define GFF_group_field_split_char  ';'
#define GFF_group_field_assign_char '='

Regexp GFF::regexp ("^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*$");

NSE::NSE()
  : seqname(), start (0), end (0)
{ }

NSE::NSE (const sstring& name, int start, int end)
  : seqname (name), start (start), end (end)
{ }

ostream& operator<< (ostream& outs, const GFF& gff)
{
  outs << gff.seqname << '\t' << gff.source << '\t' << gff.feature << '\t' << gff.start << '\t' << gff.end << '\t' << gff.score << '\t';
  outs << (gff.strand == GFF::NoStrand ? '.' : (gff.strand > 0 ? '+' : '-'));
  outs << '\t';
  if (gff.frame == GFF::NoFrame) outs << '.'; else outs << gff.frame;
  outs << '\t' << gff.group << '\n';
  return outs;
}

istream& operator>> (istream& ins, GFF& gff)
{
  sstring s;
  s.getline (ins);
  gff.initialise (s);
  return ins;
}

GFF::Strand GFF::string2strand (const sstring& strand_str)
{
  Strand strand;
  if (strand_str.size())
    switch (strand_str[0])
      {
      case '+':
	strand = PlusStrand;
	break;
      case '-':
	strand = MinusStrand;
	break;
      case '.':
	strand = NoStrand;
	break;
      default:
	CLOGERR << "Warning: illegal strand character: '" << strand_str << "'\n";
	strand = NoStrand;
	break;
      }
  else
    strand = NoStrand;
  return strand;
}

GFF::Frame GFF::string2frame (const sstring& frame_str)
{
  Frame frame;
  if (frame_str == sstring("") || frame_str == sstring("."))
    frame = NoFrame;
  else
    {
      frame = (Frame) atoi (frame_str.c_str());
      if (frame < 0 || frame > 2)
	{
	  CLOGERR << "Warning: illegal frame string: '" << frame_str << "'\n";
	  frame = NoFrame;
	}
    }
  return frame;
}

void GFF::initialise (const sstring& s)
{
  vector<sstring> f = s.split ("\t", 0, 9);
  if (f.size() < 9) THROWEXPR ("Not a GFF string");

  seqname = f[0];
  source  = f[1];
  feature = f[2];
  start   = atoi (f[3].c_str());
  end     = atoi (f[4].c_str());
  score   = atof (f[5].c_str());
  strand  = string2strand (f[6]);
  frame   = string2frame (f[7]);
  group   = f[8];
}

Regexp GFF::key_value_regexp ("^([^=]+)=(.*)$");

sstring GFF::get_value (const char* key) const
{
  const sstring GFF_group_field_split_string ((int) 1, (char) GFF_group_field_split_char);
  const sstring keystr (key);
  const vector<sstring> group_fields = group.split (GFF_group_field_split_string.c_str());
  for_const_contents (vector<sstring>, group_fields, kv)
    if (key_value_regexp.Match (kv->c_str()))
      if (key_value_regexp[1] == keystr)
	return key_value_regexp[2];
  return sstring();
}

map<sstring,sstring> GFF::get_key_value_map() const
{
  map<sstring,sstring> key_val;
  const sstring GFF_group_field_split_string ((int) 1, (char) GFF_group_field_split_char);
  const vector<sstring> group_fields = group.split (GFF_group_field_split_string.c_str());
  for_const_contents (vector<sstring>, group_fields, kv)
    if (key_value_regexp.Match (kv->c_str()))
      key_val[key_value_regexp[1]] = key_value_regexp[2];
  return key_val;
}

void GFF::set_value (const char* key, const char* val)
{
  const sstring GFF_group_field_split_string ((int) 1, (char) GFF_group_field_split_char);
  const sstring keystr (key);
  vector<sstring> group_fields = group.split (GFF_group_field_split_string.c_str());
  vector<sstring>::iterator kv  = group_fields.begin();
  while (kv < group_fields.end())
    {
      if (key_value_regexp.Match (kv->c_str()))
	if (key_value_regexp[1] == keystr)
	  {
	    kv->clear();
	    break;
	  }
      ++kv;
    }
  if (kv == group_fields.end())
    {
      group_fields.push_back (sstring());
      kv = group_fields.end() - 1;
    }
  *kv << key << GFF_group_field_assign_char << val;
  group = sstring::join (group_fields, GFF_group_field_split_string.c_str());
}

void GFF::set_values (const map<sstring,sstring>& key_val)
{
  group.clear();
  typedef map<sstring,sstring> ssmap;
  for_const_contents (ssmap, key_val, kv)
    group << kv->first << GFF_group_field_assign_char << kv->second << GFF_group_field_split_char;
  group.chomp (GFF_group_field_split_char);
}

const char* GFF::mask_tens_key = "mask";
const char* GFF::mask_units_key = "lmask";
const char* GFF::seq_key = "seq";

void GFF::get_mask (Metascore& meta, bool invert) const
{
  if ((int) meta.size() < end)
    {
      CLOGERR << "Warning --- expanding metascore vector\n";
      meta.insert (meta.end(), end - meta.size(), -InfinityScore);
    }
  sstring mask_tens_str = get_value (mask_tens_key);
  sstring mask_units_str = get_value (mask_units_key);
  if (mask_tens_str.size() == 0)
    CLOGERR << "Warning --- mask not found\n";
  else
    {
      if ((int) mask_tens_str.size() != length() || (int) mask_units_str.size() != length())
	THROWEXPR ("Length of retrieved mask attribute doesn't match length of GFF feature");
      for (int i = 0; i < length(); ++i)
	{
	  int tens = (int) (((unsigned char) mask_tens_str[i]) - '0');
	  int units = (int) (((unsigned char) mask_units_str[i]) - 'a');
	  int percent = tens * 10 + units;
	  if (tens < 0 || tens > 9 || units < 0 || units > (tens == 9 ? 10 : 9))
	    THROWEXPR ("While decoding mask attribute: bad percentage '" << mask_tens_str[i] << mask_units_str[i] << "'");
	  if (invert) percent = 100 - percent;
	  meta[start+i-1] = Prob2Score (((double) percent) / 100.0);
	}
    }
}

void GFF::set_mask (const Metascore& meta, bool invert)
{
  if ((int) meta.size() < end)
    THROWEXPR ("Supplied mask is shorter than GFF feature");
  sstring mask_tens_str (length());
  sstring mask_units_str (length());
  for (int i = 0; i < length(); ++i)
    {
      int percent = (int) (100.0 * Score2Prob (meta[start+i-1]) + .5);
      if (percent > 100) percent = 100;
      if (invert) percent = 100 - percent;
      unsigned char tens = (percent > 99 ? 9 : percent/10) + '0';
      unsigned char units = (percent > 99 ? 10 : percent%10) + 'a';
      mask_tens_str[i] = tens;
      mask_units_str[i] = units;
    }
  set_value (mask_tens_key, mask_tens_str.c_str());
  set_value (mask_units_key, mask_units_str.c_str());
}

ostream& operator<< (ostream& outs, const GFF_list& gff_list)
{
  for_const_contents (GFF_list, gff_list, gff) outs << *gff;
  return outs;
}

istream& operator>> (istream& in, GFF_list& gff_list)
{
  sstring s;
  while (in && !in.eof())
    {
      s.getline (in);
      if (GFF::regexp.Match (s.c_str()))
	{
	  GFF gff;
	  gff.initialise (s);
	  gff_list.push_back (gff);
	}
    }
  return in;
}

void GFF_list::load (const char* filename)
{
  ifstream infile (filename);
  if (!infile) THROWEXPR ("Couldn't open \"" << filename << "\" for reading");
  infile >> *this;
}

void GFF_list::save (const char* filename) const
{
  ofstream outfile (filename);
  if (!outfile) THROWEXPR ("Couldn't open \"" << filename << "\" for writing");
  outfile << *this;
}

void GFF_list::write_n (ostream& out) const
{
  out << size() << "\n";
  out << *this;
}

void GFF_list::read_n (istream& in)
{
  int n;
  in >> n;
  GFF gff;
  for (int i = 0; i < n; ++i)
    {
      in >> gff;
      push_back (gff);
    }
}

void GFF_list::acquire_mask (const Metaprob& meta, Prob min_prob, const char* seqname, const char* seq, const char* source, const char* feature, bool invert)
{
  Metascore meta_sc = Prob2ScoreVec (meta);
  // make GFFs for high-probability regions of the mask
  for (int start = 0; start < (int) meta.size(); ++start)
    {
      if (meta[start] < min_prob)
	continue;
      Prob total_pr = 0;
      int end;
      for (end = start; end < (int) meta.size(); ++end)
	if (meta[end] < min_prob)
	  break;
	else
	  total_pr += meta[end];
      
      GFF gff;
      gff.seqname = seqname;
      gff.source = source;
      gff.feature = feature;
      gff.start = start + 1;
      gff.end = end;
      gff.score = total_pr;
      gff.strand = GFF::PlusStrand;
      gff.frame = GFF::NoFrame;
      
      if (seq)
	{
	  const sstring seq_val (seq + start, (size_type) (end+1-start));
	  gff.set_value (GFF::seq_key, seq_val.c_str());
	}
      gff.set_mask (meta_sc, invert);
      push_back (gff);
      
      start = end;
    }
}

void GFF_list::apply_mask (Sequence_database_index& index, int mask_metascore_idx, bool invert) const
{
  for_const_contents (GFF_list, *this, gff)
    {
      Named_profile* np = index.name2profile (gff->seqname);
      Metascore meta (np->size(), 0);
      gff->get_mask (meta, invert);
      np->add_metascores (mask_metascore_idx, meta);
    }
}
