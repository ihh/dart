#include "hmm/pairenv.h"

Pair_envelope::Pair_envelope()
{ }

Pair_envelope::Pair_envelope (int xlen, int ylen, int init_val)
{
  initialise_full (xlen, ylen, init_val);
}

Pair_envelope& Pair_envelope::operator= (const Pair_envelope& pair_env)
{
  allow_cut = pair_env.allow_cut;
  return *this;
}

void Pair_envelope::initialise_full (int xlen, int ylen, int init_val)
{
  if (init_val)
    CTAG(6,PAIRENV) << "Initialising full " << xlen+1 << "*" << ylen+1 << " pairwise alignment envelope\n";
  else
    CTAG(5,PAIRENV) << "Creating empty " << xlen+1 << "*" << ylen+1 << " pairwise alignment envelope\n";
  // do the alloc
  allow_cut = array2d<int> (xlen+1, ylen+1, init_val);
}

void Pair_envelope::initialise_banded (int xlen, int ylen, int band_size)
{
  if (band_size < 0)
    initialise_full (xlen, ylen, 1);
  else
    {
      CTAG(6,PAIRENV) << "Initialising banded " << xlen+1 << "*" << ylen+1 << " pairwise alignment envelope, max deviation from diagonal " << band_size << "\n";
      initialise_full (xlen, ylen, 0);
      const double k = ((double) ylen) / (double) xlen;
      const double c = 0.5 - ((double) band_size) / 2.;
      for (int i = 0; i <= xlen; ++i)
	{
	  const int j_min = (int) (c + k * (double) i);
	  for (int j = max(j_min,0); j <= min(j_min+band_size, ylen); ++j)
	    allow_cut(i,j) = 1;
	}
    }
}

void Pair_envelope::initialise_from_posterior_matrix (const array2d<Prob>& post, Prob min_prob)
{
  initialise_full (post.xsize(), post.ysize(), 0);
  const double ratio = 1. / min_prob;
  for (int x = 0; x < post.xsize(); ++x)
    for (int y = 0; y < post.ysize(); ++y)
      allow_cut(x+1,y+1) = (int) (ratio * post(x,y));
  allow_cut(0,0) = 1;
  for (int x = 0; x < post.xsize(); ++x)
    allow_cut(x+1,0) = 1;
  for (int y = 0; y < post.ysize(); ++y)
    allow_cut(0,y+1) = 1;
}

Pair_envelope::Pair_envelope (const Pairwise_path& path)
  : allow_cut (path.count_steps_in_row(0) + 1, path.count_steps_in_row(1) + 1, 0)
{
  add_pairwise_path (path, 1);
}

void Pair_envelope::add_pairwise_path (const Pairwise_path& path, bool allow)
{
  if (allow)
    CTAG(4,PAIRENV) << "Adding pairwise path to alignment envelope\n";
  else
    CTAG(4,PAIRENV) << "Subtracting pairwise path from alignment envelope\n";
  int x = 0;
  int y = 0;
  int i = 0;
  while (1)
    {
      if (allow)
	++allow_cut(x,y);
      else
	allow_cut(x,y) = 0;
      if (i == path.columns())
	break;
      x += path(0,i);
      y += path(1,i);
      ++i;
    }
}

int Pair_envelope::cells() const
{
  int c = 0;
  for_const_contents (array2d<int>, allow_cut, cell)
    if (*cell)
      ++c;
  return c;
}

ostream& operator<< (ostream& out, const Pair_envelope& pair_env)
{
  for (int y = 0; y < pair_env.allow_cut.ysize(); ++y)
    {
      for (int x = 0; x < pair_env.allow_cut.xsize(); ++x)
	out << (pair_env.allow_cut(x,y) ? 'X' : '.');
      out << "\n";
    }
  return out;
}

void Pair_envelope::render_dotplot (ostream& out, const Biosequence& xseq, const Biosequence& yseq, int max_level, bool use_ansi_color) const
{
  if ((int) xseq.size() != allow_cut.xsize()-1 || (int) yseq.size() != allow_cut.ysize()-1)
    THROWEXPR ("While rendering dotplot: sequence lengths (" << xseq.size() << "," << yseq.size() << ") don't match dotplot dimensions (" << allow_cut.xsize()-1 << "," << allow_cut.ysize()-1 << ")");
  // prepare chars
  const int cols = 8;
  sstring esc (2);
  esc[0] = 27;
  esc[1] = '[';
  const double scale = (double) (cols-1) / (double) max_level;
  sstring inverse_text;
  sstring normal_text;
  vector<sstring> color_text (cols);
  vector<sstring> cell_text (cols);
  cell_text[0] = " ";
  if (use_ansi_color)
    {
      const int colseq[] = { 31, 33, 32, 34, 36, 35, 37 };  // rainbow
      for (int k = 1; k < cols; ++k)
	{
	  color_text[k] << esc << colseq[k-1] << 'm';
	  cell_text[k] << ' ';
	}
      inverse_text << esc << "7m";
      normal_text << esc << "27m";
    }
  else
    for (int k = 1; k < cols; ++k)
      cell_text[k] << k;
  // render, with axes
  out << " *" << xseq << '\n';
  for (int j = 0; j < allow_cut.ysize(); ++j)
    {
      out << (char) (j==0 ? '*' : yseq[j-1]);
      int last_col = -1;
      for (int i = 0; i < allow_cut.xsize(); ++i)
	{
	  // determine color
	  const int val = allow_cut(i,j);
	  const int col = val<=0 ? 0 : min ((int) ((scale * (double) val) + 1), (int) (cols-1));
	  // render cell
	  if (i == 0 || col != last_col)
	    {
	      if (last_col >= 0) out << color_text[last_col];
	      if (last_col <= 0 && col > 0) out << normal_text << inverse_text;
	      else if (last_col > 0 && col <= 0) out << inverse_text << normal_text;
	      out << color_text[col];
	      last_col = col;
	    }
	  out << cell_text[col];
	}
      if (last_col >= 0) out << color_text[last_col] << inverse_text << normal_text;
      out << (char) (j==0 ? '*' : yseq[j-1]);
      out << '\n';
    }
  out << " *" << xseq << '\n';
  out << "Color code: [0] " << normal_text << inverse_text;
  for (int i = 1; i < cols; ++i)
    out << color_text[i] << ' ' << color_text[i];
  out << inverse_text << normal_text << "[" << max_level << "]\n";
}
