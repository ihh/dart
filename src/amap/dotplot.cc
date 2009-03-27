#include "dotplot.h"

Dotplot::Dotplot (const Biosequence& x, const Biosequence& y)
  : array2d<Prob> (x.size(), y.size(), 0.), xseq (x), yseq (y)
{ }

void Dotplot::write_dotplot (const sstring& filename) const
{
  ofstream file (filename.c_str());
  if (!file)
    THROWEXPR ("Couldn't create dotplot file with name '" << filename << "'");
  // print horizontal sequence axis labels
  file << ".";
  for (int x = 0; x < (int) xseq.size(); ++x)
    file << ' ' << xseq[x];
  file << '\n';
  // print rows
  for (int y = 0; y < (int) yseq.size(); ++y) {
    file << yseq[y];
    for (int x = 0; x < (int) xseq.size(); ++x)
      file << ' ' << (*this) (x,y);
    file << '\n';
  }
	
  file.close();
}


void Dotplot::write_dotplot (const sstring& prefix, const sstring& seqname) const
{
  sstring filename;
  filename << prefix << '-' << seqname;
  write_dotplot (filename);
}

void Dotplot::write_dotplot (const sstring& prefix, const sstring& xseqname, const sstring& yseqname) const
{
  sstring filename;
  filename << prefix << '-' << xseqname << '-' << yseqname;
  write_dotplot(filename);
}

