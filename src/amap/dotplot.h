#ifndef DOTPLOT_INCLUDED
#define DOTPLOT_INCLUDED


#include "util/strsaver.h"
#include "util/array2d.h"
#include "seq/biosequence.h"
#include "util/score.h"

/// Base class for probability matrices.
/*
 * Both alignment and fold (for e.g. Pair SCFGs) dotplots inherit from this class.
 */
struct Dotplot : array2d<Prob> {

  /// axis labels
  Biosequence xseq;
  Biosequence yseq;

  /// constructor
  Dotplot() { }
  Dotplot (const Biosequence& x, const Biosequence& y);
  Dotplot (const Biosequence& x, const Biosequence& y, const sstring xseqname, const sstring yseqname);

  /// output methods
  /*
   * For e.g. alignment dotplots, output is formatted as:
   * <cruft>
   * .   x1   x2   x3
   * y1  .2   .2   .2
   * y2  .2   .2   .2
   */
  void write_dotplot (const sstring& filename) const;

  /// (For writing fold dotplots.)
  void write_dotplot (const sstring& prefix, const sstring& seqname) const;

  /// (For writing alignment dotplots.)
  void write_dotplot (const sstring& prefix, const sstring& xseqname, const sstring& yseqname) const;

};

#endif /*DOTPLOT_INCLUDED*/	
