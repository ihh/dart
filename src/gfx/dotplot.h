#ifndef DOTPLOT_INCLUDED
#define DOTPLOT_INCLUDED

#include <gd.h>
#include <vector>
#include "util/array2d.h"

/** Palette is a palette for dotplots
 */
struct Palette
{
  enum { MaxColors = 256 };
  vector<int> r, g, b;
  Palette (int size = MaxColors) : r(size), g(size), b(size) { }
  const int size() const { return r.size(); }
};

/** Palette_vector is a singleton with various types of palette
 */
struct Palette_vector : vector<Palette>
{
  Palette_vector();
  const Palette& get_palette (int i) const { return (*this) [i % size()]; }
};
extern Palette_vector palette_vector;

/** Dotplot is a dotplot image with dimensions (L,M) where L and M are the lengths of sequences X and Y.
 */
struct Dotplot
{
  // GD image
  gdImagePtr image;
  // save methods
  void save (const char* filename) const;  // automatically adds ".png" suffix
  void save (const sstring& prefix, const sstring& xseq_name, const sstring& yseq_name) const;  // creates unique filename
  // constructor, destructor
  Dotplot (const array2d<double>& prob, const Palette& palette);
  Dotplot();
  ~Dotplot();
  // make image methods
  void make_image (const array2d<double>& prob, const Palette& palette);
};

#endif /* DOTPLOT_INCLUDED */
