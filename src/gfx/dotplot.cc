#include <errno.h>
#include "gfx/dotplot.h"
#include "util/math_fn.h"
#include "util/sstring.h"
#include "util/logfile.h"

#define PNG_SUFFIX ".png"

Palette_vector::Palette_vector()
{
  // greyscale
  Palette greyscale;
  for (int i = 0; i < greyscale.size(); ++i)
    greyscale.r[i] = greyscale.g[i] = greyscale.b[i] = i;
  push_back (greyscale);
  // heat (red, yellow, white)
  Palette heat;
  for (int i = 0; i < heat.size(); ++i)
    {
      heat.r[i] = i;
      heat.g[i] = minmax (i*5/4-63, 0, 255);
      heat.b[i] = minmax (i*6/4-127, 0, 255);
    }
  push_back (heat);
}
Palette_vector palette_vector;

void Dotplot::save (const char* filename) const
{
  sstring name_with_suffix;
  name_with_suffix << filename << PNG_SUFFIX;
  CTAG(3,DOTPLOT) << "Saving dotplot\n";
  FILE* file = fopen (name_with_suffix.c_str(), "w");
  if (file == NULL)
    THROWEXPR ("Couldn't open file '" << name_with_suffix << "' for writing: " << strerror(errno) << "\n");
  gdImagePng (image, file);
  if (fclose (file) != 0)
    THROWEXPR ("Couldn't close file '" << name_with_suffix << "': " << strerror(errno) << "\n");
}

void Dotplot::save (const sstring& prefix, const sstring& xseq_name, const sstring& yseq_name) const
{
  sstring dotplot_name;   // create unique filename for this sequence pair
  dotplot_name << prefix << '-' << xseq_name << '-' << yseq_name;
  save (dotplot_name.c_str());
}

Dotplot::Dotplot() : image(0) { }
Dotplot::Dotplot (const array2d<double>& prob, const Palette& palette)
{
  make_image (prob, palette);
}

void Dotplot::make_image (const array2d<double>& prob, const Palette& palette)
{
  // create GD image
  image = gdImageCreate (prob.xsize(), prob.ysize());
  // allocate palette
  const int max_col = palette.size();
  vector<int> col;
  for (int i = 0; i < max_col; ++i)
    col.push_back (gdImageColorAllocate (image, palette.r[i], palette.g[i], palette.b[i]));
  // fill the dotplot
  for (int x = 0; x < prob.xsize(); ++x)
    for (int y = 0; y < prob.ysize(); ++y)
      {
	const int c = (int) (prob(x,y) * (double) max_col + .5);
	gdImageSetPixel (image, x, y, col[minmax(c,0,max_col-1)]);
      }
}

Dotplot::~Dotplot()
{
  gdImageDestroy (image);
}
