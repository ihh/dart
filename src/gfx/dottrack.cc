#include "gfx/dotplot.h"
#include "util/logfile.h"
#include "seq/biosequence.h"

// subroutine to get average pixel intensity
// (x,y) are the pixel co-ordinates on the (untiled) track
double avg_intensity (const Digitized_biosequence& dsq,
		      const array2d<Prob>& submat,
		      int x,
		      int y,
		      int bases_per_pixel)
{
  // rotate
  x -= y/2;

  // count a square array of pixels
  double psum = 0., n = 0.;
  const int seqlen = dsq.size();
  for (int i = max (x * bases_per_pixel, 0); i < (x+1) * bases_per_pixel && i < seqlen; ++i)
    {
      const int ci = dsq[i];
      for (int j = y * bases_per_pixel; j < (y+1) * bases_per_pixel; ++j)
	if (i + j < seqlen)
	  {
	    const double p = submat (ci, dsq[i + j]);
	    psum += p;
	    n++;
	  }
    }

  // average
  return n > 0 ? psum / n : 0.;
}

// main program
int main(int argc, char* argv[])
{
  // initialise the options parser
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <sequence database in FASTA format>",
		  "make a dotplot track from a sequence\n");


  int distance;
  int bases_per_pixel;
  int tile_width;
  bool complement;
  sstring tile_prefix;

  opts.print_title ("Dotplot options");
  opts.add ("d -distance", distance = 100, "max distance between bases; d/b is height of dotplot");
  opts.add ("b -bases", bases_per_pixel = 1, "bases per pixel");
  opts.add ("w -width", tile_width = 1000, "tile width in pixels");
  opts.add ("c -complement", complement = false, "look for complement, rather than identity");
  opts.add ("f -filename", tile_prefix = "TILE", "tile filename prefix");

  // parse the command line
  //
  try
    {
      if (!opts.parse()) { cerr << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      cerr << opts.short_help();
      cerr << e.what();
      exit(1);
    }

  try
    {
      // get args
      const char* seq_db_filename = opts.args[0].c_str();

      // initialise database
      const Alphabet& alph (DNA_alphabet);
      FASTA_sequence_database seq_db (seq_db_filename, &alph, Profile_flags_enum::DSQ);
      if (seq_db.size() != 1)
	THROWEXPR ("Error: expected exactly 1 sequence in database");
      const Digitized_biosequence& dsq (seq_db.get_seq(0).dsq);

      // get dimensions of track
      const int seqlen = dsq.size();

      int track_xmax = seqlen / bases_per_pixel;
      if (seqlen % bases_per_pixel)
	++track_xmax;

      int track_ymax = distance / bases_per_pixel;
      if (distance % bases_per_pixel)
	++track_ymax;
 
      int n_tiles = track_xmax / tile_width;
      if (track_xmax % tile_width)
	++n_tiles;

      // create substitution matrix
      array2d<Prob> submat (alph.size(), alph.size(), 0.);
      for (int i = 0; i < alph.size(); ++i)
	submat (i, complement ? alph.complement(i) : i) = 1.;

      // get mean and variance of average intensities
      CTAG(7,DOTTRACK) << "Computing mean and variance of pixel intensities\n";
      double isum = 0., i2sum = 0., n = 0.;
      for (int x = 0; x < track_xmax; ++x)
	for (int y = 0; y < track_ymax; ++y)
	  {
	    const double p = avg_intensity (dsq, submat, x, y, bases_per_pixel);
	    isum += p;
	    i2sum += p*p;
	    n++;
	  }
      const double imean = n > 0 ? isum / n : 0.;
      const double isd = n > 0 ? sqrt (i2sum / n - imean * imean) : 1.;
      CTAG(7,DOTTRACK) << "Mean = " << imean << ", SD = " << isd << "\n";

      // make palette
      const Palette& grayscale = palette_vector[0];

      // loop over tiles
      CTAG(7,DOTTRACK) << "Rendering tiles\n";
      for (int tile = 0; tile < n_tiles; ++tile)
	{
	  // print log message
	  CTAG(6,DOTTRACK) << "Rendering tile #" << tile+1 << " of " << n_tiles << "\n";

	  // get tile X offset
	  const int x_offset = tile * tile_width;

	  // create prob array
	  array2d<double> prob (tile_width, track_ymax);
	  for (int x = 0; x < tile_width; ++x)
	    for (int y = 0; y < track_ymax; ++y)
	      {
		double p = avg_intensity (dsq, submat, x + x_offset, y, bases_per_pixel);
		p = 0.5 + (p - imean) / (2*isd);
		p = minmax (p, 0., 1.);
		prob(x,y) = p;
	      }

	  // create and save dotplot tile
	  CTAG(3,DOTTRACK) << "Saving tile #" << tile+1 << "\n";
	  Dotplot dotplot (prob, grayscale);
	  sstring tile_filename;
	  tile_filename << tile_prefix << tile+1;
	  dotplot.save (tile_filename.c_str());
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << "ERROR: " << e.what();
      exit(1);
    }

  return 0;
}
