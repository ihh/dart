#include "handel/transducer.h"
#include "handel/rivas.h"
#include "handel/alitrans.h"
#include "handel/recorder.h"

#include "util/rnd.h"
#include "util/logfile.h"

// default prefix for auto-assigned node names
#define DEFAULT_NODE_PREFIX "Node"

int main(int argc, char* argv[])
{
  // initialise handler
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <alignment>", "estimate complexity of Indelign ANNOTATE algorithm (Kim & Sunha, 2007) for a given Stockholm alignment");

  // params
  int kmin;
  bool print_distrib;

  // add options to handler
  opts.newline();
  opts.print_title ("Block-selection parameters");

  opts.add ("kmin", kmin = 0, "print all k mutually-dependent column 'path strings' of length k >= kmin", false);
  opts.add ("dist", print_distrib = true, "print distribution of k", false);

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // display initial logging messages
      Score_fns::describe_scoring_scheme (CLOG(8));

      // read in alignment
      Sequence_database seqs;
      Stockholm stock;

      ifstream align_file (opts.args[0].c_str());
      if (!align_file) THROW String_exception ("Alignment file not found: ", opts.args[0].c_str());
      CLOG(7) << "Reading initial alignment\n";

      stock.read_Stockholm (align_file, seqs);

      // find hole boundaries.
      // a hole boundary at col means that columns (col,col+1) have different gap profiles.
      vector<int> hole_boundary;
      for (int col = 0; col < stock.columns() - 1; ++col)
	for (int row = 0; row < stock.rows(); ++row)
	  if (stock.path (row, col) != stock.path (row, col+1))
	    {
	      hole_boundary.push_back (col);
	      break;
	    }

      // record distribution of k, the number of blocks in a mutually-dependent block sequence
      // (section 3.1.2, Kim & Sinha 2007)
      // and possibly print all path strings of length k >= kmin
      typedef map<int,int> IntDistrib;
      IntDistrib k_dist;
      int k = 0;

      const bool print_path = (kmin >= 1);
      Alignment_path path (stock.rows(), (int) 0);

      // count LH edge of alignment
      if (stock.columns())
	for (int row = 0; row < stock.rows(); ++row)
	  if (!stock.path (row, 0))
	    {
	      // inc k
	      k = 1;

	      // append block 0 to path string
	      if (print_path)
		path.append_column (stock.path.get_column (0));

	      break;
	    }

      // count edges at hole boundaries
      for_const_contents (vector<int>, hole_boundary, col)
	{
	  // flag if columns (*col,*col+1) on either side of a hole boundary have:
	  //  - a gap in the same row
	  //  - a gap in the right-hand column (i.e. in *col+1)
	  bool found_mutual_gap = false, next_has_gap = false;
	  for (int row = 0; row < stock.rows(); ++row)
	    if (!stock.path (row, *col + 1))
	      {
		next_has_gap = true;
		if (!stock.path (row, *col))
		  {
		    found_mutual_gap = true;
		    break;
		  }
	      }

	  // if there's a mutual gap, increment k and append column to the path string
	  if (found_mutual_gap)
	    {
	      // inc k
	      ++k;

	      // append column (*col + 1) to path string
	      if (print_path)
		path.append_column (stock.path.get_column (*col + 1));

	    }
	  else
	    {
	      // no mutual gap => end of k mutually dependent block sequence
	      if (k)
		{
		  // store k
		  ++k_dist[k];

		  // print path string
		  if (print_path && k >= kmin)
		    path.show (cout << '\n', stock.row_name);
		}

	      // reset k, path string
	      k = 0;
	      path.erase_columns (0, path.columns());
	      if (next_has_gap)
		{
		  // inc k
		  ++k;

		  // append column (*col + 1) to path string
		  if (print_path)
		    path.append_column (stock.path.get_column (*col + 1));
		}
	    }
	}

      // test for end of path string at RH edge of alignment
      if (k)
	{
	  // store k
	  ++k_dist[k];

	  // store path string
	  if (print_path && k >= kmin)
	    path.show (cout << '\n', stock.row_name);
	}

      // print distribution of k
      if (print_distrib)
	{
	  cout << "\n;k:\tFrequency:\n";
	  for_const_contents (IntDistrib, k_dist, k_count_ptr)
	    cout << k_count_ptr->first << '\t' << k_count_ptr->second << '\n';

	  // print complexity
	  IntDistrib::const_iterator last_k_count_ptr = k_dist.end();
	  --last_k_count_ptr;
	  cout << "\n;Complexity of Indelign-ANNOTATE is O(2^(k_max)) = O(2^" << last_k_count_ptr->first << ")\n\n";
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
