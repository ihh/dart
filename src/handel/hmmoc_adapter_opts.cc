#include "handel/hmmoc_adapter_opts.h"
#include "handel/hmmoc_adapter.h"

HMMoC_adapter_options::HMMoC_adapter_options()
  : hmmoc_filename_prefix(),
    try_to_use_hmmoc_adapter (false),
    strict_banded_reversibility (true)
{ }

void HMMoC_adapter_options::init_opts_list (Opts_list& opts, bool dump_hmmoc_file)
{
#if defined(HMMOC_INCLUDED) && HMMOC_INCLUDED
  try_to_use_hmmoc_adapter = true;
#else
  CLOGERR << "Warning: without HMMoC installed, tree-sampling MCMC moves that use dynamic programming may be very slow.\n";
#endif

  opts.newline();
  opts.print_title ("HMMoC adapter options");

  if (dump_hmmoc_file)
    opts.add ("hml -hmmoc-model", hmmoc_filename_prefix, "save acyclic composite transducer as HMMoC model file", false);

  opts.add ("ha -hmmoc-adapter", try_to_use_hmmoc_adapter, "use HMMoC adapter to do DP wherever possible");
  opts.add ("hr -hmmoc-root", HMMoC_adapter::root_directory, "root directory of HMMoC installation");
  opts.add ("ht -hmmoc-tmp", HMMoC_adapter::tmp_directory, "temporary directory");
  opts.add ("hc -hmmoc-cache", HMMoC_adapter::cache_filename, "specify location of persistent compiled-binary index file", false);
  opts.add ("hk -hmmoc-keep-tmp", HMMoC_adapter::leave_debris, "don't delete temporary files");
  opts.add ("hcx -hmmoc-cc-exec", HMMoC_adapter::gcc_exec, "path to C++ compiler executable");
  opts.add ("hca -hmmoc-cc-args", HMMoC_adapter::gcc_args, "C++ compiler options");

  opts.newline();
  opts.print_title ("HMMoC adapter banding");

  opts.add ("hc2d -hmmoc-2d-banding", HMMoC_adapter::two_d_banding, "2-D DP banding width, 0 for no banding");
  opts.add ("hc3d -hmmoc-3d-banding", HMMoC_adapter::three_d_banding, "3-D DP banding diameter, 0 for no banding");
  opts.add ("hc4d -hmmoc-4d-banding", HMMoC_adapter::four_d_banding, "4-D DP banding diameter, 0 for no banding");
  opts.add ("hrev -hmmoc-strict-reversible", strict_banded_reversibility, "reject MCMC moves that start outside the band", true, "hirrev");
}

bool HMMoC_adapter_options::path_is_in_band (const Alignment_path& subpath, int width)
{
  // figure out size of subpath
  const Alignment_path::Sequence_coords len = subpath.seq_coords_end();

  // check that every cell is in-band
  bool in_band = true;
  Alignment_path::Sequence_coords pos = subpath.create_seq_coords();
  for (int col = 0; in_band && col < subpath.columns(); ++col)
    {
      // check cell
      in_band = in_band && cell_is_in_band (pos, len, width);

      // increment pos
      if (in_band)
	subpath.inc_seq_coords (pos, col);
    }

  // log
  if (in_band)
    CTAG(2,HMMOC_ADAPTER) << "Alignment path is inside upper-limit DP banding constraint\n";
  else
    CTAG(3,HMMOC_ADAPTER) << "Path contains at least one cell (" << pos << ") that is outside the upper-limit DP banding constraint\n";

  // return
  return in_band;
}

bool HMMoC_adapter_options::cell_is_in_band (const Alignment_path::Sequence_coords& pos,
					     const Alignment_path::Sequence_coords& len,
					     int width)
{
  bool in_band = true;  // return true by default
  int diag0_max, diag0_min, diag1_max, diag1_min, diag2_max, diag2_min;
  switch (len.size())
    {
      // 2D banding
    case 2:
      if (width <= 0)
	width = HMMoC_adapter::two_d_banding;

      if (width > 0)
	{
	  diag0_max = (width-1)/2 + max ( diagonal0_1(pos,len,1)-1, diagonal0_1(pos,len) );
	  diag0_min = diagonal0_1(pos,len) - width/2;

	  in_band =
	    pos[1] >= 0
	    && pos[1] <= len[1]

	    && pos[0] >= max (diag0_min, 0)
	    && pos[0] <= min (diag0_max, len[0]);
	}
      break;

      // 3D banding
    case 3:
      if (width <= 0)
	width = HMMoC_adapter::three_d_banding;

      if (width > 0)
	{
	  diag1_max = (width-1)/2 + max ( diagonal1_2(pos,len,1)-1, diagonal1_2(pos,len) );
	  diag1_min = diagonal1_2(pos,len) - width/2;

	  diag0_max = (width-1)/2 + max ( diagonal0_1(pos,len,1)-1, diagonal0_1(pos,len) );
	  diag0_min = diagonal0_1(pos,len) - width/2;

	  in_band =
	    pos[2] >= 0
	    && pos[2] <= len[2]

	    && pos[1] >= max (diag1_min, 0)
	    && pos[1] <= min (diag1_max, len[1])

	    && pos[0] >= max (diag0_min, 0)
	    && pos[0] <= min (diag0_max, len[0]);
	}
      break;

      // 4D banding
    case 4:
      if (width <= 0)
	width = HMMoC_adapter::four_d_banding;

      if (width > 0)
	{
	  diag2_max = (width-1)/2 + max ( diagonal2_3(pos,len,1)-1, diagonal2_3(pos,len) );
	  diag2_min = diagonal2_3(pos,len) - width/2;

	  diag1_max = (width-1)/2 + max ( diagonal1_2(pos,len,1)-1, diagonal1_2(pos,len) );
	  diag1_min = diagonal1_2(pos,len) - width/2;

	  diag0_max = (width-1)/2 + max ( diagonal0_1(pos,len,1)-1, diagonal0_1(pos,len) );
	  diag0_min = diagonal0_1(pos,len) - width/2;

	  in_band =
	    pos[3] >= 0
	    && pos[3] <= len[3]

	    && pos[2] >= max (diag2_min, 0)
	    && pos[2] <= min (diag2_max, len[2])

	    && pos[1] >= max (diag1_min, 0)
	    && pos[1] <= min (diag1_max, len[1])

	    && pos[0] >= max (diag0_min, 0)
	    && pos[0] <= min (diag0_max, len[0]);
	}
      break;

    default:
      break;
    }

  // return
  return in_band;
}
