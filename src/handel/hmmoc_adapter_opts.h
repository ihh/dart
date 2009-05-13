#ifndef HMMOC_ADAPTER_OPTS_INCLUDED
#define HMMOC_ADAPTER_OPTS_INCLUDED

#include "util/opts_list.h"
#include "seq/alignment.h"

// HMMoC adapter options
struct HMMoC_adapter_options
{
  // data
  sstring hmmoc_filename_prefix;
  bool try_to_use_hmmoc_adapter, strict_banded_reversibility;

  // constructor
  HMMoC_adapter_options();

  // method to hook this object up to an Opts_list
  void init_opts_list (Opts_list& opts, bool dump_hmmoc_file = true);

  // helpers: check whether a given alignment path is inside the relevant banding constraint
  static bool path_is_in_band (const Alignment_path& subpath, int width = 0);  // if width=0, then default (upper-limit) band will be used
  static bool cell_is_in_band (const Alignment_path::Sequence_coords& pos,
			       const Alignment_path::Sequence_coords& len,
			       int width = 0);  // if width=0, then default (upper-limit) band will be used

  // helpers more-or-less directly copied from relevant classes in hmmoc_banding.h
  static inline int diagonal0_1(const Alignment_path::Sequence_coords& pos,
				const Alignment_path::Sequence_coords& len,
				int inc = 0)
  { return ((pos[1] + inc) * len[0] + len[1]/2) / max (len[1], 1); } //return diagonal position in dim 0 given position in dim 1

  static inline int diagonal1_2(const Alignment_path::Sequence_coords& pos,
				const Alignment_path::Sequence_coords& len,
				int inc = 0)
  { return ((pos[2] + inc) * len[1] + len[2]/2) / max (len[2], 1); } //return diagonal position in dim 1 given position in dim 2

  static inline int diagonal2_3(const Alignment_path::Sequence_coords& pos,
				const Alignment_path::Sequence_coords& len,
				int inc = 0)
  { return ((pos[3] + inc) * len[2] + len[3]/2) / max (len[3], 1); } //return diagonal position in dim 2 given position in dim 3


};

#endif /* HMMOC_ADAPTER_OPTS_INCLUDED */
