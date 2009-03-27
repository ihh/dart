#ifndef PAIRENV_INCLUDED
#define PAIRENV_INCLUDED

#include "util/array2d.h"
#include "seq/alignment.h"

// primary sequence pair envelope
struct Pair_envelope
{
  // data
  // allow_cut(i,j)=TRUE if can make a "cut" after residues [1..i] of X and [1..j] of Y
  // NB allow_cut.xsize() == xlen+1 and allow_cut.ysize() == ylen+1
  array2d<int> allow_cut;
  // constructors
  Pair_envelope();
  Pair_envelope (int xlen, int ylen, int init_val = 1);
  Pair_envelope (const Pairwise_path& path);
  // assignment operator
  Pair_envelope& operator= (const Pair_envelope& pair_env);
  // builders
  void initialise_full (int xlen, int ylen, int init_val = 1);
  void initialise_banded (int xlen, int ylen, int band_size);  // band_size = size of alignment band (-1 to unlimit)
  void add_pairwise_path (const Pairwise_path& path, bool allow = TRUE);  // if allow=FALSE, then path is masked out
  // initialise_from_posterior_matrix() assumes that post.xsize() == xlen and post.ysize() == ylen
  void initialise_from_posterior_matrix (const array2d<Prob>& post, Prob min_prob = 0.5);
  // helper method to count allowed cells
  int cells() const;
  // output methods
  friend ostream& operator<< (ostream& out, const Pair_envelope& pair_env);
  void render_dotplot (ostream& out, const Biosequence& xseq, const Biosequence& yseq, int max_level = 1, bool use_ansi_color = TRUE) const;
};

#endif /* PAIRENV_INCLUDED */
