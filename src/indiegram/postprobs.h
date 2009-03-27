#ifndef POSTPROBS_INCLUDED
#define POSTPROBS_INCLUDED

#include "util/array2d.h"
#include "indiegram/tripletscfgdp.h"

struct Triplet_inside_outside_matrix
{
  // data
  Triplet_inside_matrix inside;
  Triplet_outside_matrix outside;

  Triplet_inside_outside_matrix (const Triplet_SCFG& scfg,
				 const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
				 const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z,
				 bool fill_now = true);

  void show (ostream& o) const;

  // Posterior score calculation methods.
  // These give the posterior probability that a parse tree contains a subtree rooted in state DEST,
  // where the subtree from DEST downwards has emitted XSUBSEQ and YSUBSEQ (including the emission from DEST).
  // (In other words, the subseq-indexing convention here follows Inside/CYK, not Outside/KYC.)
  // In the case of post_transition_sc(), it's also required that the parent node of DEST must be SRC.
  Score post_transition_sc (int src_state, int dest_state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const;
  Score post_state_sc (int dest_state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const;  // = psum_{src} post_transition_sc(src,...)
};

struct Triplet_postprobs : Triplet_inside_outside_matrix
{

  


};


struct Triplet_fold_dotplot : array2d<Prob>, SCFG_state_type_enum {

  /// axis labels
  const Biosequence& seq;

  /// constructor
  /*
   * Index: 0 for X, 1 for Y, 2 for Z
   */ 
  Triplet_fold_dotplot (const Triplet_inside_outside_matrix& in_out, int seqindex);

  /// output methods
  void write_dotplot (const sstring& filename) const;
  void write_dotplot (const sstring& prefix, const sstring& seqname) const;
};

#endif /* POSTPROBS_INCLUDED */
