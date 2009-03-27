#ifndef ALITRANS_INCLUDED
#define ALITRANS_INCLUDED

#include "handel/transducer.h"
#include "handel/handelalign.h"
#include "handel/recorder.h"
#include "handel/hmmoc_adapter_opts.h"

// default banding coefficient -- estimated using t/simband.pl
#define DEFAULT_BANDING_COEFFICIENT 10.

// Pair transducer factory: the base class for the transducer MCMC aligner
struct Pair_transducer_factory : Grammar_state_enum
{
  // dumps
  Transducer_dotfile_recorder dotfile_recorder;  /* currently unused */
  Composition_recorder composition_recorder;

  // virtual factory methods
  virtual Pair_transducer_scores prior_pair_trans_sc()
  = 0;  // returns the jointly-normalised "prior" transducer
  virtual Pair_transducer_scores branch_pair_trans_sc (double time)
  = 0;  // returns conditionally-normalised "branch" transducer

  // virtual banding method
  virtual double gap_rate() = 0;
  virtual double mean_gap_size() = 0;

  // virtual destructor
  virtual ~Pair_transducer_factory() { }
};

// Transducer-based Handel MCMC aligner & phylo-sampler
struct Transducer_alignment
  : Handel_base, Pair_transducer_factory, Transducer_state_type_enum
{
  // data
  // flag to indicate whether to use Redelings-Suchard proposal scheme
  bool use_Redelings_Suchard;

  // banding coefficient, and flag to indicate whether to use it
  // phylocomposer's banding coefficient (B) is related to this banding coefficient (C) by the following formula:
  // B = C * mean_deletion_size * deletion_opening_probability
  double banding_coefficient;
  bool use_banding_coefficient;

  // HMMoC adapter options
  HMMoC_adapter_options hmmoc_opts;

  // methods
  // constructor
  Transducer_alignment()
    : use_Redelings_Suchard (false),
      banding_coefficient (DEFAULT_BANDING_COEFFICIENT),
      use_banding_coefficient (false)
  { }

  // resample part of a multiple alignment
  // subtree node numbers are w.r.t. nodes_to_sample
  // if available_subtrees contains two trees, then the first should be the "old" tree (for Hastings ratio calculation purposes)
  void sample_subtree (const vector<int>& nodes_to_sample,
		       const vector<PHYLIP_tree>& available_subtrees,
		       double kT = 1.,
		       bool viterbi = false);

  // method to dump entire transducer composition
  void dump (ostream& out);

  // method to prepare the Handel_movement, shared by sample_subtree and dump
  Handel_movement prepare_movement (const PHYLIP_tree& tree, const vector<int>& constrained_phylip_branch, const vector<int>& topologically_identical_phylip_branch, vector<int>& etree2phylip, Alignment_path::Decomposition& decomp);

  // inherited virtual methods
  // return conditional alignment path score for a branch
  Score conditioned_branch_path_score (const Node_pair& branch) const;

  // common implementations of sample/optimise methods
  bool propose_sample_or_optimise_node (Node node, double kT, bool optimise);
  bool propose_sample_or_optimise_branch (const Undirected_pair& branch, double kT, bool optimise);

  // optimise methods
  bool propose_optimise_node (Node node);
  bool propose_optimise_branch (const Undirected_pair& branch);
  void propose_optimise_branch_length (const Undirected_pair& branch,
				       double tmax = 10., double tres = .01);  // UNIMPLEMENTED

  // sample methods
  void propose_sample_node (Node node, double kT = 1.);
  void propose_sample_branch (const Undirected_pair& branch, double kT = 1.);
  bool propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT = 1);

  // helper to return an SExpr-safe unique tape name
  static sstring safe_tape_name (const PHYLIP_tree& tree, int node);

  // align_and_infer_parent; obsolete?
  // (TO BE IMPLEMENTED, or bypassed)
  void align_and_infer_parent (const Score_profile& xseq,
			       const Score_profile& yseq,
			       Node ancestor_node,
			       Alignment_path& axy_path);  // UNIMPLEMENTED


};

#endif /* ALITRANS_INCLUDED */
