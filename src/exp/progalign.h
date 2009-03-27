#ifndef PROGALIGN_INCLUDED
#define PROGALIGN_INCLUDED

#include "tree/tree_alignment.h"
#include "hmm/pairphmm.h"
#include "seq/dirichlet.h"

// Progressive_aligner
// Abstract base class for doing progressive alignment & iterative refinement.

struct Progressive_aligner
{
  // typedefs
  typedef Alignment_path::Row           Row;
  typedef Alignment_path::Row_pair      Row_pair;
  typedef Alignment_path::Decomposition Decomposition;
  typedef Alignment_path::Row_index_set Row_index_set;
  typedef Phylogeny::Node               Node;
  // data
  PHYLIP_tree tree;
  Alignment_path::Decomposition decomp;
  // concrete methods
  Tree_alignment make_progressive_alignment (const Sequence_database_index& seqdb_index);
  // virtual methods
  // parameterisation methods
  virtual void init_default_params() = 0;
  virtual void load_params (const char* filename) = 0;
  virtual void train_params (const char* filename) = 0;
  virtual void train_params_unaligned (const Sequence_database_index& seqdb_index) = 0;
  virtual void save_params (const char* filename) = 0;
  // alignment methods
  virtual void initialise (const Sequence_database_index& seqdb_index) = 0;
  virtual void make_guide_tree() = 0;
  virtual Alignment_path align_axy (int a, int x, int y, double tx, double ty) = 0;
  virtual void refine() = 0;   // updates decomp
  // destructor
  virtual ~Progressive_aligner();
};

// Optimal_accuracy_progressive_aligner
// Implements maximum expected-SPS multiple alignment using the algorithm of Do, Brudno and Batzoglou.

struct Optimal_accuracy_progressive_aligner : Progressive_aligner
{
  // typedefs
  typedef array2d<Prob> PostMat;

  // model data
  PScores& pscore;   // taken from prior
  Pair_PHMM& phmm;
  Dirichlet_prior& prior;

  // sequence info
  vector<Named_profile*> profile;
  vector<sstring> leaf_name;
  vector<int> leaf_seqlen;
  int leaves, nodes;

  // expected accuracy info
  array2d<PostMat*> post;
  array2d<double> expected_accuracy;  // create with enough space for all nodes in guide tree

  // constructor
  Optimal_accuracy_progressive_aligner (Pair_PHMM& phmm, Dirichlet_prior& prior);
  // destructor
  ~Optimal_accuracy_progressive_aligner();
  // method to delete PostMat's
  void delete_post();

  // Progressive_aligner methods
  void init_default_params();
  void load_params (const char* filename);
  void train_params (const char* filename);
  void train_params_unaligned (const Sequence_database_index& seqdb_index);
  void save_params (const char* filename);
  // alignment methods
  void initialise (const Sequence_database_index& seqdb_index);
  void make_guide_tree();
  Alignment_path align_axy (int a, int x, int y, double tx, double ty);
  void refine();

  // internal methods
  void fill_post();
  void do_consistency_transformation();
  // method to align two subalignments, maximising expected sum-of-pairs score
  Pairwise_path align_subpaths (const Alignment_path& xpath, const Row_index_set& xrows,
				const Alignment_path& ypath, const Row_index_set& yrows);
};

#endif /* PROGALIGN_INCLUDED */
