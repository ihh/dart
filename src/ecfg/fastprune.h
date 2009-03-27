#ifndef FAST_PRUNE_INCLUDED
#define FAST_PRUNE_INCLUDED

#include "hsm/em_matrix_base.h"

struct Fast_prune
{
  // data
  // tree topology
  int root, nodes, alph_size;
  vector<pair<int,int> > branch;  // each entry is a (parent,child) node pair
  // root probabilities
  vector<vector<Prob> > pi;  // pi[n][i] = P(node n is in state i | node n is clique root)
  // edge transition probability matrices
  vector<array2d<Prob> > Q;  // Q[N](i,j) = P(node N is in state j | parent is in state i)
  // pointer to Column_matrix, for "signpost" arrays (gapped[], allowed[], and root[])
  const Column_matrix* colmat;
  // dynamic programming matrix
  vector<vector<Prob> > F;  // F[node N][state X] = P(observed seqs descended from node N | node N is in state X)

  // methods
  void prepare (const PHYLIP_tree& tree, const vector<const EM_matrix_base*>& submat, const Column_matrix& cm);
  void clear();  // zeroes the F matrix
  Prob prune();  // does pruning, returns final likelihood
};


#endif /* FAST_PRUNE_INCLUDED */
