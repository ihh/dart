// subclass of Phylogeny that builds itself using the neighbor-joining algorithm
// should get round to doing weighted neighbor-joining (Bruno) as well, one of these days...

#ifndef NJ_INCLUDED
#define NJ_INCLUDED

#include "tree/phylogeny.h"

class NJ_tree : public PHYLIP_tree {
public:
  NJ_tree();
  void build (const Node_name_vec& node_name, const array2d<double>& distance_matrix);
};

#endif
