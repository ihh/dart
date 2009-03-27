#ifndef PSTREEALIGN_INCLUDED
#define PSTREEALIGn_INCLUDED

#include "tree/tree_alignment.h"
#include "util/ps_doc.h"

struct PS_tree_alignment : PS_doc
{
  const Tree_alignment&  tree_align;

  vector <double>        node_height;
  vector <double>        node_branch_length;
  vector <double>        node_x;
  vector <double>        node_y;
  
  double y_midpoint (Phylogeny::Node node) { return node_y[node] + node_height[node] / 2; }

  double font_scale;
  double char_height;
  double char_width;
  double vertical_leaf_separation;
  double branch_length_multiplier;
  double root_stub_length;
  double line_width;

  PS_tree_alignment (const Tree_alignment& tree_align);

  virtual PS_tree& branch_font();
  virtual PS_tree& node_font();
};

#endif
