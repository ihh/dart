#ifndef PS_TREE_INCLUDED
#define PS_TREE_INCLUDED

#include "tree/phylogeny.h"
#include "util/ps_doc.h"

struct PS_tree : PS_doc
{
  const Phylogeny&       phylogeny;

  Phylogeny::Node        top;
  Phylogeny::Node        parent_of_top;

  Phylogeny::Node_vector parent_of_node;

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

  PS_tree (const Phylogeny& phylogeny, Phylogeny::Node top, Phylogeny::Node parent_of_top = -1, double hscale = 1.);
  virtual ~PS_tree() { }

  void calculate_heights ();
  void draw_branches ();
  void draw_branches (const Phylogeny::Branch_support& support, double grayscale_mul = .5);
  void label_branches (const Phylogeny::Branch_support& support);

  void label_leaf_nodes (const Phylogeny::Node_name_vec& leaf_label);
  void number_nodes ();

  void draw_scale();
  
  virtual PS_tree& branch_font();
  virtual PS_tree& node_font();
};

#endif
