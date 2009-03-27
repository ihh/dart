#include <math.h>
#include "tree/ps_tree.h"

PS_tree::PS_tree (const Phylogeny& phylogeny, Phylogeny::Node top, Phylogeny::Node parent_of_top, double hscale) :
  phylogeny (phylogeny),

  top (top),
  parent_of_top (parent_of_top),
  parent_of_node (phylogeny.find_parents (top, parent_of_top)),

  node_height (phylogeny.nodes(), (double) 0),
  node_branch_length (phylogeny.nodes(), (double) 0),
  node_x (phylogeny.nodes(), (double) 0),
  node_y (phylogeny.nodes(), (double) 0),

  font_scale (10),
  char_height (10),
  char_width (8),
  vertical_leaf_separation (5),
  branch_length_multiplier (100*hscale),
  root_stub_length (1),
  line_width (6)
{
  calculate_heights ();
  setrgbcolor (0, 0, 0);
  translate (30, 30);
}

void PS_tree::calculate_heights()
{
  node_height = vector<double> (phylogeny.nodes(), (double) 0);
  node_branch_length = vector<double> (phylogeny.nodes(), (double) 0);
  node_x = vector<double> (phylogeny.nodes(), (double) 0);
  node_y = vector<double> (phylogeny.nodes(), (double) 0);

  for_iterator (Phylogeny::Branch_iter, branch,
		phylogeny.branches_begin (top, parent_of_top, 0, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*branch).second, parent = (*branch).first;
      if (phylogeny.is_leaf (node)) node_height [node] = char_height + vertical_leaf_separation;
      if (parent != -1) node_height [parent] += node_height [node];
    }

  vector<double> next_y (phylogeny.nodes(), (double) 0);

  for_iterator (Phylogeny::Branch_iter, branch,
		phylogeny.branches_begin (top, parent_of_top, 1, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*branch).second, parent = (*branch).first;
      node_branch_length [node] = parent==-1 ? root_stub_length : max (root_stub_length, max (0., ((Phylogeny&)phylogeny).branch_length (*branch)) * branch_length_multiplier);
      node_x [node] = (parent==-1 ? 0 : node_x [parent]) + node_branch_length [node];
      next_y[node] = node_y [node] = parent==-1 ? 0 : next_y[parent];
      if (parent != -1) next_y[parent] += node_height[node];
    }
}

void PS_tree::draw_branches()
{
  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top, 1, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node parent = (*b).first, node = (*b).second;
      Phylogeny::Node first_kid = -1, last_kid = -1;
      Phylogeny::Child_iter n = phylogeny.children_begin (node, parent);
      while (n != phylogeny.children_end (node, parent)) { if (first_kid == -1) first_kid = *n; last_kid = *n++; }
      if (first_kid != -1)
	  line (node_x[node], y_midpoint(first_kid), node_x[node], y_midpoint(last_kid));
      line (node_x[node], y_midpoint(node), node_x[node] - node_branch_length[node], y_midpoint(node));
    }
}

void PS_tree::draw_branches (const Phylogeny::Branch_support& support, double grayscale_mul)
{
  vector<double> branch_width (phylogeny.nodes(), (double) 0);
  vector<double> node_gray (phylogeny.nodes(), (double) 0);

  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top, 1, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*b).second;
      double bsupp = ((Phylogeny::Branch_support&)support).hits[*b] / support.total;

      node_gray[node] = grayscale_mul * (1 - bsupp);
      branch_width[node] = (line_width-1) * bsupp + 1;
    }

  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top, 1, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node parent = (*b).first, node = (*b).second;
      gsave ();
      setgray (node_gray[node]);

      Phylogeny::Node first_kid = -1, last_kid = -1;
      if (node != -1)
	{
	  Phylogeny::Child_iter n = phylogeny.children_begin (node, parent);
	  while (n != phylogeny.children_end (node, parent)) { if (first_kid == -1) first_kid = *n; last_kid = *n++; }
	}

      if (first_kid != -1)
	rectangle (node_x[node] - branch_width[node] / 2,
		   y_midpoint(first_kid) - branch_width[first_kid] / 2,
		   node_x[node] + branch_width[node] / 2,
		   y_midpoint(last_kid) + branch_width[last_kid] / 2) . fill();

      rectangle (node_x[node] - node_branch_length[node] + (parent==-1 ? 0 : branch_width[parent] / 2),
		 y_midpoint(node) - branch_width[node] / 2,
		 node_x[node],
		 y_midpoint(node) + branch_width[node] / 2) . fill();

      grestore();
    }
}

void PS_tree::label_branches (const Phylogeny::Branch_support& support)
{
  vector<double> branch_width (phylogeny.nodes(), (double) 0);

  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top, 1, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*b).second;
      double bsupp = ((Phylogeny::Branch_support&)support).hits[*b] / support.total;
      branch_width[node] = (line_width-1) * bsupp + 1;
    }

  char _buf[100];
  branch_font().scalefont(font_scale / 3).setfont();
  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top, 1, 1),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*b).second;
      double bsupp = ((Phylogeny::Branch_support&)support).hits[*b] / support.total;
      moveto (node_x[node] - node_branch_length[node] + char_width / 3, y_midpoint(node) - char_height / 3 - branch_width[node]/2);
      sprintf (_buf, "%.2f", 100*bsupp);
      show (_buf);
    }
}

void PS_tree::label_leaf_nodes (const Phylogeny::Node_name_vec& leaf_label)
{
  node_font().scalefont(font_scale).setfont();
  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*b).second;
      if (phylogeny.is_leaf (node))
	{
	  moveto (node_x[node] + char_width / 2, y_midpoint(node) - char_height / 2);
	  show (leaf_label[node].c_str());
	  stroke();
	}
    }
}


void PS_tree::number_nodes ()
{
  char _buf[100];
  node_font().scalefont(font_scale / 4).setfont();
  for_iterator (Phylogeny::Branch_iter, b,
		phylogeny.branches_begin (top, parent_of_top),
		phylogeny.branches_end ())
    {
      Phylogeny::Node node = (*b).second;
      double chars = 1 + (int) (log((double)node) / log(10.));
      moveto (node_x[node] - chars * char_width / 4 - line_width/2, y_midpoint(node) + line_width * 2);
      sprintf (_buf, "%d", node);
      show (_buf);
      stroke();
    }
}

PS_tree& PS_tree::branch_font() { findfont("Helvetica"); return *this; }
PS_tree& PS_tree::node_font() { findfont("Helvetica"); return *this; }

void PS_tree::draw_scale()
{
  line(0,-15,branch_length_multiplier,-15);
  branch_font().scalefont(font_scale).setfont();
  moveto (branch_length_multiplier, -15 - char_height/2);
  show("1 expected substitution per site");
}
