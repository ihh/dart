#include "treeshuffler.h"

void Tree_shuffler::reset_branch_pool()
{
  branch_pool.clear();
  copy (tree->branches_begin(0), tree->branches_end(), back_inserter(branch_pool));
  for (int n = 0; n < (int) branch_pool.size(); ++n)
    swap (branch_pool[n + Rnd::rnd_int(branch_pool.size() - n)], branch_pool[n]);
}

void Tree_shuffler::reset_node_pool()
{
  node_pool.clear();
  copy (tree->internals_begin(), tree->internals_end(), back_inserter(node_pool));
  for (int n = 0; n < (int) node_pool.size(); ++n)
    swap (node_pool[n + Rnd::rnd_int(node_pool.size() - n)], node_pool[n]);
}


Phylogeny::Branch Tree_shuffler::next_branch()
{
  if (branch_pool.size() == 0)
    reset_branch_pool();
  Phylogeny::Branch b = branch_pool.back();
  branch_pool.pop_back();
  return b;
}


Phylogeny::Node Tree_shuffler::next_node()
{
  if (node_pool.size() == 0)
    reset_node_pool();
  Phylogeny::Node n = node_pool.back();
  node_pool.pop_back();
  return n;
}

Phylogeny::Branch Tree_shuffler::next_internal_branch()
{
  bool found_internal = 0;
  for_branches_pre (*tree, 0, -1, b)
    if (tree->is_internal (*b))
      {
	found_internal = 1;
	break;
      }
  if (!found_internal)
    THROWEXPR ("Tree has no internal branches")

  Phylogeny::Branch b;
  do
    b = next_branch();
  while (!tree->is_internal (b.first) || !tree->is_internal (b.second));
  return b;
}

Phylogeny::Node Tree_shuffler::next_internal_node()
{
  bool found_internal = 0;
  for_nodes_pre (*tree, 0, -1, b)
    if (tree->is_internal ((*b).second))
      {
	found_internal = 1;
	break;
      }
  if (!found_internal)
    THROWEXPR ("Tree has no internal nodes")

  Phylogeny::Node n;
  do
    n = next_node();
  while (!tree->is_internal (n));
  return n;
}



Tree_shuffler::Tree_shuffler() : tree (0)
{ }

void Tree_shuffler::set_tree (const PHYLIP_tree& tree_ref)
{
  tree = &tree_ref;
  reset();
}

void Tree_shuffler::get_next_slide (Node& grumpa, Node& dad, Node& son)
{
  // check that we can find a node to slide
  tree->assert_tree_is_binary();
  if (tree->internals() < 2)
    THROWEXPR ("Tree has too few internal nodes to do a node-slide move");
  // pick a node
  do {
    dad = next_internal_node();
  } while (tree->neighbours (dad) != 3);
  // get neighbors
  Node_vector relatives = tree->relatives (dad);
  if (relatives.size() != 3)
    THROWEXPR ("Number of relatives != 3");  // unreachable
  // pick two neighbors at random
  swap (relatives[0], relatives [Rnd::rnd_int (3)]);
  grumpa = relatives[0];
  son = relatives[1 + Rnd::rnd_int (2)];
}

void Tree_shuffler::get_next_swap (Node& aunt, Node& nephew, Node& grumpa, Node& dad)
{
  // check that we can find a node to swap
  tree->assert_tree_is_binary();
  if (tree->internals() < 2)
    THROWEXPR ("Tree has too few internal nodes to do a branch-swap move");
  // pick a branch
  const Branch b = next_internal_branch();
  grumpa = b.first;
  dad = b.second;
  // get neighbors
  const Node_vector aunts = tree->children (grumpa, tree->parent[grumpa]);
  const Node_vector nephews = tree->children (dad, grumpa);
  if (aunts.size() < 2)
    THROWEXPR ("Something stinks in this tree");
  // pick aunt & nephew
  nephew = nephews [Rnd::rnd_int (nephews.size())];
  do
    aunt = aunts [Rnd::rnd_int (aunts.size())];
  while (aunt == dad);
}

Tree_shuffler::Action Tree_shuffler::next_action()
{
  vector<double> weight (Total_moves, 0.);
  weight[0] = node_realign_rate;
  weight[1] = branch_realign_rate;
  weight[2] = node_slide_rate;
  weight[3] = branch_scale_rate;
  weight[4] = branch_swap_rate;
  weight[5] = indel_param_sampling_rate;
  weight[6] = subst_param_sampling_rate;

  switch (Rnd::choose (weight))
    {
    case 0: return Sample_node;
    case 1: return Sample_branch;
    case 2: return Slide_node;
    case 3: return Scale_branch;
    case 4: return Swap_branches;
    case 5: return Sample_indel_params;
    case 6: return Sample_subst_params;
    default: break;
    }
  return Sample_branch;  // default, should be unreachable
}

void Tree_shuffler::reset()
{
  reset_branch_pool();
  reset_node_pool();
}

void Tree_shuffler::cue (const sstring& node_name)
{
  const Node node = tree->find_node_or_die (node_name);
  node_pool.push_back (node);
  if (node != tree->root)
    branch_pool.push_back (Branch (tree->parent[node], node, *tree));
}
