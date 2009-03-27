#include <fstream>

#include "tree/phylogeny.h"
#include "tree/ps_tree.h"
#include "util/logfile.h"

#define CLOGTOP CLOG(8)

int main(int argc, char* argv[])
{
  if (argc != 4)
    {
      CLOGERR << argv[0] << " - calculate consensus branch support for a tree given a set of comparison trees\n";
      CLOGERR << "usage: " << argv[0] << " branch_length reference_tree_file comparison_trees_file\n";
      exit(1);
    }

  PHYLIP_tree ref_tree;
  Phylogeny::Branch_support support;

  double branch_length = atof (argv[1]);
  ifstream ref_file (argv[2]);
  ifstream comparison_file (argv[3]);

  try
    {
      ref_tree.read (ref_file);
      CLOGTOP << "Read reference tree\n";

      support = ref_tree.calculate_branch_support (comparison_file);
      CLOGTOP << "Read comparison trees\n";
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  PS_tree picture (ref_tree, ref_tree.root);
  picture.branch_length_multiplier = branch_length;
  picture.calculate_heights();

  picture.draw_branches (support);
  picture.label_branches (support);
  picture.label_leaf_nodes (ref_tree.node_name);
  picture.number_nodes ();
  picture.showpage ();

  cout << picture;

  return 1;
}
