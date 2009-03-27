#ifndef DRAW_PS_TREE_INCLUDED
#define DRAW_PS_TREE_INCLUDED

#include "tree/ps_tree.h"

int main (int argc, char** argv)
{
  // create Opts_list object using a macro
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <tree file>",
		  "plot a tree in Postscript");

  double hscale;
  bool drawscale;
  opts.add ("hs -hscale", hscale = 1., "horizontal scale factor");
  opts.add ("ds -draw-scale", drawscale = false, "draw scale");
  
  opts.parse_or_die();  // parse the command-line options
  try
    {
      const sstring tree_filename = opts.args[0];
      
      ifstream tree_file (tree_filename.c_str());
      if (!tree_file) THROWEXPR ("Couldn't open tree file " << tree_filename);
      PHYLIP_tree tree;
      tree.read (tree_file);

      PS_tree ps_tree (tree, tree.root, -1, hscale);
      if (drawscale)
	ps_tree.draw_scale();
      ps_tree.draw_branches();
      ps_tree.label_leaf_nodes (tree.node_name);

      cout << ps_tree;
    }
  catch (const Dart_exception& e)  // exception; bail out gracefully
    {
      CLOGERR << e.what();
      exit(1);
    }
  return 0;
}

#endif  /* DRAW_PS_TREE_INCLUDED */
