#include "tree/phylogeny.h"

int main (int argc, char** argv)
{
  // create Opts_list object using a macro
  INIT_OPTS_LIST (opts, argc, argv, 2, "[options] <tree file> <name of new root node>",
		  "re-root a Newick-format tree");

  opts.print_title ("Output options");

  int maxcols;
  opts.add ("mc -maxcols", maxcols = -1, "maximum columns for output", false);

  opts.parse_or_die();  // parse the command-line options
  try
    {
      const sstring tree_filename = opts.args[0];
      const sstring new_root_name = opts.args[1];
      
      ifstream tree_file (tree_filename.c_str());
      if (!tree_file) THROWEXPR ("Couldn't open tree file " << tree_filename);
      PHYLIP_tree tree;
      tree.read (tree_file);

      const Phylogeny::Node new_root = tree.find_node (new_root_name.c_str());
      if (new_root < 0)
	THROWEXPR ("Node '" << new_root_name << "' not found in tree");

      tree.write (cout, maxcols, new_root);
    }
  catch (const Dart_exception& e)  // exception; bail out gracefully
    {
      CLOGERR << e.what();
      exit(1);
    }
  return 0;
}
