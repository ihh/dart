#include <sstream>
#include "tree/phylogeny.h"
#include "util/piper.h"

int main (int argc, char** argv)
{
  try
    {
      //      if (argc != 2) THROWEXPR ("Usage: " << argv[0] << " <tree file>");

      const sstring tree_string = "((NP_828851.1/75-609:100.00,(Q66951/272-808:100.00,((Q84712/231-734:100.00,VGL2_CVH22/32-536:100.00):39.27,VGL2_CVCAI/259-790:100.00):37.50):34.31):18.31,(((Q64930/20-543:100.00,VGL2_IBVD2/20-538:100.00):80.50,Q82666/20-541:100.00):79.60,VGL2_IBVM/20-537:100.00):76.41);";
      istringstream tree_input (tree_string);

      PHYLIP_tree tree;
      tree.read (tree_input);

      if (tree.force_binary()) cout << "Changed to binary\n";

      vector<int> seen (tree.nodes(), 0);
      int n_seen = 0;
      for_iterator (Phylogeny::Branch_iter, bi,
		    tree.branches_begin (tree.root, -1, 1, 1),
		    tree.branches_end())
	{
	  const Phylogeny::Branch& b = *bi;
	  if (tree.parent[b.second] != b.first) THROWEXPR ("Parent/child askew");
	  if (seen[b.second]++) THROWEXPR ("Duplicate child");
	  ++n_seen;

	  vector<int> seen2 (tree.nodes(), 0);

	  if (b.first != -1)
	    {
	      ++seen2[b.first];
	      for_iterator (Phylogeny::Branch_iter, bj,
			    tree.branches_begin (b.first, b.second, 0, 0),
			    tree.branches_end())
		if (seen2[(*bj).second]++) THROWEXPR ("Duplicate child");
	    }

	  if (b.second != -1)
	    {
	      ++seen2[b.second];
	      for_iterator (Phylogeny::Branch_iter, bj,
			    tree.branches_begin (b.second, b.first, 0, 0),
			    tree.branches_end())
		if (seen2[(*bj).second]++) THROWEXPR ("Duplicate child");
	    }

	  for_const_contents (vector<int>, seen2, s)
	    if (*s == 0) THROWEXPR ("Missing child");
	}
      if (n_seen != tree.nodes()) THROWEXPR ("Missing child");

      cout << "ok\n";
    }
  catch (const Dart_exception& e)
    {
      cout << "not ok\n";
      CLOGERR << e.what();
      exit(1);
    }
  return 0;
}

