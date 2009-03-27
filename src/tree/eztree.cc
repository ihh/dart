#include "tree/phylogeny.h"
#include "seq/stockholm.h"

// main
int main (int argc, char** argv)
{
  try
  {
    // create opts
    INIT_OPTS_LIST (opts, argc, argv, 2, "[options] <Newick format tree file> <Stockholm format alignment>",
		    "convert a Newick tree & Stockholm alignment into a more easily-parsed format");

    // parse opts
    opts.parse_or_die();
    const char* tree_filename = opts.args[0].c_str();
    const char* align_filename = opts.args[1].c_str();

    // open tree file
    ifstream tree_file (tree_filename);
    if (!tree_file)
      THROWEXPR ("Can't open tree file '" << tree_filename << "'");

    // read tree
    PHYLIP_tree tree;
    tree.read (tree_file);

    // open alignment file
    ifstream align_file (align_filename);
    if (!align_file)
      THROWEXPR ("Can't open alignment file '" << align_filename << "'");

    // read tree
    Stockholm stock;
    Sequence_database seq_db;
    stock.read_Stockholm (align_file, seq_db);

    // print it out
    cout << "NODES  " << tree.nodes() << "\n";
    cout << "ROOT   " << tree.root << "\n";
    for_rooted_branches_post (tree, b)
      cout << "BRANCH FROM " << (*b).first << " TO " << (*b).second << " LENGTH " << (*b).length << "\n";
    for_iterator (Phylogeny::Node_const_iter, l, tree.leaves_begin(), tree.leaves_end())
      cout << "LEAF   " << *l << " NAME " << tree.node_specifier(*l) << "\n";
    cout << "COLS   " << stock.columns() << "\n";
    for_iterator (Phylogeny::Node_const_iter, l, tree.leaves_begin(), tree.leaves_end())
      {
	if (stock.row_index.find (tree.node_specifier(*l)) == stock.row_index.end())
	  THROWEXPR ("Can't find sequence for tree node " << *l);
	cout << "SEQ    " << *l << " ";
	const int row = stock.row_index[tree.node_specifier(*l)];
	const Biosequence& seq = stock.np[row]->seq;
	int pos = 0;
	for (int col = 0; col < stock.columns(); ++col)
	  cout << (char) (stock.path(row,col) ? seq[pos++] : '-');
	cout << "\n";
      }
  }
  catch (const Dart_exception& e)
  {
	CLOGERR << e.what();
	return 1;
  }
  return 0;
}
