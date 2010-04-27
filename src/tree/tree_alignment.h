#ifndef TREE_ALIGNMENT_INCLUDED
#define TREE_ALIGNMENT_INCLUDED

#include "tree/phylogeny.h"
#include "seq/alignment.h"
#include "seq/stockholm.h"
#include "seq/distmat.h"

struct Stockholm_tree : PHYLIP_tree
{
  // data
  vector<int> node2row;  // node-to-row mapping

  // constructor
  Stockholm_tree (const Stockholm& stock, bool die_if_tree_missing = true);
};

struct Tree_alignment
{
  PHYLIP_tree tree;
  Alignment align;

  // row<-->node index maps
  Phylogeny::Node_vector row2node;
  vector<int> node2row;

  // node_profile is a list of node-associated Score_profile's owned exclusively by this object
  // DO NOT set this directly; use set_node_profile() method instead
  typedef map<Phylogeny::Node,Score_profile*> Node_profile_map;
  Node_profile_map node_profile;   // DO NOT EDIT DIRECTLY; use set_node_profile() method
  void set_node_profile (Phylogeny::Node node, Score_profile* profile);

  // reference alignment for benchmarking & related data
  Alignment ref_align;
  bool has_ref_align;
  Alignment::Residue_pair_set ref_pair_set;

  // constructors, destructor
  Tree_alignment();
  Tree_alignment (const Tree_alignment& tree_align);
  Tree_alignment (const Stockholm& stock, bool die_if_trees_missing = true);
  virtual ~Tree_alignment();

  // observer methods
  virtual void tree_changed();             // empty observer method; mechanism for notifying derived classes that tree has changed
  virtual void align_changed();            // empty observer method; mechanism for notifying derived classes that alignment has changed

  // copy method
  virtual Tree_alignment& operator= (const Tree_alignment& t);         // virtual copy method

  // main methods
  void reset_maps();
  void build_maps_from_names();    // build row2node & node2row maps by matching up alignment row names with tree node names

  void clear_tree();
  void clear_align();
  void clear_node_profiles();
  void clear();

  void read_PHYLIP (istream& tree_stream);
  void read_MUL (istream& align_stream, Sequence_database& db);
  
  void write_Stockholm (ostream& out, const Alphabet& alphabet) const;

  void set_tree (const PHYLIP_tree& new_tree);
  void set_alignment (const Alignment& new_align);

  void check_map_sizes() const;        // throw an exception if map sizes don't match
  void update_tree_node_names_from_maps();      // copy alignment row names to tree nodes
  void update_alignment_row_names_from_maps();  // copy tree node names to alignment rows

  vector<int> unattached_rows() const;  // rows that are not in the tree

  void make_empty_alignment();        // make an empty alignment row for each tree node
  void add_empty_alignment_row_for_node (Phylogeny::Node node);  // add an empty alignment row for tree node #node
  void attach_sequences (const Sequence_database& db);  // for every sequence in db whose name is a valid tree node specifier, attach it to corresponding alignment row

  void assert_nodes_equal_rows() const;   // throw an exception if imperfect 1-to-1 map for all nodes
  void assert_leaves_equal_rows() const;  // throw an exception if imperfect 1-to-1 map for leaf nodes
  bool nodes_equal_rows() const;  // true if perfect 1-to-1 map for all nodes
  bool leaves_equal_rows() const;  // true if perfect 1-to-1 map for leaf nodes
  bool has_unattached_rows() const;  // true if any of the alignment rows are missing from the tree

  vector<Alignment_path::Row_pair>  row_pairs() const;       // list of pairs of alignment rows, each pair corresponding to a branch of the tree

  // window() returns columns (start) to (start+len-1) inclusive
  Tree_alignment window (int start, int len) const;

  // subpath() method returns pairwise subpath consistent with alignment tree
  // (i.e. no match columns unless intervening nodes have matches too)
  // group_inserts flag prevents messy gaps, e.g. _--_--__--_ becomes ------_____
  Pairwise_path subpath (Phylogeny::Node_pair pair, bool group_inserts = 0) const;

  // methods to change the alignment by realigning a branch or node - used by tkfalign
  Pairwise_path realign_pair (const Phylogeny::Node_pair& pair, const Pairwise_path& ppath);               // realign a branch
  Pairwise_path realign_node (const Phylogeny::Node node, const Pairwise_path& align_to_new_row_map);      // realign a node

  // debugging output method
  void show_decomposition (ostream& o, const Alignment_path::Decomposition& decomp) const;

  // method to estimate tree by neighbor-joining
  // NB: may break if alignment contains ancestral sequences with DART node-munged names
  //  (e.g. "A::B" for most recent common ancestor of A and B)
  // ...this shouldn't happen if this method is only called when tree is uninitialised (fingers crossed)
  void estimate_tree_by_nj (Dist_func_factory& dist_func_factory);

  // copy a tree into a Stockholm alignment (messy)
  void copy_tree_to_Stockholm (Stockholm& stock) const;

  // benchmarking methods (obsolete?)
  void read_benchmark_alignment (const char* benchmark_filename, Sequence_database& db);  // reads reference alignment in MUL format
  int  benchmark_residue_pair_overlap() const;
  void log_benchmark_results (const char* alignment_name = "Current alignment") const;  // print some stats to logfile
};

struct Tree_alignment_database
{
  // data
  Sequence_database& seq_db;
  Stockholm_database* stock_db;   // null unless initialise_from_Stockholm_database() was called
  list<Tree_alignment> tree_align_list;
  vector<Tree_alignment*> tree_align;
  vector<sstring> name;
  // size accessor
  int size() const;
  // clear, load methods
  void clear();
  void load (const char* index_filename);
  void initialise_from_Stockholm_database (const Stockholm_database& stock, bool die_if_trees_missing = TRUE);
  // helper method, ensures database has trees
  void estimate_missing_trees_by_nj (Dist_func_factory& dist_func_factory);
  // helper method, forces entire database to be binary
  void force_binary();
  // constructors
  Tree_alignment_database (Sequence_database& seq_db);
  Tree_alignment_database (Sequence_database& seq_db, const char* index_filename);  // calls load()
};

struct Insertion_database : Stockholm_database
{
  // constructor
  Insertion_database (const Stockholm& orig);
};

#endif
