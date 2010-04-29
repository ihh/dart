#ifndef PHYLOGENY_INCLUDED
#define PHYLOGENY_INCLUDED

#include <vector>
#include <set>
#include <map>
#include "util/dexception.h"
#include "util/Regexp.h"
#include "seq/biosequence.h"
#include "util/logfile.h"
#include "seq/stockholm.h"

// name of root node (auto-assigned by some methods if no other name defined)
#define DART_ROOT_IDENTIFIER "root"

// maximum branch length (used by other classes, but defined here)
#define DART_MAX_BRANCH_LENGTH 10.


// Phylogeny: an unrooted tree with a length associated with each branch
class Phylogeny
{
 public:

  // typedefs
  //
  // first, Nodes & associated variables
  //
  typedef int Node;

  typedef vector <Node>                    Node_vector;
  typedef Node_vector::iterator            Node_iter;
  typedef Node_vector::const_iterator      Node_const_iter;

  typedef set < Node, less<Node> >         Node_set;

  typedef vector <sstring>                 Node_name_vec;
  typedef map <Node, double>               Branch_length_map;

  typedef map <sstring, Node>              Node_directory;
  typedef map <Node, Node>                 Node_map;

  // Node_pair is a directed edge, implemented as a pair of Node's.
  struct Node_pair : pair <Node, Node>
  {
    Node_pair () : pair <Node, Node> (-1, -1) {}
    Node_pair (Node n1, Node n2) : pair <Node, Node> (n1, n2) {}
    Node_pair& swap() { int tmp = second; second = first; first = tmp; return *this; }
    Node_pair& undirect() { if (first > second) swap(); return *this; }
    Node_pair  transpose() const { Node_pair tmp = *this; return tmp.swap(); }
    friend ostream& operator<<(ostream&o, const Node_pair& np)
      { o << "(" << np.first << "-" << np.second << ")"; return o; }
  };

  // Undirected_pair is an undirected edge, implemented as a sorted pair of Node's
  // (inheriting from Node_pair), with first < second.
  class Undirected_pair : public Node_pair
  {
  private:
    Node_pair& swap();         // make this private to prevent public access
    Node_pair  transpose();    // make this private to prevent public access
  public:
    Undirected_pair () : Node_pair (-1, -1) {}
    Undirected_pair (Node n1, Node n2) : Node_pair (min(n1,n2), max(n1,n2)) {}
    Undirected_pair (const Node_pair& n) : Node_pair (min(n.first,n.second), max(n.first,n.second)) {}

    Undirected_pair& operator= (const Node_pair& n)
    {
      first = min (n.first, n.second);
      second = max (n.first, n.second);
      return *this;
    }
  };

  // A Branch is a Node_pair (directed edge) with a length parameter.
  struct Branch : Node_pair
  {
    double length;
    Branch () : Node_pair (), length (0) {}
    Branch (Node n1, Node n2, double l) : Node_pair (n1, n2), length (l) {}
    Branch (Node n1, Node n2, const Phylogeny& p) : Node_pair (n1, n2), length (((Phylogeny&)p).branch_length (n1, n2)) {}
    Node other(Node node) const { return node == first ? second : (node == second ? first : -1); }
    struct Shorter_or_indexed_before
    { bool operator() (const Branch& b1, const Branch& b2) const { return b1.length < b2.length || (b1.length == b2.length && (b1.first < b2.first || (b1.first == b2.first && b1.second < b2.second))); } };
  };

  // Surfing the DeathMarch
  // higher level branch-related structures
  typedef map <Node_pair, Node_set>                        Node_pair2set_map;
  typedef map <Node_set, Node_pair>                        Node_set2pair_map;
  typedef set <Branch, Branch::Shorter_or_indexed_before>  Branch_set;
  typedef vector <Branch>                                  Branch_vector;

  // iterator for nodes
  //
  struct Relative_iter : public iterator <bidirectional_iterator_tag, Node>
  {
    Branch_length_map::const_iterator  i;
    Branch_length_map::const_iterator  begin;
    Branch_length_map::const_iterator  end;
    Node                               avoid_node;
    Relative_iter (Branch_length_map::const_iterator begin, Branch_length_map::const_iterator end, Node a = -1) : i (begin), begin (begin), end (end), avoid_node (a)
      { while (i != end) if ((*i).first != avoid_node) break; else ++i; }
    const Node& operator*() const { return (*i).first; }
    bool operator== (const Relative_iter& r) const { return i == r.i && avoid_node == r.avoid_node; }
    bool operator!= (const Relative_iter& r) const { return !(*this == r); }
    Relative_iter& operator++() { while (i != end) if ((*(++i)).first != avoid_node) break; return *this; }
    Relative_iter operator++ (int) { Relative_iter r = *this; ++*this; return r; }
    Relative_iter& operator--() { while (i != begin) if ((*(--i)).first != avoid_node) break; return *this; }
    Relative_iter operator-- (int) { Relative_iter r = *this; --*this; return r; }
  };
  typedef Relative_iter Child_iter;     // very lazy hack because I can't be bothered to strip out this dead class Child_iter

  // iterator to traverse all the Branches in a tree, depth-first
  //
  // currently uses type switches to be pseudo-polymorphic... ugh
  //
  // flags are: top_down          (always visits parents before children if TRUE; vice versa if FALSE; default is FALSE)
  //            visit_grumpa      (visits a (pseudo-)branch from node #grumpa to node #dad if TRUE; default is FALSE)
  //
  // Whether top-down or bottom-up, we're always depth-first in this neighbourhood.
  //
  // Branch::first is the parent, Branch::second is the child
  //
  class Branch_iter : public iterator <forward_iterator_tag, Branch>
    {
    private:
      const Phylogeny*   phylogeny;
      Node               dad;
      Node               grumpa;
      vector<Child_iter> path;
      vector<Child_iter> end;
      Branch             branch;

      bool               top_down;         // type switch - very un-OOP, so sue me.
      bool               visit_grumpa;     // another type switch, i'm so bad
      bool               ignore_kids;      // TRUE if children of current node should be seen but not heard
      bool               finished;         // TRUE when we've finished

      Node most_recent_node() const { int s = path.size(); return s == 0 ? dad : *path[s-1]; }
      Node most_recent_parent() const { int s = path.size(); return s == 0 ? grumpa : (s == 1 ? dad : *path[s-2]); }

      bool find_children()   // returns FALSE unless this path should be extended even further (NB always FALSE if doing a bottom-up strategy)
	{
	  if (ignore_kids) { ignore_kids = 0; return 0; }
	  while (1)
	    {
	      Node n = most_recent_node(), p = most_recent_parent();
	      Child_iter c = phylogeny->children_begin (n, p);
	      Child_iter e = phylogeny->children_end (n, p);
	      if (c == e) return 0;
	      path.push_back(c);
	      end.push_back(e);
	      if (top_down) return 1;
	    }
	  return 0;
	}

      void setup_branch()      // sets up the current Branch, based on the path
	{
	  if (most_recent_parent() != -1)
	    branch = Branch (most_recent_parent(), most_recent_node(), *phylogeny);
	  else
	    branch = Branch (most_recent_parent(), most_recent_node(), 0);
	}

    public:
      const Branch& operator*() const { return branch; }

      bool operator==(const Branch_iter& i2) const
	{ return path == i2.path && finished == i2.finished; }

      bool operator!=(const Branch_iter& i2) const { return !(*this == i2); }

      Branch_iter& operator++()        // preincrement: move to the next Branch
	{
	  if (top_down)          // top-down
	    {
	      if (!find_children())
		while (path.size())
		  {
		    if (++path.back() == end.back()) { path.pop_back(); end.pop_back(); }
		    else break;
		  }
	      if (path.size() == 0) finished = 1;
	    }
	  else                   // bottom-up
	    {
	      if (path.size() == 0) finished = 1;
	      else
		{
		  if (++path.back() == end.back()) { path.pop_back(); end.pop_back(); }
		  else find_children();
		  if (path.size() == 0) finished = !visit_grumpa;
		}
	    }
	  setup_branch();
	  return *this;
	}
      
      // postincrement operator is commented out 'cos linker can't find vector<Child_iter> copy constructor
      // ... really should fix this ... best way would be to rewrite Child_iter (and Relative_iter) without using generic programming ...
      //
      //      Branch_iter operator++(int) { Branch_iter tmp = *this; (*this)++; return tmp; }

      void skip_children()              // don't explore any further down current branch. only works for top-down traversal
      {
	if (!top_down) THROW Standard_exception ("Attempt to skip children in a bottom-up tree traversal");
	ignore_kids = 1;
      }

      // default constructor returns end() iterator
      //
      Branch_iter() : phylogeny(0), finished(1) {}

      // main constructor
      //
      Branch_iter (const Phylogeny& p, Node dad, Node grumpa = -1, bool top_down = 0, bool visit_grumpa = 0) :
	phylogeny(&p),
	dad(dad),
	grumpa(grumpa),
	top_down(top_down),
	visit_grumpa(visit_grumpa),
	ignore_kids (0),
	finished (0)
	{
	  if (!(top_down && visit_grumpa))
	    if (!find_children())
	      finished = !visit_grumpa;
	  setup_branch();
	}
    };
  
  friend class Branch_iter;    // is this necessary? oh well

 protected:
  vector<Branch_length_map> _branch_length;       // length of each branch
  int                       _branch_count;        // number of branches

  Node_vector               _sorted_nodes;        // sorted list of nodes: leaves first, then internals, sorted numerically
  Node_iter                 _first_internal_node; // index of first internal node in _sorted_nodes vector
  vector<Node_iter>         _node_sort_pos;       // mapping from each node to its position in _sorted_nodes list

  Node new_node();                     // NB leaves node lists unsorted ...
  void sort_nodes();                   // ... this sorts them
  
 public:

  // constructors, assignment operator
  //
  Phylogeny();
  Phylogeny (const Phylogeny& p);

  Phylogeny& operator= (const Phylogeny& p);

  // publically accessible tree editing methods
  //
  void clear();
  Node add_node (Node neighbour = -1, double branch_length = 0);
  Node_vector delete_node (Node n);  // NB may change the index of every node in the tree; returns old-to-new-node mapping
  void move_node (Node n, Node old_neighbour, Node new_neighbour, double new_branch_length);
  void move_node (Node n, Node old_neighbour, Node new_neighbour)
    { move_node (n, old_neighbour, new_neighbour, branch_length (n, old_neighbour)); }
  void swap_nodes (Node n1, Node n2);

  void add_branch (Node node1, Node node2, double length);
  void remove_branch (Node node1, Node node2);

  // methods to get length of a given branch (fast)
  //
  double& branch_length (Node node1, Node node2)
    { return (*(_branch_length[min(node1,node2)].find(max(node1,node2)))).second; }
  double& branch_length (const pair <Node, Node>& np)
    { return (*(_branch_length[min(np.first,np.second)].find(max(np.first,np.second)))).second; }

  const double& branch_length (Node node1, Node node2) const
    { return (*(_branch_length[min(node1,node2)].find(max(node1,node2)))).second; }
  const double& branch_length (const pair <Node, Node>& np) const
    { return (*(_branch_length[min(np.first,np.second)].find(max(np.first,np.second)))).second; }

  // methods to return miscellaneous data about tree (fast)
  //
  int  nodes() const { return _sorted_nodes.size(); }                                     // number of nodes
  int  branches() const { return _branch_count; }                                         // number of branches
  int  leaves() const { return (int) (((Node_const_iter) _first_internal_node) - _sorted_nodes.begin()); }    // number of leaf nodes
  int  internals() const { return (int) (_sorted_nodes.end() - ((Node_const_iter) _first_internal_node)); }   // number of internal nodes
  int  neighbours(Node node) const { return _branch_length[node].size(); }                // number of neighbours of a node
  bool is_leaf(Node node) const { return _branch_length[node].size() <= 1; }              // TRUE if node is a leaf node
  bool is_internal(Node node) const { return !is_leaf (node); }                           // TRUE if node is an internal node
  int  leaf_index(Node node) const                                                        // number of leaf nodes before this one
    { return is_leaf(node) ? (int) (((Node_const_iter) _node_sort_pos[node]) - _sorted_nodes.begin()) : -1; }
  int  internal_index(Node node) const                                                    // number of internal nodes before this one
    { return is_internal(node) ? (int) (((Node_const_iter) _node_sort_pos[node]) - _first_internal_node) : -1; }

  bool is_internal(const Node_pair& b) const { return is_internal(b.first) && is_internal(b.second); }  // TRUE if branch b is internal

  // methods returning iterators for moving through the tree
  //
  // leaf node iterators
  //
  Node_const_iter leaves_begin() const { return _sorted_nodes.begin(); }
  Node_const_iter leaves_end() const { return _first_internal_node; }

  // internal node iterators
  //
  Node_const_iter internals_begin() const { return _first_internal_node; }
  Node_const_iter internals_end() const { return _sorted_nodes.end(); }

  // relatives of a node
  //
  Relative_iter relatives_begin(Node node) const
    { return Relative_iter (_branch_length[node].begin(), _branch_length[node].end(), -1); }

  Relative_iter relatives_end(Node node) const
    { return Relative_iter (_branch_length[node].end(), _branch_length[node].end(), -1); }

  // children of a node (caller must specify which relative is the parent, since tree is unrooted)
  //
  Child_iter children_begin(Node node, Node parent) const
    { return Relative_iter (_branch_length[node].begin(), _branch_length[node].end(), parent); }

  Child_iter children_end(Node node, Node parent) const
    { return Relative_iter (_branch_length[node].end(), _branch_length[node].end(), parent); }

  // all branches in the tree - this is also currently the only iterator for traversing all nodes in bottom-up or dad-down order
  //
  Branch_iter branches_begin (Node node = 0, Node parent = -1, bool top_down = 0, bool visit_parent = 0) const
    { return Branch_iter (*this, node, parent, top_down, visit_parent); }

  Branch_iter branches_end () const { return Branch_iter (); }

  // method to return a vector of relatives of a node
  //
  Node_vector relatives(Node node) const
    { Node_vector v; copy (relatives_begin(node), relatives_end(node), back_inserter(v)); return v; }

  // method to return a vector of children of a node (caller must specify which relative is the parent, since tree is unrooted)
  //
  Node_vector children(Node node, Node parent) const
    { Node_vector v; copy (children_begin(node,parent), children_end(node,parent), back_inserter(v)); return v; }

  // method to return a vector mapping each node to its parent, given a dad branch (just uses Branch_iter)
  //
  Node_vector find_parents (Node dad, Node grumpa = -1) const
    { Node_vector v (nodes(), -1); for (Branch_iter b = branches_begin (dad, grumpa, 1, 1); b != branches_end(); ++b) v[(*b).second] = (*b).first; return v; }

  // static methods for mapping node indices from one tree to another given node names
  //
  static Node_directory node_directory (const Node_name_vec& node_name);   // permits reverse lookup by node name
  static Node_map       node_remapper (const Node_name_vec& node1_name, const Node_directory& node2_directory);
  static Node_map       node_remapper (const Node_name_vec& node1_name, const Node_name_vec& node2_name)
    { return node_remapper (node1_name, node_directory (node2_name)); }
  static Node_set       remap_node_set (const Node_set& node_set, const Node_map& node_remapper);

  // Branch_support structure for bootstrapping etc
  //
  struct Branch_support
  {
    map <Undirected_pair, double> hits;
    double total;
    Branch_support () : hits (), total (0) {}
  };

  // method to find the most recent common ancestor, given a vector mapping each node to its parent
  //
  Node youngest_common_ancestor (Node n1, Node n2, const Node_vector& parent) const;

  // method to find path connecting two nodes
  //
  Node_vector node_path (Node n1, Node n2, const Node_vector& parent) const;
};

// PHYLIP_tree: a rooted tree with named leaf nodes.
// Intended to be compatible with the PHYLIP format.
// Horrible hack: if a branch length is missing from the PHYLIP file, this is represented as a branch length of -1.
class PHYLIP_tree : public Phylogeny
{
  // private methods, don't look too closely ;)
 private:
  // internal parser methods
  // these are pretty ugly inside
  Node read_node (istream& in, Node parent_node);       // returns index of created node; creates parents vector on the fly
  void write_node(ostream& out,                         // NB doesn't need a valid parents vector
		  Node node,
		  Node from_node,
		  int max_columns,
		  int& columns) const;

  // private method to rearrange the node indices
  void remap_nodes (const Node_vector& old_to_new_map);


  // public methods
 public:
  Node_name_vec node_name;  // name of each node
  Node_vector   parent;     // parent of each node
  Node          root;       // root node

  // constructor
  PHYLIP_tree();

  // methods to verify or maintain integrity of tree
  void assert_tree_is_binary() const;                      // throws an exception if tree isn't binary
  bool force_binary();                                     // returns TRUE if tree was changed
  void rebuild_parents() { parent = find_parents(root); }  // sets up parents vector

  // named version of add_node
  Node add_named_node (const sstring& name, Node neighbour = -1, double branch_length = 0)
  {
    node_name.push_back(name);
    parent.push_back(neighbour);
    return add_node(neighbour,branch_length);
  }

  // given a Node_pair, flip it into parent-child order
  Node_pair parent_child_pair (const Node_pair& np) const
    {
      if (np.first == -1) return Node_pair (np.first, np.second);
      if (np.first == parent[np.second]) return Node_pair (np.first, np.second);
      if (np.second != parent[np.first]) THROW Standard_exception ("Bad node pair");
      return Node_pair (np.second, np.first);
    }

  // find depth of a node, i.e. distance to root
  int depth (Node n) const { int d = 0; while (n != root) { n = parent[n]; ++d; } return d; }

  // find sum of branch lengths from root to a node
  double branch_depth (Node n) const { double d = 0; while (n != root) { d += branch_length (n, parent[n]); n = parent[n]; } return d; }

  // sets of leaves and ancestors
  // NB this can return slightly different results from the leaves_*() and internals_*() iterators in the superclass,
  // because the root may not be an internal node (and therefore would be counted as a leaf by the superclass, but an ancestor by this class)
  vector<Node> leaf_vector();
  vector<Node> ancestor_vector();
  
  // I/O
  // read() method: reads a New Hampshire format tree.
  // Parser is a bit antiquated and clunky, may not comply strictly with the NH grammar.
  // Missing branch lengths are represented as -1 (I know).
  void read (istream& in) { clear(); node_name.clear(); parent.clear(); root = read_node (in, -1); }

  // write() method: writes the tree in New Hampshire format, followed by a newline.
  // Set max_columns=0 to avoid splitting output with further newlines.
  // You can also re-root the tree for output by specifying an alternate root node.
  void write (ostream& out, int max_columns = 50, int root_for_output = -1) const
    {
      if (root_for_output < 0) root_for_output = root;
      int columns = 0;
      write_node (out, root_for_output, -1, max_columns, columns);
    }

  // write_Stockholm method
  void write_Stockholm (ostream& out) const
  {
    out << Stockholm_file_annotation << ' ';
    out << Stockholm_New_Hampshire_tag << ' ';
    write (out, 0);
  }

  // structures & functions relating to text specification of nodes & branches
  struct Node_specifier_exception : Standard_exception
  {
    Node_specifier_exception (const char* specifier) : Standard_exception ("Bad node specifier: '")
      { msg.chop(); msg.append(specifier).append("'\n"); }
  };

  // find_node methods
  Node find_node (const char* name) const        // exact matches only
  {
    const sstring namestr (name);
    Node_name_vec::const_iterator node_iter = find (node_name.begin(), node_name.end(), namestr);
    return node_iter == node_name.end() ? -1 : node_iter - node_name.begin();
  }

  Node find_node_or_die (const sstring& specifier) const
  {
    Node n = find_node(specifier.c_str());
    if (n < 0)
      THROW Node_specifier_exception(specifier.c_str());
    return n;
  }

  // node specifier methods
  // node_specifier returns the name of a node, or 'root' for the root, or 'A::B' for common ancestor of A & B
  sstring node_specifier (Node n) const;
  // specified_node decodes names in the format returned by node_specifier, or throws a Node_specifier exception
  Node specified_node (const sstring& specifier) const;
  // branch_specifier
  sstring branch_specifier (const Node_pair& b) const;
  sstring directed_branch_specifier (const Node_pair& b) const;
};

// Macros for depth-first traversal of Phylogeny

//  'for_nodes_*' macros include branch from 'PARENT' to 'NODE'
//  'for_branches_*' macros exclude this branch
// 'Branch_iter' is the iterator class in both cases.
// To use 'for_nodes_*' as a pure node iterator, take the current node as 'ITER.second'.

#define for_branches_pre(TREE,NODE,PARENT,ITER)  for_iterator (Phylogeny::Branch_iter, ITER, (TREE).branches_begin ((NODE), (PARENT), 1, 0), (TREE).branches_end())
#define for_branches_post(TREE,NODE,PARENT,ITER) for_iterator (Phylogeny::Branch_iter, ITER, (TREE).branches_begin ((NODE), (PARENT), 0, 0), (TREE).branches_end())
#define for_nodes_pre(TREE,NODE,PARENT,ITER)     for_iterator (Phylogeny::Branch_iter, ITER, (TREE).branches_begin ((NODE), (PARENT), 1, 1), (TREE).branches_end())
#define for_nodes_post(TREE,NODE,PARENT,ITER)    for_iterator (Phylogeny::Branch_iter, ITER, (TREE).branches_begin ((NODE), (PARENT), 0, 1), (TREE).branches_end())

// 'for_neighbors' iterator visits neighbours of 'NODE'

#define for_neighbors(TREE,NODE,ITER)            for_iterator (Phylogeny::Relative_iter, ITER, (TREE).relatives_begin (NODE), (TREE).relatives_end(NODE))
#define for_neighbours(TREE,NODE,ITER)           for_neighbors((TREE),(NODE),(ITER))

// 'for_children' iterator visits neighbours of 'NODE', excluding 'PARENT'

#define for_children(TREE,NODE,PARENT,ITER)      for_iterator (Phylogeny::Child_iter, ITER, (TREE).children_begin ((NODE), (PARENT)), (TREE).children_end((NODE),(PARENT)))

// Macros for traversal of rooted trees
//  'for_nodes_*' macros traverse all branches, including branch from '-1' to 'root'
//  'for_branches_*' macros ignore this parent branch
// 'Branch_iter' is the iterator class in both cases.
// To use 'for_nodes_*' as a pure node iterator, simply take the current node as 'ITER.second'.

#define for_rooted_branches_pre(TREE,ITER)  for_branches_pre ((TREE), (TREE).root, -1, ITER)
#define for_rooted_branches_post(TREE,ITER) for_branches_post ((TREE), (TREE).root, -1, ITER)
#define for_rooted_nodes_pre(TREE,ITER)     for_nodes_pre ((TREE), (TREE).root, -1, ITER)
#define for_rooted_nodes_post(TREE,ITER)    for_nodes_post ((TREE), (TREE).root, -1, ITER)
#define for_rooted_children(TREE,NODE,ITER) for_children ((TREE), (NODE), (TREE).parent[(NODE)], ITER)

#endif
