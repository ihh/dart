#include <stdlib.h>
#include <string.h>
#include "util/macros.h"
#include "tree/phylogeny.h"
#include "util/logfile.h"
#include "util/nstring.h"

#define CLOGTOP CLOG(1)

Phylogeny::Phylogeny() : _branch_count(0), _first_internal_node(_sorted_nodes.end()) {}

Phylogeny::Phylogeny (const Phylogeny& p)
{
  *this = p;
}

Phylogeny& Phylogeny::operator= (const Phylogeny& p2)
{
  Phylogeny& p ((Phylogeny&) p2);   // cast away fucking const
  _branch_length = p._branch_length;
  _branch_count = p._branch_count;
  _sorted_nodes = p._sorted_nodes;
  _first_internal_node = (p._first_internal_node  - p._sorted_nodes.begin()) + _sorted_nodes.begin();
  _node_sort_pos = p._node_sort_pos;
  return *this;
}

Phylogeny::Node Phylogeny::new_node()
{
  int n = nodes();
  _branch_length.push_back (Phylogeny::Branch_length_map());
  _sorted_nodes.push_back (n);
  _node_sort_pos.push_back (_sorted_nodes.end() - 1);
  return n;
}

void Phylogeny::add_branch (Node node1, Node node2, double length)
{
  if (node1 > node2) swap (node1, node2);

  for (Node n = nodes(); n <= node2; n++) new_node();

  _branch_length[node1][node2] = length;
  _branch_length[node2][node1] = -1;
  _branch_count++;

  sort_nodes();
}

struct IsLeafFunction : unary_function<int,bool>
{
  const Phylogeny& tree;
  IsLeafFunction (const Phylogeny& tree) : tree(tree) { }
  result_type operator() (argument_type node) const { return tree.is_leaf (node); }
};

void Phylogeny::sort_nodes()
{
  const IsLeafFunction isLeaf (*this);
  _first_internal_node = stable_partition(_sorted_nodes.begin(),
					  _sorted_nodes.end(),
					  isLeaf);
  sort(_sorted_nodes.begin(), _first_internal_node);
  sort(_first_internal_node, _sorted_nodes.end());
  for (Node_iter i = _sorted_nodes.begin(); i != _sorted_nodes.end(); i++)
    _node_sort_pos[*i] = i;
}

void Phylogeny::remove_branch (Node node1, Node node2)
{
  bool l1b = is_leaf(node1);
  bool l2b = is_leaf(node2);

  _branch_length[node1].erase(node2);
  _branch_length[node2].erase(node1);
  _branch_count--;

  bool l1a = is_leaf(node1);
  bool l2a = is_leaf(node2);

  if (l1a != l1b || l2a != l2b)
    sort_nodes();
}


void Phylogeny::clear()
{
  _branch_length.clear();
  _sorted_nodes.clear();
  _node_sort_pos.clear();
  _branch_count = 0;
  _first_internal_node = _sorted_nodes.end();
}



Phylogeny::Node_directory Phylogeny::node_directory (const Node_name_vec& node_name)
{
  Node_directory node_dir;
  for (Node node = 0; node < (int) node_name.size(); node++)
    if (node_name[node].size())
      node_dir[node_name[node]] = node;
  return node_dir;
}


Phylogeny::Node_map Phylogeny::node_remapper (const Node_name_vec&  node1_name,
					      const Node_directory& node2_directory)
{
  Node_map node1_to_node2;
  for (Node node1 = 0; node1 < (int) node1_name.size(); node1++)
      if (node1_name[node1].size())
	{
	  Node_directory::const_iterator node2_iter = node2_directory.find(node1_name[node1]);
	  if (node2_iter != node2_directory.end())
	    node1_to_node2[node1] = (*node2_iter).second;
	}
  
  return node1_to_node2;
}


Phylogeny::Node_set Phylogeny::remap_node_set (const Node_set& node_set, const Node_map& node_remapper)
{
  Node_set remapped;
  for_const_contents (Node_set, node_set, i) remapped.insert (((Node_map&) node_remapper) [*i]);
  return remapped;
}


Phylogeny::Node Phylogeny::youngest_common_ancestor (Node n1, Node n2, const Node_vector& parent) const
{
  Node_vector a1, a2;
  while (1) { a1.push_back(n1); if (n1 == -1) break; n1 = parent[n1]; }
  while (1) { a2.push_back(n2); if (n2 == -1) break; n2 = parent[n2]; }
  int i;
  for (i = 1; i <= (int) min(a1.size(),a2.size()); i++)
    if (a1[a1.size() - i] != a2[a2.size() - i])
      break;
  return a1[a1.size() + 1 - i];
}



Phylogeny::Node_vector Phylogeny::node_path (Node n1, Node n2, const Node_vector& parent) const
{
  Node_vector result (1, n1);
  Node a = youngest_common_ancestor (n1, n2, parent);
  while (n1 != a) { n1 = parent[n1]; result.push_back(n1); }
  Node_vector n2_to_a;
  while (n2 != a) { n2_to_a.push_back(n2); n2 = parent[n2]; }
  copy (n2_to_a.rbegin(), n2_to_a.rend(), back_inserter(result));
  return result;
}


void Phylogeny::swap_nodes (Node n1, Node n2)
{
  // loop through all neighbors of n1
  for_contents (Branch_length_map, _branch_length[n1], bl)
    {
      Branch_length_map& neighbour_map = _branch_length[(*bl).first];
      Branch_length_map::iterator i = neighbour_map.find (n2);
      if (i != neighbour_map.end())   // if n2 is also a neighbor of this node, then swap n1 & n2 lengths
	swap (neighbour_map[n1], (*i).second);
      else
	{
	  // n2 is not a neighbor of this node, so make it one, & remove n1
	  neighbour_map[n2] = neighbour_map[n1];
	  neighbour_map.erase(n1);
	}
    }

  // loop through all neighbors of n2
  for_contents (Branch_length_map, _branch_length[n2], bl)
    {
      Branch_length_map& neighbour_map = _branch_length[(*bl).first];
      Branch_length_map::iterator i = neighbour_map.find (n1);
      // if n1 is not a neighbor of this node, then make it one, & remove n2
      if (i == neighbour_map.end())
	{
	  neighbour_map[n1] = neighbour_map[n2];
	  neighbour_map.erase(n2);
	}
    }

  // swap n1 & n2 neighbor lists
  swap (_branch_length[n1], _branch_length[n2]);
  sort_nodes();
}

Phylogeny::Node Phylogeny::add_node (Node neighbour, double branch_length)
{
  Node n = new_node();
  if (neighbour >= 0)
    add_branch (n, neighbour, branch_length);
  sort_nodes();
  return n;
}

Phylogeny::Node_vector Phylogeny::delete_node (Node n)
{
  // calculate old-->new mapping
  const int n_nodes = nodes();
  Node_vector old_new_mapping (n_nodes);
  for (Node x = 0; x < n_nodes; ++x)
    old_new_mapping[x] = x > n ? (x-1) : x;
  old_new_mapping[n] = -1;

  // remove branches
  Node_vector rellies = relatives(n);
  if (rellies.size() == 1)
    remove_branch (n, rellies[0]);
  else if (rellies.size() == 2)
    {
      double total_length = branch_length(n,rellies[0]) + branch_length(n,rellies[1]);
      remove_branch (n, rellies[0]);
      remove_branch (n, rellies[1]);
      add_branch (rellies[0], rellies[1], total_length);
    }
  else
    THROW Standard_exception ("Attempt to delete a node with more than two neighbours");

  // update _branch_length
  _branch_length.erase (_branch_length.begin() + n);
  for_contents (vector<Branch_length_map>, _branch_length, blm)
    {
      Branch_length_map new_blm;
      for_contents (Branch_length_map, *blm, old_bl)
	{
	  const int nbr = old_bl->first;
	  const double len = old_bl->second;
	  if (nbr != n)
	    new_blm[old_new_mapping[nbr]] = len;
	}
      new_blm.swap (*blm);
    }

  // re-sort nodes
  _sorted_nodes.pop_back();
  _node_sort_pos.pop_back();
  for (Node i = 0; i < n_nodes - 1; ++i)
    _sorted_nodes[i] = i;
  sort_nodes();

  // return old-->new mapping
  return old_new_mapping;
}

void Phylogeny::move_node (Node n, Node old_neighbour, Node new_neighbour, double new_branch_length)
{
  remove_branch (n, old_neighbour);
  add_branch (n, new_neighbour, new_branch_length);
}

PHYLIP_tree::PHYLIP_tree() :
  root (-1)
{ }

Regexp name_regexp ("[ \t\n]*([^ \t\n:\\(\\)]*)[ \t\n]*");
Regexp length_regexp (":[ \t\n]*([-0123456789\\.\\+e]+)[ \t\n]*");
Regexp common_regexp ("([^:]+)::([^:]+)");
Regexp parent_regexp ("([^:]+):([0123456789]+)");

Phylogeny::Node PHYLIP_tree::read_node (istream& in,
					Node parent_node)       // returns index of created node
{
  Node node = new_node();
  while ((int) node_name.size() <= node) node_name.push_back(sstring());
  while ((int) parent.size() <= node) parent.push_back(-1);
  parent[node] = parent_node;
  
  char c;
  do { in.get(c); } while (isspace(c) && in && !in.eof());

  bool is_parent = (c == '(');
  if (is_parent)
    {
      do
	{
	  read_node (in, node);
	  if (in.eof()) THROW Format_exception(in,"Unexpected end of file, while expecting comma or closing parenthesis"); else in >> c;
	  if (c != ',' && c != ')') THROW Format_exception(in,"Expected comma or closing parenthesis");
	}
      while (c != ')');
      do { in.get(c); } while (isspace(c) && in && !in.eof());
    }

  sstring name;
  while (c != ':' && c != ';' && c != ',' && c != ')')
    {
      name.push_back(c);
      if (in.eof()) THROW Format_exception(in,"Unexpected end of file, while reading node name");
      in >> c;
    }
  in.putback(c);

  if (!(name_regexp.Match(name.c_str()) && name_regexp.SubStrings() == 1))
    THROW Format_exception(in,"Encountered something that does not look like a node name");
  
  node_name[node] = name_regexp[1];
  if (node_name[node].size() == 0 && !is_parent)
    THROWEXPR ("Leaf nodes must have a name");

  if (parent_node == -1)
    {
      do
	{
	  if (in.eof()) THROW Format_exception(in,"Unexpected end of file, while expecting semicolon");
	  in >> c;
	  if (!isspace(c) && c != ';') THROW Format_exception(in,"Expected a semicolon");
	}
      while (c != ';');
    }
  else
    {
      sstring len;
      double l = -1;  // missing branch lengths are represented as -1.... ugh
      while (true)
	{
	  if (in.eof()) THROW Format_exception(in,"Unexpected end of file, while trying to parse branch length");
	  in >> c;
	  if (c == ',' || c == ')') {
	    in.putback(c);
	    break;
	  }
	  len.push_back(c);
	}
      
      if (len.size()) {
	if (!(length_regexp.Match(len.c_str()) && length_regexp.SubStrings() == 1))
	  THROW Format_exception(in,"Encountered something that does not look like a branch length");
	l = atof(length_regexp[1].c_str());
      }

      add_branch (node, parent_node, l);
    }
  
  return node;
}

void PHYLIP_tree::write_node(ostream& out,
			     Node node,
			     Node from_node,
			     int max_columns,
			     int& columns,
			     bool newline, // Added 7/22/11 OW
			     bool numbered_branches) const
{
  Child_iter c = children_begin(node, from_node);
  Child_iter end = children_end(node, from_node);
  if (c != end)
    {
      out << '(';
      while (c != end)
	{
	  write_node (out, *c, node, max_columns, columns, newline, numbered_branches);
	  if (++c == end)
	    out << ')';
	  else
	    {
	      out << ',';
	      if (columns >= max_columns && max_columns > 0)
		{
		  out << '\n';
		  columns = 0;
		}
	    }
	}
    }

  const sstring& name = ((Node_name_vec&) node_name) [node];
  out << name;
  columns += name.size();

  if (from_node != -1)
    {
      double l = ((Phylogeny*) this) -> branch_length (from_node, node);
      if (l >= 0) {
	char tmp_buf[128];
	sprintf (tmp_buf, ":%-.10f", l);
	out << tmp_buf;
	if (numbered_branches)
	  out << "[" << node << "]"; 
	columns += strlen(tmp_buf) + 1;
      }
    }
  else
    {
      if (numbered_branches) // Added 7/22/11 OW
	out << "[" << node << "]"; 
      out << ";";
      if (newline) // Added 7/22/11 OW
	out<< "\n";
      columns = 0;
    }
}


Phylogeny::Node PHYLIP_tree::specified_node (const sstring& specifier) const
{
  Node node;
  if ((node = find_node (specifier.c_str())) == -1)
    {
      if (((Regexp&)common_regexp).Match (specifier.c_str()) && ((Regexp&)common_regexp).SubStrings() == 2)
	{
	  Node leaf1 = find_node (common_regexp[1].c_str());
	  Node leaf2 = find_node (common_regexp[2].c_str());
	  if (leaf1 == -1 || leaf2 == -1)
	    THROW Node_specifier_exception (specifier.c_str());
	  node = youngest_common_ancestor (leaf1, leaf2, parent);
	}
      else if (specifier == sstring(DART_ROOT_IDENTIFIER)) node = root;
      else if (((Regexp&)parent_regexp).Match (specifier.c_str()) && ((Regexp&)parent_regexp).SubStrings() == 2)
	{
	  node = find_node (parent_regexp[1].c_str());
	  int d = atoi (parent_regexp[2].c_str());
	  if (node == -1)
	    THROW Node_specifier_exception (specifier.c_str());
	  while (d-- > 0) node = parent[node];
	}
      else THROW Node_specifier_exception (specifier.c_str());
    }
  return node;
}

sstring PHYLIP_tree::node_specifier (Node n) const
{
  sstring specifier;
  Node_vector kids;
  copy (children_begin(n,parent[n]), children_end(n,parent[n]), back_inserter(kids));
  if (kids.size() == 0 || node_name[n].size())
    {
      specifier = node_name[n];
    }
  else if (n == root && find_node(DART_ROOT_IDENTIFIER) == -1)
    {
      specifier = DART_ROOT_IDENTIFIER;
    }
  else if (kids.size() == 1)
    {
      int d = 1;
      n = kids[0];
      while (!is_leaf(n)) { n = *children_begin(n,parent[n]); ++d; }
      char buf[10];
      sprintf (buf, "%d", d);
      specifier.append(node_name[n]).append(":").append(buf);
    }
  else
    {
      Node n1 = kids[0], n2 = kids[1];
      while (!is_leaf(n1)) n1 = *children_begin(n1,parent[n1]);
      while (!is_leaf(n2)) n2 = *children_begin(n2,parent[n2]);
      specifier.append(node_name[n1]).append("::").append(node_name[n2]);
    }
  return specifier;
}

sstring PHYLIP_tree::branch_specifier (const Node_pair& b) const
{
  sstring s;
  s << "'" << node_specifier(b.first) << "'--'" << node_specifier(b.second) << "'";
  return s;
}

sstring PHYLIP_tree::directed_branch_specifier (const Node_pair& b) const
{
  sstring s;
  s << "'" << node_specifier(b.first) << "'-->'" << node_specifier(b.second) << "'";
  return s;
}

void PHYLIP_tree::assert_tree_is_binary() const
{
  for (int n = 0; n < nodes(); n++)
    if (n != root && !is_leaf(n) && neighbours(n) != 3)
      {
	sstring ex;
	ex << "Node '" << node_specifier(n) << "' has " << neighbours(n) << " neighbours (nodes in a binary tree should have three)";
	THROW Standard_exception (ex);
      }
  if (neighbours(root) != 2)
    {
      sstring ex;
      ex << "Root node has " << neighbours(root) << " neighbours (root node in my definition of a binary tree should have two)";
      THROW Standard_exception (ex);
    }
}

void PHYLIP_tree::remap_nodes (const Node_vector& old_to_new_map)
{
  root = old_to_new_map[root];
  Node_name_vec new_node_name (nodes());
  for (int o = 0; o < (int) old_to_new_map.size(); ++o)
    {
      const int n = old_to_new_map[o];
      if (n >= 0)
	new_node_name[n] = node_name[o];
    }
  new_node_name.swap (node_name);
}

bool PHYLIP_tree::force_binary()
{
  // keep track of whether we change the tree
  bool tree_changed = false;

  // prune fake roots
  while (neighbours (root) == 1)
    {
      const Node new_root = *relatives_begin (root);
      Node_vector old_to_new_map = delete_node (root);
      remap_nodes (old_to_new_map);
      root = old_to_new_map[new_root];
      tree_changed = true;
    }

  // make lists of nodes to delete or split
  Node_vector nodes_to_split;
  Node_vector nodes_to_delete;
  for (int n = 0; n < nodes(); n++)
    {
      const int nbrs = neighbours(n);
      if (n == root)
	{
	  if (nbrs < 1) THROW Standard_exception ("Root is isolated; can't make tree into a binary tree");
	  if (nbrs > 2) nodes_to_split.push_back(n);
	}
      else
	{
	  if (nbrs == 2) nodes_to_delete.push_back(n);
	  if (nbrs > 3) nodes_to_split.push_back(n);
	}
    }

  // split marked nodes
  for (int i = 0; i < (int) nodes_to_split.size(); ++i)
    {
      const int n = nodes_to_split[i];
      do
	{
	  Node_vector rellies = relatives(n);
	  const int n2 = add_node (n, 0.);
	  node_name.push_back (sstring());
	  move_node (rellies[0], n, n2);
	  move_node (rellies[1], n, n2);
	}
      while (neighbours(n) > 3);
      tree_changed = true;
    }

  // delete marked nodes
  for (int i = 0; i < (int) nodes_to_delete.size(); ++i)
    {
      const int n = nodes_to_delete[i];
      Node_vector old_to_new_map = delete_node(n);
      for (int j = i+1; j < (int) nodes_to_delete.size(); ++j)
	nodes_to_delete[j] = old_to_new_map[nodes_to_delete[j]];
      remap_nodes (old_to_new_map);
      tree_changed = true;
    }

  // update tree vars
  rebuild_parents();
  assert_tree_is_binary();

  // log if changes
  if (tree_changed)
    CLOG(8) << "Warning: had to split/merge some nodes, to make tree strictly binary\n";

  // return
  return tree_changed;
}

vector<Phylogeny::Node> PHYLIP_tree::leaf_vector() {
  vector<Node> lvec;
  for (Phylogeny::Node_const_iter iter = leaves_begin(); iter != leaves_end(); ++iter)
    if (*iter != root)
      lvec.push_back(*iter);
  return lvec;
}

vector<Phylogeny::Node> PHYLIP_tree::ancestor_vector() {
  vector<Node> avec;
  if (is_leaf(root))
    avec.push_back(root);
  for (Phylogeny::Node_const_iter iter = internals_begin(); iter != internals_end(); ++iter)
    avec.push_back(*iter);
  return avec;
}
