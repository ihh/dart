#ifndef SUFFIX_TREE_INCLUDED
#define SUFFIX_TREE_INCLUDED

#include "seq/biosequence.h"

struct Word_count : map<vector<int>,int>
{
  const int word_len;
  Word_count (const Sequence_database_index& index, int word_len);
  int prefix_count (const int* prefix_begin, const int* prefix_end) const;
};

class Suffix_tree
{
private:
  bool allow_long_branches;
  void delete_nodes();

public:
  struct Node
  {
    int count;
    vector<int> data;
    vector<int> child;
    Node() : count(0) { }
  };
  vector<Node*> node;

  void clear();
  void swap (Suffix_tree& tree);

  void add_suffix (const int* suff_begin, const int* suff_end);
  void add_suffices (const vector<int>& seq_data, int max_len = 0);
  int find_node (const int* suff_begin, const int* suff_end);
  int suffix_count (const int* suff_begin, const int* suff_end);

  void set_null_metascores (const Sequence_database_index& index, int metasc_idx, int null_order, int alphabet_size);

  // constructors
  Suffix_tree (bool allow_long_branches = TRUE);
  Suffix_tree (const Sequence_database_index& index, int max_len = 0, int min_len = 1, int min_count = 1, bool allow_long_branches = TRUE);
  ~Suffix_tree();

  // builder
  void build (const Sequence_database_index& index, int max_len, int min_len, Word_count* word_count, int min_count);

  // postorder iterator over Node's
  class Postorder
  {
  public:
    typedef vector<int>::const_iterator Child_iter;

  private:
    const Suffix_tree* tree;
    vector<Child_iter> child_iter;

    void delve() { while ((**this).child.size()) child_iter.push_back ((**this).child.begin()); }

  public:
    int suffix_len() const
      {
	int len = 0;
	for_const_contents (vector<Child_iter>, child_iter, iter) len += tree->node[**iter]->data.size();
	return len;
      }
    vector<int> suffix() const
      {
	vector<int> s;
	s.reserve (suffix_len());
	for_const_contents (vector<Child_iter>, child_iter, iter)
	  s.insert (s.end(), tree->node[**iter]->data.begin(), tree->node[**iter]->data.end());
	return s;
      }

    int node_idx() const { return *child_iter.back(); }
    const Node& operator*() const { return *tree->node[*child_iter.back()]; }
    bool operator== (const Postorder& i2) const { return child_iter == i2.child_iter; }
    bool operator!= (const Postorder& i2) const { return !(*this == i2); }
    Postorder& operator++()
      {
	const int parent = child_iter.size() > 1 ? *child_iter[child_iter.size()-2] : 0;
	if (++child_iter.back() == (*tree).node[parent]->child.end())
	  child_iter.pop_back();
	else
	  delve();
	return *this;
      }
    //    Postorder operator++(int) { Postorder tmp = *this; ++*this; return tmp; }

    Postorder (const Suffix_tree* _tree) : tree (_tree)
      {
	if (tree)
	  if (tree->node.size())
	    if (tree->node[0]->child.size())
	      {
		child_iter.push_back (tree->node[0]->child.begin());
		delve();
	      }
      }
  };

  Postorder begin() { return Postorder (this); }
  Postorder end() { return Postorder (0); }

  // debugging output method
  void dump (ostream& o, const Alphabet* alph = 0) const;
};

#endif /*  SUFFIX_TREE_INCLUDED */
