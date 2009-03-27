#include <fstream>
#include <deque>
#include "tree/phylogeny.h"

// Tree_shuffler class offers up branches & nodes in a random order
class Tree_shuffler
{
  // typedefs
public:
  typedef Phylogeny::Branch Branch;
  typedef Phylogeny::Node Node;

  typedef Phylogeny::Branch_vector Branch_vector;
  typedef Phylogeny::Node_vector Node_vector;

  // data
private:
  // Pools sorted by Fisher-Yates shuffle
  // (thanks to Kevin Karplus for this algorithm)
  Branch_vector branch_pool;
  Node_vector   node_pool;

  const PHYLIP_tree* tree;

  void reset_branch_pool();
  void reset_node_pool();

public:
  double node_realign_rate;
  double branch_realign_rate;
  double node_slide_rate;
  double branch_swap_rate;
  double branch_scale_rate;
  double indel_param_sampling_rate;

  enum Action { Sample_node = 0, Sample_branch, Slide_node,
		Scale_branch, Swap_branches, Sample_indel_params, Total_moves };

  void set_tree (const PHYLIP_tree& tree);
  void reset();  // call after a successful branch-swapping move
  virtual Action next_action();

  Branch next_branch();
  Branch next_internal_branch();

  Node next_node();
  Node next_internal_node();

  void get_next_slide (Node& grumpa, Node& dad, Node& son);
  void get_next_swap (Node& aunt, Node& nephew, Node& grumpa, Node& dad);

  void cue (const sstring& node_name);  // adds named node to node_pool, and parent-->node branch to branch_pool
  
  Tree_shuffler();
  virtual ~Tree_shuffler() { }
};
