#ifndef ALIGNMENT_SAMPLER_H
#define ALIGNMENT_SAMPLER_H
#include <queue>
#include <vector>
#include <list>
#include <stack>

#include "protpal/utils.h"
#include "protpal/algebras.h"
#include "protpal/profile.h"
#include "tree/phylogeny.h"
#include "util/sstring.h"

typedef vector< map<node, string> > colMap; 


class AlignmentSampler
{
 public:
  AlignmentSampler(state_path, node, map<node, Profile>*, map<node, AbsorbingTransducer>*, PHYLIP_tree*);
  AlignmentSampler(void); 
  // sample all subtree alignments, filling columns along the way
  void sample_all(bool viterbi=false, bool logging_in =false );
  // Display the alignment in stockholm form
  void display(ostream&, string format="stockholm"); 
 private:
  int profile_start, profile_delete, profile_pre_end;
  bool logging, have_restored_nulls; 
  node treeNode; 
  PHYLIP_tree* tree; 
  // Path of transducer states
  state_path path; 
  // Alignment data, with alphabet characters converted back to strings.  Suitable for displaying.
  colMap columns; 
  colMap::iterator colIter; 
  // verify that all rows in the alignment have the same length
  bool is_flush(void); 
  // given a subtree alignment, verify that it is flush, then pad remaining nodes (as determined by the tree) with the appropriate gaps
  void pad(void); 
  // Add all the columns from another alignment to the current alignment, inserting them before the given position. 
  void add_cols(unsigned int, colMap); 

  map<node, Profile>* profiles; 
  map<node, AbsorbingTransducer>* AbsProfiles; 
  // sample a path through the transducer rooted at node, starting with the state at path[i][node], and ending in path[j][node]
  colMap sample_expanded_path(node, M_id, M_id, bool viterbi=false); 
  colMap sample_expanded_path(node, int, int, bool viterbi=false); 
  M_id find_last_del_state(const state_path& pi, int pathIdx, string side);

};

class IndelCounter
{
  // Some functions and variables related to tabulating indel information
 public:
  IndelCounter(Stockholm&, PHYLIP_tree*);
  void gather_indel_info(bool logging=false); 
  void display_indel_info(ostream&, bool per_branch=false); 

  double avg_insert_rate, avg_delete_rate, avg_insert_ext, avg_delete_ext; 
 private:
  PHYLIP_tree* tree; 
  map<node, string> rows; 
  map<node, string>:: iterator rowIter; 
  unsigned int L; // length of rows
  map<node, vector<string> > insertions, deletions; 
  map<node, int> matches; 

  double insert_rate(node); 
  double delete_rate(node); 

  double insert_extend(node); 
  double delete_extend(node); 
  
  void average_indel_counts(void); 
  void get_pairwise_alignment(node parent, node child, string& parentString, string& childString); 
  void all_insertions(ostream&); 
  void all_deletions(ostream&); 
};


class DP_container
{
 private:
  map< vector<int>, bfloat> data; 
 public:
  DP_container(void); 
  void set(M_id, bfloat);
  void add(M_id, bfloat);
  void multiply(M_id, bfloat);
  bfloat get(M_id);
};

bool all_parents_visited(vector<M_id> parents, vector<M_id> visited);
map<node, string>  padColumn(map<node, string>& column, PHYLIP_tree* tree );
  
  


#endif



