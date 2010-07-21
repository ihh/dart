#ifndef COMPOSITE_PATH_H
#define COMPOSITE_PATH_H

#include<iostream>
#include<vector>
#include<map>

#include "protpal/utils.h"
#include "protpal/transducer.h"
#include "protpal/profile.h"
#include "protpal/Q.h"
#include "tree/phylogeny.h"
using namespace std; 

// a map from a pair of states to the expected number of times this trans was used
// this is what we use to store expected counts: 
typedef map<pair<int, int>, double> ExpCount; 

// forward declaration
class CompositeState; 

class CompositePath
{
 public:
  PHYLIP_tree tree; 
  CompositePath(M_id, M_id, PHYLIP_tree); 
  vector<CompositeState> path; 
  vector<CompositeState>::iterator pathIter; 
  
  ExpCount get_counts(void); 
  map<node, AbsorbingTransducer>* profiles; 
 private: // most private methods are called within get_counts
  // Do the state-path-explosion operation, start to finish.  After this, we can read off the counts
  // relatively easily.  
  void explode(void); 
  
  // expand the state path at a given node starting at a given state.  Returns the number of new
  // states that were inserted after the given state
  vector<CompositeState> expand(int, node);  
  
  
};


class CompositeState
{
 public:
  PHYLIP_tree* tree; 
  map<node, M_id> profile_states; 
  map<node, state> component_states; 
  
  map<node, AbsorbingTransducer>* profiles; 
  
  // Some basic methods 
  // Set a given node to have states determined by a given M_id
  void set_node(M_id m, node n);
  
  // recursively determine profile states given a starting point w/ and M_id
  void get_profile_states(M_id start, node n);
  

  // Determine transducer states for each node, via a map from node to transducers (allowing for different transducers on each branch)
  void get_component_states(map<node, BranchTrans> transducers);
  
  // Constructors
  // Make an empty state  from a tree
  CompositeState(PHYLIP_tree& tree); 
  
  // Make a state from a tree and a previously-existing state.  State assignments (profile states only) are copied up until node n)
  CompositeState(CompositeState&, node n); 
  
  // display 
  void display(ostream& out);
 private:
  bool prof_states_known; 
  bool component_states_known; 
};
  
    
  
#endif
