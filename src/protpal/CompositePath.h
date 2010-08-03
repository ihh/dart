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
typedef map<pair<state, state>, double> ExpCount; 

// forward declaration
class CompositeState; 

class CompositePath
{
 public:
  PHYLIP_tree tree; 
  double postProb; 
  CompositePath(M_id, M_id, PHYLIP_tree, double, map<node, AbsorbingTransducer>& profiles_in, QTransducer Q_in, bool logging_in = false); 
  vector<CompositeState> path; 
  vector<CompositeState>::iterator pathIter; 
  
  void get_counts(ExpCount&); 
  ExpCount count(int, int); 

  map<node, AbsorbingTransducer>* profiles; 
 private: // most private methods are called within get_counts
  // Do the state-path-explosion operation, start to finish.  After this, we can read off the counts
  // relatively easily.  
  void explode(void); 
  bool logging; 
  // expand the state path at a given node starting at a given state.  Returns the number of new
  // states that were inserted after the given state
  vector<CompositeState> expand(int, node);  
  bool component_changed(int, int, node); 
  QTransducer Q; 
};


class CompositeState
{
 public:
  PHYLIP_tree* tree; 
  QTransducer* Q; 
  map<node, float> profile_indices; 
  map<node, M_id> profile_states; 
  map<node, state> component_states; 
  int placeholder_wait; 
  
  map<node, AbsorbingTransducer>* profiles; 
  
  // Some basic methods 
  // Set a given node to have states determined by a given M_id
  void set_node(M_id m, node n);
  
  // Find the active node - the highest (closest to root) node below which all branches are
  // synchronized
  node active_node(void); 
  bool are_synced_nodes(node, node); 
  bool are_synced_types(string, string); 
  bool synced_below(node); 
  // recursively determine profile states given a starting point w/ and M_id
  void get_profile_states(M_id start, node n, bool logging=false);
  // Access  the transducer state at a given node/branch
  state node_state(node); 

  // Access  the transducer state type at a given node/branch
  string state_type(node); 

  // Determine transducer states for each node, via a map from node to transducers (allowing for different transducers on each branch)
  void get_component_states(map<node, BranchTrans> transducers);
  
  // Constructors
  // Make an empty state  from a tree
  CompositeState(PHYLIP_tree& tree, map<node, AbsorbingTransducer>*, QTransducer*); 
  
  // Make a state from a tree and a previously-existing state.  State assignments (profile states only) are copied up until node n)
  CompositeState(CompositeState&, node n, map<node, AbsorbingTransducer>*, QTransducer*); 
  
  // display 
  void display(ostream&, QTransducer&);
 private:
  bool prof_states_known; 
  bool component_states_known; 
};
  
    
  
#endif
