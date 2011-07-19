#ifndef EXACTMATCH_H
#define EXACTMATCH_H
#include<iostream>
#include<set>
#include<map>
#include<vector>
#include<string>

#include "protpal/utils.h"
#include "protpal/algebras.h"

using namespace std;

class ExactMatch
{
 public:
  // Constructor - does everything we need to get it in absorbing form.
  ExactMatch(string &sequence, node treeNode_in, Alphabet& alphabet_in); 

  // The (leaf) node of the phylogeny this sequence is associated with.
  node treeNode; 

  // Name 
  string name; 

  //alphabet
  vector<string> alphabet;
  int alphabet_size;

  // Number of indexed del states.  Here there are just sequence.length() delete states, start, end, and Wait.  
  int num_delete_states;
  
  // explicitely keep track of start, end, pre-end, end states here...
  state start_state, pre_end_state, end_state; 
  // Incoming transitions
  MyMap <state, vector <state> > incoming;

  // Absorption matrix; one entry per state, character
  MyMap <state, vector<bfloat> > absorb; 
  
  // state-type-phylogeny. stores binary ((wait, null, end), (match, delete)) state assignments for each state
  MyMap<state, MyMap<node, string> > state_type_phylogeny;

  // leaf coordinates.  Each state maps to a pair holding the smallest and largest leaf coordinate accounted
  // for by this state. 
  MyMap<state, pair<int, int> > leaf_coords;

  // leaf sequence coordinates , for use in guide alignment-based alignment envelopes
  MyMap<state, MyMap<node, int> > leaf_seq_coords;
  
  
  // Test incoming, absorb maps of this tranducer.  
  void test(void);
};


#endif
