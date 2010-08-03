#ifndef QTransducer_H
#define QTransducer_H
#include<iostream>
#include<map>
#include<vector>
#include "transducer.h"
#include "ecfg/ecfgsexpr.h"

using namespace std;

class QTransducer
{
 public:

  // Constructor handles a variety of tasks, like building composite state space, transitions, etc.
  // Takes as input the transducers R, B_l, B_r, Upsilon, as transducer objects.  Each of these holds its own
  // transitions, weights, emission distributions, etc.  

  QTransducer(SingletTrans, BranchTrans, BranchTrans, SplittingTrans, Alphabet& );
  QTransducer(void); //dummy constructor

  // Ideally we'll only need to build QTransducer once, then reset the branch lengths with this function
  // It sets left, right branch lengths, and also re-computes emission probabilities
  void change_branch_lengths(float left_branch_length, float right_branch_length);

  // Character alphabet
  vector<string> alphabet;
  int alphabet_size;
  
  //Name
  string name; 

  // Recover the state indexes or names of the component transducers (e.g. R, B_l, B_r, Upsilon) for a given q state
  vector<int> get_state_indices(state q);
  string get_state_name(state q);

  // State type-sets for states q in Q
  state composite_start_state; 
  state composite_end_state; 

  // Show all states/transitions (by component names)
  void show_all_states(void);
  void show_all_transitions(void);

  // Retrieve certain state sets 
  vector<state> get_match_states(void);
  vector<state> get_left_emit_states(void);
  vector<state> get_right_emit_states(void);
  vector<state> get_wait_states(void);  
  
  // This one is *slightly* special:
  vector<state> get_SSSI_states(void);  


  // Connectivity information.  Return the vector of states such that q -> qPrime has nonzero weight, and
  // q has the relevant type.  
  vector<state> get_incoming(state qPrime);
  vector<state> get_incoming_match_states(state qPrime);
  vector<state> get_incoming_left_emit_states(state qPrime);
  vector<state> get_incoming_right_emit_states(state qPrime);  

  // Transition function; assume these have been cached, so this is a lookup.
  double get_transition_weight(state q, state qPrime);

  // Factor  this out later
  double get_R_transition_weight(state q, state qPrime);
  
  // query whether or not a given transition exists in the transducer    
  bool has_transition(state q, state qPrime);
  // query the state class of a given q (e.g. match, wait, left_ins, etc)
  string get_state_class(state q);

  // query the state type of a given q (e.g. IMMM, SSSI, etc)
  string get_state_type(state q);   

  //Emission function - emission of left_emit and right_emit characters
  double get_emission_weight(state q, int left_emit, int right_emit);
  double get_emission_weight(state q, string left_emit, string right_emit);  

  // Branch match weights - for computing absorption weight
  double left_match_weight(state q, int absorb, int emit);
  double right_match_weight(state q, int absorb, int emit);
  double right_delete_weight(state q, int absorb);
  double left_delete_weight(state q, int absorb);   

  // Testing functions:
  void test_emissions(void);
  void test_transitions(void);  

  //The state classes:
  string match;
  string wait;
  string null;  
  string right_del;
  string left_del;
  string left_ins;
  string right_ins;


  // adjust these weights for null states
  void marginalizeNullStates(void);

  // The map emission_weight holds these...the vector<int> is state, omega_left, omega_right
  // One of the omegas can be -1, meaning 'no character' (e.g. omega_left = -1 in a left_del state)
  map< vector<int>, double> emission_weight;

  vector<state> get_components(state); 

  // Component transducers R, B_l, B_r, Upsilon
  // Each state in q is a combination of these transducers
  SingletTrans R;
  BranchTrans B_l;
  BranchTrans B_r;
  SplittingTrans Upsilon;
  
 private:  
  // Number of states
  int num_states; 

  // This map has keys according to state 'classes': match, wait, right_del, left_del, right_ins, left_ins
  map< string, vector<state> >  state_class_set;



  // This map has keys according to state 'types': "SSSI", "MMID", etc
  map< string, vector<state> >  state_type_set; 

  // Vector for individual state classes (e.g. match, wait, right_del, etc)
  vector<string> state_class;
  
  // Vector for individual state types (e.g. IMMM, SSSI, SSIW, etc)
  vector<string> state_type;

  // Some hairy-ish functions called by the constructor: 
  // enumerate all the states of Q, via the composite state types list
  void cache_states(void);

  // enumerate transitions of Q, storing the incoming transitions according to type
  void cache_transitions(void);

  // compute/store weights of  transitions of Q
  void cache_transition_weights(void);
  

  
  // First fill these, based on the table
  void cache_in_out_info(void);
  map<state, vector<state> > outgoing;
  map<state, vector<state> > incoming;

  

  // Then fill this map with weights
  map< vector<state>, double> transition_weight;
  double compute_transition_weight(state, state, vector<int>);
  

  // compute/store weights of  emissions of Q
  void cache_emission_weights(void);  
  
  
  // State assignment s (indices) for each states in Q.   The inner vector is always 4 items long.  
  vector< vector<state> > assignments; 

  // The composite state types of Q.  Hard-coded here, though I think this is true for any Q (?)
  // This is taken directly from section 2.2.1 in transducer.tex
  map<string, vector<string> > state_type_reference;

  // The function to initialize this:
  void get_state_type_reference(void);

  // The table storing the transition information for Q
  // each entry is a vector of state-type-tuples (r,u,b_l,b_r), (r',u',b_l',b_r'), 
  // and the value is a vector of integers, telling us how many transitions each component transducer
  // must make for the composite transition to occur.  This is described in the latex document.
  // The table is taken from the latex document, and code is generated with table2cc.py
  map<vector<string>, vector<int> > transition_table;

  // The function to initialize this:
  void get_transition_table(void);  
  

};
 
#endif
