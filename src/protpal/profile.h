#ifndef PROFILE_H
#define PROFILE_H
#include<iostream>
#include<set>
#include<map>
#include<vector>
#include<queue>

#include "protpal/utils.h"
#include "protpal/Q.h"
#include "protpal/exactMatch.h"


using namespace std;

// forward declarations
class Profile; 
class M_id;


// ***** The 'mature' profile class (e.g. E_n) *****
class AbsorbingTransducer 
{
 public:

  // Constructs an absorbing transducer from a (sampled) profile object, an exact-match transducer, or void input. 
  AbsorbingTransducer(Profile *sampled_profile);
  AbsorbingTransducer(ExactMatch *EM_in);
  AbsorbingTransducer(void);   

  //Basic info about the transducer
  node treeNode; 
  vector<node> subtreeNodes; 
  string name;
  vector<string> alphabet;
  int alphabet_size;   

  // Number of delete states.  If in SNF, there will be num_delete_states + 3 states in the transducer
  int num_delete_states;

  // Keep track of these directly.  Assume integers in 0 <= i <  num_states  are delete states.  
  state start_state, pre_end_state, end_state; 
  //vector<state> pre_end_states;
  // Index to M-id map
  map<state, M_id> state2mid; 
  
  // accessor for incoming states
  vector<state> get_incoming(state);
  
  // query whether or not a given transition exists in the transducer
  bool has_transition(state e, state ePrime);

  //displaying states, transitions, absorption weights
  void display(bool states, bool transitions, bool absorptions);
  void showDOT(Profile *sampled_profile);
  map<vector<int>, int> mid2int;                                                                                                                                    


  // Transition/absorbing accessor functions
  bfloat get_transition_weight(state e, state ePrime);
  bfloat get_absorb_weight(state e, int charIndex);
  
  // state-type-phylogeny. stores binary ((wait, null, end), (match, delete)) state assignments for each state
  map<state, map<node, string> > state_type_phylogeny; 

  // between stores alignment columns in between said external states
  map<pair <state, state>, map< node, string> > between;
  map<pair <state, state>, vector<M_id> > summed_nulls;
  
  // leaf coordinates.  Each state maps to a pair holding the smallest and largest leaf coordinate accounted
  // for by this state.                                                                                 
  map<state, pair<int,int> > leaf_coords;

  //testing
  void test_transitions(void);



 private:
  // Constructor uses the following private functions to transform profile -> absorbing :
  void identify_states(Profile *sampled_profile);
  void marginalize_null_states(Profile *sampled_profile);
  void index_delete_states(Profile *sampled_profile);

  // Constructor uses the following data to transform profile -> absorbing :
  int sampled_start_state; 
  int sampled_pre_end_state;

  // vector<int> sampled_pre_end_states;  

  int sampled_end_state;   
  vector<string> sampled_state_types;
  vector<int> sampled_externals;
  map< vector<int>, vector<M_id> > sampled_externals_outgoing;
  map<pair< vector<state> , vector<state> >, bfloat > transition_weight_tmp; 
  map< vector<int>, vector<M_id> > incoming_tmp; 
  // possible sampled_state types: ("H-types")
  string h_start;
  string h_left_internal;
  string h_right_internal;
  string h_external;
  string h_pre_end;  
  string h_end;
  


  // raw data maps
  map< state, vector<state> > incoming;
  map< state, vector<double> >  absorption_weight; 
  map< pair<state, state>, double > transition_weight;
  
  //various utility functions
  bool checkState(state e);
  bool isDelState(state e);
};



// The states in Profile s are identified by M_id:

class M_id
{
 public:
  state q_state;
  state left_state;
  state right_state;
  
/*   string left_type; */
/*   string right_type; */
  int left_type;
  int right_type;

  void display(QTransducer &);
  int operator==(const M_id &);
  int operator!=(const M_id &);
  vector<int> toVector(void);
};

int index(M_id query, vector<M_id> in );
bool contains(M_id child, queue<M_id> stateQueue);
bool contains(M_id child, deque<M_id> stateQueue);
bool in(M_id query, vector<M_id> in);

// ***** The 'immature' profile class (e.g. M_n) *****
class Profile
{
 public:
  string name;  
  // Start, end, pre_end states explicitely stored as in simpler classes, but now as M_ids
  M_id start_state, pre_end_state, end_state;
  //  vector<M_id> pre_end_states; 

  // constructor: set left and right profiles, and node on the tree that this profile corresponds to
  Profile(node node_in, AbsorbingTransducer left_in, AbsorbingTransducer right_in, QTransducer Q_in);

  // Initiate and fill forward-like DP matrix
  void fill_DP(int, bool inLog=false);
  map< vector<int>, vector<M_id> > incoming; // this is filled only if the backward algorithm is requested
  map< vector<int>, vector<M_id> > outgoing; // this is filled only if the backward algorithm is requested
  map< pair<vector<int>, vector<int> >, bfloat > transition_weight; // this is filled only if the backward algorithm is requested
  
  bfloat forward_prob; // This is always  stored, for postprob calculations

  // Backward-like DP matrix methods
  void fill_backward_DP(int logging); 

  void clear_DP(void); 
  // sample set of paths from DP matrix, storing info on their transitions, transition weights, etc
  string sample_DP(int, int, bool, bool);
  void cache_state(M_id m, M_id mPrime, bfloat weight);
  void cache_path(vector<M_id>); 
  void store_summed_nulls(vector<M_id>); 

  vector<M_id> sampled_states;
  int num_sampled_externals; 
  int max_sampled_externals; 
  map< vector<int>, vector<M_id> > sampled_incoming;
  map< vector<int>, vector<M_id> > sampled_outgoing;  
  map< pair< vector<int>, vector<int> >, bfloat> sampled_transition_weight; 
  bfloat get_sampled_transition_weight(M_id m, M_id mPrime);
  bfloat get_internal_cascade_weight(M_id m);
  bfloat get_external_cascade_weight(M_id m, int charIndex);  
  bfloat compute_emission_weight(M_id m);

  bool is_start(M_id);
  bool is_external(M_id);
  bool is_right_int(M_id); 
  bool is_left_int(M_id);
  bool is_pre_end(M_id);     
  bool is_end(M_id);
  
  // after sampling, the profile is ready to be handed over the AbsorbingTransducer constructor to be
  // turned into a new absorbing child profile.  the following data is public for this purpose
  vector<string> alphabet; 
  int alphabet_size;
  vector<int> alphabet_ints; 
  vector<int>::iterator alphIter; 
  vector<int>::iterator alphIter2; 
  
  // Coordinate/tree information
  // The node on the tree where this profile lives
  node treeNode;   
  vector<node> subtreeNodes; 
  // vector mapping tree nodes to common names
  vector<string> node_names; 
  // leaf nodes below this node (at present this is all leaf nodes...possibly exists ways to get this w/phylogeny class)
  vector<Node> leaves; 

  // state-type-phylogeny. stores the state type for each state and node in the tree              
  map< vector<int>, map<node, string> > state_type_phylogeny;
  map<node, string> merge_STP(M_id);
  
  // between stores alignment columns in between said external states
  map< pair< vector<int>, vector<int> >, map<node, string> > between;
  map< pair< vector<int>, vector<int> >, vector<M_id> > summed_nulls;
  
  
  // leaf coordinates.  Each state maps to a pair holding the smallest and largest leaf coordinate accounted 
  // for by this state.
  map< vector<int>, pair<int, int> > leaf_coords;

  // Simple envelope criterion - leaf sequences must be within envelope_distance of each other
  int envelope_distance; 

  // Qtransducer object
  QTransducer Q;
  int get_profile_type(state, string); 

  // Ordering of states
  deque<M_id> backward_states;
 private:

  // Left and right (mature) profile objects
  AbsorbingTransducer left_profile;
  AbsorbingTransducer right_profile;

  // **** Fill_DP uses the following: ****
  // the actual DP matrix.  Maps states in profile (as M_id ) to double
  //  map<M_id, double> Z; //eventually a 'bigger' number class?
  map< vector<int>, bfloat> Z; //fix this up soon!
  map< vector<int>, bfloat> backward_matrix; 

  //accessor function, returns 0 if no matching entry
  bfloat get_DP_cell(M_id); 
  bfloat get_backward_DP_cell(M_id); 

  
  // Add to a DP cell - first checks if the cell exists. 
  void add_to_DP_cell(M_id, bfloat); 
  void add_to_backward_DP(M_id, bfloat); 

  // sum over all paths to the M_id cell in the forward matrix
  void sum_paths_to(M_id, bool);

  //Access functions
  // compute/access the transition weight from state M_id m  to M_id mPrime
  double compute_transition_weight(M_id m , M_id mPrime);
  double get_transition_weight(M_id m , M_id mPrime);
  
  map< pair< vector<int>, vector<int> >, bfloat> transition_weight_test; 
  
  // compute/access the emission cascade weight from state m
  
  double get_emission_weight(M_id m);
  vector<bfloat> tmpEmitVals; 
  vector<int> tmpEmitTuple; 
  vector<int> tmpMidVec;
  
  // Utility functions
  bool is_in_envelope(state left_state, state right_state); 

  void path2alignment(vector<M_id> &);
  string show_alignment(vector<M_id> &, bool);
  void show_state(M_id);


  // Utility variables
  int profile_start;
  int profile_wait;
  int profile_delete;
  int profile_pre_end;
  int profile_end;

};

void show_state_phylo(map<node, string> &);    
map<node, string> cat_STP(map<node, string>, map<node, string>);
map<node, string> pad_STP(map<node, string>, vector<node>);
bool is_flush(map<node, string> map1);
pair<int, int> merge_coords(pair<int, int>, pair<int,int>); 
#endif
