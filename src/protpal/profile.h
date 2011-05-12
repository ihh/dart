#ifndef PROFILE_H
#define PROFILE_H
#include<iostream>
#include<set>
#include<map>
#include<ext/hash_map>
#include<vector>
#include<queue>
#include<list>

#include "protpal/utils.h"
#include "protpal/Q.h"
#include "protpal/exactMatch.h"
#include "protpal/AlignmentEnvelope.h"
#include "hmm/transmat.h"
#include "util/sexpr.h"

using namespace std;
class M_id
{
 public:
  M_id(void); 
  state q_state;
  state left_state;
  state right_state;
/*   string left_type; */
/*   string right_type; */
  int left_type;
  int right_type;
  int dummy_value; 
  void display(QTransducer &);
  int operator==(const M_id &);
  int operator!=(const M_id &);
  vector<int> toVector(void);
  bool isValid(void); 
};


struct eqM_id{                                                                                                   
  bool operator()(const M_id m1, const M_id m2) const {                                                          
    return (m1.q_state == m2.q_state &&                                                                          
            m1.left_state == m2.left_state &&                                                                    
            m1.left_type == m2.left_type &&                                                                      
            m1.right_state == m2.right_state &&                                                                  
            m1.right_type == m2.right_type);                                                                     
  }                                                                                                          
};

// forward declarations
class Profile; 
typedef vector<M_id> state_path; 
//hash_map<const char*, int, hash<const char*>, eqstr> months;
/* namespace __gnu_cxx { */
/*   template<> struct hash<M_id>{ */
/*     size_t operator()(const M_id& m) const { */
/*       return size_t(m.q_state*1 + m.left_state*10 + m.left_type*100 + m.right_state*1000 + m.right_type*10000);  */
/*     } */
/*   }; */
/* }  // namespace __gnu_cxx */


namespace __gnu_cxx {

  template<> struct hash<M_id>{
    size_t operator()(const M_id& m) const {
      string s; 
      stringstream inter; 
#ifdef DART_DEBUG
      if (m.q_state == -10 || m.left_state == -10 || m.right_state == -10 || m.left_type == -10 || m.right_type == -10 )
	std::cerr<<"\nWarning: using uninitialized M_id in vector conversion!\n";
#endif
      inter << m.q_state << "," << m.left_state << "," << m.right_state << "," << m.left_type << "," << m.right_type; 
      s = inter.str();
      // SDBM hash taken from http://www.partow.net/programming/hashfunctions
      unsigned int hash = 0;
      for(size_t i = 0; i < s.length(); i++)
      {
	hash = s[i] + (hash << 6) + (hash << 16) - hash;
      }
      return size_t(hash);
    }
  };
}  // namespace __gnu_cxx


typedef __gnu_cxx::hash_map<M_id,bfloat, hash<M_id>, eqM_id> DP_hash;




// ***** The 'mature' profile class (e.g. E_n) *****
class AbsorbingTransducer 
{
 public:

  // Constructs an absorbing transducer from a (sampled) profile object, an exact-match transducer, or void input. 
  AbsorbingTransducer(Profile *sampled_profile);
  AbsorbingTransducer(ExactMatch *EM_in);
  AbsorbingTransducer(const char*, vector<string> alphabet, PHYLIP_tree& tree);   
  AbsorbingTransducer(void);   

  Profile* stored_sampled_profile; 
  bool test_equality(AbsorbingTransducer&, bool, bool); 
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
  MyMap<state, M_id> state2mid; 
  // The reverse: 
  MyMap<vector<int>, int> mid2int;
  
  // accessor for incoming states
  vector<state> get_incoming(state);
  
  // query whether or not a given transition exists in the transducer
  bool has_transition(state e, state ePrime);

  //displaying states, transitions, absorption weights
  void display(bool states, bool transitions, bool absorptions);
  void show_DOT(ostream&, string name ="MyGraph");



  // Transition/absorbing accessor functions
  bfloat get_transition_weight(state e, state ePrime);
  bfloat get_absorb_weight(state e, int charIndex);
  
  // state-type-phylogeny. stores binary ((wait, null, end), (match, delete)) state assignments for each state
  MyMap<state, MyMap<node, string> > state_type_phylogeny; 

  // between stores alignment columns in between said external states
  MyMap<pair <state, state>, MyMap< node, string> > between;
  MyMap<pair <state, state>, vector<M_id> > summed_nulls;
  
  // leaf coordinates.  Each state maps to a pair holding the smallest and largest leaf coordinate accounted
  // for by this state.                                                                                 
  MyMap<state, pair<int,int> > leaf_coords;

  // Leaf sequence coordinates - naming is a bit awkward wrt the previous data.  I originally intended the above to be what 
  // leaf_seq_coords is holding, but for a simple alignment envelope only the min/max leaf coords are necessary.  To use a guide
  // alignment we're going to want to store the actual sequence coordinates.  

  // Each state maps to a map holding the leaf coordinates of each sequence in the subtree accounted 
  // for by this state.                                                                                 
  MyMap<state, MyMap<node,int> > leaf_seq_coords;
  

  //testing
  void test_transitions(void);

  // I/O stuff
  // Output - internal
  void write_profile(ostream& out, state_path& viterbi_path);  
  void add_tag_value_pair(ostream& out, string tag, string value, bool newline=true);  
  void add_tag_value_pair(ostream& out, string tag, bfloat value, bool newline=true);  
  void add_tag_value_pair(ostream& out, string tag, int value, bool newline=true);  
  void add_basic_state(ostream& out, string type, int stateIndex);  
  void add_delete_state(ostream& out, int stateIndex);  
  void add_transition(ostream& out, int fromState, int toState, bfloat weight); 

  // Input - internal
  void read_profile(const char* profile_filename, bool logging=false); 
  void parse_state(SExpr* s); 
  void parse_transition(SExpr* s); 
  void verify(bool logging); 

 private:
  // Constructor uses the following private functions to transform profile -> absorbing :
  void identify_states(Profile *sampled_profile);
  void marginalize_null_states(Profile *sampled_profile);
  void DART_marginalize_null_states(Profile *sampled_profile);
  void BFS_marginalize_null_states(Profile *sampled_profile);

  void index_delete_states(Profile *sampled_profile);

  // Constructor uses the following data to transform profile -> absorbing :
  int sampled_start_state; 
  int sampled_pre_end_state;

  // vector<int> sampled_pre_end_states;  

  int sampled_end_state;   
  vector<string> sampled_state_types;
  vector<int> sampled_externals;
  MyMap< vector<int>, vector<M_id> > sampled_externals_outgoing;
  MyMap<pair< vector<state> , vector<state> >, bfloat > transition_weight_tmp; 
  MyMap< vector<int>, vector<M_id> > incoming_tmp; 
  // possible sampled_state types: ("H-types")
  string h_start;
  string h_left_internal;
  string h_right_internal;
  string h_external;
  string h_pre_end;  
  string h_end;
  


  // raw data maps
  MyMap< state, vector<state> > incoming;
  MyMap< state, vector<bfloat> >  absorption_weight; 
  MyMap< pair<state, state>, bfloat > transition_weight;
  
  //various utility functions
  bool checkState(state e);
  bool isDelState(state e);
};



// The states in Profile s are identified by M_id:


int index(M_id query, vector<M_id> in );
bool contains(M_id child, queue<M_id> stateQueue);
bool contains(M_id child, deque<M_id> stateQueue);
bool in(M_id query, vector<M_id> in);
state_path reversed(vector<M_id>);
state_path l2v(list<M_id>); 


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
  Profile(void); 
  // Initiate and fill forward-like DP matrix
  void fill_DP(int, bool inLog=false);
  bfloat post_prob(M_id); 
  MyMap< vector<int>, vector<M_id> > incoming; // this is filled only if the backward algorithm is requested
  MyMap< vector<int>, vector<M_id> > outgoing; // this is filled only if the backward algorithm is requested
  MyMap< pair<vector<int>, vector<int> >, bfloat > transition_weight; // this is filled only if the backward algorithm is requested
  
  bfloat forward_prob; // This is always  stored, for postprob calculations

  // Backward-like DP matrix methods
  void fill_backward_DP(int logging); 

  void clear_DP(void); 
  // sample set of paths from DP matrix, storing info on their transitions, transition weights, etc
  state_path sample_DP(int, int, bool, bool, bool viterbi=false);
  void cache_state(M_id m, M_id mPrime, bfloat weight);
  void cache_path(vector<M_id>); 
  void store_summed_nulls(vector<M_id>); 

  vector<M_id> sampled_states;
  int num_sampled_externals; 
  int max_sampled_externals; 
  MyMap< vector<int>, vector<M_id> > sampled_incoming;
  MyMap< vector<int>, vector<M_id> > sampled_outgoing;  
  MyMap< pair< vector<int>, vector<int> >, bfloat> sampled_transition_weight; 
  bfloat get_sampled_transition_weight(M_id m, M_id mPrime);
  bfloat get_internal_cascade_weight(M_id m);
  bfloat get_external_cascade_weight(M_id m, int charIndex);  
  bfloat compute_emission_weight(M_id m);

  bool is_start(M_id);
  inline bool is_external(M_id);
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
  MyMap< vector<int>, MyMap<node, string> > state_type_phylogeny;
  MyMap<node, string> merge_STP(M_id);
  
  // between stores alignment columns in between said external states
  MyMap< pair< vector<int>, vector<int> >, MyMap<node, string> > between;
  MyMap< pair< vector<int>, vector<int> >, vector<M_id> > summed_nulls;
  
  
  // leaf coordinates.  Each state maps to a pair holding the smallest and largest leaf coordinate accounted 
  // for by this state.
  MyMap< vector<int>, pair<int, int> > leaf_coords;
  // Simple envelope criterion - leaf sequences must be within envelope_distance of each other
  int envelope_distance; 
  
  // Each state maps to a map holding the leaf coordinates of each sequence in the subtree accounted 
  // for by this state.                                                                                 
  MyMap<vector<int>, MyMap<node,int> > leaf_seq_coords;

  // envelope-related variables
  MyMap<node, int> envMap; 
  MyMap<node, int>::iterator envIter1, envIter2; 
  double total;
  int num_total; 
  bool use_guide_alignment; 
  AlignmentEnvelope* envelope; 
  string name1, name2; 
  
  
  // Qtransducer object
  QTransducer Q;
  int get_profile_type(state, string); 

  // Ordering of states
  deque<M_id> backward_states;
  string show_alignment(vector<M_id> &, bool);
  MyMap<string, string> alignment_map(vector<M_id> &, bool);
  // Left and right (mature) profile objects
  AbsorbingTransducer left_profile;
  AbsorbingTransducer right_profile;
  // number of envelope-discards
  int num_discarded_states; 
  int num_zero_states; 
  int DP_size(void); 

  // Displaying the transducer
  void show_DOT(ostream&, string name="MyGraph"); 
  void depth_first_traversal(ostream& out, int sampled_start_state, int sampled_end_state);
  vector<M_id> breadth_first_traversal(M_id beginning_state, bool ordered=true); 
  bool all_parents_visited(M_id, vector<M_id>, vector<M_id>, bool logging=false); 
  string displayStates(vector<M_id>); 
  string getDOTname(M_id); 

 private:
  // **** Fill_DP uses the following: ****
  // the actual DP matrix.  Maps states in profile (as M_id ) to double
  //  MyMap< vector<int>, bfloat> Z; // change from map to hash_map ?
  DP_hash Z; 
  MyMap< vector<int>, bfloat> backward_matrix; 

  //accessor function, returns 0 if no matching entry
  inline bfloat get_DP_cell(M_id); 
  bfloat get_backward_DP_cell(M_id); 

  
  // Add to a DP cell - first checks if the cell exists. 
  inline void add_to_DP_cell(M_id, bfloat); 
  void add_to_backward_DP(M_id, bfloat); 

  // sum over all paths to the M_id cell in the forward matrix
  void sum_paths_to(M_id, bool, bool logging=false);

  //Access functions
  // compute/access the transition weight from state M_id m  to M_id mPrime
  bfloat compute_transition_weight(M_id m , M_id mPrime);
  bfloat get_transition_weight(M_id m , M_id mPrime);
  
  MyMap< pair< vector<int>, vector<int> >, bfloat> transition_weight_test; 
  
  // compute/access the emission cascade weight from state m
  
  //double get_emission_weight(M_id m);
  vector<bfloat> tmpEmitVals; 
  vector<int> tmpEmitTuple; 
  vector<int> tmpMidVec;
  //  MyMap<vector<int>, bfloat>::iterator tmpIter; // change from map to hash_map ?
  DP_hash::iterator tmpIter; 
  
  // Utility functions
  inline bool is_in_envelope(state left_state, state right_state, string checks="NA"); 

  void path2alignment(vector<M_id> &);

  void show_state(M_id);


  // Utility variables
  int profile_start;
  int profile_wait;
  int profile_delete;
  int profile_pre_end;
  int profile_end;

};

void show_state_phylo(MyMap<node, string> &);    
inline MyMap<node, string> cat_STP(MyMap<node, string>, MyMap<node, string>);
MyMap<node, string> pad_STP(MyMap<node, string>, vector<node>);
bool is_flush(MyMap<node, string> map1);
pair<int, int> merge_coords(pair<int, int>, pair<int,int>); 
MyMap<node, int> merge_seq_coords(MyMap<node, int> l , MyMap<node, int> r); 
#endif
