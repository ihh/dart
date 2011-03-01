#include<iostream>
#include<map>
#include<vector>
#include <stack>
#include <list>
#include <queue>
#include <algorithm>
#include <math.h>

#include "protpal/utils.h"
#include "protpal/profile.h"
#include "util/dexception.h"

using namespace std;

// Mini-class M_id methods
int M_id::operator==(const M_id &right)
{
  // Equality operator for M_n states.  Q, left, right states must be equal, as well as the 
  // types of left and right profiles. 
  if (this->q_state == right.q_state && this->left_state == right.left_state && this->right_state == right.right_state && \
      this->left_type == right.left_type && this->right_type == right.right_type) return 1;
  else return 0;
}

int M_id::operator!=(const M_id &right)
{
  if ( *this == right) return 0; 
  else return 1;
}

bool M_id::isValid(void)
{
  bool retVal = true; 
  if ( q_state == dummy_value)
    {
      cerr<<"Warning: uninitialized Q state.\n ";
      retVal = false; 
    }

  if ( left_state == dummy_value)
    {
      cerr<<"Warning: uninitialized l state.\n ";
      retVal = false; 
    }
  if ( right_state == dummy_value)
    {
      cerr<<"Warning: uninitialized r state state.\n ";
      retVal = false; 
    }

  if ( left_type == dummy_value)
    {
      cerr<<"Warning: uninitialized r type.\n ";
      retVal = false; 
    }

  if ( right_type == dummy_value)
    {
      cerr<<"Warning: uninitialized r type.\n ";
      retVal = false; 
    }
  return retVal; 
}

vector<int> M_id::toVector(void)
{
  // An inefficient hack: rather than storing maps from M_id to data (e.g. DP cells),
  // they're first mapped to vectors of integers.
  vector<int> out(5); 
  out[0]= q_state;
  out[1]= left_state;
  out[2]= right_state;
  out[3]= left_type;
  out[4]= right_type;	  
  if (index(dummy_value, out) != -1)
    {
      cerr<<"Warning: Uninitialized M_id used in vector conversion! The vector: "; displayVector(out); 
    }
  return out; 
}

M_id::M_id(void)
{
  dummy_value = -10; 
  q_state = dummy_value; 
  left_state = dummy_value; 
  right_state = dummy_value; 
  left_type = dummy_value; 
  right_type = dummy_value; 
}


void M_id::display(QTransducer &Q)
{
  // Display a state in M_n, showing it's Q state, left/right states and their types.
  MyMap<int, string> int2state;
  int2state[-1]="start";
  int2state[0]="wait";
  int2state[1]="delete";
  int2state[2]="pre-end";
  int2state[3]="end";  

  std::cout<<"QTransducer: \t"<<Q.get_state_name(q_state)<<endl;
  std::cout<<"\tleftProfile:\t type: "<<int2state[left_type]<<", index:"<<left_state<<endl;
  std::cout<<"\trightProfile: \t type: "<<int2state[right_type]<<", index:"<<right_state<<endl;
}


int index(M_id query, vector<M_id> in )
{
  // Find the index of an M_n state in a vector of states
  for (unsigned int i=0; i<in.size(); i++)
    {
      if (in[i] == query) return(i);
    }
  return(-1);
}

bool in(M_id query, vector<M_id> in)
{
  // Query whether or not a vector of states contains a given state.  Not super efficient.  
  if (index(query,in) == -1) return false;
  else return true; 
}

bool contains(M_id child, queue<M_id> stateQueue)
{
  //Query whether or not a queue of states contains a given state.  Not super efficient.  
  while (!stateQueue.empty())
	{
	  if (child == stateQueue.front()) return 1; 
	  else stateQueue.pop();
	}
  return 0;
}
bool contains(M_id child, deque<M_id> stateQueue)
{
  //Query whether or not a queue of states contains a given state.  Not super efficient.  
  while (!stateQueue.empty())
    {
      if (child == stateQueue.front()) return 1; 
      else stateQueue.pop_front();
    }
  return 0;
}
	  
vector<M_id> reversed(vector<M_id> in)
{
  int len = in.size(); 
  vector<M_id> out; 
  for (int i=len-1; i>=0; i--)
    out.push_back(in[i]); 
  return out; 
}
  
state_path l2v(list<M_id> in)
{
  // is this as stupid as it feels? 
  state_path out; 
  for (list<M_id>::iterator iter=in.begin(); iter!=in.end(); iter++)
    out.push_back(*iter); 
  return out; 
}

// ***** Absorbing transducer (E_n) methods  ******

void AbsorbingTransducer::test_transitions(void)
{
  // This function tests for super-normalization of transitions.  
  // Loops through states, summing the outgoing transitions for each one
  MyMap<state, bfloat> outgoing;
  MyMap<state, bfloat>::iterator stateIter;   
  vector<int>::iterator j; 
  int i; 
  for (i=0; i<num_delete_states; i++)
	outgoing[i] = 0; 
  outgoing[start_state] = 0;
  outgoing[pre_end_state] = 0; 
  outgoing[end_state] = 0;   

  for (i=0; i<num_delete_states; i++)
	{
	  for (vector<int>::iterator j=incoming[i].begin();j!=incoming[i].end();j++)
		{
		  //std::cerr<<"Adding "<< get_transition_weight(*j,i)<<" to "<<*j<<endl; 
		  outgoing[*j] += get_transition_weight(*j,i);
		}
	}
  i = pre_end_state; 
  for (j=incoming[i].begin();j!=incoming[i].end();j++)
	{
	  //std::cerr<<"Adding "<< get_transition_weight(*j,i)<<" to "<<*j<<endl; 
	  outgoing[*j] += get_transition_weight(*j,i);
	}
  
  i = end_state; 
  for (j=incoming[i].begin();j!=incoming[i].end();j++)
	{
	  //std::cerr<<"Adding "<< get_transition_weight(*j,i)<<" to "<<*j<<endl; 
	  outgoing[*j] += get_transition_weight(*j,i);
	}
  bool foundSuper = false; 
  for (stateIter = outgoing.begin(); stateIter!=outgoing.end(); stateIter++)
	{
	  if (stateIter->second>1.001)
		{
		  std::cerr<<"Error: State "<<stateIter->first<< " is super-normalized in its transitions\n";
		  std::cerr<<"Value: "<< stateIter->second<<endl; 
		  foundSuper = true;
		}
	}
  if (foundSuper)
	exit(1);
}

void AbsorbingTransducer::display( bool states, bool transitions, bool absorptions)
{
  // Display an absorbing transducer.  These can get pretty large, so there are options as to
  // which aspects to display.  
  vector<state> in;
  vector<state>::iterator inState;
  int charIndex; 
  bfloat absorbWeight, cutoff = 0.001;
 
  if (states) std::cout<<"**** Displaying  states of absorbing transducer ****\n";
  if (transitions) std::cout<<"**** Displaying  transitions of absorbing transducer ****\n";
  if (absorptions) std::cout<<"**** Displaying  absorption weights of absorbing transducer ****\n";    

  if (states) std::cout<<"\nStart state, #"<<start_state<<endl;
  if(transitions) std::cout<<"\t Has no incoming transitions\n";
  if(absorptions) std::cout<<"\t Has no absorption profile\n";  

  for (int i=0; i<num_delete_states; i++)
	{
	  std::cout<<"\n-----------------\nDelete state #"<<i<<endl;
	  if (states) 
		{
		  std::cout<<"State looks like:\n"; show_state_phylo(state_type_phylogeny[i]);
		}
	  if (transitions)
		{
		  in = get_incoming(i);
		  std::cout<<"\t Has the incoming transitions:\n";	  
		  for (inState = in.begin(); inState != in.end(); inState++)
			{
			  std::cout<<"\t\t"<<*inState<<" weight: "; 
			  if (has_transition(*inState, i))
				std::cout<<get_transition_weight(*inState, i)<<endl;			  
			  else std::cout<<"No weight, strange!\n"; 
			}
		}
	  if (absorptions)
		{
		  std::cout<<"\t Has the absorption profile (only non-small weights shown):\n";	  
		  for (charIndex = 0; charIndex < alphabet_size; charIndex++)
			{
			  absorbWeight = get_absorb_weight(i,charIndex); 
				  if (absorbWeight > cutoff) std::cout<<"\t\tCharacter "<<alphabet[charIndex]<<" is absorbed with weight "<<absorbWeight<<endl;
			}
		}
	}

  std::cout<<"\n-----------------\nPre-end state #"<<pre_end_state<<endl;
  if (transitions)
	{
	  in = get_incoming(pre_end_state);
	  std::cout<<"\tHas the incoming transitions:\n";	  
	  for (inState = in.begin(); inState != in.end(); inState++)
		{
		  std::cout<<"\t\t"<<*inState<<" weight: "<<get_transition_weight(*inState, pre_end_state)<<endl;
		}
	}
  if(absorptions)  std::cout<<"\t\tPre-end state has no absorption profile\n";
  
  std::cout<<"\n-----------------\nEnd state #"<<end_state<<", has the incoming transitions:\n";	  
  if (transitions)
	{
	  in = get_incoming(end_state);
	  for (inState = in.begin(); inState != in.end(); inState++)
		{
		  std::cout<<"\t"<<*inState<<" weight: "<<get_transition_weight(*inState, end_state)<<endl;
		}
	}
  if(absorptions)  std::cout<<"\t\tEnd state has no absorption profile\n";
}

AbsorbingTransducer::AbsorbingTransducer(void)
{
  // Placeholder constructor
}
AbsorbingTransducer::AbsorbingTransducer(Profile *sampled_profile)
{
  /* This constructor takes in a profile object that has sampled from its DP matrix,
     and has the (filled) vector sampled_states.  Via a set of transformations detailed
     in Section 3.3 of transducer.tex, we transform this null-in, null-out 
     transducer (e.g. M_n) into a non-null-in, null-out transducer (e.g. E_n).  
     
     In brief, the steps are: 
     - catalog the states as internal/external/null, etc
     - marginalize null states
     - index remaining delete states in topological order
     - calculate the absorption weight for each (delete) state, character pair

     After these steps, this object will be presentable to a Profile object constructor (as one of 
	 the left or right profiles) .  
  */
  bool logging =false;
  
  h_start = "h_start";
  h_left_internal = "left_internal";
  h_right_internal = "right_internal";
  h_external = "external";
  h_end = "h_pre_end";
  h_end = "h_end";

  /* Collect some basic information about the incoming transducer, like:
     - the tree node it is associated with
     - name, alphabet, etc

     These will be stored later, after we've moved to integer state indexing:
     - the characters of its children (MyMap<state, MyMap<node, string> > state_type_phylogeny)
     - min and max leaf_coordinates accounted for (MyMap<state, pair<int, int> > leaf_coords)

  */
  treeNode = sampled_profile->treeNode;
  subtreeNodes = sampled_profile->subtreeNodes;  
  name = sampled_profile->name;
  alphabet = sampled_profile->alphabet; 
  alphabet_size = alphabet.size();  

  /*  
      Go through the sampled states in the profile, storing its H_n type
      (e.g. ext, left_int, right_int,null, etc) as well as identifying the start
      pre_end and end states.  This is meant as a pre-processing step before
      doing the null-elimination.  store the variables num_delete_states, 
      sampled_{start, end, pre_end}_state 
  */
  if(logging) std::cerr<<"\n****** Identifying the states sampled ******\n";
  identify_states(sampled_profile);
  
  /*
    For each non-null state "e" (start state included (?)), traverse the profile DAG
    until another non-null state "ePrime" is found, adding this weight to the transition
    weight e->ePrime.  When this is finished, we will have a new transition matrix 
    that only involves non-null states (a.k.a. delete states).  
  */

  if(logging) std::cerr<<"\n****** Marginalizing null states ******\n";
  marginalize_null_states(sampled_profile);


  /*
    Traverse the delete states in topological order, assigning each an index i from 
    0 <= i < num_delete_states.  Fill the final transition matrix which is a map from
    pairs of state integers to transition weight.  
    NB. this could have been incorporating in the previous function, but I thought that
    that process was sufficiently hairy to separate the two out.  Maybe will integrate
    them later (though probably not).
   */
  if(logging) std::cerr<<"\n ****** Indexing states ******\n\n";
  index_delete_states(sampled_profile);

}

AbsorbingTransducer::AbsorbingTransducer(ExactMatch *EM_in)
{

  /* Construct an absorbing transducer from an instantiated ExactMatch transducer object
     This is fairly trivial, so all the operations (which are split up in the profile->absorbingProfile
	 constructor) are lumped together here.  
  */

  treeNode = EM_in->treeNode;
  subtreeNodes.push_back(EM_in->treeNode); 
  alphabet = EM_in->alphabet; 
  alphabet_size = alphabet.size();
  num_delete_states = EM_in->num_delete_states; 
  
  incoming = EM_in->incoming; 
  start_state = EM_in->start_state;
  // pre_end_states.push_back(EM_in->pre_end_state);
  pre_end_state = EM_in->pre_end_state;
  end_state = EM_in->end_state; 
    
  absorption_weight = EM_in->absorb;
  state_type_phylogeny = EM_in->state_type_phylogeny;
  leaf_coords = EM_in->leaf_coords;   
  // for alignment envelope
  leaf_seq_coords = EM_in->leaf_seq_coords;   

  // Set transition weights. All weight 1; we're assuming that a sequence has a single parse through the leaf/singlet grammar
  pair<state, state> transitionPair; 
  for (int i=0; i<num_delete_states; i++)
	{
	  for (vector<state>::iterator j=incoming[i].begin(); j!=incoming[i].end(); j++)
		{
		  transitionPair.first = *j; transitionPair.second = i; 
		  transition_weight[transitionPair] = 1; 
		}
	}

  transitionPair.first = num_delete_states-1; 
  transitionPair.second = pre_end_state; 
  transition_weight[transitionPair] = 1; 
  
  transitionPair.first = pre_end_state; transitionPair.second = end_state; 
  transition_weight[transitionPair] = 1; 
}


bool Profile::is_start(M_id m)
{
  return (m.q_state == Q.composite_start_state);
}

bool Profile::is_end(M_id m)
{
  return (m.q_state == Q.composite_end_state);
}

bool Profile::is_external(M_id m)
{
  // query whether an M_n state is external - that is, the 'R' transducer is emitting a character 
  // which is absorbed by the branch transducers (and matched or deleted). 
  string q_type;
  //string upsilon_type, bl_type, br_type;

  //string q_type; 
  //  char upsilon_type, bl_type, br_type; 
  int left_type, right_type; 
  q_type = Q.get_state_type(m.q_state);
//   upsilon_type = stringAt(q_type,1);
//   bl_type = stringAt(q_type,2);
//   br_type = stringAt(q_type,3);            
  right_type = m.right_type;
  left_type = m.left_type;

  //  upsilon_type = q_type[1];
  //  bl_type = q_type[2]; 
  //  br_type = q_type[3]; 



  if (q_type == "IMMM" && left_type == profile_delete && right_type == profile_delete)
    return true;
//   if (upsilon_type == "M" && bl_type == "M" && left_type == profile_delete && 
//       br_type == "M" && right_type == profile_delete)
//     return true; //   // MMDMD
      

  if (q_type == "IMMD" && left_type == profile_delete && right_type == profile_wait)
    return true;
//   else if (upsilon_type == "M" && bl_type == "M" && left_type == profile_delete &&
// 	   br_type == "D" && right_type == profile_wait)
//     return true; // ||  // MMDDW

  if (q_type == "IMDM" && left_type == profile_wait && right_type == profile_delete)
    return true;      
//   else if (upsilon_type == "M" && bl_type == "D" && left_type == profile_wait &&
// 	   br_type == "M" && right_type == profile_delete)
//     return true; //|| // MDWMD

  if (q_type == "IMDD" && left_type == profile_wait && right_type == profile_wait)
    return true;      
//   else if (upsilon_type == "M" && bl_type == "D" && left_type == profile_wait &&
//       br_type == "D" && right_type == profile_wait )  // MDWDW
//     return 1; 
  
  else 
    return false; 
}

bool Profile::is_left_int(M_id m)
{
    // query whether an M_n state is left-internal - left branch is making an insertion
  string q_type, upsilon_type, bl_type, br_type;
  int left_type, right_type; 
  q_type = Q.get_state_type(m.q_state);
  upsilon_type = stringAt(q_type,1);
  bl_type = stringAt(q_type,2);
  br_type = stringAt(q_type,3);            
  right_type = m.right_type;
  left_type = m.left_type;

  if (upsilon_type == "S" && bl_type == "I" && left_type == profile_delete && 
	  br_type == "W" && right_type == profile_wait || // SIDWW

	  upsilon_type == "M" && bl_type == "I" && left_type == profile_delete &&
	  br_type == "W" && right_type == profile_wait) // MIDWW
	{
	  return 1; 
	}
  else return 0; 
}

bool Profile::is_pre_end(M_id m)
{
    // query whether an M_n state is the pre-end state
  string q_type, upsilon_type, bl_type, br_type;
  int left_type, right_type; 
  q_type = Q.get_state_type(m.q_state);
  upsilon_type = stringAt(q_type,1);
  bl_type = stringAt(q_type,2);
  br_type = stringAt(q_type,3);            
  right_type = m.right_type;
  left_type = m.left_type;

  if (upsilon_type == "W" && bl_type == "W" && left_type == profile_pre_end && 
	  br_type == "W" && right_type == profile_pre_end) // WWWWW
	{
	  return 1; 
	}
  else return 0; 
}

bool Profile::is_right_int(M_id m)
{
    // query whether an M_n state is right-internal - right branch is making an insertion
  string q_type, upsilon_type, bl_type, br_type;
  int left_type, right_type; 
  q_type = Q.get_state_type(m.q_state);
  upsilon_type = stringAt(q_type,1);
  bl_type = stringAt(q_type,2);
  br_type = stringAt(q_type,3);            
  right_type = m.right_type;
  left_type = m.left_type;
  if ( upsilon_type == "S" && bl_type == "S" && left_type == profile_start && 
	   br_type == "I" && right_type == profile_delete ||  // SSSID
	   
	   upsilon_type == "M" && bl_type == "M" && left_type == profile_delete &&
	   br_type == "I" && right_type == profile_delete ||  // MMDID
	   
	   upsilon_type == "M" && bl_type == "D" && left_type == profile_wait &&
	   br_type == "I" && right_type == profile_delete ||  // MDWID

	   // Put this in, since in practice states like Q = IMMI, will have left-prof in wait state, not del
	   upsilon_type == "M" && bl_type == "M" && left_type == profile_wait &&
	   br_type == "I" && right_type == profile_delete)  // MMWID
	{
	  return 1; 
	}
  else return 0; 
}


void AbsorbingTransducer::identify_states(Profile *sampled_profile)
{
  vector<M_id>::iterator m; 
  // This type-system is stupid, and a consequence of M_id's storing ints rather than strings
  // to denote the component types.  I'll this this soon, I hope. 
  int profile_start = -1;
  //  int profile_wait = 0;
  //  int profile_delete = 1;
  int profile_pre_end = 2;
  int profile_end = 3;  

  //  int left_type, right_type;
  string q_type, upsilon_type, bl_type, br_type;
  int state_counter = -1;
  bool logging = false;
  
  if(logging) std::cerr<<"Number of sampled states: "<<sampled_profile->sampled_states.size()<<endl;
  for(m=sampled_profile->sampled_states.begin(); m!=sampled_profile->sampled_states.end(); m++)
    {

      state_counter++;
	  if(logging) std::cerr<<"\nExamining sampled state # "<<state_counter<<"...";
      // External states
      if(sampled_profile->is_external(*m) )
		{
		  sampled_state_types.push_back(h_external); 
		  sampled_externals.push_back(state_counter);
		  if(logging) 
			{
			  std::cerr<<"categorized as external:\n\t"; m->display(sampled_profile->Q);
			}
		}
      
      //left-internal
      else if (sampled_profile->is_left_int(*m))
		{
		  sampled_state_types.push_back(h_left_internal);
		  if(logging) 
			{
			  std::cerr<<"categorized as left-internal:\n\t"; m->display(sampled_profile->Q);
			}
		}
      
      // right-internal
      else if (sampled_profile->is_right_int(*m)) 
		{
		  sampled_state_types.push_back(h_right_internal);
		  if(logging) 
			{
			  std::cerr<<"categorized as right-internal:\n\t"; m->display(sampled_profile->Q);
			}

		}
	  
      // start
      else if (m->q_state == sampled_profile->Q.composite_start_state &&
			   m->left_type == profile_start &&
			   m->right_type == profile_start)
		{
		  sampled_state_types.push_back(h_start);
		  sampled_externals.push_back(state_counter); // consider this as an external state
		  sampled_start_state = state_counter;
		  if(logging) 
			{
			  std::cerr<<"categorized as START state:\n\t"; m->display(sampled_profile->Q);
			}
		}
      // pre-end wait state
      else if (sampled_profile->Q.get_state_class(m->q_state) == "wait" &&
			   m->left_type == profile_pre_end &&
			   m->right_type == profile_pre_end)
		{
		  sampled_state_types.push_back(h_pre_end);
		  // sampled_pre_end_states.push_back(state_counter);
		  sampled_pre_end_state = state_counter;
		  if(logging) 
			{
			  std::cerr<<"categorized as pre-end state:\n\t"; m->display(sampled_profile->Q);
			}
		}
      
      // end state
      else if (m->q_state == sampled_profile->Q.composite_end_state &&
			   m->left_type == profile_end &&
			   m->right_type == profile_end)
		{
		  sampled_state_types.push_back(h_end);
		  sampled_end_state = state_counter;
		  if(logging) 
			{
			  std::cerr<<"categorized as END state:\n\t"; m->display(sampled_profile->Q);
			}
		}
      else
		{
		  std::cerr<<"Error: encountered unknown h state type.  The individual types were:\n";
		  std::cerr<<"Q: "<<sampled_profile->Q.get_state_name(m->q_state);
		  std::cerr<<"el(idx, type): "<<m->left_state<< " " << m->left_type<<endl;
		  std::cerr<<"er(idx, type): "<<m->right_state<< " " << m->right_type<<endl;		  
		  exit(1);
		}
    }
  num_delete_states = sampled_externals.size()-1; 
}

void AbsorbingTransducer::show_DOT(ostream& out, string name)
{
  out<< "digraph " << name  << "_absorbing {\n"; 
  // show graphiviz -friendly output for the state space of an absorbing transducer.  
  vector<state> in; 
  int to; 
  for (to = 0; to<num_delete_states; to++)
    {
      in = get_incoming(to); 
      for (vector<state>::iterator fromIter=in.begin(); fromIter!=in.end(); fromIter++)
	out << *fromIter << "->" << to << "label=" << "'"<<get_transition_weight(*fromIter, to) << "']\n";
    }

  to = pre_end_state;   
  in = get_incoming(to); 
  for (vector<state>::iterator fromIter=in.begin(); fromIter!=in.end(); fromIter++)
    out << *fromIter << "->" << to << "label=" << "'"<<get_transition_weight(*fromIter, to) << "']\n";

  to = end_state; 
  in = get_incoming(to); 
  for (vector<state>::iterator fromIter=in.begin(); fromIter!=in.end(); fromIter++)
    out << *fromIter << "->" << to << "label=" << "'"<<get_transition_weight(*fromIter, to) << "']\n";
  out<<"}\n";
}

void AbsorbingTransducer::DART_marginalize_null_states(Profile *sampled_profile)
{
  // Use dart's machinery to sum out null states of the sampled profile. 
  bool logging = true, logStates = false;
  M_id e,f; 
  pair<vector<int>, vector<int> > transitionPair; 
  int eIdx, ePrimeIdx, num_states;
  vector<M_id> children;
  vector<M_id>::iterator ePrime;
  vector<int> incoming_states, nulls;

  // the main task here is filling incoming_tmp and sampled_externals_outgoing, and transition_weight_tmp
  // First we do the actual marginalization with DART's machinery...

  // Initialize DART's transition matrix                         
  num_states = sampled_profile->sampled_states.size(); 
  if (logging)
    cerr<<"Transition matrix has size: " << num_states << endl; 
  Transition_probs pre_summed(num_states);

  // Fill the transition matrix as well as the list of null states
  if (logging)
    cerr<<"Begin populating pre_summed\n"; 
  for (eIdx=0; eIdx!=num_states; eIdx++)
    {
      e = sampled_profile->sampled_states[eIdx];
      if (sampled_profile->is_external(e) || sampled_profile->is_start(e) || sampled_profile->is_pre_end(e) )
	{
	  if (logStates)
	    {
	      cerr<<"Non-null state: \n"; 
	      e.display(sampled_profile->Q); 
	    }
	}
      else
	nulls.push_back(eIdx);
      if (sampled_profile->sampled_outgoing.count(e.toVector())>0)
	{
	  if (logStates)
	    cerr<<"State has  children\n"; 
	  children = sampled_profile->sampled_outgoing[e.toVector()];
	  if (logStates)
	    cerr<<"Children pulled out\n"; 
	  for (ePrime = children.begin(); ePrime!=children.end(); ePrime++)
	    {
	      ePrimeIdx = index(*ePrime, sampled_profile->sampled_states); // hmm, this is a bit inefficient. Usually about 1K states
	      if (logStates)
		{
		  cerr<<"Index determined\n"; 
		  cerr<<"e and eprime idx: " << eIdx << " " << ePrimeIdx << endl; 
		}
	      pre_summed.transition(eIdx,ePrimeIdx) = sampled_profile->get_sampled_transition_weight(e,*ePrime);
	      if (logStates)
		cerr<<"Transition weight cataloged\n"; 
	    }
	}
    }
  if (logging)
    cerr<<"Done populating pre_summed\n"; 
  // Now DART does the heavy lifting.  Sum over null states and return a new transition matrix. 
  Transition_methods eliminator;
  Concrete_transition_probs post_summed = eliminator.eliminate(pre_summed, nulls);
  if (logging)
    cerr<<"Done eliminating states\n"; 
  // Now repopulate the data structures relevant for AbsorbingTransducer to continue its business. 
  // clear them first, to be sure.
  incoming_tmp.clear(); sampled_externals_outgoing.clear(); transition_weight_tmp.clear(); 
  for (eIdx=0; eIdx!=num_states; eIdx++)
    for (ePrimeIdx=0; ePrimeIdx!=num_states; ePrimeIdx++)
      if (post_summed.transition(eIdx,ePrimeIdx) > 0.0)
	{
	  e = sampled_profile->sampled_states[eIdx]; 
	  f = sampled_profile->sampled_states[ePrimeIdx]; 
	  
	  transitionPair.first = e.toVector();
	  transitionPair.second = f.toVector();

	  //  fill incoming_tmp and sampled_externals_outgoing, and transition_weight_tmp	  
	  incoming_tmp[f.toVector()].push_back(e); 
	  sampled_externals_outgoing[e.toVector()].push_back(f); 
	  transition_weight_tmp[transitionPair] = post_summed.transition(eIdx,ePrimeIdx); 
	}
}



string Profile::displayStates(vector<M_id> states)
{
  string out; 
  vector<M_id>::iterator e; 
  for (e= states.begin(); e!= states.end(); e++)
    out +=  " " + getDOTname(*e) + " "; 
  out += "\n"; 
  return out; 
}

bool Profile::all_parents_visited(M_id e, vector<M_id> visited, vector<M_id> subgraph, bool logging)
{
  vector<M_id>::iterator parent; 
  // verify that all the nodes upstream of e (in the relevant subgraph) have been visited 
  if (!sampled_incoming.count(e.toVector()))
    {
      if (!is_start(e))
	cerr<<"Warning: in BFS, state has no incoming transitions: " << getDOTname(e) << endl; 
      return true; 
    }

  if (logging)
    {
      cerr<<"Querying state: " << getDOTname(e) <<endl; 
      cerr<<"Subgraph:       " << displayStates(subgraph); 
      cerr<<"Parents:        " << displayStates(sampled_incoming[e.toVector()]); 
    }
  for (parent = sampled_incoming[e.toVector()].begin(); parent != sampled_incoming[e.toVector()].end(); parent++)
    {
      if (index(*parent, subgraph) == -1)
	{
	  if (logging)
	    cerr<< "Parent " << getDOTname(*parent) << " not in subgraph, continue. \n"; 
	  continue;
	}
      if (is_external(*parent) || is_start(*parent) || is_pre_end(*parent))
	{
	  if (logging)
	    cerr<< "Parent " << getDOTname(*parent) << " is external, continue. \n"; 
	  continue;
	}
      if (index(*parent, visited) == -1)
	{
	  if (logging)
	    {
	      cerr<< "Parent " << getDOTname(*parent) << " is not visited, but should be, return false. \n"; 
	      cerr<< "Its index in the subgraph is " << index(*parent, subgraph) << endl; 
	    }
	  return false;
	}
    }
  return true; 
}
vector<M_id> Profile::breadth_first_traversal(M_id beginning_state,bool ordered)
{
  // traverse the state space starting from beginning_state, terminating at external states.  
  // This is a breadth-first-traversal, but also a parents-before-children traversal, so there are some extra
  // checks in the queue loop...

  queue<M_id> stateQueue;
  M_id e, endState; 
  vector<M_id> children,visited,subgraph;
  vector<M_id>::iterator child; 
  int counter=0; 
  bool logging = false; 
  
  if (ordered)
    subgraph = breadth_first_traversal(beginning_state,false); 

  stateQueue.push(beginning_state); 
  if (logging)
    cerr << "Beginning to traverse state space\n"; 

  while (!stateQueue.empty())
    {
      e = stateQueue.front(); 
      stateQueue.pop();
      
      if (index(e,visited) != -1)
	{
	  if (logging)
	    cerr<<"Already visited state " << getDOTname(e) << endl; 
	}
      else
	{
	  if (ordered) 
	    // if we must have nodes in the correct order (parents before children), check this condition:
	    // children must have had all their parents visited before they are visited. 
	    {
	      if (! all_parents_visited(e, visited, subgraph) )
		{
		  if (logging)
		    {
		      cerr<<"Not visited becuase not all its parents have been: " << getDOTname(e) << endl; 
		      cerr<<"Parents: " << displayStates(sampled_incoming[e.toVector()]) << endl; 
		      cerr<<"Subgraph: " << displayStates(subgraph) << endl; 
		      cerr<<"Trace from all_parents_visited:\n";
		      all_parents_visited(e, visited, subgraph,true);
		    }
		  
		  continue;
		}
	    }
	  ++counter; 
	  visited.push_back(e); 
	  if (logging)
	    cerr<<"Currently visiting state: " << getDOTname(e) << endl; 
	  // we only add the children of an external state e if e is the starting point
	  if (e != beginning_state)
	    {
	      if (is_external(e) || is_pre_end(e))
		{
		  if (logging)
		    cerr<<"\tNot adding this state's children to stack, since it is external or pre-end\n";
		  continue; 
		}
	    }
	  
	  if (sampled_outgoing.count(e.toVector()) )
	    {
	      children = sampled_outgoing[e.toVector()];
	      for (child = children.begin(); child != children.end(); child++)
		{
		  if (logging)
		    cerr<<"\tAdded child to queue: "<< getDOTname(*child) << endl; 
		stateQueue.push(*child); 
		}
	    }
	  else
	    cerr<<"Warning: the following state has no outgoing transitions: " << getDOTname(e) << endl; 
	  
	}
    }
  if (!visited.size())
    {
      cerr<<"Warning: BFS of size zero returned.  The state " << getDOTname(beginning_state) << " had outgoing transitions: ";
      for (child = sampled_outgoing[e.toVector()].begin(); child != sampled_outgoing[e.toVector()].end(); child++)
	cerr<<" " << getDOTname(*child)<< " ";
      cerr<<endl; 
      exit(0); 
    }
  return visited; 
}


void AbsorbingTransducer::BFS_marginalize_null_states(Profile *sampled_profile)
{
  // This function has a lot going on.  Lots of comments follow. 
  vector<M_id> children, BFS, parents;
  vector<int>::iterator sourceIdx;
  vector<M_id>::iterator e, parent;
  MyMap<vector<int>, bfloat> weights;
  MyMap<vector<int>, bfloat> outgoingSum;   
  MyMap<vector<int>, int> outgoingCount;   
  M_id source; 
  pair< vector<int> , vector<int> > transitionPair;
  bool logging = false;
  bfloat toAdd, small=1e-500;
  int counter; 
  
  
  // Loop over delete/external states.  We'd like to know the transition weight
  // of 'source'-> ePrime, where ePrime is another non-null state. 
  // the main task here is filling incoming_tmp and sampled_externals_outgoing, and transition_weight_tmp
  for (sourceIdx=sampled_externals.begin();sourceIdx!=sampled_externals.end();sourceIdx++)
    {
      source = sampled_profile->sampled_states[*sourceIdx];
      BFS = sampled_profile->breadth_first_traversal(source);     
      if (logging)
	cerr<<"\nSource state: " << sampled_profile->getDOTname(source) << endl; 
      if (logging)
	{
	  cerr<<"\tLength of BFS list: " << BFS.size(); 
	  for (e= BFS.begin(); e!= BFS.end(); e++)
	    cerr<< " " << sampled_profile->getDOTname(*e) << " "; 
	  cerr<<endl; 
	}
      

      weights.clear(); 
      weights[source.toVector()]=1.0;
      counter = 0; 

      for (e = BFS.begin(); e!= BFS.end(); e++)
	{
	  if (logging)
	    cerr<<"New state from BFS: " << sampled_profile->getDOTname(*e) <<endl; 
	  if (*e == source)
	    continue; 

	  // The possible parents are the null states within the (external-terminating) subgraph rooted at "source"
	  parents.clear();
	  for (parent = sampled_profile->sampled_incoming[e->toVector()].begin(); 
	       parent != sampled_profile->sampled_incoming[e->toVector()].end(); parent++)
	    {
	      if (*parent == source)
		parents.push_back(*parent);
	      else if (sampled_profile->is_external(*parent))
		continue;
	      else if ( index(*parent, BFS) != -1)
		parents.push_back(*parent); 
	    }

	  if (sampled_profile->is_external(*e) || sampled_profile->is_pre_end(*e) )
	    {
	      // We've found the next delete ( or pre_end) state that is connected to 
	      // source
	      if (logging)
		cerr<<"New destination state: " << sampled_profile->getDOTname(*e) << endl; 
	      transitionPair.first = source.toVector();
	      transitionPair.second = e->toVector();
	      // store the connectivity between source and e:
	      sampled_externals_outgoing[source.toVector()].push_back(*e);
	      incoming_tmp[e->toVector()].push_back(source);
	      // Set the transition weight from source -> e  based on parents of e
	      // This trans matrix is temporary, becuase states haven't been indexed
	      // yet, so we're still storing data via M_id vectors.
	      toAdd = 0.0; 
	      for (parent = parents.begin(); parent != parents.end(); parent++)
		{
		  if (! weights.count(parent->toVector()))
		    {
		      cerr<< "Parent state had no weights entry.  State: " << sampled_profile->getDOTname(*e) <<endl; 
		      cerr<< "                                   Parent: " << sampled_profile->getDOTname(*parent) <<endl; 
		      exit(0); 
		    }
		  toAdd += weights.safe_get(parent->toVector()) * sampled_profile->get_sampled_transition_weight(*parent, *e); 
		}
	      
	      if(logging)
		cerr<<"New transition added, value: "<<toAdd<<endl; 

	      if (toAdd>1.001)
		{
		  cerr<<"Error: cached transition weight is greater than 1 after summing over null states: " << toAdd <<endl; 
		  exit(1);
		}
	      
	      if (toAdd < small)
		{
		  cerr<<"Summed transition weight is excessivelysmall.  Previous weight: "<<weights[parent->toVector()]<<endl;
		  cerr<<"Current trans weight: "<< sampled_profile->get_sampled_transition_weight(*parent, *e)<<endl;
		  
		  cerr<<"Parent state: " << sampled_profile->getDOTname(*parent) << endl; 
		  exit(1);
		}
	      if (transition_weight_tmp.count(transitionPair)==1)
		cerr<<"Warning: filling a cell that was already filled.  The state: " << sampled_profile->getDOTname(*e); 
	      else 
		transition_weight_tmp[transitionPair] = toAdd; 
	      
	    }
	  else // The new state is a null state - fill its "weight" cell and continue.  
	    { 
	      if (logging)
		cerr<<"\n\n***Bypassing non-external state: " << sampled_profile->getDOTname(*e) << endl; 

	      toAdd = 0.0; 

	      for (parent = parents.begin(); parent != parents.end(); parent++)
		{
		  if (! weights.count(parent->toVector()))
		    {
		      cerr<< "Parent state had no weights entry.  State: " << sampled_profile->getDOTname(*e) <<endl; 
		      cerr<< "                                   Parent: " << sampled_profile->getDOTname(*parent) <<endl; 
		      exit(0); 
		    }

		toAdd += weights.safe_get(parent->toVector()) * sampled_profile->get_sampled_transition_weight(*parent, *e); 
	      //toAdd *= sampled_profile->get_internal_cascade_weight(*e); 	      
		}

	      if (toAdd < small)
		{
		  cerr<<"Stored transition weight is very small.  Previous weight: "<<weights.safe_get(parent->toVector())<<endl;
		  cerr<<"Current trans weight: "<< sampled_profile->get_sampled_transition_weight(*parent, *e) << endl; 
		  cerr<<"Internal emission weight: "<< sampled_profile->get_internal_cascade_weight(*e) << endl; 
		  exit(1);
		}
	      else if(toAdd >1.001)
		{
		  cerr<<"Error: Stored transition weight is greater than 1:" << toAdd << endl;
		  exit(1);
		}
	      if(logging) 
		cerr<<"Caching  state: " << sampled_profile->getDOTname(*e) << "\n          weight: " << toAdd << endl; 
	      if (weights.count(e->toVector()))
		cerr<<"Warning: attempting to fill weight cell  that's already been filled!\n";
	      weights[e->toVector()] = toAdd;
	      if (weights[e->toVector()] < small)
		cerr<<"A small entry made it into the 'weights' map: " << toAdd << " " << 
		  weights[e->toVector()]<< endl; 
	    }
	}
    }
}

void AbsorbingTransducer::marginalize_null_states(Profile *sampled_profile)
{
  // another option, use DART:
  //  DART_marginalize_null_states(sampled_profile); 
  // return; 

  // yet another option - use a slimmer BFS/DP version
//   ofstream fh;
//   fh.open("saved_profile.dot");
//   sampled_profile->show_DOT(fh); 
//   cerr<<"Profile saved as saved_profile.dot\n";
  
  BFS_marginalize_null_states(sampled_profile); 
  return; 
  
  // This function has a lot going on.  Lots of comments follow. 
  // I've tried to follow the notation in transducer.tex as close as I can
  vector<int>::iterator sourceIdx;
  vector<M_id>::iterator ePrime;
  stack<M_id> stateStack;
  MyMap<vector<int>, bfloat> weights;
  MyMap<vector<int>, bfloat> outgoingSum;   
  MyMap<vector<int>, int> outgoingCount;   
  M_id e, source; 
  vector<M_id> children, BFS;
  pair< vector<int> , vector<int> > transitionPair;
  bool logging = false;
  bfloat toAdd, small=1e-500;
  int counter; 
  
  BFS = sampled_profile->breadth_first_traversal(sampled_profile->sampled_states[sampled_start_state]); 
  cerr<<"Size of breadth-first search: " << BFS.size() << endl; 
  // Loop over delete/external states.  We'd like to know the transition weight
  // of 'source'-> ePrime, where ePrime is another non-null state. 
  // the main task here is filling incoming_tmp and sampled_externals_outgoing, and transition_weight_tmp
  for (sourceIdx=sampled_externals.begin();sourceIdx!=sampled_externals.end();sourceIdx++)
    {
      source = sampled_profile->sampled_states[*sourceIdx];
      if (logging)
	cerr<<"\nSource state: " << sampled_profile->getDOTname(source) << endl; 
      
//       if (logging)
// 	{
// 	  std::cerr<<"\nSource state: \n\t"; source.display(sampled_profile->Q);
// 	}
      if (stateStack.size() != 0)
	cerr<<"Warning: Stack not empty when it ought to be. \n"; 
      stateStack.push(source);
      weights.clear(); 
      weights[source.toVector()]=1.0;
      outgoingSum[source.toVector()] =0; 
      counter = 0; 
      // search depth-first for next delete/external state
      while (!stateStack.empty())
		{
		  counter++; 
		  if (counter >= 10000000)
		    {
		      // Error - we're likely caught in an infinite loop
		      ofstream fh;
		      fh.open("saved_profile.dot");
		      sampled_profile->show_DOT(fh); 
		      cerr<<"Infinite loop in null-elimination detected.  Profile saved as saved_profile.dot\n";

		      fh.close(); 
		      fh.open("saved_DFS.txt"); 
		      sampled_profile->depth_first_traversal(fh, sampled_start_state, sampled_end_state); 
		      cerr<<"Profile DFS pattern  saved as saved_DFS.txt\n";
		      exit(0);
		    }
		  if (logging) 
		    cerr<<"Stack size is currently: " << stateStack.size() << endl; 
		  
		  // Pop off next state, inspect its children
		  e = stateStack.top();
		  stateStack.pop();
		  if (logging)
		    { 
		      std::cerr<<"New state from stack: \n\t"; 
		      e.display(sampled_profile->Q);
		    }
		  if (sampled_profile->sampled_outgoing.count(e.toVector())>0)
			children = sampled_profile->sampled_outgoing[e.toVector()];

		  if(logging) 
		    {
		      cerr<<"This state has "<< children.size()<< " children\n";
		      cerr<<" ****** Begin children *****\n"; 
		      for (ePrime = children.begin(); ePrime!=children.end(); ePrime++)
			ePrime->display(sampled_profile->Q); 
		      cerr<<" ****** End children *****\n"; 
		    }
		  
		  for (ePrime = children.begin(); ePrime!=children.end(); ePrime++)
			{
			  if (sampled_profile->is_external(*ePrime) || sampled_profile->is_pre_end(*ePrime) )
				{
				  if (logging){ std::cerr<<"New destination state: \n\t"; ePrime->display(sampled_profile->Q);}
				  // We've found the next delete ( or pre_end) state that is connected to 
				  // source
				  transitionPair.first = source.toVector();
				  transitionPair.second = ePrime->toVector();
				  // store the connectivity between source and ePrime:
				  sampled_externals_outgoing[source.toVector()].push_back(*ePrime);
				  incoming_tmp[ePrime->toVector()].push_back(source);
				  // add to, or initialize the transition weight from source -> ePrime
				  // this trans matrix is temporary, becuase states haven't been indexed
				  // yet, so we're still storing data via M_id vectors.

				  toAdd = weights.safe_get(e.toVector()) * sampled_profile->get_sampled_transition_weight(e, *ePrime); 
				  if(logging)
					{
					  cerr<<"New transition added, value: "<<toAdd<<endl; 
					}
				  if (toAdd>1.001)
					{
					  cerr<<"Error: cached transition weight is greater than 1 after summing over null states: " << toAdd <<endl; 
					  exit(1);
					}
				  
				  if (toAdd < small)
					{
					  cerr<<"Summed transition weight small.  Previous weight: "<<weights[e.toVector()]<<endl;
					  cerr<<"Current trans weight: "<< sampled_profile->get_sampled_transition_weight(e, *ePrime)<<endl;

					  cerr<<"Parent state: \n\t"; 
					  e.display(sampled_profile->Q); 
					  exit(1);
					}
				  if (transition_weight_tmp.count(transitionPair)==1)
					transition_weight_tmp[transitionPair] += toAdd; 
				  else 
				    transition_weight_tmp[transitionPair] = toAdd; 

				}
			  else // child ePrime is a null state - add to stack
				{ 
				  stateStack.push(*ePrime); 
				  
				  if (logging)
					{ 
					  std::cerr<<"\n\n***Bypassing non-external state (e.g. pushed to stack):*** \n\t"; 
					  ePrime->display(sampled_profile->Q);
					  std::cerr<<"\n\n";
					}
				  
				  toAdd = weights.safe_get(e.toVector()) *
				    sampled_profile->get_sampled_transition_weight(e, *ePrime);
				    //					*sampled_profile->get_internal_cascade_weight(*ePrime);

				  if (toAdd < small)
					{
					  std::cerr<<"Stored transition weight is very small.  Previous weight: "<<weights.safe_get(e.toVector())<<endl;
					  std::cerr<<"Current trans weight: "<< sampled_profile->get_sampled_transition_weight(e, *ePrime) << endl; 
					  std::cerr<<"Internal emission weight: "<< sampled_profile->get_internal_cascade_weight(*ePrime) << endl; 
					  exit(1);
					}
				  else if(toAdd >1.001)
					{
					  std::cerr<<"Error: Stored transition weight is greater than 1.\n";
					  exit(1);
					}
				  if(logging) 
					{
					  std::cerr<<"Caching state's weight: " << toAdd << "\n\t"; 
					  ePrime->display(sampled_profile->Q); 
					}
				  weights[ePrime->toVector()] = toAdd;
				  if (weights[ePrime->toVector()] < small)
				    std::cerr<<"A small entry made it into the 'weights' map: " << toAdd << " " << 
				      weights[ePrime->toVector()]<< endl; 
				}
			}
		}
    }
  // Check for super-normalized transitions after null-state-marginalization.  
  // Things may be sub-normalized, but super shouldn't occur...
  // NB this is now implemented as test_transitions function - I'm leaving this here
  // because it might be good for real-time debugging of the above function 
  /*
  for (MyMap< vector<int>, bfloat>::iterator src=outgoingSum.begin(); src!=outgoingSum.end(); src++)
	{
	  if (src->second > 1.0)
		{
		  std::cerr<<"Warning, outgoing sum greater than 1: "<< src->second<<endl; 
		  exit(1);
		}
	}
  */
}

void AbsorbingTransducer::index_delete_states(Profile *sampled_profile)
{
  // Basically this is just to move from M_id -based bookkeeping to state/integer based
  // Still, it's a pretty long function
  // Objectives:
  //  - index states in topological order
  //  - store transitions/connectivity
  //  - fill leaf_coords and state_type_phylogeny 
  bfloat outSum, small=pow(.1,10000); 
  queue<M_id> stateQueue;
  vector<int>::iterator sourceIdx;
  int charIndex, state_counter = -2; //we want the start state to have index -1.  
  vector<M_id> children;
  vector<M_id>::iterator child,dad;
  vector<int>::iterator i; 
  M_id parent, source;
  pair< vector<int> , vector<int> > transitionPair;
  bool logging = false, allParentsVisited;

  stateQueue.push(sampled_profile->sampled_states[sampled_start_state]);
  start_state = -1;
  // Traverse external states breadth first.  
  while (!stateQueue.empty())
    {
	  // get the next parent
      parent = stateQueue.front();
      stateQueue.pop();

	  // figure out if we ought to continue. 
	  // The cases:
	  // We're at the pre-end state.  We want this to have a special index, so we save it until the end
	  if (parent == sampled_profile->sampled_states[sampled_pre_end_state] ) continue;

	  // if not all of a node's parents have been visited. We want the final order to be parents-before-children.
	  allParentsVisited = true; 
	  for (dad=incoming_tmp[parent.toVector()].begin(); dad!=incoming_tmp[parent.toVector()].end(); dad++)
	   {
		 if (mid2int.count(dad->toVector() )<1) allParentsVisited = false; 
	   }
	  if (!allParentsVisited && logging) 
		{std::cerr<<"State not yet ready- not all parents have been visited\n"; continue; }	  
	  if (!allParentsVisited)
		continue; 
	  if (mid2int.count(parent.toVector()) >0) continue; 

      state_counter++;
      mid2int[parent.toVector()] = state_counter;
      state2mid[state_counter] = parent; 


	  if(logging) { std::cerr<<"\nThe following state was indexed as "<< state_counter <<":\n\t"; parent.display(sampled_profile->Q);}

	  // **** Store these important data containers, for non-(start or pre_end) states. ****
	  // Absorption weight associated with this state index 
	  if (parent != sampled_profile->sampled_states[sampled_start_state] && 
	      parent != sampled_profile->sampled_states[sampled_pre_end_state] )
	    {
	      for (charIndex = 0; charIndex < alphabet_size; charIndex++)
		absorption_weight[state_counter].push_back(sampled_profile->get_external_cascade_weight(parent, charIndex));
		// underflow test
		//absorption_weight[state_counter].push_back(1.0)


		  // The leaf  coordinates accounted for by this state. 
		  leaf_coords[state_counter] = sampled_profile->leaf_coords[parent.toVector()];
		  
		  // for alignment envelope
		  leaf_seq_coords[state_counter] = sampled_profile->leaf_seq_coords[parent.toVector()]; 
		  // The alignment column accounted for by this state.
		  state_type_phylogeny[state_counter] = sampled_profile->state_type_phylogeny[parent.toVector()];
		}

	  // If parent is the start state, we give it a "dummy" STP map
	  if (parent == sampled_profile->sampled_states[sampled_start_state] ) 
	    {
	      state_type_phylogeny[state_counter] = sampled_profile->merge_STP(parent); 	  
	      // for alignment envelope
	      leaf_seq_coords[state_counter] = sampled_profile->leaf_seq_coords[parent.toVector()]; 
	    }

	  if ( sampled_externals_outgoing.count(parent.toVector()) < 1 )
	    {
	      std::cerr<<"Error: in traversal, the following state has no outgoing transitions!\n\t";
	      parent.display(sampled_profile->Q);
	      exit(1);
	    }
	  children = sampled_externals_outgoing[parent.toVector()];
	  for (child = children.begin(); child != children.end(); child++)
	    {
	      transitionPair.first = parent.toVector();
	      transitionPair.second = child->toVector();
	      
	      //		  if ( !contains( *child, stateQueue) && mid2int.count(child->toVector()) == 0) stateQueue.push(*child);
	      stateQueue.push(*child);
	    }
    }    

  // give the pre-end state an index
  state_counter++;
  pre_end_state = state_counter;
  parent = sampled_profile->sampled_states[sampled_pre_end_state];
  mid2int[parent.toVector()] = pre_end_state; 
  state2mid[pre_end_state] = parent; 
  state_type_phylogeny[state_counter] = sampled_profile->merge_STP(parent); 
  // for alignment envelope
  leaf_seq_coords[pre_end_state] = sampled_profile->leaf_seq_coords[parent.toVector()]; 

  if (num_delete_states != pre_end_state)
	{
	  std::cerr<<"Pre-end state set as "<<state_counter<<".  There were "<<num_delete_states<<" external states.\n";
	  std::cerr<<"This indicates a problem with num_delete_states and/or pre_end_state 's values\n";
	}
  if(logging) 
	{ 
	  std::cerr<<"\nIndexing pre_end state...\n";
	  std::cerr<<"The following (pre_end) state was indexed as "<< state_counter <<":\n\t"; parent.display(sampled_profile->Q);}


  // give the end state an index
  state_counter++;
  mid2int[sampled_profile->sampled_states[sampled_end_state].toVector()] = state_counter;
  state2mid[state_counter] = sampled_profile->sampled_states[sampled_end_state];
  end_state = state_counter; 

  // Now that all state indices are known (and we can map M_id -> index), 
  // loop over states again, catalog incoming transitions 
  // and store transition weights. These will go into the 'final' data maps.

  pair<int, int> transitionIntPair;
  pair< vector<int>, vector<int> > transitionVecPair ;  
  for (sourceIdx=sampled_externals.begin(); sourceIdx!=sampled_externals.end(); sourceIdx++)
    {
      source = sampled_profile->sampled_states[*sourceIdx];      
	  // if in(*sourceIdx, sampled_pre_end_states) continue;
	  if (*sourceIdx == sampled_pre_end_state) continue;
	  if (sampled_externals_outgoing.count(source.toVector()) < 1)
		{
		  std::cerr<<"Error: in mid->int loop, the following state has no outgoing transitions!\n\t";
		  source.display(sampled_profile->Q);
		  exit(1);
		}

      children = sampled_externals_outgoing[source.toVector()];
      for (child = children.begin(); child!=children.end(); child++)
		{
		  // if in (*child, sampled_profile->sampled_states[sampled_pre_end_state])
		  if ( *child == sampled_profile->sampled_states[sampled_pre_end_state])
			{
 			  if(logging) std::cerr<<"Pre-end state considered as child to some state\n";
// 			  std::cerr<<"index of parent:"<<mid2int[source.toVector()]<<endl;
// 			  std::cerr<<"index of child:"<<mid2int[child->toVector()]<<endl;			  
// 			  std::cerr<<"The parent:\n\t";
//			  source.display(sampled_profile->Q);
			}
		  // store incoming transition info among integers
		  incoming[mid2int[child->toVector()]].push_back( mid2int[source.toVector()] );
		  //		  std::cerr<<"Incoming to child: "<<mid2int[child->toVector()]<<"\n"; displayVector(incoming[mid2int[child->toVector()]]);

		  // store transition weight among integers.  We already know the transition
		  // weight, but it's stored in terms of a pair of M_id vectors
		  transitionVecPair.first = source.toVector();
		  transitionVecPair.second = child->toVector();	  
		  transitionIntPair.first = mid2int[source.toVector()];
		  transitionIntPair.second = mid2int[child->toVector()];	  
		  if (transition_weight_tmp[transitionVecPair] > 0.0)
		    transition_weight[transitionIntPair] = transition_weight_tmp[transitionVecPair];
		  else
		    std::cerr<<"Adding a zero-weight transition to new profile's transition matrix...\n";
		    // underflow test
		  //max(transition_weight_tmp[transitionVecPair], bfloat(1));
		}
    }
  // pre-end state
  source = sampled_profile->sampled_states[sampled_pre_end_state];      

  if(logging) std::cerr<<"\nIndexing end state...\n";
  M_id end = sampled_profile->sampled_states[sampled_end_state];
  if(logging) 
	{
	  std::cerr<<"The following (end) state was indexed as "<<mid2int[end.toVector()]<<":\n\t";
	  end.display(sampled_profile->Q);
	}
  
  incoming[mid2int[end.toVector()]].push_back( mid2int[source.toVector()] );

  transitionIntPair.first = mid2int[source.toVector()];
  transitionIntPair.second = mid2int[end.toVector()];	  
  transition_weight[transitionIntPair] = 1;  

  // transfer the 'between' map from sampled_profile.    
  MyMap<pair<vector<int>, vector<int> >, MyMap<node, string> >::iterator nodeState; 
  for (nodeState=sampled_profile->between.begin(); nodeState!=sampled_profile->between.end(); nodeState++)
	{
	  transitionIntPair.first = mid2int[nodeState->first.first];
	  transitionIntPair.second = mid2int[nodeState->first.second]; 	  
	  between[transitionIntPair] = nodeState->second; 
	}
}


bfloat AbsorbingTransducer::get_absorb_weight(state e, int charIndex)
{
  // Accessor function for absorption weight.  
  // Safety checks that the state exists and is an absorbing state, etc
  #ifdef DART_DEBUG
  if (! isDelState(e))
	{
	  std::cerr<<"Error: the state "  << e << " is not a delete state\n";
	  std::cerr<<"The offending call was: get_absorb_weight in transducer: "<<name<<" state: "<<e<<" character index: "<<charIndex<<endl;
	  exit(1); 
	}
  else if ( charIndex >=0 && charIndex < alphabet_size)
	{
	  if (absorption_weight.count(e)<1)
		{
		  std::cerr<<"Error: the state "<<e<<" appears OK, but have no entry in the emission matrix.\n";
		  exit(1);
		}
	  else if (charIndex >= absorption_weight[e].size())
		{
		  std::cerr<<" The character index "<< charIndex<<" appears too big for the state "<<e<<" which has emission vector of size: "<< absorption_weight[e].size()<<endl;
		  exit(1);
		}
	}
  else
	{
	  std::cerr<<"The alphabet index "<<charIndex<<" is not valid\n"; 
	  std::cerr<<"The offending call was: get_absorb_weight in transducer: "<<name<<" state: "<<e<<" character index: "<<charIndex<<endl;
	  exit(1);
	}
  #endif
  bfloat toRet = absorption_weight[e][charIndex];
  bool isReal=false;
  if (toRet >0.0)
    isReal = true;
  if (isReal)
    return toRet; 
  else
    {
      std::cerr<<"Bad absorption value\n";
      exit(0); 
    }

}

bfloat AbsorbingTransducer::get_transition_weight(state e, state ePrime)
{
  // Transition accessor - checks validity of transition between the queried states
  if (! checkState(e))
	{
	  std::cerr<<"Error: the state "  << e << " is not valid\n";
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<e<<endl;
	  std::cerr<<"There are "<<num_delete_states<<" delete states in this transducer\n"; 
	  std::cerr<<"Start, pre-end, end states are: "<<start_state<<","<<pre_end_state<<","<<end_state<<endl;
	  exit(1);
	}
  if (! checkState(ePrime))
	{
	  std::cerr<<"Error: the state "  << ePrime << " is not valid\n";
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" destination state: "<<ePrime<<endl;
	  std::cerr<<"There are "<<num_delete_states<<" delete states in this transducer\n"; 
	  std::cerr<<"Start, pre-end, end states are: "<<start_state<<","<<pre_end_state<<","<<end_state<<endl;
	  exit(1);
	}
  if (ePrime == start_state) 
	{
	  std::cerr<<"Error: calling incoming transition to start state\n"; 
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<e<<endl;
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" destin state: "<<ePrime<<endl;
	  exit(1);
	}
  else if(index(e, incoming[ePrime]) == -1 )
	{
	  std::cerr<<"Error: state "<<ePrime<<" does not have "<<e<< " as one of its incoming state transitions\n";
	  std::cerr<<"States incoming to "<<ePrime<<" are: "; displayVector(incoming[ePrime]);
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<ePrime<<endl;
	  exit(1);
	}
  else 
	{
	  pair<state, state> transitionPair; transitionPair.first=e, transitionPair.second=ePrime;
	  if (transition_weight.count(transitionPair)) 
	    {
	      bool isReal = false;
	      if (transition_weight[transitionPair] > 0.0)
		isReal = true; 
	      if (!isReal)
		{
		  std::cerr<< "Error: the transition between states "<<e<<" , "<<ePrime<<" is valid, but it is zero\n"; 
		  exit(1);
		}
	      else
		return transition_weight[transitionPair];
	    }
	  else
	    {
	      std::cerr<<"Error: the transition between states "<<e<<" , "<<ePrime<<" ought to be valid, but there is no entry in the transition matrix for it!\n";
	      exit(1);
	    }
	  
	}
}

vector<state> AbsorbingTransducer::get_incoming(state e)
{
  // Access the incoming states of a given state - checks that the state is valid and not 
  // the start state.
  if (e == start_state) 
	{
	  std::cerr<<"Error: in get_incoming from transducer "<<name<<" , start state has no incoming transitions\n";
	  exit(1);
	}

  else if (incoming.count(e) >0) return incoming[e]; 
  else
	{
	  std::cerr<<"State "<<e<< " appears not to have any incoming transitions, though it is not the start state.  Very peculiar. \n";
	  exit(1);
	}
}

// Some boolean 'checking/testing' functions
bool AbsorbingTransducer::has_transition(state e, state ePrime)
{
  //Query whether a transducer has a given transition between states
  if (ePrime == start_state) 
	{
	  return 0; 
	}
  
  if (index(e, incoming[ePrime]) == -1) return 0;
  else return 1;
}

bool AbsorbingTransducer::checkState(state e)
 {
   // query whether a state is valid or not
   if (e == start_state || e == end_state || e == pre_end_state) return 1; 
   else if(e >=0 && e<num_delete_states) return 1; 
   else return 0; 
 }

bool AbsorbingTransducer::isDelState(state e)
{
  // query whether or not a state is a delete/absorbing state.  This is easy to do given our 
  // indexing convention (delete states are numbered 0 to num_delete_states). 
  return (e >=0 && e<num_delete_states);
}



// Profile transducer methods 
Profile::Profile(void)
{
  //placeholder
}

Profile::Profile(node node_in, AbsorbingTransducer left_in, AbsorbingTransducer right_in, QTransducer Q_in)
{
  //constructor - input left, right absorbers and a Q object on top, as well as a phylogenetic node identifier
  Q = Q_in; 
  left_profile = left_in;
  right_profile = right_in; 
  if (! left_profile.alphabet_size == left_profile.alphabet_size == Q.alphabet_size) 
	{
	  std::cerr<<"Alphabet sizes do not match up.  Check construction of Q, left_profile, right profile\n";
	}
  else
	{
	  alphabet = Q.alphabet; 
	  alphabet_size = Q.alphabet_size; 
	  for(int i=0; i<alphabet_size; i++)
	    alphabet_ints.push_back(i); 
	}
  treeNode = node_in;   
  vector<node>::iterator n; 
  
  // The subtree that the profile accounts for is the union of the profiles' subtreeNodes and the 'top' node
  for (n=left_in.subtreeNodes.begin(); n!=left_in.subtreeNodes.end(); n++)
	{ subtreeNodes.push_back(*n); }
  for (n=right_in.subtreeNodes.begin(); n!=right_in.subtreeNodes.end(); n++)
	{ subtreeNodes.push_back(*n); }  
  subtreeNodes.push_back(treeNode); 

  //  envelope_distance = 200; //hard-coded envelope distance
  num_sampled_externals = 0; 
  
  // This is a bit awkward...owing to the quirk that M_ids store types as integers
  // rather than strings.  If I store all types now as variables it gets a bit less confusing.
  profile_start = -1;
  profile_wait = 0;
  profile_delete = 1;
  profile_pre_end = 2;
  profile_end = 3;
  
  num_discarded_states = 0; 

}
// Top-level public methods
state_path Profile::sample_DP(int num_paths, int logging, bool showAlignments, bool leaves_only, bool viterbi_in)
{
  // sample num_paths alignments from the DP matrix, showing the alignments according to 
  // showAlignments.  If showing alignments, they are printed in stockholm format, with
  // their bitscore and posterior probability as #=GF annotations.  
  
  bool testing = false, viterbi=false, nextBreak=false;   //  bool fromStart;
  vector<M_id> states_in; // states for sampling
  vector<bfloat> weights; // weights for each of these states
  bfloat weight, pathWeight, emissionWeight, small=0.00001, left_transition_weight; 
  vector<M_id> pi, pi_previous; // sampled state path
  string S = "S";
  string I = "I";   
  int sampleIdx, loopIdx;
  vector<state> left_incoming, right_incoming, q_incoming;
  vector<state>::iterator e_l, e_r, q; 
  pair<vector<int>, vector<int> > transitionPair; 
  // these variables are the 'from' and 'to' M_n states, and the M_n start state, which tells us when we're finished
  M_id m, mPrime, bigStart;
  string alignString; 

  // Set start state here
  bigStart.q_state = Q.composite_start_state;
  bigStart.left_state = left_profile.start_state;
  bigStart.right_state = right_profile.start_state;
  bigStart.right_type = -1;
  bigStart.left_type = -1;
  if (num_paths == 1 && viterbi==false && treeNode != 0) // assume root=0
    std::cerr<<"Warning: sampling only one path with viterbi option off.  Consider more paths or turning on viterbi..\n"; 
  for (int pathIdx = 0; pathIdx < num_paths; pathIdx++)
	{
	  if (logging >=1)
	    if (pathIdx % 1000 == 0 and pathIdx!=0)
	      cerr<<" " << pathIdx << " "; 
	  
	  if (logging >=1 and nextBreak )
	    {
	      std::cerr<<"(sampling truncated after caching " << num_sampled_externals << " states via ";
	      std::cerr<< pathIdx << " paths)...";
	      break; 
	    }
	  
	  if (pathIdx == num_paths-1)
	    viterbi = viterbi_in; 
	  else
	    viterbi =false; 
	  if (num_sampled_externals > max_sampled_externals)
	    {
	      nextBreak = true; 
	      viterbi = true; 
	    }

	  pathWeight = 1.0; 
	  // initialize pi with the M_n end state:
	  mPrime.q_state = Q.composite_end_state;
	  mPrime.left_state = left_profile.end_state;
	  mPrime.right_state = right_profile.end_state;
	  mPrime.right_type = 3;
	  mPrime.left_type = 3;
  
	  pi.clear(); 
	  pi.push_back(mPrime);
	  if (index(mPrime, sampled_states) == -1) 
		{
		  sampled_states.push_back(mPrime); 
		}
	  // Keep sampling states until we reach M_n's start state
	  // Keep a counter, in case things get out of control. 
	  loopIdx = 0;

	  while( mPrime != bigStart ) 
		{
		  try
		    {
		  loopIdx++; 
		  if (loopIdx > 100000) {std::cerr<<"Possibly caught in  while loop. Exiting...\n";exit(1);}
		  if(logging>=2) 
			{ std::cerr<<"\nCurrent state:\n\t"; mPrime.display(Q); }

		  states_in.clear(); weights.clear();
	
		  // if we're at the END state, possibly move to the pre-end state
		  if (mPrime.q_state == Q.composite_end_state)
			{
			  q_incoming = Q.get_wait_states();
			  for (q=q_incoming.begin(); q!=q_incoming.end(); q++)
				{
				  // for left_pre_end in left_profile.pre_end_states
				  // for right_pre_end in right_profile.pre_end_states				  
				  if (!Q.has_transition(*q, Q.composite_end_state)) 
				    continue;
				  m.q_state = *q;
				  m.left_state = left_profile.pre_end_state;
				  m.right_state = right_profile.pre_end_state; 			  
				  m.left_type = 2; 
				  m.right_type = 2;
				  weight = get_DP_cell(m)*Q.get_transition_weight(*q, mPrime.q_state) *
					left_profile.get_transition_weight(left_profile.pre_end_state,left_profile.end_state) *
					right_profile.get_transition_weight(right_profile.pre_end_state,right_profile.end_state); 
				  
				  // for test
				  if (weight > 0.0 && testing)
					{
					  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
					  if (transition_weight_test.count(transitionPair)<1)
					    {std::cerr<<"Using a novel transition in traceback!\n"; 
					      std::cerr<<"Source:"; m.display(Q);
					      std::cerr<<"Dest:"; mPrime.display(Q); 					  
					      exit(1);}
					  
					  else if( absoluted(transition_weight_test[transitionPair] + small - (weight/get_DP_cell(m)+small))>small )
					    {std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
					      std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
					      exit(1);}
					} //end testing loop
				  if (weight > 0.0)
					{
					  if (m.isValid())
					    {
					      states_in.push_back(m);
					      weights.push_back( weight );
					    }
					}
				}
			  if (weights.size() == 0) 
				{
				  std::cerr<<"Error: end state determined to have no incoming states! the state was:\n";
				  std::cerr<<"Pre-end forward value: "<<get_DP_cell(m)<<endl;
				  mPrime.display(Q);
				}
			  
			}

		  else
			{
			  // we could possibly move to the start state here. 
			  if(left_profile.has_transition(left_profile.start_state, mPrime.left_state) || 
				 mPrime.left_state == left_profile.start_state)
				{
				  if (Q.has_transition(Q.composite_start_state, mPrime.q_state)  && \
					  right_profile.has_transition(right_profile.start_state, mPrime.right_state))
					{
					  m.q_state = Q.composite_start_state; 
					  m.left_state = left_profile.start_state;
					  m.right_state = right_profile.start_state;
					  m.right_type = -1;// start type
					  m.left_type = -1; // start type
					  // NB DP cell of start = 1 .  Left in here for clarity
					  weight = get_DP_cell(m) *							\
						Q.get_transition_weight(Q.composite_start_state, mPrime.q_state)* 						
						right_profile.get_transition_weight(right_profile.start_state, mPrime.right_state);

					  if (mPrime.left_state != left_profile.start_state)
						{ weight*=left_profile.get_transition_weight(left_profile.start_state, mPrime.left_state);}
					  
				  // for test
				  if (weight > 0.0 && testing)
				    {
				      transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
				      if (transition_weight_test.count(transitionPair)<1)
					{
					  std::cerr<<"Using a novel transition in traceback!\n"; 
					  std::cerr<<"Source:"; m.display(Q);
					  std::cerr<<"Dest:"; mPrime.display(Q); 					  
					  exit(1);}
					  
					  else if( absoluted(transition_weight_test[transitionPair] - weight/get_DP_cell(m))>small )
						{std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
						  std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
						  std::cerr<<"Source:"; m.display(Q);
						  std::cerr<<"Dest:"; mPrime.display(Q); 					  
						  exit(1);}
				    }//end test
				  
				  if (weight > 0.0)
				    {
					  if (m.isValid())
					    {
					      states_in.push_back(m);
					      weights.push_back( weight );
					    }
				    }
					}
				}
			  // incoming states where Q is match
			  // both left and right profiles must not be in their start state
			  if (mPrime.right_state != right_profile.start_state && mPrime.left_state != left_profile.start_state)
				{
				  q_incoming = Q.get_incoming_match_states(mPrime.q_state);
				  left_incoming = left_profile.get_incoming(mPrime.left_state);
				  right_incoming = right_profile.get_incoming(mPrime.right_state);  
			  
				  for (e_l = left_incoming.begin(); e_l != left_incoming.end(); e_l++)
					{
					  m.left_state = *e_l;
					  if (m.left_state == left_profile.start_state) 
					    m.left_type = -1;
					  else 
					    m.left_type = 1;

					  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
						{
						  //						  if (! is_in_envelope(*e_l, *e_r))
						  //							continue;
						  m.right_state = *e_r;
						  if (m.right_state == right_profile.start_state) 
						    m.right_type = -1;
						  else 
						    m.right_type = 1;

						  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
							{

							  m.q_state = *q;
							  weight = get_DP_cell(m); 
							  if (weight>0.0)
							    weight *=
							      Q.get_transition_weight(*q, mPrime.q_state) * 
							      left_profile.get_transition_weight(*e_l, mPrime.left_state) *
							      right_profile.get_transition_weight(*e_r, mPrime.right_state) ;
							  else
							    continue;


							  // for test
							  if (weight > 0.0 && testing)
								{
								  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
								  if (transition_weight_test.count(transitionPair)<1)
								    {std::cerr<<"Using a novel transition in traceback!\n"; 
								      std::cerr<<"Source:"; m.display(Q);
								      std::cerr<<"Dest:"; mPrime.display(Q); 					  
								      exit(1);
								    }
							  
							  else if( absoluted(transition_weight_test[transitionPair] - weight/get_DP_cell(m))>small )
								{std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
								  std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
								  std::cerr<<"Source:"; m.display(Q);
								  std::cerr<<"Dest:"; mPrime.display(Q); 			
								  std::cerr<<"Q transition: "<<Q.get_transition_weight(*q, mPrime.q_state)<<endl; 
								  std::cerr<<"left transition: "<<left_profile.get_transition_weight(*e_l, mPrime.left_state)<<endl; 
								  std::cerr<<"right transition: "<<right_profile.get_transition_weight(*e_r, mPrime.right_state)<<endl; 
								  exit(1);}
								} // end test

							  if (weight > 0.0)
							    {
							      if (m.isValid())
								{
								  states_in.push_back(m);
								  weights.push_back( weight );
								}
							    }
							}
						}
					}
				}
			  // incoming states where Q is left-emit. Left type is 1, right type 0. 
			  // left state must not be in start state
			  if (mPrime.left_state != left_profile.start_state)				
				{
				  q_incoming = Q.get_incoming_left_emit_states(mPrime.q_state);
				  left_incoming = left_profile.get_incoming(mPrime.left_state);
  
				  for (e_l = left_incoming.begin(); e_l != left_incoming.end(); e_l++)
					{
					  //					  if (! is_in_envelope(*e_l, mPrime.right_state))
					  //						continue; 

					  m.left_state = *e_l;
					  if (m.left_state == left_profile.start_state) 
					    m.left_type = -1;
					  else 
					    m.left_type = 1;
					  
					  m.right_state = mPrime.right_state;
					  m.right_type = 0; 
					  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
						{
						  m.q_state = *q;
						  weight = get_DP_cell(m);
						  if (weight >0.0)
						    weight *=
						      Q.get_transition_weight(*q, mPrime.q_state)*
						      left_profile.get_transition_weight(*e_l, mPrime.left_state);
						  else
						    continue;
						  // right profile transition weight is implicitely 1 here 

						  // for test
						  if (weight > 0.0 && testing)
							{
							  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
							  if (transition_weight_test.count(transitionPair)<1)
							    {std::cerr<<"Using a novel transition in traceback!\n"; 
							      std::cerr<<"Source:"; m.display(Q);
							      std::cerr<<"Dest:"; mPrime.display(Q); 					  
							      exit(1);
							    }
							  
							  else if( absoluted(transition_weight_test[transitionPair] - weight/get_DP_cell(m)) > small )
							    {
							      std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
							      std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
							      std::cerr<<"Source:"; m.display(Q);
							      std::cerr<<"Dest:"; mPrime.display(Q); 					  
							      exit(1);}
							} // end test
						  if (weight > 0.0)
						    {
						      if (m.isValid())
							{
							  states_in.push_back(m);
							  weights.push_back( weight );
							}
						    }
						}
					}
				}
			  
			  // incoming states where Q is right-emit
			  // right profile must not be in start state, left can be either wait or 'start' (type AND index)
			  if (mPrime.right_state != right_profile.start_state) 				
				{
				  q_incoming = Q.get_incoming_right_emit_states(mPrime.q_state);
				  right_incoming = right_profile.get_incoming(mPrime.right_state);
  
				  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
					{
					  //					  if (! is_in_envelope(mPrime.left_state , *e_r)) 
					  //						continue;

					  m.right_state = *e_r;
					  if (m.right_state == right_profile.start_state) 
					    continue; //m.right_type = -1;
					  else
					    m.right_type = 1;
					  
   					  

					  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
						{
						  m.q_state = *q;
						  m.left_type = get_profile_type(*q, "left");
						  // begin careful stuff
						  if (m.left_type == profile_start)
						    m.left_state = left_profile.start_state;
						  else
						    m.left_state = mPrime.left_state; // note this
						  
						  if (m.left_state != mPrime.left_state)
						    {
						      if (left_profile.has_transition(m.left_state, mPrime.left_state))
							left_transition_weight  = 
							  left_profile.get_transition_weight(m.left_state, mPrime.left_state);
						      else
							continue;
						    }
						  else
						    left_transition_weight  = 1.0;
						  // end careful stuff



// 						  if ( stringAt(Q.get_state_type(*q),0) == S &&
// 							   stringAt(Q.get_state_type(*q),1) == S &&
// 							   stringAt(Q.get_state_type(*q),2) == S &&
// 							   stringAt(Q.get_state_type(*q),3) == I) 
// 							{
// 							  if( mPrime.left_state != left_profile.start_state ) continue;
// 							  else 
// 								{
// 								  if( m.left_type != -1){std::cerr<<"Left profile incorrectly typed!\n"; exit(1);}
// 								}
// 							}// end if SSSI

						  weight = get_DP_cell(m);
						  if (weight>0.0)
						    weight *=
						      Q.get_transition_weight(*q, mPrime.q_state)*
						      right_profile.get_transition_weight(*e_r, mPrime.right_state) * 
						      left_transition_weight; 
						  else
						    continue;

						  // for test
						  if (weight > 0.0 && testing)
							{
							  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
							  if (transition_weight_test.count(transitionPair)<1)
							    {std::cerr<<"Using a novel transition in traceback!\n"; 
							      std::cerr<<"Source:"; m.display(Q);
							      std::cerr<<"Dest:"; mPrime.display(Q); 					  
							      exit(1);
							    }
							  
							  else if( absoluted(transition_weight_test[transitionPair] + small - (weight/get_DP_cell(m)+small))>small )
							    {
							      std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
							      std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
							      std::cerr<<"Source:"; m.display(Q);
							      std::cerr<<"Dest:"; mPrime.display(Q); 					  
							      exit(1);
							    }
							} //end test
						  if (weight > 0.0)
							{
							  if (m.isValid())
							    {
							      states_in.push_back(m);
							      weights.push_back( weight );
							    }
							}
						}
					}
				}
			}
		    }
		  catch(char * str)
		    {
		      cerr<<"Error in computing weights\n"; 
		      exit(1); 
		    }
		  try
		    {
		      if (weights.size() > 1)
			{
			  if (logging>=2) 
			    {
			      /*
				std::cerr<<"\nWeights for incoming states: "; 
				displayVector(weights);
				std::cerr<<"\nThese states were:\n";
				for (vector<M_id>::iterator in=states_in.begin(); in!=states_in.end(); in++)
				{
				MyMap<node, string> tmpMap = merge_STP(*in); 
				show_state_phylo(tmpMap);
				std::cerr<<endl;
				} */
			    }

			  if (viterbi) 
			    sampleIdx = maxIndex(weights);
			  else 
			    sampleIdx = sample(weights);
			  if (sampleIdx >= states_in.size())
			    cerr<<"ERROR: Sampled a state beyond end of vector of possible states!\n";
			  //pathWeight *= (weights[sampleIdx] / get_DP_cell(states_in[sampleIdx]))*
			  //				compute_emission_weight(states_in[sampleIdx]);

			  if (weights.size() != states_in.size())
			    {
			      cerr<<"Error: candidate states and their relative weights are not size-matched!\n"; 
			      exit(1); 
			    }
			  if (logging>=2) std::cerr<<"Sampled state number: "<<sampleIdx<<endl;
			  if (! states_in.at(sampleIdx).isValid())
			    cerr<<"WARNING: the sampled state is not 'valid':\n" << getDOTname(states_in[sampleIdx]) << endl; 
			  pi.push_back(states_in.at(sampleIdx));
			  cache_state(states_in.at(sampleIdx), mPrime, weights.at(sampleIdx) );
			  state_type_phylogeny[states_in.at(sampleIdx).toVector()] = merge_STP(states_in.at(sampleIdx));
			  if(logging>=2)
				{
				  std::cerr<<"Next state:"<<endl;
				  states_in[sampleIdx].display(Q);
				  show_state(states_in[sampleIdx]);
				}
		    }
		  else if (weights.size() == 1)
			{
			  if(logging>=2)std::cerr<<"(Deterministic sample: only one incoming state)\n";
			  if (!states_in[0].isValid())
			    cerr<<"WARNING: the sampled state is not 'valid':\n" << getDOTname(states_in[0]) << endl; 
			  pi.push_back(states_in[0]);
			  cache_state(states_in[0], mPrime, weights[0] );
			  state_type_phylogeny[states_in[0].toVector()] = merge_STP(states_in[0]);
			  //pathWeight *= (weights[0] / get_DP_cell(states_in[0]))*
			  //				compute_emission_weight(states_in[0]);

			  if(logging>=2)
				{
				  std::cerr<<"Next state:"<<endl;
				  states_in[0].display(Q);
				  show_state(states_in[0]);
				  std::cerr<<"Should be same as:\n";
				  show_state_phylo(state_type_phylogeny[states_in[0].toVector()]);				  
				}


			}
		  else //weights has no elements!
			std::cerr<<"Error: the current state has no incoming transitions.\n";
		  if (pi.back() == mPrime) 
			{
			  std::cerr<<"Error: Sampled state is same as current state.\nCurrent state: ";
			  mPrime.display(Q); 
			  exit(1);
			}
		  mPrime = pi.back();
		    }
		  catch(char * str)
		    {
		      cerr<<"Error in caching state\n";
		    }
		
		  if(logging>=2) 
		    { std::cerr<<"\nFinal state:\n\t"; mPrime.display(Q); }

		  // Ideally this would be handled with a nice compression/fast representation of the sampled path
		  // This accounts just for the simple cases of sampling the same path over and over again. 
		  bool samePath = true; 
		  if (pi_previous.size() == pi.size())
		    {
		      for (int pathCount=0; pathCount < pi.size(); pathCount++)
			{
			  if (pi_previous[pathCount] != pi[pathCount])
			    {
			      samePath =false;
			      break;
			    }
			}
		    }
		  else
		    samePath = false;
		  if (!samePath)
		    {
		      // Store the null states' accounted characters:
		      cache_path(pi); 
		      pi_previous = pi; 
		    }
		  // this operation is now handled within cache_path:
		  //  store_summed_nulls(pi); 
		  
		  if (showAlignments) 
		    {
		      std::cout<<"#=GF bit_score "<< -log(pathWeight)/log(2)<<endl; 
		      std::cout<<  "#=GF post_prob "<< pathWeight/forward_prob<<endl; 
		      std::cout<< show_alignment(pi, leaves_only); 
		      std::cerr<<"\n";
		    }
		  
		}
	}
  return pi; 
}
void Profile::cache_state(M_id m, M_id mPrime, bfloat weight)
{
  // upon sampling a move from m to mPrime, store the info about the state
  // m and the transition m -> mPrime

  pair< vector<int>, vector<int> > transitionPair;
  transitionPair.first=m.toVector();   transitionPair.second=mPrime.toVector();
  if (sampled_transition_weight.count(transitionPair))
    return; 
  // m is the new state
  if (index(m, sampled_states) == -1 ) 
	{
	  sampled_states.push_back(m);
	  if ( is_external(m) )
		  num_sampled_externals += 1; 
	}
  if (sampled_outgoing.count(m.toVector()))
    {
      if (index(mPrime, sampled_outgoing[m.toVector()]) == -1) 
	sampled_outgoing[m.toVector()].push_back(mPrime);
    }
  else
    sampled_outgoing[m.toVector()].push_back(mPrime);
      
  if (sampled_incoming.count(mPrime.toVector()))
    {
      if (index(m, sampled_incoming[mPrime.toVector()]) == -1)
	sampled_incoming[mPrime.toVector()].push_back(m); 
    }
  else
    sampled_incoming[mPrime.toVector()].push_back(m); 

  // for alignment envelope - merge a pair of leaf sequence coords, similar to above. 
  if (is_external(m) )
    {
      string side; 
      if (index(m.q_state, Q.get_match_states()) != -1)
	side = "both"; 
      else if (index(m.q_state, Q.get_left_emit_states()) != -1)
	side = "left";
      else if (index(m.q_state, Q.get_right_emit_states()) != -1)
	side = "right";
      else
	cerr<<"Unknown external q type: " << Q.get_state_type(m.q_state) << endl; 


      if (leaf_coords.count(m.toVector())<1)
	{
	  if ( m.left_state  == left_profile.end_state || m.left_state  == left_profile.pre_end_state)
	    {
	      if ( m.right_state  == right_profile.end_state || m.right_state  == right_profile.pre_end_state)
		// do not store coordinates
		{}
	      else 
		if (side == "right")
		  leaf_coords[m.toVector()] = right_profile.leaf_coords[m.right_state];
	    }
	  
	  else if ( m.right_state  == right_profile.end_state || m.right_state  == right_profile.pre_end_state)
	    {
	      if (side == "left")
		leaf_coords[m.toVector()] = left_profile.leaf_coords[m.left_state];
	    }
	  
	  else
	    {
	      if (side == "both")
		leaf_coords[m.toVector()] = merge_coords(left_profile.leaf_coords[m.left_state],
							      right_profile.leaf_coords[m.right_state]);
	      else if (side == "left")
		leaf_coords[m.toVector()] = left_profile.leaf_coords[m.left_state];
	      else if (side == "right")
		leaf_coords[m.toVector()] = right_profile.leaf_coords[m.right_state];
	    }
	}

      if (leaf_seq_coords.count(m.toVector())<1)
	{
	  if ( m.left_state  == left_profile.end_state || m.left_state  == left_profile.pre_end_state || 
	       m.left_state == left_profile.start_state)
	    {
	      if ( m.right_state  == right_profile.end_state || m.right_state  == right_profile.pre_end_state ||
		   m.right_state  == right_profile.start_state)
		// do not store coordinates
		{}
	      //	  if ( m.left_state  == left_profile.end_state || m.right_state  == right_profile.end_state )
	      //	    {}
	      else 
		{
		  if (side == "right")
		    leaf_seq_coords[m.toVector()] = right_profile.leaf_seq_coords[m.right_state];
		}
	    }
 	  else if ( m.right_state  == right_profile.end_state || m.right_state  == right_profile.pre_end_state || 
		    m.right_state  == right_profile.start_state)
	    {
	      if (side=="left")
		leaf_seq_coords[m.toVector()] = left_profile.leaf_seq_coords[m.left_state];
	    }
	  
	  else
	    {
	      if (side=="both")
		leaf_seq_coords[m.toVector()] = merge_seq_coords(left_profile.leaf_seq_coords[m.left_state],
								 right_profile.leaf_seq_coords[m.right_state]);
	      // added afterwards
	      // 	      else if (side=="left")
//  		leaf_seq_coords[m.toVector()] = left_profile.leaf_seq_coords[m.left_state];
// 	      else if (side=="right")
// 		leaf_seq_coords[m.toVector()] = right_profile.leaf_seq_coords[m.right_state];
	    }

	}
    }
  

  // if this is a novel transition, catalog it, possibly factoring out R-transition weights
  if (sampled_transition_weight.count(transitionPair)<1)
	{
	  if (is_external(mPrime) || is_pre_end(mPrime) )
		// we need only factor out R's transition weight if it is an external
		// transition.  Otherwise, R has remained the same. 
		sampled_transition_weight[transitionPair] = 
		  weight / (get_DP_cell(m) * Q.get_R_transition_weight(m.q_state, mPrime.q_state));
	  else
		sampled_transition_weight[transitionPair] = 
		  weight / (get_DP_cell(m));

	  if (sampled_transition_weight[transitionPair] > 1.001)
		{
		  //		  cerr<< "Warning: sampled transition weight was greater than 1.0\n"; 
		  //		  cerr<< "Large trans:" << sampled_transition_weight[transitionPair] << endl; 
		  sampled_transition_weight[transitionPair] = 0.99999999; //for now..

// 		  std::cerr<<"Warning: transition weight was cached as greater than one\nSome info:\n";
// 		  std::cerr<<"Source state:\n\t"; m.display(Q);
// 		  std::cerr<<"Destination state:\n\t"; mPrime.display(Q); 		  
// 		  std::cerr<<"Sampling weight: "<<weight<<endl; 
// 		  std::cerr<<"DP cell for source state: "<<get_DP_cell(m)<<endl; 
// 		  if(is_external(mPrime))
// 			{
// 			  std::cerr<<"External state,\n";
// 			  std::cerr<<"Emission weight for destination state: "<<compute_emission_weight(mPrime)<<endl; 
// 			  std::cerr<<"R's transition weight: "<<Q.get_R_transition_weight(m.q_state, mPrime.q_state)<<endl;
			  
// 			}
// 		  else 
// 			{
// 			  std::cerr<<"External state,\n";
// 			  std::cerr<<"Emission weight for destination state: "<<compute_emission_weight(mPrime)<<endl; 
// 			}
		}
	}
}

void Profile::store_summed_nulls(vector<M_id> path)
{
  M_id m; 
  int i; 
  pair<vector<int>, vector<int> > externalPair; 
  vector<M_id> nulls; 
  // the path is reverse (e.g. start state is at the end), so the indexing is reversed...
  for (i=path.size()-1; i>-1; i--)
    {
      m = path[i];
      if (is_start(m) || is_external(m) || is_pre_end(m) || is_end(m) )
	{
	  if (nulls.empty())
	    externalPair.first = m.toVector(); 
	  else
	    {
	      externalPair.second = m.toVector(); 
	      summed_nulls[externalPair] = nulls; 
	      //std::cerr<<"Nulls added at node" << treeNode << " : "<<summed_nulls.size() <<endl; 
	      nulls.clear(); 
	    }
	}
      else
	nulls.push_back(m);
    }
}
     
      

void Profile::cache_path(vector<M_id> path)
{
//   store the characters in-between two external states in a sampled path.  
//   it is much easier to do this on a path than on a DAG
//   The basic idea is that we fix a unique sequence of null states between each pair of 
//   external states.  Here we store that mapping, and call it up later when displaying multiple 
//   alignments.  
  bool logging = false; 
  int i,j; 
  M_id ext_start, right;
  MyMap<node, string> btwn_tmp, left_tmp, right_tmp; 
  pair< vector<int> , vector<int> > externalPair; 
  pair< state, state > profileExternalPair; 
  vector<node> toPad; // we'll need to pad one side with gaps.  
  int latestLeft, latestRight; 
  // the path is reverse (e.g. start state is at the end), so the indexing is reversed...
  for (i=path.size()-1; i>-1; i--)
	{
	  ext_start = path[i]; 
	  latestLeft = ext_start.left_state;
	  latestRight = ext_start.right_state; 	  

	  if (!is_external(ext_start) && i!=path.size()-1) continue;
	  btwn_tmp.clear(); right_tmp.clear(); left_tmp.clear(); 
	  
	  for (j=i-1; j>-1; j--)
		{
		  // right refers to the right state in the path, not the right profile.  
		  right=path[j];
		  // we'll add an empty map or two here, but I think that's OK
		  // if external        or   pre_end_state   or  end state
		  if (is_external(right) || is_pre_end(right) || j==0)
			{
			  profileExternalPair.first = latestLeft; 
			  profileExternalPair.second = right.left_state; 
			  if(left_profile.between.count(profileExternalPair) > 0)
				{
				  if(logging) std::cerr<<"concatenating left w/ left\n";
				  left_tmp = cat_STP(left_tmp, left_profile.between[profileExternalPair]);
				}
			  if (!is_flush(left_tmp)) std::cerr<<"STP not flush after a ext-ext joining, left\n";

			  profileExternalPair.first = latestRight;
			  profileExternalPair.second = right.right_state; 
			  if(right_profile.between.count(profileExternalPair) > 0)
				{
				  if(logging) std::cerr<<"concatenating right w/ right\n";
				  right_tmp = cat_STP(right_tmp, right_profile.between[profileExternalPair]);
				}
			  
			  if (!left_tmp.empty())
				{
				  toPad = right_profile.subtreeNodes; toPad.push_back(treeNode);
				  left_tmp = pad_STP(left_tmp, toPad);
				  if(logging) std::cerr<<"concatenating left w/ btwn \n";
				  btwn_tmp = cat_STP(btwn_tmp, left_tmp);
				}			  

			  if (!right_tmp.empty())
				{
				  toPad = left_profile.subtreeNodes; toPad.push_back(treeNode);
				  right_tmp = pad_STP(right_tmp, toPad); 
				  if(logging) std::cerr<<"concatenating right w/ btwn \n";
				  btwn_tmp = cat_STP(btwn_tmp, right_tmp);			  
				}
			  
			  if (!is_flush(btwn_tmp)) std::cerr<<"STP not flush after a ext-ext joining, right. \n";

			  if (btwn_tmp.size() > 0)
				{
				  externalPair.first = ext_start.toVector();
				  externalPair.second = right.toVector(); 			  
				  between[externalPair] = btwn_tmp; 
				  if(logging)
					{
					  std::cerr<<"Stored the in-between characters for the two states:\n";
					  std::cerr<<"Start:\n";
					  ext_start.display(Q); 
					  std::cerr<<"End:\n";
					  right.display(Q); 
					  std::cerr<<"State:\n";
					  show_state_phylo(btwn_tmp);
					  std::cerr<<"flush?:"<<is_flush(btwn_tmp)<<endl;
					}
				}
			  btwn_tmp.clear(); right_tmp.clear(); left_tmp.clear(); 
			  break;
			}
		  else if(is_left_int(right))
			{
			  if(logging) std::cerr<<"Left internal...\n";
			  profileExternalPair.first = latestLeft;
			  profileExternalPair.second = right.left_state; 
			  latestLeft = right.left_state; 
			  if(left_profile.between.count(profileExternalPair) > 0)
				{
				  if(logging) std::cerr<<"concatenating left w/ left, internal \n";
				  left_tmp = cat_STP(left_tmp, left_profile.between[profileExternalPair]);
				  toPad = right_profile.subtreeNodes; toPad.push_back(treeNode);
				  left_tmp = pad_STP(left_tmp, toPad);
				}
			  left_tmp = cat_STP(left_tmp, state_type_phylogeny[right.toVector()]); 
			  btwn_tmp = cat_STP(btwn_tmp, left_tmp); 
			  left_tmp.clear(); 
			  
			  if (!is_flush(left_tmp)) std::cerr<<"STP not flush after a left-int joining\n";
			  
// 			  if (left_tmp.size() < 1 ) std::cerr<<"Warning: tmp is still empty, after additions...\n";
// 			  else if(logging) {std::cerr<<"Current btwn: \n"; show_state_phylo(left_tmp); }
			}

		  else if(is_right_int(right))
			{
			  if(logging) std::cerr<<"Right internal...\n";
			  profileExternalPair.first = latestRight;
			  profileExternalPair.second = right.right_state; 
			  latestRight = right.right_state; 
			  if(right_profile.between.count(profileExternalPair) > 0)
				{
				  if(logging) std::cerr<<"concatenating right w/ right, internal \n";
				  right_tmp = cat_STP(right_tmp, right_profile.between[profileExternalPair]);
				  toPad = left_profile.subtreeNodes; toPad.push_back(treeNode);
				  right_tmp = pad_STP(right_tmp, toPad);
				}
			  right_tmp = cat_STP(right_tmp, state_type_phylogeny[right.toVector()]); 
			  btwn_tmp = cat_STP(btwn_tmp, right_tmp); 
			  right_tmp.clear(); 

			  if (!is_flush(right_tmp)) std::cerr<<"STP not flush after a right-int joining\n";

	// 		  if (right_tmp.size() <1) std::cerr<<"Warning: tmp is still empty, after additions...\n";
// 			  else if(logging) {std::cerr<<"Current btwn: \n"; show_state_phylo(right_tmp); }
			}
		}
	}
}
	
  

void Profile::sum_paths_to(M_id mPrime, bool inLog, bool logging)
{
  // The core DP function - sum over paths into the composite state mPrime.  
  #ifdef DART_DEBUG
  if (get_DP_cell(mPrime) != 0.0){
	std::cerr<<"Asked to fill cell which already has entry!\n"; exit(1);
  }
  #endif
  
  // This function fills the cell at mPrime in Z
  bool testing = false;
  bool fromStart = 0;
  bfloat toAdd, finalSum = 0.0, left_transition_weight;
  bfloat small = 1e-500; 
  vector<int> mPrimeVec = mPrime.toVector(); 
  vector<state> left_incoming;
  vector<state> right_incoming;   
  vector<state> q_incoming; 
  pair< vector<int>, vector<int> > transitionPair; 
  vector<state>::iterator e_l, e_r, q, qPrime; 
  M_id m; 
  vector<M_id> pi; 
  bool showState = 0; 
  bfloat emissionWeight = compute_emission_weight(mPrime);
  // underflow test
  //if (!emissionWeight>0)
  //return;

  if (showState)
	{
	  state_type_phylogeny[mPrime.toVector()] = merge_STP(mPrime);
	  pi.clear(); 
	  pi.push_back(mPrime);
	  std::cerr<<"The state summed to, emissionweight: "<< emissionWeight<<"\n";
	  mPrime.display(Q);
	  show_alignment(pi, false); 
	}
	
  if(logging) std::cerr<<"Emission weight for this state: "<<emissionWeight<<endl;
  #ifdef DART_DEBUG
  if (emissionWeight < small) { std::cerr<<"Warning: very small emission weight for state: \n";mPrime.display(Q);}
  #endif

  //  Z[mPrime.toVector()] = 0;
  //  add_to_DP_cell(mPrime, 0.0 ); 

  if (Q.has_transition(Q.composite_start_state, mPrime.q_state) &&  \
	  left_profile.has_transition(left_profile.start_state, mPrime.left_state) && \
	  right_profile.has_transition(right_profile.start_state, mPrime.right_state))
	{ 
	  fromStart = 1; 
	  toAdd = Q.get_transition_weight(Q.composite_start_state, mPrime.q_state);
	  toAdd *= left_profile.get_transition_weight(left_profile.start_state, mPrime.left_state);
	  toAdd *= right_profile.get_transition_weight(right_profile.start_state, mPrime.right_state)*\
		emissionWeight;

	  // for test
	  if (testing)
		{
		  m.q_state = Q.composite_start_state;
		  m.left_state = left_profile.start_state;
		  m.right_state = right_profile.start_state; 
		  m.left_type = -1;
		  m.right_type = -1; 	 
		  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
		  transition_weight_test[transitionPair] = toAdd/(emissionWeight*get_DP_cell(m)); 
		}	  
	  
	}

  else if ( Q.has_transition(Q.composite_start_state, mPrime.q_state) && \
	   mPrime.left_state == left_profile.start_state &&					\
	   right_profile.has_transition(right_profile.start_state, mPrime.right_state))
	// Implicitly set start->start transition to 1 here. 
	{ 
	  fromStart =1;
	  toAdd = \
		Q.get_transition_weight(Q.composite_start_state, mPrime.q_state)* \
		right_profile.get_transition_weight(right_profile.start_state, mPrime.right_state)*	\
		emissionWeight;

	  // for test
	  if(testing)
		{
		  m.q_state = Q.composite_start_state;
		  m.left_state = left_profile.start_state;
		  m.right_state = right_profile.start_state; 
		  m.left_type = -1;
		  m.right_type = -1; 	 
		  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
		  transition_weight_test[transitionPair] = toAdd/(emissionWeight*get_DP_cell(m)); 
		}
	}
  else if (logging) std::cerr<<"No incoming transitions from start detected!\n";

  if (fromStart)
	{
	  // inLog
	  // underflow test
	  // if (true)
	  if (toAdd > 0.0) 
	    {
	      if (inLog)
		{
		  M_id start_state;
		  start_state.q_state = Q.composite_start_state; 
		  start_state.left_state = left_profile.start_state; 
		  start_state.right_state = right_profile.start_state; 
		  start_state.left_type = -1;
		  start_state.right_type = -1; 

		  incoming[mPrimeVec].push_back(start_state);
		  outgoing[start_state.toVector()].push_back(mPrime); 
		  transitionPair.first = start_state.toVector(); transitionPair.second = mPrimeVec; 
		  transition_weight[transitionPair] = toAdd / emissionWeight; 
		  if (1) std::cerr<<"Adding contribution from start state as source: " << toAdd<<endl;
		}


	      finalSum +=toAdd; 
	    }
	}


  // incoming states where Q is match
  //  if (mPrime.left_state !=0 && mPrime.right_state !=0 && mPrime.left_state != left_profile.start_state)
  if (mPrime.right_state != right_profile.start_state && mPrime.left_state != left_profile.start_state)	
	{
	  q_incoming = Q.get_incoming_match_states(mPrime.q_state);
	  left_incoming = left_profile.get_incoming(mPrime.left_state);
	  right_incoming = right_profile.get_incoming(mPrime.right_state);  

	  for (e_l = left_incoming.begin(); e_l != left_incoming.end(); e_l++)
		{
	  
		  m.left_state = *e_l;
		  if (*e_l == left_profile.start_state) 
		    m.left_type = -1; 	  
		  else
		    m.left_type = 1;

		  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
			{
			  m.right_state = *e_r;
			  if (*e_r == right_profile.start_state) 
			    m.right_type = -1; 	  
			  else 
			    m.right_type = 1;
		  
			  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
				{
				  m.q_state = *q;
				  toAdd = get_DP_cell(m); 
				  // underflow test
				  // if (true) 
				  if (toAdd > 0.0)
				    {
				      toAdd *=
					Q.get_transition_weight(*q, mPrime.q_state)*
					left_profile.get_transition_weight(*e_l, mPrime.left_state)*
					right_profile.get_transition_weight(*e_r, mPrime.right_state)*
					emissionWeight;
				      finalSum += toAdd; 
				      // inLog
				      if (inLog && toAdd >0.0)
					{
					  incoming[mPrimeVec].push_back(m);
					  outgoing[m.toVector()].push_back(mPrime); 
					  transitionPair.first = m.toVector(); transitionPair.second = mPrimeVec; 
					  transition_weight[transitionPair] = toAdd /  ( get_DP_cell(m) * emissionWeight); 
					}
				    }
				  // for test
				  if (testing)
					{
					  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
					  transition_weight_test[transitionPair] = toAdd/(emissionWeight*get_DP_cell(m)); 
					}
				  
				  if (logging) 
				    {
				      std::cerr<<"Adding contribution from source state having Q: "<<Q.get_state_name(m.q_state)<<" "<<toAdd<<endl;
				      std::cerr<<"The full state: forward value: "<<get_DP_cell(m)<<" \n\t"; m.display(Q);
				    }
				  
				}
			}
		}
	}
  //  if (mPrime.left_state !=0  && mPrime.left_state != left_profile.start_state)
  if (mPrime.left_state != left_profile.start_state)	
	{
	  // incoming states where Q is left-emit
	  q_incoming = Q.get_incoming_left_emit_states(mPrime.q_state);
	  left_incoming = left_profile.get_incoming(mPrime.left_state);  
	  for (e_l = left_incoming.begin(); e_l != left_incoming.end(); e_l++)
		{
		  m.left_state = *e_l;
		  if (*e_l == left_profile.start_state) continue; //m.left_type = -1; 	  
		  else m.left_type = 1;

		  m.right_state = mPrime.right_state; // note this
		  m.right_type = 0; // this is always wait, since right side "winds" back first. 
		  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
			{
			  m.q_state = *q;
			  toAdd = get_DP_cell(m); 
			  // underflow test
			  //if (true) //
			   if (toAdd > 0.0)
			    {
			      toAdd *= Q.get_transition_weight(*q, mPrime.q_state)*
				left_profile.get_transition_weight(*e_l, mPrime.left_state)*
				emissionWeight;
			  // right profile transition weight is implicitely 1 here 
			      finalSum += toAdd; 
			      // inLog
			      if (inLog && toAdd > 0.0)
				{
				  incoming[mPrimeVec].push_back(m);
				  outgoing[m.toVector()].push_back(mPrime); 
				  transitionPair.first = m.toVector(); transitionPair.second = mPrimeVec; 
				  transition_weight[transitionPair] = toAdd /  ( get_DP_cell(m) * emissionWeight); 
				}
			    }
			  // for test
			  if (testing)
				{
				  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
				  transition_weight_test[transitionPair] = toAdd/(emissionWeight*get_DP_cell(m)); 
				}

			  if (logging) 
			    {
			      std::cerr<<"Adding contribution from source state having Q: "<<Q.get_state_name(m.q_state)<<" "<<toAdd<<endl;
			      std::cerr<<"The full state: forward value: "<<get_DP_cell(m)<<" \n\t"; m.display(Q);
			    }
			  
			}
		}
	}
  //  if (mPrime.right_state !=0)
  if (mPrime.right_state != right_profile.start_state)	
	{
	  // incoming states where Q is right-emit
	  q_incoming = Q.get_incoming_right_emit_states(mPrime.q_state);
	  right_incoming = right_profile.get_incoming(mPrime.right_state);
  
	  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
		{

		  m.right_state = *e_r;
		  if (*e_r == right_profile.start_state) 
		    {
		      continue;
		      m.right_type = -1; 	  
		    }
		  else 
		    m.right_type = 1;

		  

		  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
			{
			  m.q_state = *q;
			  m.left_type = get_profile_type(*q, "left");
			  // begin careful stuff
			  if (m.left_type == profile_start)
 			    m.left_state = left_profile.start_state; // note this
			  else
			    m.left_state = mPrime.left_state; // note this

			  if (m.left_state != mPrime.left_state)
			    {
			      if (left_profile.has_transition(m.left_state, mPrime.left_state))
				left_transition_weight  = left_profile.get_transition_weight(m.left_state, mPrime.left_state); 
			      else
				continue; 
			    }
			  else
			    left_transition_weight  = 1.0; 
			  // end careful stuff
			  toAdd = get_DP_cell(m); 
			  // underflow test
			  // (true) // 

			  if (toAdd >0.0)
			    {
			      toAdd *= Q.get_transition_weight(*q, mPrime.q_state)*
				right_profile.get_transition_weight(*e_r, mPrime.right_state)* left_transition_weight *
				emissionWeight;
			      // left profile transition weight is implicitely 1 here 
			      finalSum +=toAdd; 
			      // inLog
			      if (inLog && toAdd > 0.0)
				{
				  incoming[mPrimeVec].push_back(m);
				  outgoing[m.toVector()].push_back(mPrime); 
				  transitionPair.first = m.toVector(); transitionPair.second = mPrimeVec; 
				  transition_weight[transitionPair] = toAdd /  ( get_DP_cell(m) * emissionWeight); 
				}
			    }
			  //for test
			  if(testing)
				{
				  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
				  transition_weight_test[transitionPair] = toAdd/(emissionWeight*get_DP_cell(m)); 
				}

			  if (logging) 
				{
				  std::cerr<<"Adding contribution from source state having Q: "<<Q.get_state_name(m.q_state)<<" "<<toAdd<<endl;
				  std::cerr<<"The full state: forward value: "<<get_DP_cell(m)<<" \n\t"; m.display(Q);
				}
		  
			  
			}
		}
	}
  if (finalSum > 0.0)
    add_to_DP_cell(mPrime, finalSum);
}
  

void Profile::clear_DP(void)
{
  //  MyMap< vector<int>, bfloat> Z_tmp;
  //  DP_hash Z_tmp; 
  //  Z_tmp.swap(Z); 
}
void Profile::fill_backward_DP(int logging)
{
  bool debug = false; 
  M_id m; //, mPrime; 
  vector<M_id>::iterator m_iter; 
  bfloat toAdd; 
  pair<vector<int>, vector<int> > transitionPair; 
  
  // Initialize recursion w/ end state 
  m.q_state = Q.composite_end_state; 
  m.left_state = left_profile.end_state; 
  m.right_state = right_profile.end_state; 
  m.left_type = 3; 
  m.right_type = 3; 
  add_to_backward_DP(m, 1.0);

  for (int backward_idx = backward_states.size()-1; backward_idx>=0; backward_idx--)
    {
      toAdd = 0.0;
      m = backward_states[backward_idx]; 
      if ( m.q_state == Q.composite_end_state)
	continue; 

      if (debug)
	{
	  std::cout<<"Size of state queue: " << backward_states.size()<<endl ;
	  std::cout<<"Current state:";
	  m.display(Q); 
	}
      
      for (m_iter = outgoing[m.toVector()].begin(); m_iter != outgoing[m.toVector()].end(); m_iter++)
	{
	  transitionPair.first = m.toVector(); transitionPair.second = m_iter->toVector(); 
	  if (get_backward_DP_cell(*m_iter) > 0.0 && transition_weight.count(transitionPair) > 0)
	    {
	      toAdd += compute_emission_weight(m) * get_backward_DP_cell(*m_iter) * transition_weight[transitionPair]; 
	      if (!toAdd > 0.0)
		{
		  std::cerr<< "Taking the product of: "<< compute_emission_weight(m) << " " << get_backward_DP_cell(*m_iter) << " " << transition_weight[transitionPair] <<endl; 
		  std::cerr<<"their product: " << toAdd <<endl; 
		  m.display(Q); 	      
		  m_iter->display(Q); 
		  exit(1);
		}
	    }
	}
      if(toAdd >0.0)
	add_to_backward_DP(m, toAdd);
      if (debug)
	std::cout<<"The current state has DP bits:  " << -log(get_backward_DP_cell(m))/log(2) << endl; 
    }

  m.q_state = Q.composite_start_state; 
  m.left_state = left_profile.start_state; 
  m.right_state = right_profile.start_state; 
  m.left_type = -1; 
  m.right_type = -1; 
  if ( abs(-log(get_backward_DP_cell(m))/log(2) - -log(forward_prob)/log(2)) > 0.0001)
    {
      std::cerr<< "The forward and backward recursions have arrived at sum-over-alignment bit scores: \n";
      std::cerr<<"\tBackward value of start state: " << -log(get_backward_DP_cell(m))/log(2) <<endl; 
      std::cerr<<"\tForward value of end state: " << -log(forward_prob)/log(2) <<endl; 
      std::cerr<<" This is problem, exiting...\n";
      exit(1); 
    }
}

void Profile::fill_DP(int logging, bool inLog)
{
  inLog = false; 
  // top -level DP function - traverse over all states in M_n (in the proper order) and call
  // up sum_paths_to for each of them. 
  M_id mPrime, m;
  vector<state>::iterator qPrime, e_l, e_r, q;
  vector<state> qStates, left_incoming, right_incoming;
  bfloat toAdd;
  state el_Prime, er_Prime;
  // bool testing = false; 
  tmpEmitTuple.resize(3);
  tmpEmitVals.resize(alphabet_size);

  // initialize Z at start state M_id
  mPrime.q_state = Q.composite_start_state;
  mPrime.left_state = left_profile.start_state;
  mPrime.right_state = right_profile.start_state;  
  mPrime.left_type = -1;
  mPrime.right_type = -1; 
  if(logging>=2)
	{
	  std::cerr<<"\nThe following state is being filled:\n";
	  mPrime.display(Q);
	}
  add_to_DP_cell(mPrime, 1);  
  // bLog
  //  backward_states.push_back(mPrime); 
  if (logging>=2) 
	{
	  std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
	}

  qStates = Q.get_SSSI_states();
  for (er_Prime=0; er_Prime < right_profile.num_delete_states; er_Prime++)
	{
	  //	  if (!is_in_envelope( left_profile.start_state, er_Prime )) // always returns true
	  //	   		  continue; 
	  for (qPrime = qStates.begin();qPrime != qStates.end(); qPrime++)
		{
		  mPrime.q_state = *qPrime; 
		  mPrime.left_state = left_profile.start_state;
		  mPrime.left_type = -1;		  
		  mPrime.right_state = er_Prime;
		  mPrime.right_type = 1; 		  
		  
		  if(logging>=2)
			{
			  std::cerr<<"\nThe following state is being filled:\n";
			  mPrime.display(Q);
			}
		  sum_paths_to(mPrime, inLog, logging>=2);
		  // bLog
		  //		  backward_states.push_back(mPrime); 
		  if (logging>=2) 
			{
			  std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
			}
		}
	}

  for (el_Prime=0; el_Prime < left_profile.num_delete_states; el_Prime++)
	{
	  for (er_Prime=0; er_Prime < right_profile.num_delete_states; er_Prime++)
		{
		  if (is_in_envelope(el_Prime, er_Prime, "both"))
		    {
		      //  q is match state, both right and left types are "1"
		      qStates = Q.get_match_states();
		      for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
			{
			  mPrime.q_state = *qPrime; 
			  mPrime.left_state = el_Prime; 
			  mPrime.left_type = 1;
			  mPrime.right_state = er_Prime; 
			  mPrime.right_type = 1;
			  if(logging>=2)
			    {
			      std::cerr<<"\nThe following state is being filled:\n";
			      mPrime.display(Q);
			    }
			  
			  sum_paths_to(mPrime, inLog, logging>=2);
			  // bLog
			  //			  backward_states.push_back(mPrime); 
			  if (logging>=2) 
			    {
			      std::cerr<<"\tForward value of state matching: "<< el_Prime << " and " << er_Prime << " : " << get_DP_cell(mPrime)<<endl;
			    }
			}
		    }
		  if (is_in_envelope(el_Prime, er_Prime, "left"))
		    {
		      //q is a left-emit state, left type is 1, right is 0
		      qStates = Q.get_left_emit_states();
		      for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
			{
			  mPrime.q_state = *qPrime; 
			  mPrime.left_state = el_Prime; 
			  mPrime.left_type = 1; 
			  mPrime.right_state = er_Prime; 
			  mPrime.right_type = 0; //get_profile_type(*qPrime, "right");
			  if(logging>=2)
			    {
			      std::cerr<<"\nThe following state is being filled:\n";
			      mPrime.display(Q);
			    }
			  
			  sum_paths_to(mPrime, inLog, logging>=2);
			  // bLog
			  //			  backward_states.push_back(mPrime); 
			  if (logging>=2) 
			    {
			      std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
			    }
			}
		    }
		  if (is_in_envelope(el_Prime, er_Prime, "right"))
		    {
		      //q is a right-emit state, right type is 1, left is *something* determined by get_profile_type(q, left)
		      qStates = Q.get_right_emit_states();
		      string S = "S";
		      string I = "I";			  
		      for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
			{
			  // avoid SSSI state, we fill these in a separate loop above
			  if ( stringAt(Q.get_state_type(*qPrime),0) == S &&
			       stringAt(Q.get_state_type(*qPrime),1) == S &&
			       stringAt(Q.get_state_type(*qPrime),2) == S &&
			       stringAt(Q.get_state_type(*qPrime),3) == I) continue; 
			  
			  mPrime.q_state = *qPrime; 
			  mPrime.left_state = el_Prime; 
			  mPrime.left_type = get_profile_type(*qPrime, "left");
			  mPrime.right_state = er_Prime; 
			  mPrime.right_type = 1; 
			  if(logging>=2)
			    {
			      std::cerr<<"\nThe following state is being filled:\n";
			      mPrime.display(Q);
			    }
			  
			  sum_paths_to(mPrime, inLog, logging>=2);
			  // bLog
			  //			  backward_states.push_back(mPrime); 
			  if (logging>=2) 
			    {
			      std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
			    }
			}
		    }
		}
	}

  // fill left-emit and right-emit where one side is pre_end_state, wait
  // somewhat skittish about this section...
  if (logging >=1)
    std::cerr<<" filling boundary cases..."; 
  // Q is left-emit 
  for (el_Prime=0; el_Prime < left_profile.num_delete_states; el_Prime++)
	{
	  //	  if (! is_in_envelope(el_Prime, right_profile.pre_end_state, "left" ))
	  //     continue; 
	  //q is a left-emit state, left type is 1, right is 0.  right index is right_profile.pre_end_state
	  qStates = Q.get_left_emit_states();
	  string S = "S";
	  string I = "I";			  
	  for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
		{
		  // for right_pre_end in right_profile.pre_end_states
		  mPrime.q_state = *qPrime; 
		  mPrime.left_state = el_Prime; 
		  mPrime.left_type = 1;
		  mPrime.right_state = right_profile.pre_end_state; 
		  mPrime.right_type = 0; //get_profile_type(*qPrime, "right");
		  if(logging>=2)
			{
			  std::cerr<<"\nThe following state is being filled:\n";
			  mPrime.display(Q);
			}

		  sum_paths_to(mPrime, inLog, logging>=2);
		  // bLog
		  //		  backward_states.push_back(mPrime); 
		  if (logging>=2) 
			{
			  std::cerr<<"\tForward value (workaround): "<< get_DP_cell(mPrime)<<endl;
			}
		}
	}

  // Q is right-emit 
  if (logging >=2)
    std::cerr<<"Beginning to fill right emit/left pre-end state..\n"; 
  for (er_Prime=0; er_Prime < right_profile.num_delete_states; er_Prime++)
	{
	  //	  if (! is_in_envelope(left_profile.pre_end_state, er_Prime, "right" ))
	  //	      continue; 
	  //q is a right-emit state, left type is 0, right is 1.  left index is left_profile.pre_end_state
	  qStates = Q.get_right_emit_states();
	  string S = "S";
	  string I = "I";			  
	  for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
		{
		  // avoid SSSI state, we fill these in a separate loop above
		  if ( stringAt(Q.get_state_type(*qPrime),0) == S &&
			   stringAt(Q.get_state_type(*qPrime),1) == S &&
			   stringAt(Q.get_state_type(*qPrime),2) == S &&
			   stringAt(Q.get_state_type(*qPrime),3) == I) continue; 
		  // for left_pre_end in left_profile.pre_end_states
		  mPrime.q_state = *qPrime; 
		  mPrime.left_state = left_profile.pre_end_state; 
		  mPrime.left_type = get_profile_type(*qPrime, "left");
		  mPrime.right_state = er_Prime;
		  mPrime.right_type = 1;
		  if(logging>=2)
			{
			  std::cerr<<"\nThe following state is being filled:\n";
			  mPrime.display(Q);
			}

		  sum_paths_to(mPrime, inLog, logging>=2);
		  // bLog
		  //		  backward_states.push_back(mPrime); 
		  if (logging>=2) 
			{
			  std::cerr<<"\tForward value (workaround): "<< get_DP_cell(mPrime)<<endl;
			}
		}
	}

  
  // Fill the WWWW,pre-end, pre-end state.  
  qStates = Q.get_wait_states();
  for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
	{
	  if (!Q.has_transition(*qPrime, Q.composite_end_state)) 
	    continue; // don't fill states which don't connect to end state
	  mPrime.q_state = *qPrime; 
	  mPrime.left_state = left_profile.pre_end_state; 
	  mPrime.left_type = 2;
	  mPrime.right_state = right_profile.pre_end_state; 
	  mPrime.right_type = 2;
	  if(logging>=2) // log for now
		{
		  std::cerr<<"\nThe following state is being filled:\n";
		  mPrime.display(Q);
		}

	  sum_paths_to(mPrime, inLog,logging>=2); //log for now, debugging
	  // bLog
	  //	  backward_states.push_back(mPrime); 
	  if (logging>=2) 
		{
		  std::cerr<<"\tForward value of wait/pre-end state: "<< get_DP_cell(mPrime)<<endl;
		}

	}

  // Initialize the end state
  mPrime.q_state = Q.composite_end_state;
  mPrime.left_state = left_profile.end_state;
  mPrime.right_state = right_profile.end_state;
  mPrime.right_type = 3;
  mPrime.left_type = 3; 
  
  //  Z[mPrime.toVector()] = 0;
  add_to_DP_cell(mPrime, 0.0 ); 
  // bLog
  //  backward_states.push_back(mPrime); 
  
  // Fill the end state cell
  // for left_pre_end in left_profile.pre_end_states:
    // for right_pre_end in right_profile.pre_end_states:

  //for test
  pair< vector<int>, vector<int> > transitionPair; 


  m.left_state = left_profile.pre_end_state;
  m.right_state = right_profile.pre_end_state;
  m.right_type = 2;
  m.left_type = 2;
  toAdd = 0.0; 
  for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
	{	  
	  if (!Q.has_transition(*qPrime, Q.composite_end_state)) continue;
	  m.q_state = *qPrime;
	  toAdd += get_DP_cell(m) * Q.get_transition_weight(*qPrime, Q.composite_end_state);
	} 
  
  if (true) //(toAdd > 0.0)
    {
      add_to_DP_cell(mPrime, toAdd); 
      // inLog
      if (inLog)
	{
	  incoming[mPrime.toVector()].push_back(m);
	  outgoing[m.toVector()].push_back(mPrime); 
	  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
	  transition_weight[transitionPair] = toAdd /  get_DP_cell(m) ;
	  if (! transition_weight[transitionPair] > 0 || transition_weight[transitionPair] > 1.0)
	    {
	      std::cerr<<"Faulty transition weight stored: " << transition_weight[transitionPair] << endl; 
	      m.display(Q);
	      mPrime.display(Q);
	    }
// 	  else
// 	    {
// 	      std::cerr<<"Correct transition weight stored: " << transition_weight[transitionPair] << endl; 
// 	      m.display(Q); 
// 	      mPrime.display(Q); 
// 	    }
	  
	}
    }

  if(logging>=2)
	{
	  std::cerr<<"\nThe following state has been filled:\n";
	  mPrime.display(Q);
	}
  if (!Z.count(mPrime))
    cerr<<"End state was never filled! Thus it is zero\n";
  forward_prob = get_DP_cell(mPrime);
  bool isReal = false;
  if (forward_prob>0.0)
    isReal = true; 
  if (!isReal)
    {
      std::cerr<<"Error: sum-over-alignments likelihood is zero: "<< forward_prob <<"  This is not good.  Check your input alignment and tree for irregularities\n";
      std::cerr<<"If there are no such irregularities, consider submitting a bug report, see www.biowiki.org/ProtPal for details.\n";
      exit(1);
    }
	
  if(logging>=2) std::cerr<<"forward-bits: "<<-log(get_DP_cell(mPrime))/log(2.0)<<endl;
  
}


// Medium-weight computation functions (transitions, emissions)
bfloat Profile::compute_emission_weight(M_id m)
{
  // Compute the emission weight for cascades within m.
  // A somewhat hairy function, though luckily the relevant defs are in 2.3.2  of transducer.tex

  // left, right, character indices, respectively
  // unsigned int char1, char2; 
  // NB we use '-1' as a proxy for the null character in the non-match states below. 

  bfloat weight = 0.0;
  string qClass = Q.get_state_class(m.q_state);
  // q is a match state
  if ( qClass == "match" ) 
    {
      for (alphIter=alphabet_ints.begin(); alphIter!=alphabet_ints.end(); alphIter++)
	tmpEmitVals[*alphIter] = 0.0;
      tmpEmitTuple[0] = m.q_state; 
      for (alphIter=alphabet_ints.begin(); alphIter!=alphabet_ints.end(); alphIter++)
	{
	  tmpEmitTuple[1] = *alphIter; 
	  for (alphIter2=alphabet_ints.begin(); alphIter2!=alphabet_ints.end(); alphIter2++)
	    {
	      tmpEmitTuple[2] = *alphIter2; 
	      tmpEmitVals[*alphIter2] += Q.emission_weight[tmpEmitTuple]*left_profile.get_absorb_weight(m.left_state, *alphIter);
	    }
	}
      for (alphIter=alphabet_ints.begin(); alphIter!=alphabet_ints.end(); alphIter++)
	weight += right_profile.get_absorb_weight(m.right_state, *alphIter)*tmpEmitVals[*alphIter];
      return weight;
    }
	  
  //older, slow way:

//   if ( qClass == "match")
// 	{
// 	  for (char1 = 0; char1 < alphabet_size; char1++)
// 	    {
// 	      for (char2 = 0; char2 < alphabet_size; char2++)		
// 		{
// 		  weight +=
// 		    Q.get_emission_weight(m.q_state, char1, char2)*
// 		    left_profile.get_absorb_weight(m.left_state, char1)*
// 		    right_profile.get_absorb_weight(m.right_state, char2);
// 		}
// 	    }
// 	  return weight;
// 	}
  
  // if q is left-emit
  if ( qClass == "left_ins" || qClass == "right_del")
	{
	  for (alphIter=alphabet_ints.begin(); alphIter!=alphabet_ints.end(); alphIter++)
		{
		  weight += \
			Q.get_emission_weight(m.q_state, *alphIter, -1)* \
			left_profile.get_absorb_weight(m.left_state, *alphIter);
		}
	  return weight;
	}

  // if q is right-emit
  if ( qClass == "right_ins" || qClass == "left_del")
	{
	  for (alphIter=alphabet_ints.begin(); alphIter!=alphabet_ints.end(); alphIter++)
	    {
	      weight +=						  
		Q.get_emission_weight(m.q_state, -1, *alphIter)*   
		right_profile.get_absorb_weight(m.right_state, *alphIter);
	    }
	  return weight;
	}
  return 1.0; 
}
		  

bfloat Profile::get_sampled_transition_weight(M_id e, M_id ePrime)
{
  // Look up the transition weight between two sampled states.  This ought to be in 
  // the table sampled_transition_weight, but becuase of the hairy-ness of operations
  // around this step, check first before returning

  pair< vector<int>, vector<int> > transitionPair;
  transitionPair.first = e.toVector(); 
  transitionPair.second = ePrime.toVector(); 
  
  if (sampled_transition_weight.count(transitionPair)<1 )
    {
      std::cerr<<"Error: requested a non-existant sampled transition weight.\n The offending call was:\n";
      std::cerr<<"Source state:\n\t";e.display(Q);
      std::cerr<<"Destination state:\n\t";ePrime.display(Q);      
      exit(1); 
    }
  return sampled_transition_weight[transitionPair]; 
}

bfloat Profile::get_internal_cascade_weight(M_id e)
{
  // for a non-external state e, return the weight of its internal cascades.  
  if (Q.get_state_class(e.q_state) == Q.left_ins || Q.get_state_class(e.q_state) == Q.right_ins)
    {
      return compute_emission_weight(e);
    }
  else
    {
      std::cerr<<"Error: calling internal cascade weight on a state that is not right-int or left-int.  The state was:\n\t";
      e.display(Q); 
      exit(1);
    }
}

bfloat Profile::get_external_cascade_weight(M_id e, int charIndex)
{
  // Return the weight of state e absorbing character alphabet[charIndex]
  string q_type = Q.get_state_class(e.q_state);
  int omegaPrime;
  bfloat d, out;
  out = 1;
  if (q_type == Q.match)
    {
      // sum both
      d=0;
      for (omegaPrime=0; omegaPrime<alphabet_size; omegaPrime++)
		{
		  d += Q.left_match_weight(e.q_state, charIndex, omegaPrime)*
			left_profile.get_absorb_weight(e.left_state, omegaPrime);
		}
      out *= d;
	  
      d=0;
      for (omegaPrime=0; omegaPrime<alphabet_size; omegaPrime++)
		{
		  d += Q.right_match_weight(e.q_state, charIndex, omegaPrime)*
			right_profile.get_absorb_weight(e.right_state, omegaPrime);
		}
      out =  out*d;
    }
      
  else if (q_type == Q.left_del)
    {
      // only sum D_r
      d=0;
      for (omegaPrime=0; omegaPrime<alphabet_size; omegaPrime++)
		{
		  d += Q.right_match_weight(e.q_state, charIndex, omegaPrime)*
			right_profile.get_absorb_weight(e.right_state, omegaPrime);
		}
      out = d*Q.left_delete_weight(e.q_state, charIndex);
      
    }
  
  else if (q_type == Q.right_del)
    {
      // only sum D_l
      d=0;
      for (omegaPrime=0; omegaPrime<alphabet_size; omegaPrime++)
		{
		  d += Q.left_match_weight(e.q_state, charIndex, omegaPrime)*
			left_profile.get_absorb_weight(e.left_state, omegaPrime);
		}
      out = d*Q.right_delete_weight(e.q_state, charIndex);
    }
  else
	{
	  std::cerr<<"Error: called external delete weight on a non-recognized q type:"<< q_type<<endl;
	  std::cerr<<Q.get_state_name(e.q_state)<<endl;
	  exit(1);
	}
  if (out>1.001)
    {
      std::cerr<<"Warning: absorb weight is greater than one, not reasonable\n";
      exit(1); 
    }
  return out; 
}

//Utility functions for DP containers - forward is the 'unnamed', and the backward is so explicitely named
inline void Profile::add_to_DP_cell(M_id m , bfloat toAdd)
{
  // Add the value toAdd to the DP cell for state m, if it exists.  otherwise initialize this value. 
  //  tmpIter = Z.find(tmpMidVec);
  tmpIter = Z.find(m); 
  if (tmpIter == Z.end())
    Z[m] = toAdd; 
  else
    tmpIter->second += toAdd; 
}
void Profile::add_to_backward_DP(M_id m , bfloat toAdd)
{
  tmpMidVec = m.toVector(); 
  // Add the value toAdd to the DP cell for state m
  if (backward_matrix.count(tmpMidVec) <1) backward_matrix[tmpMidVec] = toAdd; 
  else backward_matrix[tmpMidVec] += toAdd; 
}

bfloat Profile::get_backward_DP_cell(M_id m)
{
  tmpMidVec = m.toVector(); 
  // Access the value of DP cell for state m.  If this is nonexistant, return 0.  
  if (backward_matrix.count(tmpMidVec) <1) return 0; 
  else return backward_matrix[tmpMidVec];
} 



inline bfloat Profile::get_DP_cell(M_id m)
{
  tmpIter = Z.find(m); 
  if (tmpIter == Z.end())
    return 0.0; 
  else 
    return tmpIter->second; 
} 

inline bool Profile::is_in_envelope(state left_state, state right_state, string checks)
{
  //debugging logging
  bool  logging = false; 

  if ( use_guide_alignment ) // using guide alignment
    {  
      if (
	  left_state == left_profile.end_state || right_state == right_profile.end_state || 
	  left_state == left_profile.pre_end_state || right_state == right_profile.pre_end_state || 
	  left_state == left_profile.start_state || right_state == right_profile.start_state 
	  )
      	return true; 
      // for alignment envelope
      bool bad = false; 
      // leaf_seq_coords is a map from leaf nodes to integer.  Each absorbing profile (e.g. left_profile, right_profile) has one of these for each of its
      // absorbing states.  leaf_seq_coords[s][n] = a  means that the state s accounts for the a'th position of the sequence at 
      // leaf node n
      for ( envIter1 = left_profile.leaf_seq_coords[left_state].begin(); envIter1 != left_profile.leaf_seq_coords[left_state].end(); envIter1++)
	{
	  for ( envIter2 = right_profile.leaf_seq_coords[right_state].begin(); envIter2 != right_profile.leaf_seq_coords[right_state].end(); envIter2++)
	    {
	      name1 = node_names[envIter1->first]; 
	      name2 = node_names[envIter2->first]; 
	      //  int leftMin, leftMax, rightMin, rightMax; 	      
// 	      if (envelope->coordinates.count(name1))
// 		if (envelope->coordinates[name1].count(envIter1->second))
// 		  if (envelope->coordinates[name1][envIter1->second].count(name2))
// 		    {
// 		      leftMin = envelope->coordinates[name1][envIter1->second][name2].first; 
// 		      leftMax = envelope->coordinates[name1][envIter1->second][name2].second; 
// 		    }
// 		  else
// 		    cerr<<"Error\n";
// 		else
// 		  cerr<<"Error\n";
// 	      else
// 		cerr<<"Error\n";
	      
// 	      if (envelope->coordinates.count(name2))
// 		if (envelope->coordinates[name2].count(envIter2->second))
// 		  if (envelope->coordinates[name2][envIter2->second].count(name1))
// 		    {
// 		      rightMin = envelope->coordinates[name2][envIter2->second][name1].first;
// 		      rightMax = envelope->coordinates[name2][envIter2->second][name1].second;
// 		    }
// 		  else
// 		    cerr<<"Error\n";
// 		else
// 		  cerr<<"Error\n";
// 	      else
// 		cerr<<"Error\n";
	      
	      if ( checks == "left" )
		if ( envelope->coordinates[name1][envIter1->second][name2].first  > envIter2->second || // min is greater than R coordinate
		     envelope->coordinates[name1][envIter1->second][name2].second +1  < envIter2->second )// max+1 is less than R coordinate
		  {
		    ++num_discarded_states;
		    return false; 
		    //		    bad = true;
		  }

	      //old:
	      // the following line reads like this:
	      // is  coordinate of R character  within guide_sausage of the index at which the L character is  aligned 
	      // in the guide alignment ?
	      /*
		if ( abs( envIter2->second - envelope->coordinates[name1][envIter1->second][name2] ) > guide_sausage )
		  {
		    ++num_discarded_states;
		    bad  = true; 
		  }
	      */
	      if ( checks == "right" )
		if ( envelope->coordinates[name2][envIter2->second][name1].first > envIter1->second ||// min is greater than L coordinate
		     envelope->coordinates[name2][envIter2->second][name1].second +1 < envIter1->second )// max+1 is less than L coordinate
		  {
		    ++num_discarded_states;
		    return false; 
		    //		    bad = true;
		  }

	      /*
		// Reworded - 
		// is the left coordinate within g_s of the right character's 'true' aligned coordinate?

		if ( abs( envIter1->second - envelope->coordinates[name2][envIter2->second][name1] ) > guide_sausage )
		  {
		    ++num_discarded_states;
		    bad = true; 
		  }
	      */
	      if ( checks == "both" )
		if ( envelope->coordinates[name2][envIter2->second][name1].first  > envIter1->second || // min is greater than R coordinate
		     envelope->coordinates[name2][envIter2->second][name1].second < envIter1->second  ||// max is less than R coordinate

		     envelope->coordinates[name1][envIter1->second][name2].first  > envIter2->second ||// min is greater than R coordinate
		     envelope->coordinates[name1][envIter1->second][name2].second   < envIter2->second )// max is less than R coordinate
		  {
		    ++num_discarded_states;
		    return false; 
		    //		    bad = true;
		  }

		  
	      // Warning -don't turn on logging unless you're working on a very small (e.g. 3-5 characters) alignment

		if (logging)
		{
		  if (bad)
		    std::cerr<<"\nThe alignment of these two columns was NOT in the envelope: (checktype: " << checks << ") \n"; 
		  else
		    std::cerr<<"\nThe alignment of these two columns WAS in the envelope: (checktype: " << checks << ") \n"; 
		  
		  //		  std::cerr<<"The distance from guide alignment was: " << abs( envIter2->second - envelope->lookup(name1, envIter1->second, name2) ) << " " <<
		  //		    abs( envIter1->second - envelope->lookup(name1, envIter2->second, name1) ) << endl; 
// 		  std::cerr<<"Left state:\n"; 
// 		  show_state_phylo(left_profile.state_type_phylogeny[left_state]);
// 		  std::cerr<<"\nRight state:\n"; 
// 		  show_state_phylo(right_profile.state_type_phylogeny[right_state]);
	      
		  std::cerr<<"\nThe coordinates, min, max of the left state are: ";
		  for (MyMap<node,int>::iterator envIter = left_profile.leaf_seq_coords[left_state].begin(); envIter !=left_profile.leaf_seq_coords[left_state].end(); envIter++)
		    std::cerr<< envIter->second << " " << envelope->coordinates[name1][envIter->second][name2].first<< " " <<
		      envelope->coordinates[name1][envIter->second][name2].second << " "; 
	      
		  std::cerr<<"\nThe coordinates, min, max of the right state are: ";
		  for (MyMap<node,int>::iterator envIter = right_profile.leaf_seq_coords[right_state].begin(); envIter !=right_profile.leaf_seq_coords[right_state].end(); envIter++)
		    std::cerr<< envIter->second << " " <<  envelope->coordinates[name2][envIter->second][name1].first << " " << 
		      envelope->coordinates[name2][envIter->second][name1].second << " ";
		  
		  std::cerr<<"\n"; 
		}
		
		
		if (bad)
		  return false; 
	    }
	}
      return true; 
    }
  else //  a more  primitive  envelope
    {
      if (left_state == left_profile.pre_end_state || left_state == left_profile.end_state )
	return true;
      
      if (right_state == right_profile.pre_end_state || right_state == right_profile.end_state )
	return true; 

      pair<int,int> left_coords = left_profile.leaf_coords[left_state];
      pair<int,int> right_coords = right_profile.leaf_coords[right_state];

      // alignment envelope function.  Simple check for distance between most distant leaf coordinates>
      if (max(left_coords.second, right_coords.second) - min(left_coords.first, right_coords.first) < envelope_distance )
	return true;
      else 
	//std::cerr<<"Not in envelope!\n";
	return false; 
    }
}

bool is_flush(MyMap<node, string> map1)
{  
  // Query whether an STP has all its rows the same length.  
  MyMap<node, string>::iterator nodeState;  
  int alignLength = -1; 

  // check to make sure that all alignment rows are the same length
  for (nodeState=map1.begin(); nodeState!=map1.end(); nodeState++)
	{
	  if (alignLength == -1) alignLength = nodeState->second.size(); 
	  else if (alignLength != nodeState->second.size()) return false; 
	}
  return true; 
}
MyMap<node, string> pad_STP(MyMap<node, string> map1, vector<node> toPad)
{
  // Construct a new STP where entries for nodes in toPad have gaps of appropriate length.
  MyMap<node, string> outMap;
  MyMap<node, string>::iterator nodeState;  
  int alignLength = -1; 

  // check to make sure that all alignment rows are the same length
  for (nodeState=map1.begin(); nodeState!=map1.end(); nodeState++)
	{
	  if (alignLength == -1) alignLength = nodeState->second.size(); 
	  else if (alignLength != nodeState->second.size()) 
		{
		  std::cerr<<"Error: not all alignment rows are the same length; cannot pad with gaps\n";
		  exit(1);
		}
	}
  outMap = map1; 
  bool rePad=false; 
  string pad = "";
  for (int i=0; i <alignLength ; i++) pad += "-";
  for (vector<node>::iterator n=toPad.begin(); n!=toPad.end(); n++)
	{
	  //if (outMap.count(*n) != 0) { std::cerr<<"Error: padding a node which already exists!\n"; rePad=true;}
	  outMap[*n] = pad; 
	}
  if (rePad) 
	{
	  std::cerr<<"You asked to pad:\n";
	  show_state_phylo(map1); 
	  std::cerr<<"With the nodes: \n";
	  displayVector(toPad); 
	}
  return outMap; 
}

MyMap<node, string> cat_STP(MyMap<node, string> map1, MyMap<node, string> map2)
{
  /*  construct a new STP where all entries of m are added to those of mPrime.  
  for example: 
  if mPrime =       m = 
   l   a              t   
   r   a              -
   Q   -              -
  then cat_STP(m, mPrime)=
  l: at
  r: a-
  Q: --				 
  
  This is used when summing over null states.  
  */

  // if one is null, just return the other
  if (map1.size() == 0) return map2;
  if (map2.size() == 0) return map1;
   
  MyMap<node, string> outMap;
  MyMap<node, string>::iterator nodeState;  
  vector<node> map1_nodes, map2_nodes; 

  for (nodeState = map1.begin(); nodeState != map1.end(); nodeState++) map1_nodes.push_back(nodeState->first);
  for (nodeState = map2.begin(); nodeState != map2.end(); nodeState++) map2_nodes.push_back(nodeState->first);  

  sort( map1_nodes.begin(), map1_nodes.end());
  sort( map2_nodes.begin(), map2_nodes.end());  
 
  if (map1_nodes != map2_nodes)
	{
	  std::cerr<<"Error in concatenating STPs.  Tree nodes of the two states differ, cannot concatenate.\n";
	  std::cerr<<"Left side tree nodes:\n"; displayVector(map1_nodes);
	  std::cerr<<"\nRight side tree nodes:\n"; displayVector(map2_nodes); 	  
	  
	  std::cerr<<"\nLeft STP:\n";
	  show_state_phylo(map1);

	  std::cerr<<"Right STP:\n";
	  show_state_phylo(map2);	  
	  
	  exit(1);
	}
  
  for (nodeState = map1.begin(); nodeState != map1.end(); nodeState++)  outMap[nodeState->first] = nodeState->second; 
  for (nodeState = map2.begin(); nodeState != map2.end(); nodeState++)  outMap[nodeState->first] += nodeState->second; 

  return outMap; 
}

	  
	  
	  
  
MyMap<node, string> Profile::merge_STP(M_id m)
{
  // construct a new state_type_phylogeny map based on the input state m
  // This may involve padding states with gaps.  for instance, if one side is delete and the other
  // wait, the wait side will need gaps for each of its nodes of appropriate length, etc
  MyMap<node, string> outMap;
  MyMap<node, string>::iterator nodeState;  
  string qClass  = Q.get_state_class(m.q_state);
  bool logging = 0; 
  int left_length=1, right_length=1, difference=0;
  string gap = "-";
  if(logging) 
	{
	  std::cerr<<"Merging the state type phylogeny for:\n\t"; m.display(Q);
	  std::cerr<<"Size of left_profile STP:"<< left_profile.state_type_phylogeny[m.left_state].size()<<endl;
	  std::cerr<<"Size of right_profile STP:"<< right_profile.state_type_phylogeny[m.right_state].size()<<endl;
	}
  if (m.q_state == Q.composite_end_state) return outMap; //dummy return value	  

  if (qClass == Q.wait || m.q_state == Q.composite_start_state )
	{
	  outMap[treeNode]="";
	  if(logging) std::cerr<<"In left profile...\n";
	  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
		   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = "";
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  if(logging) std::cerr<<"In right profile...\n";
	  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
		   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = ""; 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}
	}


  if ( qClass == Q.match ) 
	{
	  // Both sides get to keep their node mappings, top gets *
	  if(logging) std::cerr<<"In left profile...\n";
	  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
		   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = nodeState->second;
		  left_length = nodeState->second.size(); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  if(logging) std::cerr<<"In right profile...\n";
	  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
		   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = nodeState->second; 
		  right_length = nodeState->second.size(); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  difference = left_length - right_length; 
	  if (difference !=0) {std::cerr<<"Joining states of unequal size, problem!\n"; exit(1);}
	  if (difference>0)
		{
		  // left side is longer, pad right side with difference number of gaps
		  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
			   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
			{ 
			  outMap[nodeState->first] = rep(difference, gap) + outMap[nodeState->first];
			}
		}
	  else if (difference < 0 )
		{
		  // right side is longer, pad left side with difference number of gaps
		  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
			   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
			{ 
			  outMap[nodeState->first] = rep(-difference, gap) + outMap[nodeState->first];
			}
		}
	  outMap[treeNode]=rep(max(left_length, right_length)-1, "-")+"*";
	}

   else if ( qClass == Q.right_del ) 
	{
	  // Left side gets to keep its node mappings, top gets *
	  if(logging) std::cerr<<"In left profile...\n";
	  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
		   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = nodeState->second;
		  left_length = nodeState->second.size(); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  if(logging) std::cerr<<"In right profile...\n";
	  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
		   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = rep(left_length, "-"); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}
	  outMap[treeNode]=rep(left_length-1, "-") + "*";
	}

   else if ( qClass == Q.left_del ) 
	{
	  // Right side gets to keep its node mappings, top gets *


	  if(logging) std::cerr<<"In right profile...\n";
	  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
		   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = nodeState->second; 
		  right_length = nodeState->second.size(); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  if(logging) std::cerr<<"In left profile...\n";
	  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
		   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
		{ outMap[nodeState->first] = rep(right_length, "-"); }

	  outMap[treeNode]=rep(right_length-1, "-") + "*";
	}
  
   else if ( qClass == Q.right_ins ) 
	{
	  // Right side gets to keep its node mappings, top gap
	  if(logging) std::cerr<<"In right profile...\n";
	  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
		   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = nodeState->second;
		  right_length = nodeState->second.size();
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  if(logging) std::cerr<<"In left profile...\n";
	  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
		   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = rep(right_length, "-"); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  outMap[treeNode]=rep(right_length-1, "-") + gap;
	}

   else if ( qClass == Q.left_ins ) 
	{
	  // Left side gets to keep its node mappings, top gap
	  if(logging) std::cerr<<"In left profile...\n";
	  for (nodeState = left_profile.state_type_phylogeny[m.left_state].begin();
		   nodeState != left_profile.state_type_phylogeny[m.left_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = nodeState->second; 
		  left_length = nodeState->second.size();
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}

	  if(logging) std::cerr<<"In right profile...\n";
	  for (nodeState = right_profile.state_type_phylogeny[m.right_state].begin();
		   nodeState != right_profile.state_type_phylogeny[m.right_state].end(); nodeState++)
		{ 
		  outMap[nodeState->first] = rep( left_length, "-"); 
		  if (logging) std::cerr<<"Mapped node:  "<<nodeState->first<< " to "<<nodeState->second<<endl;
		}
	  outMap[treeNode]=rep(left_length-1, "-") + gap;
	}
  return outMap;
}

void Profile::show_state(M_id m)
{
  // show a single state in M_n, by making an alignment and viewing it.  This is largely obsolete
  // by show_state_phylo, which doesn't require the input to be a state, just a map from nodes->strings
  vector<M_id> pi; 
  pi.push_back(m); 
  std::cout<< show_alignment(pi, false);
}


void show_state_phylo(MyMap<node, string>& stp)
{
  // visualize a single mapping from nodes to strings (e.g. an alignment).
  MyMap<node, string>::iterator nodeState;
  MyMap<node, string> subAlignment;    

  for (nodeState = stp.begin(); nodeState != stp.end(); nodeState++)
	{
	  if (subAlignment.count(nodeState->first) < 1) subAlignment[nodeState->first] = nodeState->second;
	  else subAlignment[nodeState->first] += nodeState->second;
		}

  for (nodeState = subAlignment.begin(); nodeState!= subAlignment.end(); nodeState++)
	{
	  std::cerr<< nodeState->first <<"\t"<<nodeState->second<<endl;
	}
}

int Profile::DP_size(void)
{
  return Z.size();
}

string Profile::show_alignment(vector<M_id> &pi, bool leaves_only)
{
  // Visualize a path through a DP matrix as a multiple sequence alignment.  All null states that 
  // were summed over in previous iterations are shown according to bool showNulls (I'm not
  // sure why we'd ever set this to false, but possibly if the root ancestor was the ONLY
  // sequence of interest, we'd like to trim out the extra alignment info).  
  bool logging = false, showNulls=true,foundJ; 
  string out; 
  int i,j; 
  M_id m,n; 
  MyMap<node, string>::iterator nodeState;
  MyMap<node, string> subAlignment;    
  pair<vector<int>, vector<int> > statePair;
  MyMap<node, string>::iterator b_iter;   
  if (!showNulls)
	{
	  for(i=pi.size()-1; i>-1; i--)
		{
		  m =pi[i]; 
		  if(logging)
			{
			  std::cerr<<"Adding state to alignment:\n";
			  m.display(Q); 
			  std::cerr<<"Alignment characters:\n";
			  show_state_phylo(state_type_phylogeny[m.toVector()]);
			}
	  
		  if (! is_flush(state_type_phylogeny[m.toVector()])) 
			{
			  std::cerr<<"STP is not flush!\n";
			  exit(1);
			}
		  for (nodeState = state_type_phylogeny[m.toVector()].begin(); 
			   nodeState != state_type_phylogeny[m.toVector()].end(); nodeState++)
			{
			  if (subAlignment.count(nodeState->first) < 1) subAlignment[nodeState->first] = nodeState->second;
			  else subAlignment[nodeState->first] += nodeState->second;
			}
		}
	}
	  // add the 'between' characters
  if(showNulls) // do it completely differnetly...use the cached path
	{
	  for(i=pi.size()-1; i>-1; i--)
		{
		  if (is_external(pi[i]) || pi[i].q_state == Q.composite_start_state)
			{
			  foundJ = false; 
			  if(logging)
				{
				  std::cerr<<"Left external state:"<<i<<endl;
				  show_state_phylo(state_type_phylogeny[m.toVector()]);
				  std::cerr<<"\n";
				}
			  m=pi[i];
			  for (nodeState = state_type_phylogeny[m.toVector()].begin(); 
				   nodeState != state_type_phylogeny[m.toVector()].end(); nodeState++)
				{
				  if (subAlignment.count(nodeState->first) < 1) subAlignment[nodeState->first] = nodeState->second;
				  else subAlignment[nodeState->first] += nodeState->second;
				}
			  for(j=i-1; j>-1; j--)
				{
				  if (foundJ) break; 
				  if(logging) std::cerr<<"Possible j: "<< j<<endl;
				  if (is_external(pi[j]) || is_pre_end(pi[j]))
					{
					  foundJ = true; 
					  if(logging) std::cerr<<"External j: "<< j<<endl;
					  n=pi[j];
					  statePair.first=m.toVector(); 
					  statePair.second=n.toVector(); 
					  if (between.count(statePair)>0)
						{
						  if(between[statePair].size() > 0)
							{
							  if(logging)
								{
								  std::cerr<<"Adding btw characters for states start: "<<i<<"\n";
								  m.display(Q); 
								  std::cerr<<"end: "<<j<<"\n";
								  n.display(Q);
								  std::cerr<<"state:\n";
								  show_state_phylo(between[statePair]);
								}
							  if (subAlignment.empty()) subAlignment=between[statePair]; 
							  else{
								for (b_iter=between[statePair].begin(); b_iter!=between[statePair].end(); b_iter++)
								  {
									if(subAlignment.count(b_iter->first) >0) 
									  subAlignment[b_iter->first] += b_iter->second; 
									else {
									  std::cerr<<"Error: in-between characters' tree nodes do not match external nodes!\n";
									  std::cerr<<"Tried to add\n";
									  show_state_phylo(between[statePair]); 
									  std::cerr<<"To:\n";
									  show_state_phylo(subAlignment); 
									  exit(1);
									}
								  }
							  }

							}
						}
					}
				}
			}
		}
	}
  // this is a late add-on, to make alignments a bit more readable...
  int maxNameLength = 0; 
  for (nodeState = subAlignment.begin(); nodeState!= subAlignment.end(); nodeState++)
	{
	  int nameSize = node_names[nodeState->first].size(); 
	  maxNameLength = max(maxNameLength, nameSize );
	}

  for (nodeState = subAlignment.begin(); nodeState!= subAlignment.end(); nodeState++)	  
	{
	  if (leaves_only && (!in(nodeState->first, leaves)))
		continue;
	  out +=node_names[nodeState->first] + rep(maxNameLength-node_names[nodeState->first].size()+4, " ") + nodeState->second + "\n"; 
	}
  return out; 
}

MyMap<string, string> Profile::alignment_map(vector<M_id> &pi, bool leaves_only)
{
  // Return a multiple sequence alignment as a string-string map (names to seqs).  All null states that 
  // were summed over in previous iterations are shown according to bool showNulls (I'm not
  // sure why we'd ever set this to false, but possibly if the root ancestor was the ONLY
  // sequence of interest, we'd like to trim out the extra alignment info).  
  bool logging = false, foundJ; 
  string out; 
  int i,j; 
  M_id m,n; 
  MyMap<node, string>::iterator nodeState;
  MyMap<node, string> subAlignment;    
  MyMap<string, string> alignment; 
  pair<vector<int>, vector<int> > statePair;
  MyMap<node, string>::iterator b_iter;   
  for(i=pi.size()-1; i>-1; i--)
    {
      if (is_external(pi[i]) || pi[i].q_state == Q.composite_start_state)
	{
	  foundJ = false; 
	  if(logging)
	    {
	      std::cerr<<"Left external state:"<<i<<endl;
	      show_state_phylo(state_type_phylogeny[m.toVector()]);
	      std::cerr<<"\n";
	    }
	  m=pi[i];
	  for (nodeState = state_type_phylogeny[m.toVector()].begin(); 
	       nodeState != state_type_phylogeny[m.toVector()].end(); nodeState++)
	    {
	      if (subAlignment.count(nodeState->first) < 1) subAlignment[nodeState->first] = nodeState->second;
	      else subAlignment[nodeState->first] += nodeState->second;
	    }
	  for(j=i-1; j>-1; j--)
	    {
	      if (foundJ) break; 
	      if(logging) std::cerr<<"Possible j: "<< j<<endl;
	      if (is_external(pi[j]) || is_pre_end(pi[j]))
		{
		  foundJ = true; 
		  if(logging) std::cerr<<"External j: "<< j<<endl;
		  n=pi[j];
		  statePair.first=m.toVector(); 
		  statePair.second=n.toVector(); 
		  if (between.count(statePair)>0)
		    {
		      if(between[statePair].size() > 0)
			{
			  if(logging)
			    {
			      std::cerr<<"Adding btw characters for states start: "<<i<<"\n";
			      m.display(Q); 
			      std::cerr<<"end: "<<j<<"\n";
			      n.display(Q);
			      std::cerr<<"state:\n";
			      show_state_phylo(between[statePair]);
			    }
			  if (subAlignment.empty()) subAlignment=between[statePair]; 
			  else{
			    for (b_iter=between[statePair].begin(); b_iter!=between[statePair].end(); b_iter++)
			      {
				if(subAlignment.count(b_iter->first) >0) 
				  subAlignment[b_iter->first] += b_iter->second; 
				else {
				  std::cerr<<"Error: in-between characters' tree nodes do not match external nodes!\n";
				  std::cerr<<"Tried to add\n";
				  show_state_phylo(between[statePair]); 
				  std::cerr<<"To:\n";
				  show_state_phylo(subAlignment); 
				  exit(1);
				}
			      }
			  }
			  
			}
		    }
		}
	    }
	}
    }
  for (nodeState = subAlignment.begin(); nodeState!= subAlignment.end(); nodeState++)	  
	{
	  if (leaves_only && (!in(nodeState->first, leaves)))
		continue;
	  alignment[node_names[nodeState->first]]  = nodeState->second; 
	}
  return alignment; 
}


void Profile::show_DOT(ostream& out, string name)
{
  string from, to; 
  out << "digraph " << name << "{" <<endl; 
  pair< vector<int>, vector<int> > transitionPair; 
  for (vector<M_id>::iterator fromIter = sampled_states.begin(); fromIter!=sampled_states.end(); fromIter++)
    {
      from = getDOTname(*fromIter); 
      for (vector<M_id>::iterator toIter = sampled_outgoing[fromIter->toVector()].begin(); 
	   toIter!= sampled_outgoing[fromIter->toVector()].end(); toIter++)
	{
	  transitionPair.first = fromIter->toVector(); transitionPair.second = toIter->toVector(); 
	  to = getDOTname(*toIter); 
	  out << '"'<<from<<'"' << " -> " << '"' << to <<'"' << ";" <<endl; // " [label='" << sampled_transition_weight[transitionPair] << "']\n";
	}
    }
  out<< "}";
}

string Profile::getDOTname(M_id m)
{
  stringstream out; 
  if (is_start(m))
    out << "S";
  else if (is_end(m))
    out << "E";
  else if (is_pre_end(m))
    out << "P";
  else if (is_external(m))
    out << "D";
  else
    out << "I";

  out << "_";
  out << Q.get_state_type(m.q_state); 
  out << "_" << m.left_state; 
  out << "_" << m.right_state; 
  
  return string(out.str()); 
}


int Profile::get_profile_type(state qState, string side)
{
  // deduce the type of profile on the 'side' side, given the state of QTransducer is qState
  // This is handy when dealing with states like SSSI, where the left profile state is in START
  // rather than WAIT
  string qClass = Q.get_state_class(qState);
  pair<int, int> profileTypes;
  if (qClass == Q.match) 
	{
	  profileTypes.first = profile_delete;
	  profileTypes.second = profile_delete;
	}

  else if (qClass == Q.right_del) 
	{
	  profileTypes.first = profile_delete;
	  profileTypes.second = profile_wait;
	}

  else if (qClass == Q.left_del) 
	{
	  profileTypes.first = profile_wait; 
	  profileTypes.second = profile_delete;
	}

  else if (qClass == Q.null) 
	{
	  profileTypes.first = profile_wait; 
	  profileTypes.second = profile_wait; 
	}
  else if (qClass == Q.wait) 
	{
	  profileTypes.first = profile_wait; 
	  profileTypes.second = profile_wait;
	}
  else
	{
	  if (qClass == Q.left_ins ) 
		{
		  profileTypes.first = profile_delete; 
		  profileTypes.second = profile_wait;
		}

	  else if (qClass == Q.right_ins) 
		{

		  string qType = Q.get_state_type(qState);
		  profileTypes.second = profile_delete; 		  
		  if (qType == "SSSI") profileTypes.first = profile_start;
		  else if (qType == "IMMI") profileTypes.first = profile_delete;
		  else if (qType == "IMDI") profileTypes.first = profile_wait;		  		  
		}
	  else
		{
		  std::cerr<<"Unknown Q type, unable to determine profile types!"; 
		  std::cerr<<Q.get_state_name(qState); 
		  exit(1);
		}
	}
  //  std::cerr<<"Returning state types: "<<profileTypes.first<<" "<<profileTypes.second<<endl;
  if (side == "left") return profileTypes.first;
  else if (side =="right") return profileTypes.second;
  else 
    {
      std::cerr<<"Unrecognized side in get_profile_type: "<<side<<endl;
      exit(1); 
    }
  return 0; 
}

void Profile::depth_first_traversal(ostream& out, int sampled_start_state, int sampled_end_state)
{
    stack<M_id> stateStack;
    queue<M_id> stateQueue;
    M_id e, endState; 
    vector<M_id> children,visited;
    vector<M_id>::iterator child; 
    int counter=0; 
    
    endState = sampled_states[sampled_end_state];

    stateStack.push(sampled_states[sampled_start_state]); 
    //    stateQueue.push(sampled_states[sampled_start_state]); 
    out << "Beginning to traverse state space\n"; 
    while (!stateStack.empty())
      //while (!stateQueue.empty())
      {
	if (false)//counter > 100000000)
	  {
	    cerr<<"Infinite loop detected, exiting!\n";
	    out<<"Infinite loop detected, exiting!\n";
	    //	    exit(0); 
	  }
	//cerr<<"Current stack size:" << stateStack.size() << endl; 
	e = stateStack.top(); 
	stateStack.pop();
	// e = stateQueue.front(); 
	// stateQueue.pop();
	if (false) //index(e,visited) != -1)
	  cerr<<"Already visited state " << getDOTname(e) << endl; 
	else
	  {
	    ++counter; 
	    //visited.push_back(e); 
	    cerr<<"Current state off of stack: " << getDOTname(e) << endl; 
	    out<<getDOTname(e) << endl; 
	    if (sampled_outgoing.count(e.toVector()) )
	      {
		children = sampled_outgoing[e.toVector()];
		if (e == endState)
		  {
		    cerr<<"End state has " << children.size() << " children" << endl; 
		    cerr<<"Count in sampled_outgoing: " << sampled_outgoing.count(e.toVector()) << endl; 
		  }
		for (child = children.begin(); child != children.end(); child++)
		  {
		    stateStack.push(*child); 
		    //stateQueue.push(*child); 
		    if (e == endState)
		      cerr<<"pushing end's child: " << getDOTname(*child) << endl; 
		  }
	      }
	    else
	      cerr<<"The following state has no outgoing transitions: " << getDOTname(e) << endl; 
	  }
      }
    cerr << "Reached end of state space\n"; 
    out << "Reached end of state space\n"; 
}



	  
pair<int, int> merge_coords(pair<int, int> l, pair<int,int> r)
{
//   std::cerr<<"Merging coordinates:\n";
//   std::cerr<<"Left: "<< l.first<<" "<<l.second<<endl;
//   std::cerr<<"Right: "<< r.first<<" "<<r.second<<endl;   
  pair<int,int> out; 
  out.first=min(l.first, r.first);
  out.second=max(l.second, r.second);  
  //  std::cerr<<"Result: "<< out.first<<" "<<out.second<<endl;   
  return out; 
}


MyMap<node, int> merge_seq_coords( MyMap<node, int> l, MyMap<node, int> r)
{
  MyMap<node,int> out; 
  out.insert(l.begin(), l.end()); 
  out.insert(r.begin(), r.end()); 
  return out; 
}
