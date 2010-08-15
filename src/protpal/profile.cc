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
  return out; 
}

void M_id::display(QTransducer &Q)
{
  // Display a state in M_n, showing it's Q state, left/right states and their types.
  map<int, string> int2state;
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
  map<state, bfloat> outgoing;
  map<state, bfloat>::iterator stateIter;   
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
     - the characters of its children (map<state, map<node, string> > state_type_phylogeny)
     - min and max leaf_coordinates accounted for (map<state, pair<int, int> > leaf_coords)

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

  if(logging) std::cerr<<"\n ****** Marginalizing null states ******\n";
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
  string q_type, upsilon_type, bl_type, br_type;
  int left_type, right_type; 
  q_type = Q.get_state_type(m.q_state);
  upsilon_type = stringAt(q_type,1);
  bl_type = stringAt(q_type,2);
  br_type = stringAt(q_type,3);            
  right_type = m.right_type;
  left_type = m.left_type;
  if (upsilon_type == "M" && bl_type == "M" && left_type == profile_delete && 
      br_type == "M" && right_type == profile_delete ||  // MMDMD
      
      upsilon_type == "M" && bl_type == "M" && left_type == profile_delete &&
      br_type == "D" && right_type == profile_wait ||  // MMDDW
      
      upsilon_type == "M" && bl_type == "D" && left_type == profile_wait &&
      br_type == "M" && right_type == profile_delete || // MDWMD
      
      upsilon_type == "M" && bl_type == "D" && left_type == profile_wait &&
      br_type == "D" && right_type == profile_wait )  // MDWDW
	{
	  return 1; 
	}
  else return 0;
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

void AbsorbingTransducer::showDOT(Profile *sampled_profile)
{
  // show graphiviz -friendly output for the state space of an absorbing transducer.  
  // Warning: I haven't used this a whole lot, and so it's not very clean.  
  string name, nameDest; 
  vector<M_id> children;
  vector<M_id>::iterator child;
  unsigned int i; 
  for (i=0; i<sampled_profile->sampled_states.size(); i++)
	{
	  if (i == sampled_start_state) name="start";
	  else if (i == sampled_end_state) name="end";
	  //	  else if in(i, sampled_pre_end_state) name="pre_end";	  	  
	  else if (i == sampled_pre_end_state) name="pre_end";	  
	  else
		{
		  if (sampled_profile->is_external(sampled_profile->sampled_states[i]) ) name = "D";
		  else name = "N";
		}
	  if( sampled_profile->sampled_outgoing.count(sampled_profile->sampled_states[i].toVector()) >0)
		children = sampled_profile->sampled_outgoing[sampled_profile->sampled_states[i].toVector()];
	  else continue;
	  
	  for (child = children.begin(); child != children.end(); child++)
		{
		  if (*child == sampled_profile->sampled_states[sampled_start_state]) nameDest="start";
		  else if (*child == sampled_profile->sampled_states[sampled_end_state]) nameDest="end";
		  //else if in(*child, sampled_profile->sampled_states[sampled_pre_end_states]) nameDest="pre_end";	  		  
		  else if (*child == sampled_profile->sampled_states[sampled_pre_end_state]) nameDest="pre_end";	  
		  else
			{
			  if (sampled_profile->is_external(*child)) nameDest ="D"; 
			  else nameDest = "N";
			}

// 		  std::cerr<<"'"; displayVector(sampled_profile->sampled_states[i].toVector());
// 		  std::cerr<<name<<"'-> '"; displayVector(child->toVector());
// 		  std::cerr<<nameDest<<"';\n";
		}
	}
}
void AbsorbingTransducer::marginalize_null_states(Profile *sampled_profile)
{
  // This function has a lot going on.  Lots of comments follow. 
  // I've tried to follow the notation in transducer.tex as close as I can
  vector<int>::iterator sourceIdx;
  vector<M_id>::iterator ePrime;
  stack<M_id> stateStack;
  map<vector<int>, bfloat> weights;
  map<vector<int>, bfloat> outgoingSum;   
  map<vector<int>, int> outgoingCount;   
  M_id e, source; 
  vector<M_id> children;
  pair< vector<int> , vector<int> > transitionPair;
  bool logging = false;
  bfloat toAdd, small=pow(.1,1000); 

  // Loop over delete/external states.  We'd like to know the transition weight
  // of 'source'-> ePrime, where ePrime is another non-null state. 

  for (sourceIdx=sampled_externals.begin();sourceIdx!=sampled_externals.end();sourceIdx++)
    {
      source = sampled_profile->sampled_states[*sourceIdx];
	  if (logging){ std::cerr<<"\nSource state: \n\t"; source.display(sampled_profile->Q);}
      weights.clear();
      stateStack.push(source);
	  weights.clear(); 
      weights[source.toVector()]=1.0;
	  outgoingSum[source.toVector()] =0; 

      // search depth-first for next delete/external state
      while (!stateStack.empty())
		{
		  e = stateStack.top();
		  stateStack.pop();
		  if (sampled_profile->sampled_outgoing.count(e.toVector())>0)
			children = sampled_profile->sampled_outgoing[e.toVector()];
		  else continue;
		  if(logging) std::cerr<<"This state has "<< children.size()<< " children\n";

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

				  toAdd = weights[e.toVector()] * sampled_profile->get_sampled_transition_weight(e, *ePrime); 
				  if(logging)
					{
					  std::cerr<<"New transition added, value: "<<toAdd<<endl; 
					}
				  if (toAdd>1.001)
					{
					  std::cerr<<"Error: cached transition weight is greater than 1, after summing over null states\n";
					  exit(1);
					}
				  
				  // for test.  This is now done as a test_transitions function
				  /*
				  outgoingSum[source.toVector()] += toAdd;
				  outgoingCount[source.toVector()] += 1; 

				  if(outgoingSum[source.toVector()] >1.0)
					{
					  std::cerr<<"Error: outgoing sum of delete state has passed 1.  The last transition added: "<<
						toAdd<<endl; 
					  std::cerr<<"State has "<<outgoingCount[source.toVector()]<<" outgoing transitions\n";
					}
				  */
				  if (toAdd < small)
					{
					  std::cerr<<"Summed transition weight is very small.  Previous weight: "<<weights[e.toVector()]<<endl;
					  std::cerr<<"Current trans weight: "<< sampled_profile->get_sampled_transition_weight(e, *ePrime)<<endl;

					  std::cerr<<"Parent state: \n\t"; 
					  e.display(sampled_profile->Q); 
					  exit(1);
					}
				  if (transition_weight_tmp.count(transitionPair)==1)
					transition_weight_tmp[transitionPair] += toAdd; 
				  else transition_weight_tmp[transitionPair] = toAdd; 
				}
			  else // child ePrime is a null state - add to stack
				{ 
				  stateStack.push(*ePrime); 
				  
				  if (logging)
					{ 
					  std::cerr<<"\n\n***Bypassing non-external state:*** \n\t"; 
					  ePrime->display(sampled_profile->Q);
					  std::cerr<<"\n\n";
					}
				  
				  toAdd = weights[e.toVector()] *
					sampled_profile->get_sampled_transition_weight(e, *ePrime) *
					sampled_profile->get_internal_cascade_weight(*ePrime);

				  if (toAdd < small)
					{
					  std::cerr<<"Stored transition weight is very small.  Previous weight: "<<weights[e.toVector()]<<endl;
					  std::cerr<<"Current trans weight: "<< sampled_profile->get_sampled_transition_weight(e, *ePrime);
					  std::cerr<<"Internal emission weight: "<< sampled_profile->get_internal_cascade_weight(*ePrime);
					  exit(1);
					}
				  else if(toAdd >1.001)
					{
					  std::cerr<<"Error: Stored transition weight is greater than 1.\n";
					  exit(1);
					}
				  if(logging) 
					{
					  std::cerr<<"Caching state's weight: \n\t"; 
					  ePrime->display(sampled_profile->Q); 
					}
				  weights[ePrime->toVector()] = toAdd;
				}
			}
		}
    }
  // Check for super-normalized transitions after null-state-marginalization.  
  // Things may be sub-normalized, but super shouldn't occur...
  // NB this is now implemented as test_transitions function - I'm leaving this here
  // because it might be good for real-time debugging of the above function 
  /*
  for (map< vector<int>, bfloat>::iterator src=outgoingSum.begin(); src!=outgoingSum.end(); src++)
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
	      //double absorbSum= 0.0; 
	      for (charIndex = 0; charIndex < alphabet_size; charIndex++)
		absorption_weight[state_counter].push_back(
							   sampled_profile->get_external_cascade_weight(parent, charIndex));
	      
		  // The leaf  coordinates accounted for by this state. 
		  leaf_coords[state_counter] = sampled_profile->leaf_coords[parent.toVector()];
		  // The alignment column accounted for by this state.
		  state_type_phylogeny[state_counter] = sampled_profile->state_type_phylogeny[parent.toVector()];
		}

	  // If parent is the start state, we give it a "dummy" STP map
	  if (parent == sampled_profile->sampled_states[sampled_start_state] ) 
	    state_type_phylogeny[state_counter] = sampled_profile->merge_STP(parent); 	  

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
		  transition_weight[transitionIntPair] = 
			transition_weight_tmp[transitionVecPair];
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

//   std::cerr<<"\n Transferring 'between'\n"; 
//   // transfer the 'between' map from sampled_profile.    
//   map<pair<vector<int>, vector<int> >, map<node, string> >::iterator nodeState; 
//   for (nodeState=sampled_profile->between.begin(); nodeState!=sampled_profile->between.end(); nodeState++)
// 	{
// 	  transitionIntPair.first = mid2int[nodeState->first.first];
// 	  transitionIntPair.second = mid2int[nodeState->first.second]; 	  
// 	  between[transitionIntPair] = nodeState->second; 
// 	}
//   std::cerr<<"\n Transferring 'summed nulls'\n"; 
//   // transfer the 'summed_nulls' map from sampled_profile.    
//   map<pair<vector<int>, vector<int> >, vector<M_id> >::iterator nullState; 
//   for (nullState=sampled_profile->summed_nulls.begin(); nullState!=sampled_profile->summed_nulls.end(); nullState++)
// 	{
// 	  transitionIntPair.first = mid2int[nullState->first.first];
// 	  transitionIntPair.second = mid2int[nodeState->first.second]; 	  
// 	  summed_nulls[transitionIntPair] = nullState->second; 
// 	}


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
  return absorption_weight[e][charIndex];
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
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<ePrime<<endl;
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
	  if (transition_weight.count(transitionPair)) return transition_weight[transitionPair];
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

  envelope_distance = 200; //hard-coded envelope distance
  num_sampled_externals = 0; 

  // This is a bit awkward...owing to the quirk that M_ids store types as integers
  // rather than strings.  If I store all types now as variables it gets a bit less confusing.
  profile_start = -1;
  profile_wait = 0;
  profile_delete = 1;
  profile_pre_end = 2;
  profile_end = 3;
  

}
// Top-level public methods
state_path Profile::sample_DP(int num_paths, int logging, bool showAlignments, bool leaves_only, bool viterbi)
{
  // sample num_paths alignments from the DP matrix, showing the alignments according to 
  // showAlignments.  If showing alignments, they are printed in stockholm format, with
  // their bitscore and posterior probability as #=GF annotations.  
  
  bool testing = false;   //  bool fromStart;
  vector<M_id> states_in; // states for sampling
  vector<bfloat> weights; // weights for each of these states
  bfloat weight, pathWeight, emissionWeight, small=0.00001; 
  vector<M_id> pi; // sampled state path
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
	  if (num_sampled_externals > max_sampled_externals)
	    {
		  if (logging >=1 )
			{
			  std::cerr<<"(sampling truncated after caching " << num_sampled_externals << " states via ";
			  std::cerr<< pathIdx << " paths)...";
			}
		  break;
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
				  if (!Q.has_transition(*q, Q.composite_end_state)) continue;
				  m.q_state = *q;
				  m.left_state = left_profile.pre_end_state;
				  m.right_state = right_profile.pre_end_state; 			  
				  m.left_type = 2; m.right_type = 2;
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
					}
				  if (weight > 0.0)
					{
					  states_in.push_back(m);
					  weights.push_back( weight );
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
						{std::cerr<<"Using a novel transition in traceback!\n"; 
					  std::cerr<<"Source:"; m.display(Q);
					  std::cerr<<"Dest:"; mPrime.display(Q); 					  
					  exit(1);}
					  
					  else if( absoluted(transition_weight_test[transitionPair] - weight/get_DP_cell(m))>small )
						{std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
						  std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
						  std::cerr<<"Source:"; m.display(Q);
						  std::cerr<<"Dest:"; mPrime.display(Q); 					  
						  exit(1);}
					}

					  if (weight > 0.0)
						{
						  states_in.push_back(m);
						  weights.push_back( weight );
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
					  if (m.left_state == left_profile.start_state) m.left_type = -1;
					  else m.left_type = 1;

					  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
						{
						  if (! is_in_envelope(*e_l, *e_r))
							continue;
						  m.right_state = *e_r;
						  if (m.right_state == right_profile.start_state) m.right_type = -1;
						  else m.right_type = 1;

						  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
							{

							  m.q_state = *q;
							  weight = 
								get_DP_cell(m) *
								Q.get_transition_weight(*q, mPrime.q_state) * 
								left_profile.get_transition_weight(*e_l, mPrime.left_state) *
								right_profile.get_transition_weight(*e_r, mPrime.right_state) ;


							  // for test
							  if (weight > 0.0 && testing)
								{
							  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
							  if (transition_weight_test.count(transitionPair)<1)
								{std::cerr<<"Using a novel transition in traceback!\n"; 
								  std::cerr<<"Source:"; m.display(Q);
								  std::cerr<<"Dest:"; mPrime.display(Q); 					  
								  exit(1);}
							  
							  else if( absoluted(transition_weight_test[transitionPair] - weight/get_DP_cell(m))>small )
								{std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
								  std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
								  std::cerr<<"Source:"; m.display(Q);
								  std::cerr<<"Dest:"; mPrime.display(Q); 			
								  std::cerr<<"Q transition: "<<Q.get_transition_weight(*q, mPrime.q_state)<<endl; 
								  std::cerr<<"left transition: "<<left_profile.get_transition_weight(*e_l, mPrime.left_state)<<endl; 
								  std::cerr<<"right transition: "<<right_profile.get_transition_weight(*e_r, mPrime.right_state)<<endl; 
								  exit(1);}
								}

							  if (weight > 0.0)
								{
								  states_in.push_back(m);
								  weights.push_back( weight );
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
					  if (! is_in_envelope(*e_l, mPrime.right_state))
						continue; 

					  m.left_state = *e_l;
					  if (m.left_state == left_profile.start_state) m.left_type = -1;
					  else m.left_type = 1;
					  
					  m.right_state = mPrime.right_state;
					  
					  m.right_type = 0; 
					  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
						{
						  m.q_state = *q;
						  weight =
							get_DP_cell(m)*
							Q.get_transition_weight(*q, mPrime.q_state)*
							left_profile.get_transition_weight(*e_l, mPrime.left_state);
						  // right profile transition weight is implicitely 1 here 

						  // for test
						  if (weight > 0.0 && testing)
							{
						  transitionPair.first = m.toVector(); transitionPair.second = mPrime.toVector(); 
						  if (transition_weight_test.count(transitionPair)<1)
							{std::cerr<<"Using a novel transition in traceback!\n"; 
							  std::cerr<<"Source:"; m.display(Q);
							  std::cerr<<"Dest:"; mPrime.display(Q); 					  
							  exit(1);}
						  
						  else if( absoluted(transition_weight_test[transitionPair] - weight/get_DP_cell(m)) > small )
							{std::cerr<<"Transition in traceback does not match same computed in forward!\n"; 
							  std::cerr<<"DP, trace: "<<transition_weight_test[transitionPair]<<" "<<weight/get_DP_cell(m)<<endl; 
							  std::cerr<<"Source:"; m.display(Q);
							  std::cerr<<"Dest:"; mPrime.display(Q); 					  
							  exit(1);}
							}
						  if (weight > 0.0)
							{
							  states_in.push_back(m);
							  weights.push_back( weight );
							}
						}
					}
				}

			  // incoming states where Q is right-emit
			  // right profile must not be in start state
			  if (mPrime.right_state != right_profile.start_state) 				
				{
				  q_incoming = Q.get_incoming_right_emit_states(mPrime.q_state);
				  right_incoming = right_profile.get_incoming(mPrime.right_state);
  
				  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
					{
					  if (! is_in_envelope(mPrime.left_state , *e_r)) 
						continue;

					  m.right_state = *e_r;
					  if (m.right_state == right_profile.start_state) m.right_type = -1;
					  else m.right_type = 1;
					  
   					  m.left_state = mPrime.left_state; // note this

					  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
						{
						  m.left_type = get_profile_type(*q, "left");
						  m.q_state = *q;

						  if ( stringAt(Q.get_state_type(*q),0) == S &&
							   stringAt(Q.get_state_type(*q),1) == S &&
							   stringAt(Q.get_state_type(*q),2) == S &&
							   stringAt(Q.get_state_type(*q),3) == I) 
							{
							  if( mPrime.left_state != left_profile.start_state ) continue;
							  else 
								{
								  if( m.left_type != -1){std::cerr<<"Left profile incorrectly typed!\n"; exit(1);}
								}
							}// end if SSSI

						  weight =
							get_DP_cell(m)*		
							Q.get_transition_weight(*q, mPrime.q_state)*
							right_profile.get_transition_weight(*e_r, mPrime.right_state);
						  // left profile transition weight is implicitely 1 here 			  


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
							  std::cerr<<"Source:"; m.display(Q);
							  std::cerr<<"Dest:"; mPrime.display(Q); 					  
							  exit(1);}
							}
						  if (weight > 0.0)
							{
							  states_in.push_back(m);
							  weights.push_back( weight );
							}
						}
					}
				}
			}

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
					  map<node, string> tmpMap = merge_STP(*in); 
					  show_state_phylo(tmpMap);
					  std::cerr<<endl;
					  } */
				}

			  if (viterbi) sampleIdx = maxIndex(weights);
			  else sampleIdx = sample(weights);
			  pathWeight *= (weights[sampleIdx] / get_DP_cell(states_in[sampleIdx]))*
				compute_emission_weight(states_in[sampleIdx]);

			  if (logging>=2) std::cerr<<"Sampled state number: "<<sampleIdx<<endl;
			  pi.push_back(states_in[sampleIdx]);
			  cache_state(states_in[sampleIdx], mPrime, weights[sampleIdx] );
			  state_type_phylogeny[states_in[sampleIdx].toVector()] = merge_STP(states_in[sampleIdx]);
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
			  pi.push_back(states_in[0]);
			  cache_state(states_in[0], mPrime, weights[0] );
			  state_type_phylogeny[states_in[0].toVector()] = merge_STP(states_in[0]);
			  pathWeight *= (weights[0] / get_DP_cell(states_in[0]))*
				compute_emission_weight(states_in[0]);

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
  if(logging>=2) 
	{ std::cerr<<"\nFinal state:\n\t"; mPrime.display(Q); }

  // Store the null states' accounted characters:
  cache_path(pi); 
  //  store_summed_nulls(pi); 
  if(pathWeight/forward_prob > 1.001)
	{
	  std::cerr<<"\nError: an alignment's posterior probability was calculated as greater than 1.  This is not reasonable, and represents a calculation error. \n";
	  std::cout<<"#Sampled alignment, bit-score: "<< -log(pathWeight)/log(2)<<":\n";	  
	  std::cout<<  "       posterior probability: "<< pathWeight/forward_prob<<":\n";
	  show_alignment(pi, false); 
	  std::cerr<<"\n";
	  exit(1); 
	}
  alignString = show_alignment(pi, leaves_only); 
  if (showAlignments && !viterbi) 
	{
	  std::cout<<"#=GF bit_score "<< -log(pathWeight)/log(2)<<endl; 
	  std::cout<<  "#=GF post_prob "<< pathWeight/forward_prob<<endl; 
	  std::cout<<alignString; 
	  std::cerr<<"\n";
	}

	}	
  return pi; 
}
void Profile::cache_state(M_id m, M_id mPrime, bfloat weight)
{
  // upon sampling a move from m to mPrime, store the info about the state
  // m and the transition m -> mPrime
  
  // m is the new state
  if (index(m, sampled_states) == -1 ) 
	{
	  sampled_states.push_back(m);
	  if ( is_external(m) )
		  num_sampled_externals += 1; 
	}
  if (index(mPrime, sampled_outgoing[m.toVector()]) == -1) 
    { sampled_outgoing[m.toVector()].push_back(mPrime); }
  if (index(m, sampled_incoming[mPrime.toVector()]) == -1) 
    { sampled_incoming[mPrime.toVector()].push_back(m);   }
  


  if (leaf_coords.count(m.toVector())<1)
	{
	  if ( m.left_state  == left_profile.end_state || m.left_state  == left_profile.pre_end_state)
		if ( m.right_state  == right_profile.end_state || m.right_state  == right_profile.pre_end_state)
		  // do not store coordinates
		  {}
		else 
		  leaf_coords[m.toVector()] = right_profile.leaf_coords[m.right_state];

	  else if ( m.right_state  == right_profile.end_state || m.right_state  == right_profile.pre_end_state)
		leaf_coords[m.toVector()] = left_profile.leaf_coords[m.left_state];

	  else
		leaf_coords[m.toVector()] = merge_coords(left_profile.leaf_coords[m.left_state],
											 right_profile.leaf_coords[m.right_state]);
	}
  pair< vector<int>, vector<int> > transitionPair;
  transitionPair.first=m.toVector();   transitionPair.second=mPrime.toVector();
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

	  if (sampled_transition_weight[transitionPair] > 1.0)
		{
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
  // store the characters in-between two external states in a sampled path.  
  // it is much easier to do this on a path than on a DAG
  // The basic idea is that there is a unique sequence of null states between each pair of 
  // external states.  Here we store that mapping, and call it up later when displaying multiple 
  // alignments.  
  bool logging = false; 
  unsigned int i,j; 
  M_id ext_start, right;
  map<node, string> btwn_tmp, left_tmp, right_tmp; 
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
	
  

void Profile::sum_paths_to(M_id mPrime, bool inLog)
{
  // The core DP function - sum over paths into the composite state mPrime.  
  #ifdef DART_DEBUG
  if (get_DP_cell(mPrime) != 0.0){
	std::cerr<<"Asked to fill cell which already has entry!\n"; exit(1);
  }
  #endif
  
  // This function fills the cell at mPrime in Z
  bool logging = false, testing = false;
  bool fromStart = 0;
  bfloat toAdd, finalSum = 0.0;
  bfloat small = pow(.1,100);
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

  Z[mPrime.toVector()] = 0;

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
		}

	      if (logging) std::cerr<<"Adding contribution from start state as source: " << toAdd<<endl;
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
		  if (*e_l == left_profile.start_state) m.left_type = -1; 	  
		  else m.left_type = 1;

		  for (e_r = right_incoming.begin(); e_r != right_incoming.end(); e_r++)
			{
			  m.right_state = *e_r;
			  if (*e_r == right_profile.start_state) m.right_type = -1; 	  
			  else m.right_type = 1;
		  
			  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
				{
				  m.q_state = *q;
				  toAdd = get_DP_cell(m); 
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
				  
				  if (logging) std::cerr<<"Adding contribution from source state having Q: "<<Q.get_state_name(m.q_state)<<" "<<toAdd<<endl;
				   
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
		  if (*e_l == left_profile.start_state) m.left_type = -1; 	  
		  else m.left_type = 1;

		  m.right_state = mPrime.right_state; // note this
		  m.right_type = 0; // this is always wait, since right side "winds" back first. 
		  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
			{
			  m.q_state = *q;
			  toAdd = get_DP_cell(m); 
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

			  if (logging) std::cerr<<"Adding contribution from source state having Q: "<<Q.get_state_name(m.q_state)<<" "<<toAdd<<endl;
			  
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
		  if (*e_r == right_profile.start_state) m.right_type = -1; 	  
		  else m.right_type = 1;

		  m.left_state = mPrime.left_state; // note this

		  for (q = q_incoming.begin(); q != q_incoming.end(); q++)
			{
			  m.q_state = *q;
			  m.left_type = get_profile_type(*q, "left");
			  toAdd = get_DP_cell(m); 
			  if (toAdd >0.0)
			    {
			      toAdd *= Q.get_transition_weight(*q, mPrime.q_state)*
				right_profile.get_transition_weight(*e_r, mPrime.right_state)*
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
  map< vector<int>, bfloat> Z_tmp;
  Z_tmp.swap(Z); 
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
  // top -level DP function - traverse over all states in M_n (in the proper order) and call
  // up sum_paths_to for each of them. 
  M_id mPrime, m;
  vector<state>::iterator qPrime, e_l, e_r, q;
  vector<state> qStates, left_incoming, right_incoming;
  bfloat toAdd;
  state el_Prime, er_Prime;
  // bool testing = false; 
  inLog = true; // used to be that we only logged incoming info for param estimation...but currently we need it for recursive sampling. 
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
  backward_states.push_back(mPrime); 
  if (logging>=2) 
	{
	  std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
	}

  qStates = Q.get_SSSI_states();
  for (er_Prime=0; er_Prime < right_profile.num_delete_states; er_Prime++)
	{
	  if (!is_in_envelope( left_profile.start_state, er_Prime )) // always returns true
		  continue; 
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
		  sum_paths_to(mPrime, inLog);
		  // bLog
		  backward_states.push_back(mPrime); 
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
		  // check for envelope compatibility.  
		  if (! is_in_envelope( el_Prime, er_Prime ) )
			continue;
		  else
			{
			  //if (logging>=2) std::cerr<<"Examining state pair: "<< el_Prime <<" , " <<er_Prime<<endl;

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
				  sum_paths_to(mPrime, inLog);
				  // bLog
				  backward_states.push_back(mPrime); 
				  if (logging>=2) 
					{
					  std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
					}
				}

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
				  sum_paths_to(mPrime, inLog);
				  // bLog
				  backward_states.push_back(mPrime); 
				  if (logging>=2) 
					{
					  std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
					}
				}

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
				  sum_paths_to(mPrime, inLog);
				  // bLog
				  backward_states.push_back(mPrime); 
				  if (logging>=2) 
					{
					  std::cerr<<"\tForward value: "<< get_DP_cell(mPrime)<<endl;
					}
				}
	

			  // I don't think we actually need this set of states...
			  //q is a wait state, right, left types are 0, both are set to 0
			  if (0) //el_Prime == left_profile.num_delete_states-1 && er_Prime == right_profile.num_delete_states-1)
				{
				  qStates = Q.get_wait_states();
				  for (qPrime = qStates.begin(); qPrime != qStates.end(); qPrime++)
					{
					  mPrime.q_state = *qPrime; 
					  mPrime.left_state = el_Prime; 
					  mPrime.left_type = 0;
					  mPrime.right_state = er_Prime; 
					  mPrime.right_type = 0;
					  if(logging>=2)
						{
						  std::cerr<<"\nThe following state is being filled:\n";
						  mPrime.display(Q);
						}
					  sum_paths_to(mPrime, inLog);
					  // bLog
					  backward_states.push_back(mPrime); 
					}

				  if (logging>=2) 
					{
					  std::cerr<<"\tForward value of wait state: "<< get_DP_cell(mPrime)<<endl;
					}
				  if (get_DP_cell(mPrime)==0.0) {std::cerr<<"The following wait state has zero forward weight:"; mPrime.display(Q);}
				}
			}
		}
	}

  // fill left-emit and right-emit where one side is pre_end_state, wait
  // somewhat skittish about this section...

  // Q is left-emit 
  for (el_Prime=0; el_Prime < left_profile.num_delete_states; el_Prime++)
	{
	  if (! is_in_envelope(el_Prime, right_profile.pre_end_state ))//always true
		continue; 
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
		  sum_paths_to(mPrime, inLog);
		  // bLog
		  backward_states.push_back(mPrime); 
		  if (logging>=2) 
			{
			  std::cerr<<"\tForward value (workaround): "<< get_DP_cell(mPrime)<<endl;
			}
		}
	}

  // Q is right-emit 
  for (er_Prime=0; er_Prime < right_profile.num_delete_states; er_Prime++)
	{
	  if (! is_in_envelope(left_profile.pre_end_state, er_Prime ))
		continue; 
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
		  sum_paths_to(mPrime, inLog);
		  // bLog
		  backward_states.push_back(mPrime); 
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
	  // for left_pre_end in left_profile.pre_end_states:
	  // for right_pre_end in right_profile.pre_end_states:
	  mPrime.q_state = *qPrime; 
	  mPrime.left_state = left_profile.pre_end_state; 
	  mPrime.left_type = 2;
	  mPrime.right_state = right_profile.pre_end_state; 
	  mPrime.right_type = 2;
	  if(logging>=2)
		{
		  std::cerr<<"\nThe following state is being filled:\n";
		  mPrime.display(Q);
		}
	  sum_paths_to(mPrime, inLog);
	  // bLog
	  backward_states.push_back(mPrime); 
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
  
  Z[mPrime.toVector()] = 0;
  // bLog
  backward_states.push_back(mPrime); 
  
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
  
  if (toAdd > 0.0)
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
void Profile::add_to_DP_cell(M_id m , bfloat toAdd)
{
  tmpMidVec = m.toVector(); 
  // Add the value toAdd to the DP cell for state m
  if (Z.count(tmpMidVec) <1) Z[tmpMidVec] = toAdd; 
  else Z[tmpMidVec] += toAdd; 
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



bfloat Profile::get_DP_cell(M_id m)
{
  tmpMidVec = m.toVector(); 
  // Access the value of DP cell for state m.  If this is nonexistant, return 0.  
  if (Z.count(tmpMidVec) <1) return 0; 
  else return Z[tmpMidVec];
} 

bool Profile::is_in_envelope(state left_state, state right_state)
{
  if (left_state == left_profile.pre_end_state || left_state == left_profile.end_state )
	return true;

  if (right_state == right_profile.pre_end_state || right_state == right_profile.end_state )
	return true; 

  pair<int,int> left_coords = left_profile.leaf_coords[left_state];
  pair<int,int> right_coords = right_profile.leaf_coords[right_state];

  // alignment envelope function.  Simple check for distance between most distant leaf coordinates>
  if (max(left_coords.second, right_coords.second) - min(left_coords.first, right_coords.first) < envelope_distance )
	{
	  return true;
	}
  else 
	{
	  //std::cerr<<"Not in envelope!\n";
	  return false; 
	}
}

bool is_flush(map<node, string> map1)
{
  // Query whether an STP has all its rows the same length.  
  map<node, string>::iterator nodeState;  
  int alignLength = -1; 

  // check to make sure that all alignment rows are the same length
  for (nodeState=map1.begin(); nodeState!=map1.end(); nodeState++)
	{
	  if (alignLength == -1) alignLength = nodeState->second.size(); 
	  else if (alignLength != nodeState->second.size()) return false; 
	}
  return true; 
}
map<node, string> pad_STP(map<node, string> map1, vector<node> toPad)
{
  // Construct a new STP where entries for nodes in toPad have gaps of appropriate length.
  map<node, string> outMap;
  map<node, string>::iterator nodeState;  
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

map<node, string> cat_STP(map<node, string> map1, map<node, string> map2)
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
   
  map<node, string> outMap;
  map<node, string>::iterator nodeState;  
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

	  
	  
	  
  
map<node, string> Profile::merge_STP(M_id m)
{
  // construct a new state_type_phylogeny map based on the input state m
  // This may involve padding states with gaps.  for instance, if one side is delete and the other
  // wait, the wait side will need gaps for each of its nodes of appropriate length, etc
  map<node, string> outMap;
  map<node, string>::iterator nodeState;  
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
  show_alignment(pi, false);
}


void show_state_phylo(map<node, string>& stp)
{
  // visualize a single mapping from nodes to strings (e.g. an alignment).
  map<node, string>::iterator nodeState;
  map<node, string> subAlignment;    

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
  map<node, string>::iterator nodeState;
  map<node, string> subAlignment;    
  pair<vector<int>, vector<int> > statePair;
  map<node, string>::iterator b_iter;   
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



