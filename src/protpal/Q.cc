#include<iostream>
#include<math.h>
#include "Q.h"
#include "transducer.h"
#include "algebras.h"
#include "ecfg/ecfgsexpr.h"

QTransducer::QTransducer(void)
{
}

QTransducer::QTransducer(SingletTrans R_in, BranchTrans B_l_in, BranchTrans B_r_in, SplittingTrans Upsilon_in, Alphabet& alphabet_in) 
{
  // QTransducer's constructor - takes in pre-made singlet, branch, and upsilon transducer objects.  
  // From this, enumerates the states, transitions, and computes weights of transitions. 
  // Also, sets corresponding to key types are stored (e.g. match states, incoming_match_states, etc)


  bool logging = 0;
  bool testing = 0;
  if (logging) std::cerr<<"In QTransducer constructor, about to assign objects...\n";

  // Set a few basic variables
  alphabet =  string(alphabet_in.nondegenerate_chars());
  
  alphabet_size = alphabet.size();
  name = "QTransducer";
  
  // composite state classes
  match =  "match";
  wait =  "wait";
  null =  "null";  
  right_del =  "right_del";
  left_del =  "left_del";
  left_ins = "left_ins";
  right_ins = "right_ins";

  R = R_in; 
  B_l = B_l_in;
  B_r = B_r_in; 
  Upsilon = Upsilon_in; 

  get_state_type_reference();
  if (logging) std::cerr<<"Enumerating states of Q...\n";
  cache_states(); 
  if (logging) std::cerr<<"Building transition 'table' of Q...\n";

  get_transition_table();
  
  if (logging) show_all_states();

  if (logging) std::cerr<<"Enumerating/computing transitions of Q...\n";
  cache_transition_weights();

  if (logging) show_all_transitions();

  if (logging) std::cerr<<"Computing emissions of Q...\n";
  cache_emission_weights();
  
  if (testing) test_emissions();
  if (testing) test_transitions();  
	}



// ***** PUBLIC MEMBER FUNCTIONS *****

vector<state> QTransducer::get_state_indices(state q)
{
  //TODO
}


vector<state> QTransducer::get_wait_states(void)
{
  return state_class_set[wait];
}

vector<state> QTransducer::get_match_states(void)
{
  return state_class_set[match]; 
}
  
vector<state> QTransducer::get_left_emit_states(void)
{
  vector<state> out =  state_class_set[left_ins];
  vector<state> rdel =  state_class_set[right_del];  
  
  for (int i=0; i<rdel.size(); i++) out.push_back(rdel[i]);
  return out; 
}

vector<state> QTransducer::get_SSSI_states(void)
{
  return state_type_set["SSSI"];
}

vector<state> QTransducer::get_right_emit_states(void)
{
  vector<state> out =  state_class_set[right_ins];
  vector<state> ldel =  state_class_set[left_del];  
  
  for (int i=0; i<ldel.size(); i++) out.push_back(ldel[i]);
  return out; 
}


vector<state> QTransducer::get_incoming_match_states(state qPrime)
{
  vector<state> out;
  if (incoming.count(qPrime) <1)
	{
	  std::cerr<<"State " <<qPrime<< " has no incoming transitions!\n";
	  std::cerr<<"The offending call was: get_incoming_match_states, in transducer: "<<name<<" state: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	}
  else 
	{
	  for (int i=0; i<incoming[qPrime].size(); i++)
		{
		  if (state_class[incoming[qPrime][i]] == match) out.push_back(incoming[qPrime][i]);
		}
	}
  return out;
}
vector<state> QTransducer::get_incoming(state qPrime)
{
  if (incoming.count(qPrime) <1)
	{
	  std::cerr<<"State " <<qPrime<< " has no incoming transitions!\n";
	  std::cerr<<"The offending call was: get_incoming_match_states, in transducer: "<<name<<" state: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	}
  else return incoming[qPrime];

}
			

	  
vector<state> QTransducer::get_incoming_left_emit_states(state qPrime)
{
  vector<state> out;
  if (incoming.count(qPrime) <1)
	{
	  std::cerr<<"State " <<qPrime<< " has no incoming transitions!\n";
	  std::cerr<<"The offending call was: get_incoming_left_emit_states, in transducer: "<<name<<" state: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	}
  else 
	{
	  for (int i=0; i<incoming[qPrime].size(); i++)
		{
		  if (state_class[incoming[qPrime][i]] == left_ins || state_class[incoming[qPrime][i]] == right_del ) out.push_back(incoming[qPrime][i]);
		}
	}
  return out;
}
vector<state> QTransducer::get_incoming_right_emit_states(state qPrime)
{
  vector<state> out;
  if (incoming.count(qPrime) <1)
	{
	  std::cerr<<"State " <<qPrime<< " has no incoming transitions!\n";
	  std::cerr<<"The offending call was:get_incoming_right_emit_states in transducer:  "<<name<<" state: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	}
  else 
	{
	  for (int i=0; i<incoming[qPrime].size(); i++)
		{
		  if (state_class[incoming[qPrime][i]] == left_del || state_class[incoming[qPrime][i]] == right_ins ) out.push_back(incoming[qPrime][i]);
		}
	}
  return out;
}
double QTransducer::get_R_transition_weight(state q, state qPrime)
{
  return R.get_transition_weight(assignments[q][0], assignments[qPrime][0]);
}

double QTransducer::get_transition_weight(state q, state qPrime)
{
  // check state existence, connectivity, then lookup in transition_weight map
  if (q >= num_states || qPrime >= num_states)
	{
	  std::cerr<<"One of the states " <<q<<" "<<qPrime<< " is not valid.\n";
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<q<<" named: "<<get_state_name(q)<<" Destination: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	  exit(1);
	}
  else if(outgoing.count(q)<1) 
	{
	  std::cerr<<"The state "<<q<< " has no outgoing transitions..\n";
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<q<<" named: "<<get_state_name(q)<<" Destination: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	  exit(1);
	}
  else if (index(qPrime, outgoing[q]) == -1)
	{
	  std::cerr<<"The states "<<q<<" "<<qPrime<< " do not have a valid transition between them.\n";
	  std::cerr<<"The offending call was: get_transition_weight in transducer: "<<name<<" source state: "<<q<<" named: "<<get_state_name(q)<<" Destination: "<<qPrime<<" named: "<<get_state_name(qPrime)<<endl;
	  exit(1);
	}
  else 
	{
	  vector<state> transitionPair; 
	  transitionPair.push_back(q);transitionPair.push_back(qPrime);
	  return transition_weight[transitionPair];
	}
  
}

bool QTransducer::has_transition(state q, state qPrime)
{
  if (index(q, incoming[qPrime]) != -1) return 1;
  else return 0;
}

string QTransducer::get_state_class(state q)
{
  return state_class[q];
}

string QTransducer::get_state_type(state q)
{
  return state_type[q];
}

double QTransducer::left_delete_weight(state q, int absorb)
{
  // set at 1 for now...
  return 1; 
}

double QTransducer::right_delete_weight(state q, int absorb)
{
  // set at 1 for now...
  return 1; 
}


double QTransducer::left_match_weight(state q, int absorb, int emit)
{
  state branchState = assignments[q][2]; //
  return B_l.get_match_weight(branchState, absorb, emit);
}

double QTransducer::right_match_weight(state q, int absorb, int emit)
{
  state branchState = assignments[q][3]; //
  return B_r.get_match_weight(branchState, absorb, emit);
}

double QTransducer::get_emission_weight(state q, int left_emit, int right_emit)
{
  // Check state existence, character validity, then lookup in emission_weight matrix
  vector<int> emission_tuple;
  emission_tuple.push_back(q);  emission_tuple.push_back(left_emit);  emission_tuple.push_back(right_emit);
  if (q >= num_states)
	{
	  std::cerr<<"The state " <<q<<" "<<q<< " is not valid.\n";
	  std::cerr<<"The offending call was: get_emission_weight in transducer: "<<name<<" state: "<<q<<" named: "<<get_state_name(q)<<endl;
	  exit(1);
	} 
  if (emission_weight.count(emission_tuple)<1)
	{
	  std::cerr<<"One of the character indices " <<left_emit<<" "<<right_emit<< " is not valid.\n";
	  std::cerr<<"The offending call was: get_emission_weight in transducer: "<<name<<" state: "<<q<<" named: "<<get_state_name(q)<<" left-symbol index: "<<left_emit<<" right-symbol index "<<right_emit<<endl;
	  exit(1);
	}
  else
	{  

	  return emission_weight[emission_tuple];
	}
}


// ***** TESTING FUNCTIONS *****
// Make sure things are summing as they should, etc...

void QTransducer::test_emissions(void)
{
  std::cerr<<"\n****** Testing emission distributions of all states ******\n";
  double sum; 
  int left_emit;
  int right_emit;
  for (state q=0 ; q < num_states; q++)
	{
	  if (state_class[q] == match)
		{
		  sum= 0;
		  for (left_emit =0; left_emit < alphabet_size; left_emit++)
			{
			  for (right_emit =0; right_emit < alphabet_size; right_emit++)
				{			  
				  sum += get_emission_weight(q, left_emit, right_emit);
					}
			}
		  std::cerr<<"Sum over all emissions from state "<<get_state_name(q)<<" (class: "<<state_class[q]<<") :\t " <<sum<<endl;
		}
	  if (state_class[q] == left_ins)
		{
		  sum = 0;
		  for (left_emit =0; left_emit < alphabet_size; left_emit++)
			{
			  sum += get_emission_weight(q, left_emit, -1);
			}
		  std::cerr<<"Sum over all emissions from state "<<get_state_name(q)<<" (class: "<<state_class[q]<<") :\t " <<sum<<endl;		
		}
	  
	  if (state_class[q] == right_ins)
		{
		  sum = 0;
		  for (right_emit =0; right_emit < alphabet_size; right_emit++)
			{
			  sum += get_emission_weight(q, -1, right_emit);
			}
		  std::cerr<<"Sum over all emissions from state "<<get_state_name(q)<<" (class: "<<state_class[q]<<") :\t " <<sum<<endl;		
		}



	  if (state_class[q] == right_del)
		{
		  sum = 0;
		  for (left_emit =0; left_emit < alphabet_size; left_emit++)
			{
			  sum += get_emission_weight(q, left_emit, -1);
			}
		  std::cerr<<"Sum over all emissions from state "<<get_state_name(q)<<" (class: "<<state_class[q]<<") :\t " <<sum<<endl;		
		}
	  if (state_class[q] == left_del)
		{
		  sum = 0;
		  for (right_emit =0; right_emit < alphabet_size; right_emit++)
			{
			  sum += get_emission_weight(q, -1, right_emit);
			}
		  std::cerr<<"Sum over all emissions from state "<<get_state_name(q)<<" (class: "<<state_class[q]<<") :\t " <<sum<<endl;		
	   }
	  
	}
}

void QTransducer::test_transitions(void)
{
  std::cerr<<"\n****** Testing transitions of all states ******\n";
  double sum=0;
  for (state q=0; q<num_states; q++)
	{
	  sum=0;
	  if (outgoing.count(q) <1) std::cerr<<"State "<< get_state_name(q)<<" has no outgoing transitions\n";
	  else
		{	
		  std::cerr<<"Examining state "<<get_state_name(q)<<endl;
		  //		  std::cerr<<"\tOutgoing transitions: ";
		  //displayVector(outgoing[q]);
		  
		  for (vector<int>::iterator qPrime=outgoing[q].begin();qPrime!=outgoing[q].end();qPrime++)
			{
			  sum+= get_transition_weight(q, *qPrime);
			}
		  std::cerr<<"\tSum of outgoing transitions : "<<sum<<endl;
		}
	}
}
			   


// ***** MISC. *****


// Not using this anymore...branch transducers are not set up to re-set their parameters.  
void QTransducer::change_branch_lengths(float left_branch_length, float right_branch_length)
{
  B_l.branch_length = left_branch_length;
  B_r.branch_length = right_branch_length;  

  //transitions weights are (typically) branch-length dependent (e.g. ins_open_rate, etc), so we must re-compute them:
  cache_transition_weights();
  //emission weights are (typically) branch-length dependent, so we must re-compute them:
  cache_emission_weights();
}

/*
void QTransducer::change_branch_lengths(float left_branch_length, float right_branch_length)
{
  B_l.branch_length = left_branch_length;
  B_r.branch_length = right_branch_length;  

  B_l.conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(left_branch_length);
  B_r.conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(right_branch_length);   

  cache_transition_weights();
  cache_emission_weights(); 
}
*/


// ***** PRINTING/SHOWING FUNCTIONS *****

void QTransducer::show_all_transitions(void)
{
  for (state q=0; q<num_states; q++)
	{
	  std::cout<<"\n\nSource state: "<<q<<"  "<<get_state_name(q)<<endl;
	  if (q == composite_end_state) continue;
	  if (outgoing.count(q) < 1)
		{
		  std::cout<<"The state "<<get_state_name(q)<<" has no outgoing transitions.  Possibly problematic...\n";
		  continue;
		}
	  for (vector<int>::iterator qPrime = outgoing[q].begin(); qPrime != outgoing[q].end(); qPrime++)
		{
		  std::cout<<"  There is a transition to state: " <<get_state_name(*qPrime)<<", weight: "<<-log(get_transition_weight(q,*qPrime))/log(2)<<" bits"<<endl;
		}
	}
}

void QTransducer::show_all_states(void)
{
  std::cout<<"There are this many states in total: "<<num_states<<endl;
  for (int i=0; i<num_states; i++)
	{
	  std::cout<<"State "<<i<< " is : "<<get_state_name(i)<<", class:"<<get_state_class(i)<<endl;
	}	
}


string QTransducer::get_state_name(state q)
{
  string out; 
  vector<state> stateIndices = assignments[q]; 
  out = R.get_state_name(stateIndices[0])+" ";
  out += Upsilon.get_state_name(stateIndices[1])+" ";
  out += B_l.get_state_name(stateIndices[2]) + " ";
  out += B_r.get_state_name(stateIndices[3])+" "; 
  return out; 
}


// ****** CACHING FUNCTIONS FOR TRANSITIONS, EMISSIONS, STATES, ETC ******
  
// enumerates, computes, and stores transitions of Q
void QTransducer::cache_transition_weights(void)
{
  // Go through transition_table, record incoming, outgoing connections for each state.  
  // Each key in the table is a vector of  two strings, each representing the state types of the
  // component transducers.  
  // Here, we loop through the entries in the table, then call up all the states that have that type-string,
  // which is stored in the map state_type_set 

  map<vector<string>, vector<int> >::iterator it;
  vector<state>::iterator i;   vector<state>::iterator j; 
  string to; string from; vector<int> componentChanges;
  vector<state> transitionPair;
  bfloat weight;
  for (it=transition_table.begin(); it!=transition_table.end(); it++)
	{
	  from = (*it).first[0]; 
	  to = (*it).first[1]; 
	  componentChanges = (*it).second; 
	  for (i=state_type_set[from].begin(); i!=state_type_set[from].end(); i++)
		{
		  for (j=state_type_set[to].begin(); j!=state_type_set[to].end(); j++)
			{
			  transitionPair.clear();
			  transitionPair.push_back(*i); 
			  transitionPair.push_back(*j); 
			  // std::cerr<<"Computing weight between " << *i << " and " << *j<<endl;
			   weight = compute_transition_weight(*i, *j, componentChanges);
			   if (weight > 0.0)
				 {
				   transition_weight[transitionPair] = weight; 
				   incoming[*j].push_back(*i);outgoing[*i].push_back(*j);
				 }
 			}
		}
	}
}

double QTransducer::compute_transition_weight(state q, state qPrime, vector<int> componentChanges)
{
  // Using the componentChanges vector of ints, sum over intermediate states in the proper way.    
  bfloat weight = 1;
  bfloat sum;
  state waitIndex;
  vector<state> q_state_assignments = assignments[q];
  vector<state> qPrime_state_assignments = assignments[qPrime];
  vector<state> wait_indices;
  
  // Weight of R's transitions, possibly summing over wait states
  wait_indices = R.get_state_type_set("W");
  if (componentChanges[0] == 0 && q_state_assignments[0] != qPrime_state_assignments[0]) return 0; 
  if (componentChanges[0] == 1) weight *= R.get_transition_weight(q_state_assignments[0],qPrime_state_assignments[0]);
  else if(componentChanges[0]==2) // the intermediate is always *a* wait state, but there might be several of these!
	{
	  sum = 0;
	  for (int i=0; i<wait_indices.size(); i++)
		{
		  waitIndex = wait_indices[i];
		  if (R.has_transition(q_state_assignments[0], waitIndex) && 
			  R.has_transition(waitIndex, qPrime_state_assignments[0]))
			{
			  sum +=													\
				R.get_transition_weight(q_state_assignments[0], waitIndex)*
				R.get_transition_weight(waitIndex, qPrime_state_assignments[0]);
			}
		}
	  weight *= sum;
	}
  if (weight > 1.0) std::cerr<<"Error: it seems that weight has become larger than 1, after adding R's transitions\n";
  // Weight of U's transitions are always 1, no need to include them.

  if (componentChanges[2] == 0 && q_state_assignments[2] != qPrime_state_assignments[2]) return 0; 
  // Weight of B_l's transitions, possibly summing over wait states
  if (componentChanges[2] == 1) 
	{
	  if (! B_l.has_transition(q_state_assignments[2],qPrime_state_assignments[2]) ) return 0; 
	  else weight *= B_l.get_transition_weight(q_state_assignments[2],qPrime_state_assignments[2]);
	}

  else if(componentChanges[2]==2)
	{
	  wait_indices = B_l.get_state_type_set("W");
	  sum = 0;
	  for (int i=0; i<wait_indices.size(); i++)
		{
		  waitIndex = wait_indices[i];
		  if (B_l.has_transition(q_state_assignments[2], waitIndex) && 
			  B_l.has_transition( waitIndex, qPrime_state_assignments[2]))
			{
			  sum +=
				B_l.get_transition_weight(q_state_assignments[2], waitIndex)*
				B_l.get_transition_weight(waitIndex, qPrime_state_assignments[2]);
			}
		}
	  weight *= sum;
	}
  if (weight > 1.0) std::cerr<<"Error: it seems that weight has become larger than 1, after adding Bl's transitions\n";

  if (componentChanges[3] == 0 && q_state_assignments[3] != qPrime_state_assignments[3]) return 0; 
  // Weight of B_r's transitions, possibly summing over wait states
  if (componentChanges[3] == 1) 
	{
	  if (! B_r.has_transition(q_state_assignments[3],qPrime_state_assignments[3]) ) return 0; 
	  else weight *= B_r.get_transition_weight(q_state_assignments[3],qPrime_state_assignments[3]);
	}
  else if(componentChanges[3]==2)
	{
	  wait_indices = B_r.get_state_type_set("W");
	  sum = 0;
	  for (int i=0; i<wait_indices.size(); i++)
		{
		  waitIndex = wait_indices[i];
		  if (B_r.has_transition(q_state_assignments[3], waitIndex) && 
			  B_r.has_transition( waitIndex, qPrime_state_assignments[3]))
			{
			  sum +=
				B_r.get_transition_weight(q_state_assignments[3], waitIndex)*
				B_r.get_transition_weight(waitIndex, qPrime_state_assignments[3]);
			}
		}
	  weight *= sum;
	}
  if (weight > 1.0) std::cerr<<"Error: it seems that weight has become larger than 1, after adding Br's transitions\n";
  return weight;

}


//  enumerates the composite states of Q
void QTransducer::cache_states(void)
{
  bool logging = 0;

  int state_count = 0;
  map<string, vector<string> >::iterator it;
  vector<string> type_strings; 
  vector<state> state_indices;
  string compositeClass, type_R, type_B_l, type_B_r, type_Upsilon;
  
  for ( it = state_type_reference.begin(); it!= state_type_reference.end(); it++)
	{
	  compositeClass = (*it).first; 
	  type_strings = (*it).second;

	  if(logging) std::cerr<<"The composite class is: " << compositeClass<<endl;	  
	  //	  state_set[compositeType] = vector<
	  for (int i=0; i<type_strings.size(); i++)
		{
		  type_R = stringAt(type_strings[i], 0);
		  type_Upsilon = stringAt(type_strings[i], 1); 
		  type_B_l = stringAt(type_strings[i], 2); 
		  type_B_r = stringAt(type_strings[i], 3);
		  if (logging)
			{
			  std::cerr<<"\tThe individual types are (r,u,bl, br): " << type_R<<" "<<type_Upsilon<<" "<<type_B_l<<" "<<type_B_r<<" \n";
			}
		  // Loop over states of R
		  for (int r=0; r<R.states.size(); r++)
			{
			  if (R.get_state_type(r) != type_R ) continue;

			  // Loop over states of Upsilon
			  for (int u=0; u<Upsilon.states.size(); u++)
				{
				  if (Upsilon.get_state_type(u) != type_Upsilon ) continue;

				  // Loop over states of B_l
				  for (int bl=0; bl<B_l.states.size(); bl++)
					{
					  if (B_l.get_state_type(bl) != type_B_l ) continue;

					  // Loop over states of B_r
					  for (int br=0; br<B_r.states.size(); br++)
						{
						  if (B_r.get_state_type(br) != type_B_r ) continue;
						  
						  else 
							{
							  if (logging)
								{
								  std::cerr<<"Adding the following state as number "<<state_count<<endl;
								  std::cerr<<"\t"<<R.get_state_name(r)<<" "<<Upsilon.get_state_name(u)<<" "<<B_l.get_state_name(bl)<<" "<<B_r.get_state_name(br)<<endl;
								}
						
							  // If it is the start or end composite state, note this 
							  string end = "end";
							  if (compositeClass == end ) 
								{
								  composite_end_state = state_count; 
								  if (logging) std::cerr<<"Recording the end state as state: "<< composite_end_state<<endl;
								}
							  string start = "start";
							  if (compositeClass == start)
								{
								  composite_start_state = state_count; 
								  if (logging) std::cerr<<"Recording the start state as state: "<< composite_start_state<<endl;
								}

							  
							  //  add to list state_class_set
							  state_class_set[compositeClass].push_back(state_count);
							  state_type_set[type_strings[i]].push_back(state_count);

							  //  add to map state_class
							  state_class.push_back(compositeClass);
							  //  add to map state_type
							  state_type.push_back(type_strings[i]);							  

							  //  add to the matrix  assignments.  
							  //  assignments[m][s] tells what state the single transducer (R,U,B_l,B_r)[s] is in for 
							  //  composite state number m. 

							  state_indices.clear();
							  state_indices.push_back(r);
							  state_indices.push_back(u);
							  state_indices.push_back(bl);
							  state_indices.push_back(br);											
							  assignments.push_back(state_indices);

							  // increment Q's state counter
							  state_count++;
							}
						}
					}
				}
			}
		}
	}
  num_states = state_count;
}


void QTransducer::cache_emission_weights( void )
{
  double weight;
  vector<int> state_l_r; 
  bool logging = 0;
  //  store emission weight of state q, omega_l, omega_r
  // This will get stored in the map emission_weight, which maps a vector state, left-emit, right-emit to a double
  for (state q = 0; q< num_states; q++)
	{
	  if (logging)std::cerr<<"Caching emissions for: "<<q<<endl;
	  if (logging)std::cerr<<"\tComponents: "<<get_state_name(q)<<endl;
	  int state_R, state_B_l,state_B_r; 
	  int omega_left, omega_right, omega_R;

	  // Assignments vector holds [R, Upsilon, B_l, B_r]
	  state_R = assignments[q][0];
	  state_B_l = assignments[q][2]; 
	  state_B_r = assignments[q][3]; 

	  string qClass = state_class[q]; 

	  if (logging)std::cerr<<"\tClass: "<<qClass<<endl;
     
	  // q is a match state
	  if (qClass == match)
		{
		  for (omega_left = 0; omega_left<alphabet_size; omega_left++)
			{
			  for (omega_right = 0; omega_right<alphabet_size; omega_right++)
				{
				  state_l_r.clear(); 
				  state_l_r.push_back(q); state_l_r.push_back(omega_left); state_l_r.push_back(omega_right); 
				  weight = 0; 
				  // omega_R is the character at R
				  for (omega_R=0; omega_R<alphabet_size; omega_R++)
					{
					  weight += R.get_emission_weight(state_R, omega_R)*\
						B_l.get_match_weight(state_B_l, omega_R, omega_left)*\
						B_r.get_match_weight(state_B_r, omega_R, omega_right);
					}
				  emission_weight[state_l_r] = weight;
				  if(logging) 
					{
					  std::cerr<<"Caching emission weight for "<<get_state_name(q)<<" chars: ";
					  std::cerr<<alphabet[omega_right]<<" "<<alphabet[omega_left]<<endl;
					  std::cerr<<" : "<<weight<<endl;
					} 
				  
				}
			}
		}
			  

	  //left-del - omega_left is -1
	  else if (qClass == left_del)
		{
		  omega_left = -1; 
		  for (omega_right = 0; omega_right<alphabet_size; omega_right++)
			{
			  state_l_r.clear(); 
			  state_l_r.push_back(q); state_l_r.push_back(omega_left); state_l_r.push_back(omega_right); 
			  weight = 0; 
			  // omega_R is the character at R
			  for (omega_R=0; omega_R<alphabet_size; omega_R++)
				{
				  weight += R.get_emission_weight(state_R, omega_R)*	\
					B_r.get_match_weight(state_B_r, omega_R, omega_right);
					}
			  emission_weight[state_l_r] = weight;
			}
		}
	  //right-del - omega_right is -1
	  else if (qClass == right_del)
		{
		  omega_right = -1; 
		  for (omega_left = 0; omega_left<alphabet_size; omega_left++)
			{
			  state_l_r.clear(); 
			  state_l_r.push_back(q); state_l_r.push_back(omega_left); state_l_r.push_back(omega_right); 
			  weight = 0; 
			  // omega_R is the character at R
			  for (omega_R=0; omega_R<alphabet_size; omega_R++)
				{
				  weight += R.get_emission_weight(state_R, omega_R)*	\
					B_l.get_match_weight(state_B_l, omega_R, omega_left);
					}
			  emission_weight[state_l_r] = weight;
			}
		}
	  
	  //left-ins - omega_right is -1
	  else if (qClass == left_ins)
		{
		  omega_right = -1;
		  for (omega_left = 0; omega_left<alphabet_size; omega_left++)
			{
			  state_l_r.clear(); 
			  state_l_r.push_back(q); state_l_r.push_back(omega_left); state_l_r.push_back(omega_right); 

			  weight = B_l.get_emission_weight(state_B_l, omega_left);
			  emission_weight[state_l_r] = weight;
			  
			}
		}
		  
	  //right-ins - omega_left is -1
	  else if (qClass == right_ins)
		{
		  omega_left = -1;
		  for (omega_right = 0; omega_right<alphabet_size; omega_right++)
			{
			  state_l_r.clear(); 
			  state_l_r.push_back(q); state_l_r.push_back(omega_left); state_l_r.push_back(omega_right); 
			  
			  weight = B_r.get_emission_weight(state_B_r, omega_right);
			  emission_weight[state_l_r] = weight;
			}
		}
	  
	  
	}
}


void QTransducer::marginalizeNullStates(void)
{
  // Marginalize over null states in Q, those having type (I,M,D,D)
  // NB this is for the simplified case where there is only one such state...
  // This means that in all cases we're eliminating null 'B' state from a set like this:
  // A -> B -> C where B is possibly self-looping 
  bool logging = 0;
  state a;
  vector<state>::iterator b,c;
  vector<state> b_states, c_states;
  bool selfLoop=0; 
  vector<state> transitionPair;
  
  for (a=0; a<num_states; a++)
	{
	  if (outgoing.count(a) <1 || state_class[a] == "null") continue;
	  b_states = outgoing[a];
	  for (b=b_states.begin(); b!=b_states.end(); b++)
		{
		  if (state_class[*b] == "null")
			{
			  if (has_transition(*b,*b)) selfLoop = 1;
			  else selfLoop = 0; 

			  c_states = outgoing[*b];
			  for (c=c_states.begin(); c!=c_states.end(); c++)
				{
				  if (*c == *b) continue;
				  else if (state_class[*c] == "null")
					{
					  std::cerr<<"Error: multiple null states in Q.  Can't marginalize with this simple method\n";
					  exit(1);
					}
				  else
					{
					  if(logging) std::cerr<<"connecting states:"<<get_state_name(a)<<" and "<<get_state_name(*c)<<endl;

					  if (!has_transition(a,*c))
						{
						  incoming[*c].push_back(a); outgoing[a].push_back(*c);
						}
						
					  transitionPair.clear();
					  transitionPair.push_back(a); transitionPair.push_back(*c); 
					  //sum their transition until the difference is very small
					  double toAdd = 100; 
					  int counter = 1;
					  bfloat ab, bc, bb; 
					  ab = get_transition_weight(a, *b); bc = get_transition_weight( *b, *c); 
					  if (!selfLoop) transition_weight[transitionPair] = ab*bc;
					  else
						{
						  bb=get_transition_weight(*b,*b);
						  while (toAdd> .0001)
							{
							  toAdd = ab*pow(bb, counter)*bc;
							  transition_weight[transitionPair] += toAdd;
							  counter += 1;
							}
						}
					}
				}
			}
		}
	}
}
			  
	  

  
// The following functions initialize hard-coded aspects of the Qtransducer.  Most of this was auto-generated
// via emacs or python code.  Probably best to check it over via other methods rather than looking at this directly



void QTransducer::get_state_type_reference(void)
{
  vector<string> types;
  //start
  types.push_back("SSSS");  
  state_type_reference["start"] = types; types.clear(); 


  //match
  types.push_back("IMMM");
  state_type_reference["match"] = types; types.clear(); 
  
  
  //right del
  types.push_back("IMMD");
  state_type_reference["right_del"] = types; types.clear(); 

  //left del
  types.push_back("IMDM");
  state_type_reference["left_del"] = types; types.clear(); 

  //right ins
  types.push_back("SSSI");    
  types.push_back("IMMI");
  types.push_back("IMDI");
  state_type_reference["right_ins"] = types; types.clear(); 
  
  //left ins
  types.push_back("SSIW");
  types.push_back("IMIW");
  state_type_reference["left_ins"] = types; types.clear(); 

  //null
  types.push_back("IMDD");
  state_type_reference["null"] = types; types.clear(); 


  //wait
  types.push_back("WWWW");  
  state_type_reference["wait"] = types; types.clear();
  
  //end
  types.push_back("EEEE");  
  state_type_reference["end"] = types; types.clear(); 

}
 
void QTransducer::get_transition_table(void)
{
  string from; string to; vector<int> changes; vector<string> states;

  // NB the following code was generated with cat table.txt | python table2cc.py 
  // use e.g. insert-file Qtransitions.txt here
  from="SSSS"; to="SSSI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSSS"; to="SSIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSSS"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSS"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSS"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSS"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSS"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSSI"; to="SSSI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSSI"; to="SSIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSSI"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSI"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSI"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSI"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="SSSI"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSIW"; to="SSIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(0);
  transition_table[states]=changes;


  from="SSIW"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSIW"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSIW"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSIW"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="SSIW"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(0);
  transition_table[states]=changes;


  from="IMMI"; to="IMMI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMI"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMI"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMI"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMI"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMI"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMI"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMD"; to="IMMI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMD"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMD"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMD"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMD"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMD"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMD"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMM"; to="IMMI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMM"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMMM"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMM"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMM"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMM"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMMM"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDI"; to="IMDI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDI"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDI"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDI"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDI"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDI"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDI"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDD"; to="IMDI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDD"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDD"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDD"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDD"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDD"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDD"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDM"; to="IMDI"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDM"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMDM"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDM"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDM"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDM"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(2);
  transition_table[states]=changes;


  from="IMDM"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMIW"; to="IMIW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(0);
  changes.push_back(0);
  changes.push_back(1);
  changes.push_back(0);
  transition_table[states]=changes;


  from="IMIW"; to="IMMM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMIW"; to="IMMD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMIW"; to="IMDM"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMIW"; to="IMDD"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(2);
  changes.push_back(2);
  changes.push_back(1);
  transition_table[states]=changes;


  from="IMIW"; to="WWWW"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(0);
  transition_table[states]=changes;


  from="WWWW"; to="EEEE"; states.clear();
  states.push_back(from); states.push_back(to);
  changes.clear();
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  changes.push_back(1);
  transition_table[states]=changes;



}
