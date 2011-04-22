#include<iostream>

#include "protpal/utils.h"
#include "protpal/ReadProfileScore.h"

void Read::pad(void)
{
  list<sstring>::iterator it = this->begin();
  // this->insert(it, "start"); 
  this->push_back("pre-end"); 
}

void Read::set(sstring in)
{
  this->clear(); 
  for (sstring::iterator strIter=in.begin(); strIter!=in.end(); strIter++)
    this->push_back(*strIter); 
  pad(); 
}

ReadProfileScore::ReadProfileScore(AbsorbingTransducer* prof_in)
{
  profile = prof_in;
  // Collect the various types of profile states
  profile_states.clear();
  for (int i=0; i<profile->num_delete_states; i++)
    profile_states.push_back(i); 
  profile_states.push_back(profile->pre_end_state); 
  DP_ID.resize(3); 
}

void ReadProfileScore::score_and_print(const Read& read, ostream& out, bool viterbi)
{
  bfloat value;
  readSize = read.size(); 
  if ( viterbi )
    {
      pointer_ID.resize(3); 
      fill_DP_matrix(read, true, true); 
      value = get_viterbi_value(); 
    }
  else
    {
      fill_DP_matrix(read, false, false); 
      value = get_forward_value(); 
    }
  write_read_info(out, read, value); 
  clear_DP_matrix(); 
}
void ReadProfileScore::clear_DP_matrix(void)
{
  map< vector<int>, bfloat> DP_tmp; 
  DP_tmp.swap(DP_matrix); 
}

void ReadProfileScore::write_read_info(ostream& out, const Read& read, bfloat value)
{
  out << read.identifier << "\t" << value << "\n";
}

inline list<state> ReadProfileScore::get_possible_HMM_states(int i, state j)
{
  list<state> toReturn; 
  if (i == readSize-1) // pre-end "state"
    {
      if (j == profile->pre_end_state)
	{
	  toReturn.push_back(pairHMM.end_state); 
	  return toReturn;
	}
      else
	{
	  toReturn.push_back(pairHMM.state_index["post_read_ins"]); 
	  return toReturn;
	}
    }
  else
    {
      return pairHMM.non_special_states; // not start, end
    }
}

void ReadProfileScore::fill_DP_matrix(const Read& read, bool backPointers, bool logging)
{
  int i, iPrime; 
  list<int>::iterator kPrime,k; 
  vector<state>::iterator j, jPrime; 
  list<state> incoming_HMM_states, possible_HMM_states;
  vector<state> incoming_profile_states; 
  bfloat toAdd;
  if (DP_matrix.size() != 0)
    {
      cerr<<"Filling DP matrix that is not empty! (This should not happen)\n";
      exit(1);
    }
  if (logging)
    cerr<<"Beginning to fill DP matrix.\nRead size: " << readSize << ".  \nProfile size: " << profile_states.size() << "\nHMM size: " << pairHMM.num_states<<endl; 

  if (logging)
    cerr<<"Seting initial state of DP matrix\n"; 
  set_DP_cell(-1,profile->start_state, pairHMM.start_state, 1.0); 
  
  // Range is 35 to ~500
  for (i=0; i<readSize; ++i) 
    {
      if (logging)
	cerr<<"\n\nProcessing position " << i  <<  " of read...\n"; 
      // Range is 300 - 1000s
      for (j=profile_states.begin(); j<profile_states.end(); j++)
	{
	  if (logging)
	    cerr<<"\nProcessing state " << *j  <<  " of profile...\n"; 
	  // The HMM states compatible with a given (read, profile) state pair is determined:
	  possible_HMM_states = get_possible_HMM_states(i, *j); 
	  if (logging)
	    {
	      cerr<<"The possible HMM states are: "; 
	      for (k=possible_HMM_states.begin(); k!=possible_HMM_states.end(); k++)
		cerr<< pairHMM.state_name[*k] << " "; 
	      cerr<<"\n"; 
	    }
	  // Range is 5-10
	  for (k=possible_HMM_states.begin(); k!=possible_HMM_states.end(); k++)
	    {
	      if (logging)
		cerr<<"Processing state " << pairHMM.state_name[*k]  <<  " of HMM...\n"; 
	      incoming_HMM_states = pairHMM.incoming[*k]; 
	      if (logging)
		cerr<<"Summing over " << incoming_HMM_states.size() << " incoming HMM states\n"; 
	      for (kPrime = incoming_HMM_states.begin(); kPrime != incoming_HMM_states.end(); ++kPrime)
		{
		  iPrime = get_incoming_read_state(i, *kPrime); 
		  // We can be sure this DP cell will be -inf, so we continue
		  if (iPrime == -1 && pairHMM.state_type[*kPrime] != "start")
		    continue;
		  hmm_transition_weight = pairHMM.get_transition_weight(*kPrime, *k); 
		  incoming_profile_states = get_incoming_profile_states(*j, *k); 
		  for (jPrime = incoming_profile_states.begin(); jPrime != incoming_profile_states.end(); ++jPrime)
		    {
		      // We can be sure these DP cells will be -inf, so we continue
		      if (*jPrime == profile->start_state || iPrime == -1)
			{
			  if (pairHMM.state_type[*kPrime] != "start")
			    continue;
			  if (*jPrime != profile->start_state || iPrime != -1)
			    continue; 
			}
		      toAdd = get_DP_cell(iPrime, *jPrime, *kPrime); 
		      if (! toAdd > 0.0)
			continue;
		      if (*j != *jPrime)
			toAdd *= profile->get_transition_weight(*jPrime, *j); 
		      if (pairHMM.state_type[*k] != "start" && pairHMM.state_type[*k] != "end")
			toAdd *= pairHMM.get_emission_weight(i,*j,*k); 
		      toAdd *= hmm_transition_weight; 
		      // Set the value
		      add_to_DP_cell(i,*j,*k, toAdd); 
		      if ( backPointers ) // backpointers might be requested for viterbi traceback
			add_backPointer(i,*j,*k, iPrime, *jPrime, *kPrime); 
		      if (logging)
			{
			  cerr<<"Added to DP matrix:: \n";
			  cerr<<"\tHMM state " << pairHMM.state_name[*k] << endl; 
			  cerr << "\tRead position " << i <<"\n\tProfile state " << *j << "\n\tValue " << toAdd << "\n\n"; 
			}
		    }
		}
	    }
	}
    }
}

inline void ReadProfileScore::add_backPointer(int toRead, state toProfile , state toHMM, int fromRead , state fromProfile , state fromHMM )
{
  DP_ID[0] = toRead; 
  DP_ID[1] = toProfile; 
  DP_ID[2] = toHMM; 

  pointer_ID[0] = fromRead; 
  pointer_ID[1] = fromProfile; 
  pointer_ID[2] = fromHMM; 

  backPointers[DP_ID] = pointer_ID; 
}



inline void ReadProfileScore::add_to_DP_cell(int i, state j, state k, bfloat toAdd)
{
  DP_ID[0] = i; 
  DP_ID[1] = j; 
  DP_ID[2] = k; 
  
  if (DP_matrix.count(DP_ID))
    DP_matrix[DP_ID]  = DP_matrix[DP_ID] + toAdd; 
  else
    DP_matrix[DP_ID] = toAdd; 
}

inline void ReadProfileScore::set_DP_cell(int i, state j, state k, bfloat toAdd)
{
  DP_ID[0] = i; 
  DP_ID[1] = j; 
  DP_ID[2] = k; 
  
  if (DP_matrix.count(DP_ID))
    {
      cerr<<"ERROR: Attempting to 'set' DP cell that has already a value!\n";
      cerr<<"The states were: " << i  << " " << j << " " << pairHMM.state_name[k] << endl; 
      exit(1); 
    }

  else
    DP_matrix[DP_ID] = toAdd; 
}  

inline bfloat ReadProfileScore::get_DP_cell(int i, state j, state k)
{
  DP_ID[0] = i; 
  DP_ID[1] = j; 
  DP_ID[2] = k; 
  
  if (!DP_matrix.count(DP_ID))
    {
      // cerr<<"Warning: Attempting to 'get'  a DP cell that has no value!\n";
      // cerr<<"The states were: " << i  << " " << j << " " << pairHMM.state_name[k] << endl; 
      return 0.0; 
    }
  else
    return DP_matrix[DP_ID];
}  



inline bfloat ReadProfileScore::get_forward_value(void)
{
  return get_DP_cell(readSize-1, profile->pre_end_state, pairHMM.end_state); 
}      

bfloat ReadProfileScore::get_viterbi_value(void)
{
  cerr<<"Warning: Viterbi score not yet implemented!\n";
  return 1.0; 
}      
  
inline int ReadProfileScore::get_incoming_read_state(int readIndex, state hmm_state)
{
  if (pairHMM.state_type[hmm_state] == "insert")
    return readIndex; 
  else
    return readIndex-1; 
}

inline vector<state> ReadProfileScore::get_incoming_profile_states(state profileState, state hmm_state)
{
  vector<state> toReturn; 
  if (pairHMM.state_type[hmm_state] == "delete")
    {
      toReturn.push_back(profileState); 
      return toReturn; 
    }
  else
    return profile->get_incoming(profileState); 
}



// ReadProfileModel methods.  Basic pairHMM to 'align' a read to a profile
ReadProfileModel::ReadProfileModel(void)
{
  num_states=-1; 
  // A basic model - here for now as a placeholder
  add_state("start","start");
  add_state("insert","pre_read_ins");
  add_state("insert","post_read_ins");
  add_state("insert","read_ins");
  add_state("delete","delete");
  add_state("match","match");
  add_state("end","end");

  // from start
  add_transition("start", "pre_read_ins", .9);
  add_transition("start", "match", .05);
  add_transition("start", "delete", .05);

  // from pre_read_ins
  add_transition("pre_read_ins", "pre_read_ins", .99); 
  add_transition("pre_read_ins", "match", .005); 
  add_transition("pre_read_ins", "delete", .005); 

  // from match
  add_transition("match", "match", 0.75); 
  add_transition("match", "read_ins", 0.1); 
  add_transition("match", "delete", 0.1); 
  add_transition("match", "post_read_ins", 0.05); 


  // from read_ins
  add_transition("read_ins", "read_ins", 0.1); 
  add_transition("read_ins", "match", 0.75); 
  add_transition("read_ins", "delete", 0.1); 
  add_transition("read_ins", "post_read_ins", 0.05); 

  // from delete
  add_transition("delete", "read_ins", 0.1); 
  add_transition("delete", "match", 0.75); 
  add_transition("delete", "delete", 0.1); 
  add_transition("delete", "post_read_ins", 0.05); 

  // from post_read_ins
  add_transition("post_read_ins", "post_read_ins", 0.99); 
  add_transition("post_read_ins", "end", 0.01); 
}

void ReadProfileModel::add_transition(string from, string to, bfloat weight)
{
  if (!state_index.count(from) || !state_index.count(to))
    {
      cerr<<"Defining transition between unknown states: " << to << " and " << from << endl;
      exit(1); 
    }
  // Add transition weight info
  pair<state, state> transitionPair; 
  transitionPair.first=state_index[from]; 
  transitionPair.second=state_index[to]; 
  transition[transitionPair] = weight; 
  // Add connectivity information
  incoming[state_index[to]].push_back(state_index[from]); 
}
  
  
void ReadProfileModel::add_state(string type, string name)
{
  // Add a state to a read-prof HMM
  for (vector<string>::iterator stateIter=state_name.begin(); stateIter != state_name.end(); stateIter++)
    {
      if (*stateIter == name)
	{
	  cerr<<"State name " << name  << " already exists in HMM!\n";
	  exit(1); 
	}
    }
  state_name.push_back(name);
  num_states++; 
  state_index[name] = num_states; 
  state_type[num_states] = type; 
  if (type == "start")
    start_state = num_states;
  else 
    {
      if (type == "end")
	end_state = num_states;
      else
	non_special_states.push_back(num_states); 
    }
}


inline bfloat ReadProfileModel::get_transition_weight(state from, state to)
{
  transitionPair.first = from; 
  transitionPair.second = to; 
  
  tmpIter = transition.find(transitionPair);
  if (tmpIter == transition.end())
    {
      cerr<<"Requested an uninitialized transition: " << from <<  " to " << to << endl; 
      exit(1); 
    }
  else
    return tmpIter->second;
}

bfloat ReadProfileModel::get_emission_weight(int readIndex, state profileState, state hmm_state)
{
  // cerr<<"Warning: emission weight is not yet meaningful!\n"; 
  return 1.0; 
}
