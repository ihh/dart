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
  this->sstringRep = in;  
  this->clear(); 
  for (sstring::iterator strIter=in.begin(); strIter!=in.end(); strIter++)
    this->push_back(*strIter); 
  pad(); 
}



ReadProfileScore::ReadProfileScore(AbsorbingTransducer* prof_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix_in)
{
  profile = prof_in;
  if (profile->name == "")
    {
      cerr<<"ERROR! Input profile for matching read has no name, cannot continue\n"; 
      exit(1); 
    }
  else
    name = profile->name; 

  // Collect the various types of profile states
  profile_states.clear();
  for (int i=0; i<profile->num_delete_states; i++)
    profile_states.push_back(i); 
  profile_states.push_back(profile->pre_end_state); 
  DP_ID.resize(3); 
  pairHMM.set_substitution_model(alphabet_in, rate_matrix_in, profile);
}

bfloat ReadProfileScore::get_score(const Read& read, bool viterbi, bool logging)
{
  bfloat value;
  ofstream dummy_ostream; 
  readSize = read.size(); 
  if ( viterbi )
    {
      pointer_ID.resize(3); 
      fill_DP_matrix(read, 
		     dummy_ostream, // hmmoc filestream
		     false, // only write hmmoc file, don't do actual DP
		     true,  // keep backPointers - for finding viterbi traceback
		     logging); // display logging messages
      value = get_viterbi_value(); 
    }
  else
    {
      fill_DP_matrix(read, 
		     dummy_ostream, // hmmoc filestream
		     false, // only write hmmoc file, don't do actual DP
		     false,  // keep backPointers - for finding viterbi traceback
		     logging); // display logging messages
      value = get_forward_value(); 
    }
  clear_DP_matrix(); 
  return value; 
}

void ReadProfileScore::score_and_print(const Read& read, ostream& out, bool viterbi)
{
  bfloat value = get_score(read, viterbi); 
  write_read_info(out, read, value); 
}

void ReadProfileScore::score_and_store(const Read& read, ScoreMap& scores, bool viterbi )
{
  bfloat value = get_score(read, viterbi); 
  scores[read.identifier][name] = value;
}


void ReadProfileScore::clear_DP_matrix(void)
{
  map< vector<int>, bfloat> DP_tmp; 
  DP_tmp.swap(DP_matrix); 
}

void ReadProfileScore::write_read_info(ostream& out, const Read& read, bfloat value)
{
  out << read.identifier << "\t" << name << "\t" << value << endl; 
}
//                                                                  i                   j           
inline list<state> ReadProfileScore::get_possible_HMM_states(int readPosition, state profileState)
{
  //  cerr<<"Error: must add start state to this function\n"; 
  list<state> toReturn; 
  if (readPosition == readSize-1) // pre-end "state"
    {
      if (profileState == profile->pre_end_state)
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
      if (profileState == profile->pre_end_state)
	{
	  toReturn.push_back(pairHMM.end_state);
	  return toReturn; 
	}
      else
	return pairHMM.non_special_states; // not start, end
    }
}

void ReadProfileScore::fill_DP_matrix(const Read& read, ostream& hmmoc, bool hmmoc_only, bool backPointers, bool logging)
{
  // First, store the read as a Weight_profile in the pairHMM - the object which actually
  // does the probability calculations
  pairHMM.read_profile = pairHMM.sub_alphabet.new_seq2weight(read.sstringRep); 
    
  // The main DP recursion for scoring a read to a profile via an arbitrary pairHMM
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
  
  // "practical" range is 35 to ~500

  for (i=0; i<readSize; ++i) 
    {
      if (logging)
	cerr<<"\n\nProcessing position " << i  <<  " of read...\n"; 

      // practical range is 300 - 1000s
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
		cerr<<"\n\nProcessing state " << pairHMM.state_name[*k]  <<  " of HMM...\n"; 
	      // the pairHMM "changes" state at every  step, so this is a trivial lookup
	      incoming_HMM_states = pairHMM.incoming[*k]; 
	      if (logging)
		{
		  cerr<<"Summing over " << incoming_HMM_states.size() << " incoming HMM states:\t"; 
		  for (kPrime = incoming_HMM_states.begin(); kPrime != incoming_HMM_states.end(); ++kPrime)
		    cerr<< pairHMM.state_name[*kPrime] << "  "; 
		  cerr<<endl; 
		}
	      for (kPrime = incoming_HMM_states.begin(); kPrime != incoming_HMM_states.end(); ++kPrime)
		{
		  if (logging)
		    cerr<<"\tInvestigating transition from : "<< pairHMM.state_name[*kPrime] <<endl; 
		  iPrime = get_incoming_read_state(i, *kPrime); 
		  // We can be sure this DP cell will be -inf, so we continue
		  if (iPrime == -1 && pairHMM.state_type[*kPrime] != "start")
		    continue;
		  hmm_transition_weight = pairHMM.get_transition_weight(*kPrime, *k); 
		  if (logging)
		    cerr<<"\tTransition weight between HMM states: " << hmm_transition_weight << endl; 
		  //                                                 (prof  hmm)
		  incoming_profile_states = get_incoming_profile_states(*j, *kPrime); 
		  if (logging)
		    {
		      cerr<<"\tDetermined that there were " << incoming_profile_states.size() << " incoming profile states:\t"; 
		      for (vector<state>::iterator pIter=incoming_profile_states.begin(); pIter!=incoming_profile_states.end(); pIter++)
			cerr<<*pIter<<" "; 
		      cerr<<endl; 
		    }
			
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
		      if ( ! bfloat_is_nonzero( toAdd ) )
			{
			  //cerr<<"\tThis incoming trans had 1e-inf Forward value\n";
			  continue;
			}
		      if (*j != *jPrime)
			toAdd *= profile->get_transition_weight(*jPrime, *j); 
		      if (pairHMM.state_type[*k] != "start" && pairHMM.state_type[*k] != "end")
			toAdd *= pairHMM.get_emission_weight(i,*j,*k); 
		      toAdd *= hmm_transition_weight; 
		      // if reqested, add this composite state's info (emission, transition) to the HMMoC stream
		      if (hmmoc_only)
			add_state_to_HMMoC(hmmoc, *j,*k,*jPrime, *kPrime);
		      
		      // Set the value
		      add_to_DP_cell(i,*j,*k, toAdd); 
		      if (backPointers) // backpointers might be requested for viterbi traceback
			add_backPointer(i,*j,*k, iPrime, *jPrime, *kPrime); 
		      if (logging)
			{
			  cerr<<"Added to DP matrix: \n";
			  cerr<<"\tHMM state " << pairHMM.state_name[*kPrime] << endl; 
			  cerr << "\tRead position " << iPrime <<"\n\tProfile state " << *jPrime << "\n\tValue " << toAdd << "\n"; 
			  cerr<<"Forward value for this state: " << get_DP_cell(iPrime,
										*jPrime, *kPrime) << "\n";
			  cerr<<"Emission value for the dest. state: " << pairHMM.get_emission_weight(i,*j,*k) << endl; 
			  cerr<< "HMM transition weight: " << hmm_transition_weight  << "\n\n"; 
			}
		    }
		}
	      if (logging)
		cerr<<"\n\n Final forward value for this state: " << get_DP_cell(i,*j,*k)  << "\n\n"; 
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
  //  if (readIndex == 0)
    //    cerr<<"Warning: incoming on start of read called!\n";
  
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
ReadProfileModel::ReadProfileModel(void) //Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix_in)
{
  // eventually we'd like to optimize this value!
  branch_length = 0.01; 
  
  // Build up the HMM 
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
  //  add_transition("start", "read_ins", .01);

  // from pre_read_ins
  add_transition("pre_read_ins", "pre_read_ins", .99); 
  add_transition("pre_read_ins", "match", .005); 
  add_transition("pre_read_ins", "delete", .005); 

  // from match
  add_transition("match", "match", 0.999); 
  add_transition("match", "read_ins", 0.00025); 
  add_transition("match", "delete", 0.00025); 
  add_transition("match", "post_read_ins", 0.00025); 
  add_transition("match", "end", 0.00025); 


  // from read_ins
  add_transition("read_ins", "read_ins", 0.00025); 
  add_transition("read_ins", "match", .999); 
  add_transition("read_ins", "delete", 0.00025); 
  add_transition("read_ins", "post_read_ins", 0.00025); 
  add_transition("read_ins", "end", 0.00025); 

  // from delete
  add_transition("delete", "read_ins", 0.00025); 
  add_transition("delete", "match", 0.999); 
  add_transition("delete", "delete", 0.00025); 
  add_transition("delete", "post_read_ins", 0.00025); 
  add_transition("delete", "end", 0.00025); 

  // from post_read_ins
  add_transition("post_read_ins", "post_read_ins", 0.999); 
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
    {
      if (! tmpIter->second > 0)
	{
	  cerr <<"ERROR: transition weight was zero!\n"; 
	  exit(1); 
	}
      return tmpIter->second;
    }
}
inline bfloat ReadProfileModel::get_emission_weight_by_alphabet_index(int readIndex, state profileState, state hmm_state)
{
  // Warning: not yet meaningful!!
  return 1.0; 
}

inline bfloat ReadProfileModel::get_emission_weight(int readIndex, state profileState, state hmm_state)
{
  // emissionWeight is an attribute to avoide initializing it each time this fn is called
  emissionWeight = 0.0;
  //  cerr<< "Size of read profile: " << read_profile.size() << endl; 
  //  cerr<< "Index queried: " << readIndex << endl; 
  vector<sstring> toks = sub_alphabet.tokens();
  
  // Match state: emit to both read and profile
  if (state_type[hmm_state]  == "match")
    {
      for (symbolIter = read_profile[readIndex].begin(); 
	   symbolIter != read_profile[readIndex].end(); symbolIter++)
	for (alphIdx = 0; alphIdx < alphabet_size; alphIdx++)
	  {
	    emissionWeight += bfloat(equilibrium_dist[alphIdx])* // prior on this character
	      profile->get_absorb_weight(profileState, alphIdx)* //absorption by profile
	      bfloat(conditional_sub_matrix(alphIdx,symbolIter->first))* // substitute to read's character
	      bfloat(symbolIter->second); // read's affinity for the given character
	  }
      //cerr<<"\n\nThe emission weight was: " << emissionWeight << endl; 
    }

  // Insert state: emit only to profile - product of equilibrium prob and absorb weight
  else if ( state_type[hmm_state] == "insert" )
    for (alphIdx = 0; alphIdx < alphabet_size; alphIdx++)
      emissionWeight += profile->get_absorb_weight(profileState, alphIdx) * // absorption by profile
	bfloat(equilibrium_dist[alphIdx]); // prior on this character

  // Delete state: emit only to read - product of equilibrium prob and affinity for character
  else if ( state_type[hmm_state] == "delete" )
      for (symbolIter = read_profile[readIndex].begin(); 
	   symbolIter != read_profile[readIndex].end(); symbolIter++)
	{
	  // cerr<<toks[*symbolIter] << " has weight: " << equilibrium_dist[alphIdx] <<endl; 
	  emissionWeight += bfloat(equilibrium_dist[symbolIter->first])*// prior on character
	    bfloat(symbolIter->second); // read's affinity for given character
	}
  else
    emissionWeight = bfloat(1.0);

  if (! emissionWeight > 0)
    cerr <<"\n\nWARNING: emission weight is zero!\n\n"; 
  //  cerr<<"Done\t";
  return emissionWeight; 
}

void ReadProfileModel::set_substitution_model(Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix_in, AbsorbingTransducer *profile_in)
{
  profile = profile_in; 
  sub_alphabet = alphabet_in; 
  alphabet_size = alphabet_in.size();
  conditional_sub_matrix = rate_matrix_in.create_conditional_substitution_matrix(branch_length); 
  equilibrium_dist = rate_matrix_in.create_prior(); 
}

// HMMoC Adapter - an option to the ReadProfileScore class.  Not really functional just yet...
void ReadProfileScore::HMMoC_adapter(const char* filename, bool precompute)
{
  start_clique.clear(); mainClique.clear(); end_clique.clear(); 

  stringstream name; 
  name << "Node_"<<profile->treeNode; 

  ofstream out; 
  out.open(filename); 
  // Export a Profile / pairHMM composition as a HMMoC xml file
  // Initialize the HML
  out << "<hml debug='true'>\n\n";
  // Author - Westesson-Holmes
  out<< "<author>Oscar Westesson, Ian Holmes</author>\n\n";

  // Define alphabet - get this from profile
  //  out << "<alphabet id='profileAlphabet'>\n <range from='0' to='"<< profile->alphabet_size << "'/>\n\n"; 
  alphabet = "ACGT";   // for testing
  out << "<alphabet id='profileAlphabet'>\nACGT\n"; 
  profile->alphabet_size = 4;   // for testing
  out << "</alphabet>"; 

  // Define output tape - a single sequence read
  out <<"<output id='sequenceRead'>\n";
  out <<" <alphabet idref='profileAlphabet'/>\n"; 
  out <<" <identifier type='length' value='iLen'/>\n"; 
  out <<" <identifier type='sequence' value='aSeq'/>\n"; 
  out <<" <code type='parameter' value='char *aSeq'/>\n";
  out <<" <code type='parameter' value='int iLen'/>\n"; 
  out <<"</output>\n\n";


  
  // Define, describe the HMM
  out <<"<hmm id='" << name.str() << "_HMM'>\n";
  out <<"<description> Profile for node " << name.str() <<endl; 
  out <<"States in the HMM are built of a profile state and a pairHMM state, denoted like this:\nn~S\twhere n is the nth absorbing state of the profile, and S is a particular state in the pairHMM </description>\n\n"; 

  // The HMM outputs only to the above-defined output "sequenceRead"
  out <<"<outputs id='mainOut'>\n"; 
  out <<"<output idref='sequenceRead'/>\n"; 
  out <<"</outputs>\n\n"; 
  
  // A special null emission for start and end states
  out << "<emission id='NEempty'>\n"; 
  out << "<probability>\n"; 
  out << "<code> 1.0 </code>\n"; 
  out << "</probability>\n"; 
  out << "</emission>\n"; 

  // call the DP routine with a special option. 
  Read dummyRead; dummyRead.set("AAAAA"); dummyRead.identifier="dummyRead"; 
  readSize = 5; 
  fill_DP_matrix(dummyRead, // read
		 out, // hmmoc outstream
		 true, // write to hmmoc
		 false, // keep backpointers
		 false); // debugging logging
  
  // Define the 3 basic cliques based on pairHMM states
  vector<string>::iterator cliqueIter; 
  out << "<clique id='start' >\n"; 
  for ( cliqueIter = start_clique.begin(); cliqueIter != start_clique.end(); cliqueIter++)
    out << "<state idref='" << *cliqueIter << "' />\n"; 
  out << "</clique>\n\n"; 

  out << "<clique id='mainClique' >\n"; 
  for ( cliqueIter = mainClique.begin(); cliqueIter != mainClique.end(); cliqueIter++)
    out << "<state idref='" << *cliqueIter << "' />\n"; 
  out << "</clique>\n\n"; 

  out << "<clique id='end' >\n"; 
  for ( cliqueIter = end_clique.begin(); cliqueIter != end_clique.end(); cliqueIter++)
    out << "<state idref='" << *cliqueIter << "' />\n"; 
  out << "</clique>\n\n"; 
  
  
  // Refernce the 3 basic cliques in the 'graph' 
  out << "<graph>\n"; 
  out << "<clique idref='start' />\n";
  out << "<clique idref='mainClique' />\n";
  out << "<clique idref='end' />\n";
  out << "</graph>\n\n"; 

  // Write the transitions within the hmm
  out <<"<transitions>\n"; 
  out << transitionStream.str();  
  out <<"</transitions>\n\n"; 

  // Finished defining HMM
  out <<"</hmm>\n\n"; 

  //   Request Forward algorithm for our HMM
  out <<"<forward  outputTable='yes' name='Forward' id='"<<name.str()<<"_forward'>\n";
  out <<" <hmm idref='"<<name.str()<<"_HMM'/>\n";
  out <<"</forward>\n\n"; 

  // Ask for the actual code to be made:
  out << "<codeGeneration realtype='bfloat' file='" <<name.str()<<".cc' header='" << name.str() << ".h' language='C++'>\n";
  out << "<forward idref='" << name.str()<<"_forward' />\n"; 
  out << "</codeGeneration>\n\n"; 
  
  out << "</hml>\n\n"; 
  out.close();
}

void ReadProfileScore::add_state_to_HMMoC(ostream& out, state toProfile, state toHMM, state fromProfile, state fromHMM, bool logging )
{
  bool have_added; 
  compositeState to, from;
  pair<compositeState, compositeState> transition; 
  
  to.first = toProfile; to.second = toHMM; 
  from.first = fromProfile; from.second = fromHMM; 
  transition.first = from; transition.second=to; 

  vector<compositeState> statesToCheck; statesToCheck.push_back(to); statesToCheck.push_back(from); 
  for (vector<compositeState>::iterator stateToAdd=statesToCheck.begin(); stateToAdd != statesToCheck.end(); stateToAdd++)
    {
      have_added = false; 
      for (vector<compositeState>::iterator addedIter = states_added.begin(); addedIter!=states_added.end(); addedIter++)
	if (*addedIter == *stateToAdd)
	  {
	    have_added = true; 
	    break;
	  }
      if ( !have_added )
	{
	  // Add it to the record of states
	  states_added.push_back(*stateToAdd); 
	  if ( stateToAdd->first != profile->start_state && stateToAdd->first != profile->pre_end_state )
	    {
	      // Create the emission that this state uses - this is taken fairly closely from the HMMoC doc example....
	      out << "<emission id='" << composite2string(*stateToAdd)<<"_emit'>\n";
	      out << "<output idref='sequenceRead' />\n"; 
	      out << "<probability>\n"; 
	      out << "<code type='statement'>\n"; 
	      out << "<identifier output='sequenceRead' value='iSymbol'/>\n";
	      out << "<identifier type='result' value='iProb'/>\n";

	      out << "<![CDATA[\n"; 
		// the way suggested in HMMoC mnl - doesn't seem to compile:
	      //	      out << "<![CDATA[\n   switch (iSymbol) { \n"; 
		//		out << "      case '" << alphabet << "': iProb = " << pairHMM.get_emission_weight_by_alphabet_index(i, stateToAdd->first, stateToAdd->second)  << "; break;\n"; 
	      for (int i=0; i< profile->alphabet_size; ++i)
		{
		  continue; // for testing
		  if (i)
		    out <<"else "; 
		  out << "if ( iSymbol == " << i << ") iProb = " << pairHMM.get_emission_weight_by_alphabet_index(i, stateToAdd->first, stateToAdd->second)  << ";\n"; 
		}
	      out << "iProb=1.0;\n\n"; 
	      // Need for switch statement way:
	      //	      out << "   }\n"; 
	      out << "]]>\n"; 
	      out << "</code>\n"; 
	      out << "</probability>\n"; 
	      out << "</emission>\n\n"; 

	      // Add state info - basically just the emission distribution along with some other basic info
	      out << "<state id = '" << composite2string(*stateToAdd) << "' emission='" << composite2string(*stateToAdd) << "_emit' />\n\n";
	    }	  
	  else
	    out << "<state id = '" << composite2string(*stateToAdd) << "' emission='NEempty' />\n"; 
	  if (stateToAdd->second == pairHMM.start_state && stateToAdd->first == profile->start_state)
	    start_clique.push_back(composite2string(*stateToAdd)); 
	  else if (stateToAdd->second == pairHMM.end_state && stateToAdd->first == profile->pre_end_state)
	    end_clique.push_back(composite2string(*stateToAdd)); 
	  else
	    mainClique.push_back(composite2string(*stateToAdd)); 				   
	}
    }
  // add the transition, if necessary
  have_added = false;
  for (vector<pair<compositeState, compositeState > >::iterator transIter = transitions_added.begin(); transIter!=transitions_added.end(); transIter++)
    if (*transIter == transition)
      {
	have_added = true; 
	break;
      }
  if ( !have_added )
    {
      // Add it to the record of transitions
      transitions_added.push_back(transition);
      
      // The transition weight is a product of the HMM transition weight and the profile's transition weight, but only if the states differ
      bfloat transWeight = 1.0;
      if (toProfile != fromProfile)
	transWeight *= profile->get_transition_weight(fromProfile, toProfile);
      if (toHMM != fromHMM)
	transWeight *= pairHMM.get_transition_weight(fromHMM, toHMM);
      // for testing
      transWeight = 1.0; 
      if ( index(transWeight, probValues ) == -1 )
	{
	  probValues.push_back(transWeight); 
	  transitionStream << "<probability id='" << transWeight <<"'>\n"; 
	  transitionStream << "<code value='" << transWeight <<"' />\n";
	  transitionStream << "</probability>\n";
	}
      transitionStream << "<transition from='" << composite2string(from) << "' to='" << composite2string(to)<<"' probability='" << transWeight <<"' />\n\n"; 
    }

}

string ReadProfileScore::composite2string(compositeState in)
{
  stringstream s; 
  s << in.first << "~" << pairHMM.state_name[in.second]; 
  return string(s.str()); 
}
