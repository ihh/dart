#include<iostream>
#include<math.h>

#include "protpal/transducer.h"
#include "ecfg/ecfgsexpr.h" // for array2d class
#include "util/sstring.h"

//General transducer operations
string Transducer::get_state_name(state m)
{
  if (m >= int(states.size()))
	{ 
	  std::cerr<<"Error: state "<<m<<"has no name, since it is not in the transducer " << name << " which has "<<states.size()<<" states \n";
	  exit(1);
	}
  else 
	{
	  return state_names[m];
	}
}

vector<state> Transducer::get_state_index(string state_name)
{
  vector<state> out;
  for (unsigned int i=0; i<state_names.size(); i++)
	{ 
	  if(state_names[i] == state_name) out.push_back(i);
	}
  
  if (out.empty()) 
	{
	  std::cerr<<"Error: state "<<state_name<<" has no index, since it is not in the transducer  " << name << " which has the following states  ";
  displayVector(state_names);
  exit(1);
	}
  else return out;

}

string Transducer::get_state_type(state m)
{
  if (m >= int(states.size()))
	{
	  std::cerr<<"Error: state "<<m<<"has no type, since it is not in the transducer  " << name << " which has "<<states.size()<<" states \n";
	  exit(1);
	}
  else return(state_types[m]);
}

vector<state> Transducer::get_state_type_set(string type)
{
  vector<state> out; 
  for (unsigned int i=0; i<state_types.size(); i++)
	{
	  if (state_types[i] == type) out.push_back(i); 
	}
  return out; 
}

vector<state> Transducer::get_outgoing(state m)
{
  if (outgoing.count(m)>0) return outgoing[m];
  else
	{
	  std::cerr<<"Error: state "<<m<<" has no outgoing connections\n";
	  exit(1);
	}
}

bool Transducer::has_transition(state m, state mPrime)
{
  if (outgoing.count(m)<0) return false; 
  else
	{
	  if (index(mPrime, outgoing[m]) == -1) return false;
	  else return true;
	}
}

void Transducer::test_transitions(void)
{
  std::cerr<<"Testing transducer "<<name<<"'s transitions...\n";
  double sum; 
  for (unsigned int i=0; i<state_names.size(); i++)
	{
	  sum = 0; 
	  std::cerr<<"State "<<i<<": "<<state_names[i]<<" has outgoing transitions to:\n";
	  if (state_names[i]=="end") std::cerr<<"\tNone\n";
	  for (vector<int>::iterator j=outgoing[i].begin(); j!=outgoing[i].end(); j++)
		{
		  std::cerr<<"\tState "<<*j<<": "<<state_names[*j]<<" with weight "<<get_transition_weight(i,*j)<<endl;
		  sum += get_transition_weight(i,*j); 		  
		}
	  std::cerr<<"\tTotal outgoing trans weight for this state: "<<sum<<endl; 
	}
}
  
double Transducer::get_transition_weight(state m, state mPrime)
{
  vector<state> transitionPair;
  transitionPair.push_back(m);  transitionPair.push_back(mPrime);
  if (transition_weight.count(transitionPair) >0) return transition_weight[transitionPair];
  else 
	{
	  std::cerr<<"Error: state "<<m<<" has no nonzero-weight transition to "<<mPrime<<" in transducer  " << name <<endl;
	  std::cerr<<"The outgoing transitions of "<<m<<" are: ";
	  displayVector(get_outgoing(m));
	  exit(1);
	}
}

vector<double> Transducer::get_emission_distribution(state b)
{
  if (get_state_type(b) == "I")
	return emission_weight_matrix[b];
  else
	{
	  std::cerr<<"Can't get emission distribution of non-emit state!\n";
	  exit(1);
	}
}

double Transducer::get_emission_weight(state b, int charIndex)
{
  string insert = "I";
  if (get_state_type(b) != insert) 
	{
	  std::cerr<<"Calling get_emission_weight on a non-insert state is not allowed!\n";
	  std::cerr<<"The offending call was: transducer "<<name<<" state "<<b<<" named: "<<get_state_name(b)<<endl;
	  exit(1);
	}
  else if (charIndex >= alphabet_size || charIndex <0 )
	{	
	  std::cerr<<"Calling get_emission_weight on a non-alphabet index!\n";
	  std::cerr<<"The offending call was: transducer "<<name<<" character index: "<<charIndex<<endl;
	  exit(1);
	}
  else
	{
	  return(emission_weight_matrix[b][charIndex]);	  
	}
}



// Operations for singlet transducers (constructor w/ hard-coded stuff, etc)
SingletTrans::SingletTrans(void)
{
}

SingletTrans::SingletTrans(Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix, double insertProb)
{
  name = "Singlet";
  
  vector<sstring> toks = alphabet_in.tokens(); 
  for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
	alphabet.push_back(string(a->c_str()));

  alphabet_size = alphabet.size();
  
/*  states:  */
/*   0 = start */
/*   1 = insert */
/*   2 = wait */
/*   3 = end */  
  state_names.push_back("start");  
  state_names.push_back("insert");
  state_names.push_back("wait");
  state_names.push_back("end");

  state_types.push_back("S");  
  state_types.push_back("I");
  state_types.push_back("W");
  state_types.push_back("E");

  for (unsigned int i=0; i<state_names.size(); i++) states.push_back(i); 


  // outgoing transitions
  vector<state> out;
  out.push_back(1);out.push_back(2) ;outgoing[0] = out;  // start -> (insert|wait)
  out.clear(); 
  out.push_back(2);out.push_back(1); outgoing[1] = out; // insert -> wait|insert
  out.clear(); 
  out.push_back(3); outgoing[2] = out; // wait -> end

  // Transition weights - hard coded here.  
  double start2ins = insertProb;
  double start2wait = 1.0 - insertProb;
  double ins2ins = insertProb; 
  double ins2wait = 1.0 - insertProb; 

  vector<state> transitionPair; 
  transitionPair.push_back(0); transitionPair.push_back(1); 
  transition_weight[transitionPair] = start2ins; 

  transitionPair.clear();
  transitionPair.push_back(0); transitionPair.push_back(2); 
  transition_weight[transitionPair] = start2wait; 

  transitionPair.clear();
  transitionPair.push_back(1); transitionPair.push_back(1); 
  transition_weight[transitionPair] = ins2ins; 

  transitionPair.clear();
  transitionPair.push_back(1); transitionPair.push_back(2); 
  transition_weight[transitionPair] = ins2wait; 

  transitionPair.clear();
  transitionPair.push_back(2); transitionPair.push_back(3); 
  transition_weight[transitionPair] = 1; 

  // Emission weights
  // The states of type I are the only states which have an emission distribution
  vector<double> equilibrium = rate_matrix.create_prior(); 
  for (unsigned int i=0; i<states.size(); i++)
	{
	  if (state_types[i] == "I") emission_weight_matrix[i] = equilibrium; 
	}
}


// Splitting/forking transducer constructor
// In retrospect, this is a bit unnecessary...not sure why I put this in :-)
SplittingTrans::SplittingTrans(void)
{
  name = "Upsilon";

  state_names.push_back("start");  
  state_names.push_back("match");
  state_names.push_back("wait");
  state_names.push_back("end");


  state_types.push_back("S");  
  state_types.push_back("M");
  state_types.push_back("W");
  state_types.push_back("E");

  /*  states:  */
  /*   0 = start */
  /*   1 = match */
  /*   2 = wait */
  /*   3 = end */  
  for (unsigned int i=0; i<state_names.size(); i++) states.push_back(i);


  // outgoing transitions
  vector<state> out;
  out.push_back(1); out.push_back(2); outgoing[0] = out;  // start(0) -> match(1)
  out.clear(); 
  out.push_back(2); outgoing[1] = out; // match(1) -> wait(2)
  out.clear(); 
  out.push_back(1); out.push_back(3); outgoing[2] = out; // wait(2) -> (match(1) | end(3))

  // Transition weights are all 1. This tranducer is dependent on the others around it to make transitions.

  vector<state> transitionPair;  
  transitionPair.push_back(0); transitionPair.push_back(1); 
  transition_weight[transitionPair] = 1; 

  transitionPair.clear();
  transitionPair.push_back(0); transitionPair.push_back(2); 
  transition_weight[transitionPair] = 1; 
  
  transitionPair.clear();
  transitionPair.push_back(1); transitionPair.push_back(2); 
  transition_weight[transitionPair] = 1;

  transitionPair.clear();
  transitionPair.push_back(2); transitionPair.push_back(1); 
  transition_weight[transitionPair] = 1; 

  transitionPair.clear();
  transitionPair.push_back(2); transitionPair.push_back(3); 
  transition_weight[transitionPair] = 1; 

 }

// Branch transducer constructor

BranchTrans::BranchTrans(void)
{
  // Placeholder constructor
  //  std::cerr<<"Initiating a placeholder branch transducer object\n";
}


double BranchTrans::get_match_weight(state b, int incoming_character, int outgoing_character )
{
  // Match weight function.  Branch-length dependent.  Eventually, we'd like to have this use matrix
  // exponentiation to allow for a general character-evo model.  For now, a uniform model is coded
  // here, for testing purposes only!!

  string match = "M"; 
  #ifdef DART_DEBUG
  if (get_state_type(b) != match) 
	{
	  std::cerr<<"Calling get_match_weight on a non-match state!\n";
	  std::cerr<<"The offending call was: transducer "<<name<<" state "<<b<<" named: "<<get_state_name(b)<<endl;
	  exit(1);
	}
  else if (incoming_character >= alphabet_size || outgoing_character >= alphabet_size )
	{
	  
  	  std::cerr<<"Calling get_match_weight on a non-alphabet index!\n";
	  std::cerr<<"The offending call was: transducer "<<name<<" character indices (in,out): "<<incoming_character<<" "<<outgoing_character<<endl;
	  exit(1);
	}
#endif
  double result = conditional_sub_matrix(incoming_character, outgoing_character); 

  #ifdef DART_DEBUG
  if (result <= 0)
    cerr<<"Warning: match weight is not positive.  This may indicate a problem with the substitution matrix\n";
  #endif
  return result; 
  // 	  float n = alphabet_size;
  // 	  // Uniform match distribution, from early testing days
  // 	  if (incoming_character == outgoing_character) return  1/n + ((n-1)/n)*exp(-branch_length*n); 
  // 	  else return 1/n - (1/n)*exp(-branch_length*n); 
}  


// affine gap branch transducer
// BranchTrans::BranchTrans(double branch_length_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix, 
// 			 double ins_open_rate, double del_open_rate, double gap_extend, double subRate)
// {
//   branch_length = branch_length_in;
//   name = "Branch";

//   vector<sstring> toks = alphabet_in.tokens(); 
//   for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
// 	alphabet.push_back(string(a->c_str()));

//   alphabet_size = alphabet.size();
  
//   conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(branch_length*subRate); 

//   // Initialize state names
//   state_names.push_back("start");  
//   state_names.push_back("match");
//   state_names.push_back("wait");
//   state_names.push_back("delete_wait");  
//   state_names.push_back("delete");
//   state_names.push_back("insert");
//   state_names.push_back("end");

//   // Initialize state names
//   state_types.push_back("S");  
//   state_types.push_back("M");
//   state_types.push_back("W");
//   state_types.push_back("W");  
//   state_types.push_back("D");
//   state_types.push_back("I");
//   state_types.push_back("E");

//   /*  states:  */
//   /*   0 = start */
//   /*   1 = match */
//   /*   2 = wait */
//   /*   3 = delete_wait */
//   /*   4 = delete */
//   /*   5 = insert */  
//   /*   6 = end */    

//   int start=0, match=1, wait=2, delete_wait=3, del=4, insert=5, end=6; 

//   for (int i=0; i<state_names.size(); i++) states.push_back(i); 

//   // Set outgoing transitions
//   vector<state> out;

//   out.push_back(insert); out.push_back(wait); outgoing[start] = out; // start -> ( wait() | insert() )

//   out.clear(); 
//   out.push_back(insert); out.push_back(wait); outgoing[match] = out; // match(1) -> (insert() | wait() )

//   out.clear(); 
//   out.push_back(del); out.push_back(match); out.push_back(end); outgoing[wait] = out; // wait -> D|M|E

//   out.clear(); 
//   // out.push_back(end);
//   out.push_back(del); outgoing[delete_wait] = out; // delete_wait -> delete
  
//   out.clear(); 
//   out.push_back(delete_wait); out.push_back(insert); out.push_back(wait); outgoing[del] = out; // del -> delete_wait|I|W
  
//   out.clear(); 
//   out.push_back(wait); out.push_back(insert);  outgoing[insert] = out; // insert -> I|W
  

//   // Transition weights are determined by the insertion/deletion rates and the branch_length (set on input)

//   double ins_extend = gap_extend;
//   double del_extend = gap_extend;
//   double wait2end = 1; // This is a forced transition...should have weight 1.

//   // Probability/weight of an insertion 
//   double ins_open = 1-exp(-ins_open_rate*branch_length);
//   double del_open = 1-exp(-del_open_rate*branch_length);
  
//   // Define their complements, for cleanliness 
//   double no_ins_open = 1-ins_open;
//   double no_del_open = 1-del_open; 
  
//   double no_ins_extend = 1-ins_extend;
//   double no_del_extend = 1-del_extend;
  

//   vector<state> transitionPair;  

//   // Transitions from start 
//   // start -> wait = no insertion
//   transitionPair.push_back(start); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_open;   transitionPair.clear();
//   // start -> insert =  insertion
//   transitionPair.push_back(start); transitionPair.push_back(insert); 
//   transition_weight[transitionPair] =ins_open;   transitionPair.clear();

//   // Transitions from match 
//   // match -> wait = no insertion
//   transitionPair.push_back(match); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_open;   transitionPair.clear();
//   // match -> insert =  insertion
//   transitionPair.push_back(match); transitionPair.push_back(insert); 
//   transition_weight[transitionPair] = ins_open;   transitionPair.clear();

//   // Transitions from wait 
//   // wait -> match = no deletion
//   transitionPair.push_back(wait); transitionPair.push_back(match); 
//   transition_weight[transitionPair] = no_del_open;   transitionPair.clear();
//   // wait -> delete =  deletion
//   transitionPair.push_back(wait); transitionPair.push_back(del); 
//   transition_weight[transitionPair] = del_open;   transitionPair.clear();
//   // wait -> end =  forced transition
//   transitionPair.push_back(wait); transitionPair.push_back(end); 
//   transition_weight[transitionPair] = wait2end;   transitionPair.clear();


//   // Transitions from delete 
//   // delete -> wait = no insertion, no delete-extend
//   transitionPair.push_back(del); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_del_extend*no_ins_open;   transitionPair.clear();
//   // delete -> delete_wait = delete-extend
//   transitionPair.push_back(del); transitionPair.push_back(delete_wait); 
//   transition_weight[transitionPair] = del_extend;   transitionPair.clear();
//   // delete -> insert = yes insertion, no delete-extend
//   transitionPair.push_back(del); transitionPair.push_back(insert); 
//   transition_weight[transitionPair] = no_del_extend*ins_open;   transitionPair.clear();

//   //Transitions from delete_wait
//   // delete_wait -> del
//   transitionPair.push_back(delete_wait); transitionPair.push_back(del); 
//   transition_weight[transitionPair] = 1;   transitionPair.clear();


//   // NB I've removed this transition, as it has caused problems later on.  I don't think it's really
//   // necessary, though I'm somewhat upset about having to remove it. -OW
//   // delete_wait -> end 
//   //  transitionPair.push_back(delete_wait); transitionPair.push_back(end); 
//   //  transition_weight[transitionPair] = 1;   transitionPair.clear();

//   // Transitions from insert 
//   // insert -> wait = no extend insertion
//   transitionPair.push_back(insert); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_extend;   transitionPair.clear();
//   // insert -> insert =  extend insertion
//   transitionPair.push_back(insert); transitionPair.push_back(insert); 
//   transition_weight[transitionPair] = ins_extend;   transitionPair.clear();

//   // Emission weights
//   // Equilibrium over alphabet characters
//   // The states of type I are the only states which have an emission distribution

//   vector<double> equilibrium = rate_matrix.create_prior(); 
//   for (int i=0; i<states.size(); i++)
// 	{
// 	  if (state_types[i] == "I") emission_weight_matrix[i] = equilibrium;
// 	}
// }


// affine gap branch transducer - "ian" style (e.g. wait before action)
BranchTrans::BranchTrans(double branch_length_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix, 
			 double ins_open_rate, double del_open_rate, double gap_extend, double subRate)
{
  branch_length = branch_length_in;
  name = "Branch";

  vector<sstring> toks = alphabet_in.tokens(); 
  for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
	alphabet.push_back(string(a->c_str()));

  alphabet_size = alphabet.size();
  
  conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(branch_length*subRate); 

  // Initialize state names
  state_names.push_back("start");  
  state_names.push_back("match");
  state_names.push_back("match_wait");
  state_names.push_back("delete_wait");  
  state_names.push_back("delete");
  state_names.push_back("insert");
  state_names.push_back("end");
  state_names.push_back("end_wait");

  // Initialize state types
  state_types.push_back("S");  
  state_types.push_back("M");
  state_types.push_back("W");
  state_types.push_back("W");  
  state_types.push_back("D");
  state_types.push_back("I");
  state_types.push_back("E");
  state_types.push_back("W");

  /*  states:  */
  /*   0 = start */
  /*   1 = match */
  /*   2 = match_wait */
  /*   3 = delete_wait */
  /*   4 = delete */
  /*   5 = insert */  
  /*   6 = end */    
  //   7 = end_wait

  int start=0, match=1, match_wait=2, delete_wait=3, del=4, insert=5, end=6, end_wait = 7; 

  for (unsigned int i=0; i<state_names.size(); i++) states.push_back(i); 

  // Set outgoing transitions
  vector<state> out;

  // out from start
  out.push_back(insert); out.push_back(end_wait); out.push_back(match_wait); out.push_back(delete_wait); 
  outgoing[start] = out; // start -> ( match_wait() | delete_wait | insert() | end_wait

  // out from match
  out.clear(); 
  out.push_back(end_wait);   out.push_back(insert); out.push_back(match_wait); out.push_back(delete_wait); 
  outgoing[match] = out; // match(1) -> (insert() | delete_wait | match_wait() | end_wait

  // out from match_wait
  out.clear(); 
  out.push_back(match);
  outgoing[match_wait] = out; // match_wait -> M

  // out from delete
  out.clear(); 
  out.push_back(end_wait); out.push_back(delete_wait); out.push_back(insert); out.push_back(match_wait); 
  outgoing[del] = out; // del -> delete_wait | insert | match_wait | end_wait

  // out from delete_wait
  out.clear(); 
  out.push_back(del);
  outgoing[delete_wait] = out; // delete_wait -> delete

  // out from insert 
  out.clear(); 
  out.push_back(end_wait); out.push_back(match_wait); out.push_back(insert);  out.push_back(delete_wait); 
  outgoing[insert] = out; // insert -> insert | match_wait | delete_wait | end_wait
 
  // out from end_wait
  out.clear(); 
  out.push_back(end); 
  outgoing[end_wait] = out; // end_wait -> end

  // Transition weights are determined by the insertion/deletion rates and the branch_length (set on input)

  double ins_extend = gap_extend;
  double del_extend = gap_extend;
  //double wait2end = 1; // This is a forced transition...should have weight 1.

  // Probability/weight of an insertion 
  double ins_open = 1-exp(-ins_open_rate*branch_length);
  double del_open = 1-exp(-del_open_rate*branch_length);
  
  // Define their complements, for cleanliness 
  double no_ins_open = 1-ins_open;
  double no_del_open = 1-del_open; 
  
  double no_ins_extend = 1-ins_extend;
  double no_del_extend = 1-del_extend;
  

  vector<state> transitionPair;  

  // Transitions from start 
  // start -> match_wait = no insertion, no deletion
  transitionPair.push_back(start); transitionPair.push_back(match_wait); 
  transition_weight[transitionPair] = no_ins_open*no_del_open;   transitionPair.clear();
  // start -> insert =  insertion
  transitionPair.push_back(start); transitionPair.push_back(insert); 
  transition_weight[transitionPair] =ins_open;   transitionPair.clear();
  // start -> delete_wait =  no insertion, deletion
  transitionPair.push_back(start); transitionPair.push_back(delete_wait); 
  transition_weight[transitionPair] = no_ins_open*del_open;   transitionPair.clear();
  // start -> end_wait =  prob 1
  transitionPair.push_back(start); transitionPair.push_back(end_wait); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();


  // Transitions from match 
  // match -> match_wait = no insertion, no deletion
  transitionPair.push_back(match); transitionPair.push_back(match_wait); 
  transition_weight[transitionPair] = no_ins_open*no_del_open;   transitionPair.clear();
  // match -> insert =  insertion
  transitionPair.push_back(match); transitionPair.push_back(insert); 
  transition_weight[transitionPair] = ins_open;   transitionPair.clear();
  // match -> delete_wait =  deletion
  transitionPair.push_back(match); transitionPair.push_back(delete_wait); 
  transition_weight[transitionPair] = no_ins_open*del_open;   transitionPair.clear();
  // match -> end_wait =  prob 1
  transitionPair.push_back(match); transitionPair.push_back(end_wait); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();


  // Transitions from match_wait 
  // match_wait -> match = always happens
  transitionPair.push_back(match_wait); transitionPair.push_back(match); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();

  // Transitions from delete_wait 
  // delete_wait -> delete =  always happens
  transitionPair.push_back(delete_wait); transitionPair.push_back(del); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();

  // Transitions from end_wait
  // end_wait -> end =  always happens
  transitionPair.push_back(end_wait); transitionPair.push_back(end); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();


  // Transitions from delete 
  // delete -> match_wait = no delete-extend, no insertion, 
  transitionPair.push_back(del); transitionPair.push_back(match_wait); 
  transition_weight[transitionPair] = no_del_extend*no_ins_open;   transitionPair.clear();
  // delete -> delete_wait = delete-extend
  transitionPair.push_back(del); transitionPair.push_back(delete_wait); 
  transition_weight[transitionPair] = del_extend;   transitionPair.clear();
  // delete -> insert = yes insertion, no delete-extend
  transitionPair.push_back(del); transitionPair.push_back(insert); 
  transition_weight[transitionPair] = no_del_extend*ins_open;   transitionPair.clear();
  // delete -> end_wait = prob 1
  transitionPair.push_back(del); transitionPair.push_back(end_wait); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();

  //Transitions from delete_wait
  // delete_wait -> del
  transitionPair.push_back(delete_wait); transitionPair.push_back(del); 
  transition_weight[transitionPair] = 1.0;   transitionPair.clear();


  // NB I've removed this transition, as it has caused problems later on.  I don't think it's really
  // necessary, though I'm somewhat upset about having to remove it. -OW
  // delete_wait -> end 
  //  transitionPair.push_back(delete_wait); transitionPair.push_back(end); 
  //  transition_weight[transitionPair] = 1;   transitionPair.clear();

  // Transitions from insert 
  // insert -> matchwait = no extend insertion, no deletion
  transitionPair.push_back(insert); transitionPair.push_back(match_wait); 
  transition_weight[transitionPair] = no_ins_extend*no_del_open;   transitionPair.clear();
  // insert -> delete_wait = no extend insertion, yes deletion
  transitionPair.push_back(insert); transitionPair.push_back(delete_wait); 
  transition_weight[transitionPair] = no_ins_extend*del_open;   transitionPair.clear();
  // insert -> insert =  extend insertion
  transitionPair.push_back(insert); transitionPair.push_back(insert); 
  transition_weight[transitionPair] = ins_extend;   transitionPair.clear();
  // insert -> end_wait =  prob 1
  transitionPair.push_back(insert); transitionPair.push_back(end_wait); 
  transition_weight[transitionPair] = 1;
  transitionPair.clear();


  // Emission weights
  // Equilibrium over alphabet characters
  // The states of type I are the only states which have an emission distribution

  vector<double> equilibrium = rate_matrix.create_prior(); 
  for (unsigned int i=0; i<states.size(); i++)
    if (state_types[i] == "I") 
      emission_weight_matrix[i] = equilibrium;
}






BranchTrans::BranchTrans(double branch_length_in, bool linear)
{
  std::cerr<<"Using linear gap branch transducer...this should no longer be used!\n";
  exit(1);
  branch_length = branch_length_in;
  name = "Branch";
  
  alphabet.push_back("arndcqeghilkmfpstwyv"); //dummy thing
  alphabet_size = alphabet.size();
  
  // Initialize state names
  state_names.push_back("start");  
  state_names.push_back("match");
  state_names.push_back("wait");
  state_names.push_back("delete");
  state_names.push_back("insert");
  state_names.push_back("end");

  // Initialize state names
  state_types.push_back("S");  
  state_types.push_back("M");
  state_types.push_back("W");
  state_types.push_back("D");
  state_types.push_back("I");
  state_types.push_back("E");

  /*  states:  */
  /*   0 = start */
  /*   1 = match */
  /*   2 = wait */
  /*   3 = delete */
  /*   4 = insert */  
  /*   5 = end */    

  for (unsigned int i=0; i<state_names.size(); i++) states.push_back(i); 

  // Set outgoing transitions
  vector<state> out;

  out.push_back(2); out.push_back(4); outgoing[0] = out; // start(0) -> ( wait(2) | insert(4) )

  out.clear(); 
  out.push_back(2); out.push_back(4); outgoing[1] = out; // match(1) -> (insert(4) | wait(2) )

  out.clear(); 
  out.push_back(1); out.push_back(3); out.push_back(5); outgoing[2] = out; // wait(2) -> (match(1) | end(5) | delete(3) )

  out.clear(); 
  out.push_back(2); out.push_back(4); outgoing[3] = out; // delete(3) -> (insert(4) | wait(2))

  out.clear(); 
  out.push_back(2); out.push_back(4); outgoing[4] = out; // insert(4) -> (insert(4) | wait(2))

  // Transition weights are determined by the insertion/deletion rates (set here) and the branch_length (on input)
  double insRate = .1;
  double delRate = .1;
  double wait2end = 1; // This is a forced transition...should have weight 1.

  // Probability/weight of an insertion 
  double insertion_weight = 1-exp(-insRate*branch_length);
  double deletion_weight = 1-exp(-delRate*branch_length);
  
  // Define their complements, for cleanliness 
  double no_insertion_weight = 1-insertion_weight;
  double no_deletion_weight = 1-deletion_weight;
  

  vector<state> transitionPair;  

  // Transitions from start 
  // start -> wait = no insertion
  transitionPair.push_back(0); transitionPair.push_back(2); 
  transition_weight[transitionPair] = no_insertion_weight;   transitionPair.clear();
  // start -> insert =  insertion
  transitionPair.push_back(0); transitionPair.push_back(4); 
  transition_weight[transitionPair] =insertion_weight;   transitionPair.clear();

  // Transitions from match 
  // match -> wait = no insertion
  transitionPair.push_back(1); transitionPair.push_back(2); 
  transition_weight[transitionPair] = no_insertion_weight;   transitionPair.clear();
  // match -> insert =  insertion
  transitionPair.push_back(1); transitionPair.push_back(4); 
  transition_weight[transitionPair] = insertion_weight;   transitionPair.clear();

  // Transitions from wait 
  // wait -> match = no deletion
  transitionPair.push_back(2); transitionPair.push_back(1); 
  transition_weight[transitionPair] = no_deletion_weight;   transitionPair.clear();
  // wait -> delete =  deletion
  transitionPair.push_back(2); transitionPair.push_back(3); 
  transition_weight[transitionPair] = deletion_weight;   transitionPair.clear();

  // wait -> end =  forced transition
  transitionPair.push_back(2); transitionPair.push_back(5); 
  transition_weight[transitionPair] = wait2end;   transitionPair.clear();


  // Transitions from delete 
  // delete -> wait = no insertion
  transitionPair.push_back(3); transitionPair.push_back(2); 
  transition_weight[transitionPair] = no_insertion_weight;   transitionPair.clear();
  // match -> insert =  insertion
  transitionPair.push_back(3); transitionPair.push_back(4); 
  transition_weight[transitionPair] = insertion_weight;   transitionPair.clear();

  // Transitions from insert 
  // insert -> wait = no insertion
  transitionPair.push_back(4); transitionPair.push_back(2); 
  transition_weight[transitionPair] = no_insertion_weight;   transitionPair.clear();
  // insert -> insert =  insertion
  transitionPair.push_back(4); transitionPair.push_back(4); 
  transition_weight[transitionPair] = insertion_weight;   transitionPair.clear();

  // Emission weights
  //  MyMap< state, vector<double> > emission_weight_matrix;     
  MyMap< string, double> nullDist; 
  // NB These are taken from xrate's nullprot.eg grammar's equil distribution.  
  nullDist["a"]=0.0962551;
//   nullDist["r"]=0.0373506;
//   nullDist["n"]=0.0422718;
//   nullDist["d"]=0.0405318;
//   nullDist["c"]=0.0266675;
//   nullDist["q"]=0.0336235;
//   nullDist["e"]=0.0466521;
//   nullDist["g"]=0.0716804;
//   nullDist["h"]=0.0308094;
//   nullDist["i"]=0.0301717;
//   nullDist["l"]=0.114685;
//   nullDist["k"]=0.0651995;
//   nullDist["m"]=0.0179707;
//   nullDist["f"]=0.0477261;
//   nullDist["p"]=0.0691653;
//   nullDist["s"]=0.0696891;
//   nullDist["t"]=0.0602472;
//   nullDist["w"]=0.00452625;
//   nullDist["y"]=0.0430247;
//   nullDist["v"]=0.0517526;
  
  // Order these appropriately for the emission weight matrix
  vector<double> nullDistVector; 
  string omega; 
  for (int i=0; i<alphabet_size; i++)	
	{
	  omega = alphabet[i];
	  nullDistVector.push_back(nullDist[omega]);
	}
  // The states of type I are the only states which have an emission distribution
  string insert="I";
  for (unsigned int i=0; i<states.size(); i++)
	{
	  if (state_types[i] == insert) emission_weight_matrix[i] = nullDistVector;
	}

}



// mixture-of-affine gaps branch transducer - manually inputted
// BranchTrans::BranchTrans(double branch_length_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix, 
// 			 double ins_open_rate, double del_open_rate, double gap_extend, double gap_extend2, double mixPrior, string name)
// {
//   branch_length = branch_length_in;
//   name = "Mixture of affine gaps branch transducers";

//   vector<sstring> toks = alphabet_in.tokens(); 
//   for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
// 	alphabet.push_back(string(a->c_str()));

//   alphabet_size = alphabet.size();
  
//   conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(branch_length); 

//   // Initialize state names
//   state_names.push_back("start");  
//   state_names.push_back("match");
//   state_names.push_back("wait");
//   state_names.push_back("delete_wait_1");  
//   state_names.push_back("delete_1");
//   state_names.push_back("delete_wait_2");  
//   state_names.push_back("delete_2");
//   state_names.push_back("insert_1");
//   state_names.push_back("insert_2");
//   state_names.push_back("end");

//   // Initialize state names
//   state_types.push_back("S");  
//   state_types.push_back("M");
//   state_types.push_back("W");
//   state_types.push_back("W");  
//   state_types.push_back("D");
//   state_types.push_back("W");  
//   state_types.push_back("D");
//   state_types.push_back("I");
//   state_types.push_back("I");
//   state_types.push_back("E");

//   /*
//    states:  
//     0 = start 
//     1 = match
//     2 = wait 
//     3 = delete_wait_1 
//     4 = delete_1
//     5 = delete_wait_2
//     6 = delete_2
//     7 = insert_1
//     8 = insert_2
//     9 = end    
//   */
//   int start=0, match=1, wait=2, delete_wait_1=3, delete_1=4, delete_wait_2=5, delete_2=6, insert_1=7, insert_2=8, end=9; 

//   for (int i=0; i<state_names.size(); i++) states.push_back(i); 

//   // Set outgoing transitions
//   vector<state> out;

//   out.push_back(insert_1); out.push_back(wait); out.push_back(insert_2); outgoing[start] = out; // start -> ( wait | insert_1 | insert_2 )

//   out.clear(); 
//   out.push_back(delete_1); out.push_back(delete_2); out.push_back(match); out.push_back(end); outgoing[wait] = out; // wait -> D1|D2|M|E

//   out.clear(); 
//   out.push_back(insert_1); out.push_back(insert_2); out.push_back(wait); outgoing[match] = out; // match(1) -> I1|I2|W

//   out.clear(); 
//   out.push_back(delete_1); outgoing[delete_wait_1] = out; // delete_wait_1 -> delete_1

//   out.clear(); 
//   out.push_back(delete_2); outgoing[delete_wait_2] = out; // delete_wait_2 -> delete_2
  
//   out.clear(); 
//   out.push_back(delete_wait_1); out.push_back(insert_1); out.push_back(insert_2); out.push_back(wait); outgoing[delete_1] = out; // D1 -> delete_wait1|I1|I2|W

//   out.clear(); 
//   out.push_back(delete_wait_2); out.push_back(insert_1); out.push_back(insert_2); out.push_back(wait); outgoing[delete_2] = out; // D1 -> delete_wait1|I1|I2|W
  
//   out.clear(); 
//   out.push_back(wait); out.push_back(insert_1);  outgoing[insert_1] = out; // I1 -> I1|W

//   out.clear(); 
//   out.push_back(wait); out.push_back(insert_2);  outgoing[insert_2] = out; // I2 -> I2|W
  

//   // Transition weights are determined by the insertion/deletion rates and the branch_length (set on input)

//   double ins_extend_1 = gap_extend;
//   double del_extend_1 = gap_extend;

//   double ins_extend_2 = gap_extend2;
//   double del_extend_2 = gap_extend2;

//   // double mixPrior is prob of entering a insert2 or del2 state
//   double class1 = 1.0-mixPrior;
//   double class2 = mixPrior;
  
//   double wait2end = 1; // This is a forced transition...should have weight 1.

//   // Probability/weight of an insertion 
//   double ins_open = 1-exp(-ins_open_rate*branch_length);
//   double del_open = 1-exp(-del_open_rate*branch_length);
  
//   // Define their complements, for cleanliness 
//   double no_ins_open = 1-ins_open;
//   double no_del_open = 1-del_open; 
  
//   double no_ins_extend_1 = 1-ins_extend_1;
//   double no_del_extend_1 = 1-del_extend_1;

//   double no_ins_extend_2 = 1-ins_extend_2;
//   double no_del_extend_2 = 1-del_extend_2;
  

//   vector<state> transitionPair;  

//   // Transitions from start 
//   // start -> wait = no insertion
//   transitionPair.push_back(start); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_open;   transitionPair.clear();
//   // start -> insert1 =  insertion, choose 1
//   transitionPair.push_back(start); transitionPair.push_back(insert_1); 
//   transition_weight[transitionPair] =ins_open*class1;   transitionPair.clear();
//   // start -> insert2 =  insertion, choose 2
//   transitionPair.push_back(start); transitionPair.push_back(insert_2); 
//   transition_weight[transitionPair] = ins_open*class2;   transitionPair.clear();


//   // Transitions from match 
//   // match -> wait = no insertion
//   transitionPair.push_back(match); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_open;   transitionPair.clear();
//   // match -> insert1 =  insertion, choose 1
//   transitionPair.push_back(match); transitionPair.push_back(insert_1); 
//   transition_weight[transitionPair] = ins_open*class1;   transitionPair.clear();
//   // match -> insert1 =  insertion, choose 2
//   transitionPair.push_back(match); transitionPair.push_back(insert_2); 
//   transition_weight[transitionPair] = ins_open*class2;   transitionPair.clear();

//   // Transitions from wait 
//   // wait -> match = no deletion
//   transitionPair.push_back(wait); transitionPair.push_back(match); 
//   transition_weight[transitionPair] = no_del_open;   transitionPair.clear();
//   // wait -> delete1 =  deletion, choose 1
//   transitionPair.push_back(wait); transitionPair.push_back(delete_1); 
//   transition_weight[transitionPair] = del_open*class1;   transitionPair.clear();
//   // wait -> delete2 =  deletion, choose 2
//   transitionPair.push_back(wait); transitionPair.push_back(delete_2); 
//   transition_weight[transitionPair] = del_open*class2;   transitionPair.clear();
//   // wait -> end =  forced transition
//   transitionPair.push_back(wait); transitionPair.push_back(end); 
//   transition_weight[transitionPair] = wait2end;   transitionPair.clear();


//   // Transitions from delete_1
//   // delete1 -> wait = no insertion, no delete-extend_1
//   transitionPair.push_back(delete_1); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_del_extend_1*no_ins_open;   transitionPair.clear();
//   // delete1 -> delete_wait1 = delete-extend1
//   transitionPair.push_back(delete_1); transitionPair.push_back(delete_wait_1); 
//   transition_weight[transitionPair] = del_extend_1;   transitionPair.clear();
//   // delete1 -> insert1 = yes insertion, no delete-extend, choose 1
//   transitionPair.push_back(delete_1); transitionPair.push_back(insert_1); 
//   transition_weight[transitionPair] = no_del_extend_1*ins_open*class1;   transitionPair.clear();
//   // delete1 -> insert2 = yes insertion, no delete-extend, choose 2
//   transitionPair.push_back(delete_1); transitionPair.push_back(insert_2); 
//   transition_weight[transitionPair] = no_del_extend_1*ins_open*class2;   transitionPair.clear();

//   // Transitions from delete_2
//   // delete2 -> wait = no insertion, no delete-extend_2
//   transitionPair.push_back(delete_2); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_del_extend_2*no_ins_open;   transitionPair.clear();
//   // delete2 -> delete_wait2 = delete-extend2
//   transitionPair.push_back(delete_2); transitionPair.push_back(delete_wait_2); 
//   transition_weight[transitionPair] = del_extend_2;   transitionPair.clear();
//   // delete2 -> insert1 = yes insertion, no delete-extend, choose 1
//   transitionPair.push_back(delete_2); transitionPair.push_back(insert_1); 
//   transition_weight[transitionPair] = no_del_extend_2*ins_open*class1;   transitionPair.clear();
//   // delete1 -> insert2 = yes insertion, no delete-extend, choose 2
//   transitionPair.push_back(delete_2); transitionPair.push_back(insert_2); 
//   transition_weight[transitionPair] = no_del_extend_2*ins_open*class2;   transitionPair.clear();


//   //Transitions from delete_wait_1
//   // delete_wait_1 -> delete_1
//   transitionPair.push_back(delete_wait_1); transitionPair.push_back(delete_1); 
//   transition_weight[transitionPair] = 1;   transitionPair.clear();

//   //Transitions from delete_wait_2
//   // delete_wait_2 -> delete_2
//   transitionPair.push_back(delete_wait_2); transitionPair.push_back(delete_2); 
//   transition_weight[transitionPair] = 1;   transitionPair.clear();


//   // Transitions from insert1
//   // insert -> wait = no extend insertion
//   transitionPair.push_back(insert_1); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_extend_1;   transitionPair.clear();
//   // insert -> insert =  extend insertion
//   transitionPair.push_back(insert_1); transitionPair.push_back(insert_1); 
//   transition_weight[transitionPair] = ins_extend_1;   transitionPair.clear();

//   // Transitions from insert2
//   // insert -> wait = no extend insertion
//   transitionPair.push_back(insert_2); transitionPair.push_back(wait); 
//   transition_weight[transitionPair] = no_ins_extend_2;   transitionPair.clear();
//   // insert -> insert =  extend insertion
//   transitionPair.push_back(insert_2); transitionPair.push_back(insert_2); 
//   transition_weight[transitionPair] = ins_extend_2;   transitionPair.clear();


//   // Emission weights
//   // Equilibrium over alphabet characters
//   // The states of type I are the only states which have an emission distribution

//   vector<double> equilibrium = rate_matrix.create_prior(); 
//   for (int i=0; i<states.size(); i++)
// 	{
// 	  if (state_types[i] == "I") emission_weight_matrix[i] = equilibrium;
// 	}
// }


// 2 - mixture-of-affine gaps branch transducer                                                                      
BranchTrans::BranchTrans(double branch_length_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix,
                         double ins_open_rate, double del_open_rate, double gap_extend_0, double mix_prior_0, double gap_extend_1, double mix_prior_1)

{
branch_length = branch_length_in;
name = "mixture";

vector<sstring> toks = alphabet_in.tokens();
for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
   alphabet.push_back(string(a->c_str()));

alphabet_size = alphabet.size();

double ins_open = 1-exp(-ins_open_rate*branch_length);
double del_open = 1-exp(-del_open_rate*branch_length);

conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(branch_length);

// State identifiers
int start = 0;
int end = 1;
int wait_end = 2;
int match = 3;
int wait_match = 4;
int delete_0 = 5;
int wait_delete_0 = 6;
int insert_0 = 7;
int delete_1 = 8;
int wait_delete_1 = 9;
int insert_1 = 10;

// State list
states.push_back(0);
states.push_back(1);
states.push_back(2);
states.push_back(3);
states.push_back(4);
states.push_back(5);
states.push_back(6);
states.push_back(7);
states.push_back(8);
states.push_back(9);
states.push_back(10);

// State names
state_names.push_back("start");
state_names.push_back("end");
state_names.push_back("wait_end");
state_names.push_back("match");
state_names.push_back("wait_match");
state_names.push_back("delete_0");
state_names.push_back("wait_delete_0");
state_names.push_back("insert_0");
state_names.push_back("delete_1");
state_names.push_back("wait_delete_1");
state_names.push_back("insert_1");

// State types
state_types.push_back("S");
state_types.push_back("E");
state_types.push_back("W");
state_types.push_back("M");
state_types.push_back("W");
state_types.push_back("D");
state_types.push_back("W");
state_types.push_back("I");
state_types.push_back("D");
state_types.push_back("W");
state_types.push_back("I");

// State transitions
vector<state> out;

// Transitions out of delete_0
out.clear();
out.push_back(insert_0);
out.push_back(insert_1);
out.push_back(wait_match);
out.push_back(wait_delete_0);
out.push_back(wait_end);
outgoing[delete_0] = out;

// Transitions out of delete_1
out.clear();
out.push_back(insert_0);
out.push_back(insert_1);
out.push_back(wait_match);
out.push_back(wait_delete_1);
out.push_back(wait_end);
outgoing[delete_1] = out;

// Transitions out of insert_0
out.clear();
out.push_back(insert_0);
out.push_back(wait_delete_0);
out.push_back(wait_delete_1);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[insert_0] = out;

// Transitions out of insert_1
out.clear();
out.push_back(insert_1);
out.push_back(wait_delete_0);
out.push_back(wait_delete_1);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[insert_1] = out;

// Transitions out of start
out.clear();
out.push_back(insert_0);
out.push_back(wait_delete_0);
out.push_back(insert_1);
out.push_back(wait_delete_1);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[start] = out;

// Transitions out of wait_delete_1
out.clear();
out.push_back(delete_1);
outgoing[wait_delete_1] = out;

// Transitions out of wait_delete_0
out.clear();
out.push_back(delete_0);
outgoing[wait_delete_0] = out;

// Transitions out of wait_match
out.clear();
out.push_back(match);
outgoing[wait_match] = out;

// Transitions out of match
out.clear();
out.push_back(wait_delete_0);
out.push_back(insert_0);
out.push_back(wait_delete_1);
out.push_back(insert_1);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[match] = out;

// Transitions out of wait_end
out.clear();
out.push_back(end);
outgoing[wait_end] = out;

//Transition weights
vector<state> transitionPair;

// Transition weights out of delete_0
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = (1-gap_extend_0)*ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = (1-gap_extend_0)*ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_0)*(1-ins_open);
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = gap_extend_0;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of delete_1
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = (1-gap_extend_1)*ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = (1-gap_extend_1)*ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_1)*(1-ins_open);
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = gap_extend_1;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of insert_0
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = gap_extend_0;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-gap_extend_0)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-gap_extend_0)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_0)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of insert_1
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = gap_extend_1;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-gap_extend_1)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-gap_extend_1)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_1)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of start
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-ins_open)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-ins_open)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-ins_open)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_delete_1
transitionPair.clear();
transitionPair.push_back(wait_delete_1); transitionPair.push_back(delete_1);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_delete_0
transitionPair.clear();
transitionPair.push_back(wait_delete_0); transitionPair.push_back(delete_0);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_match
transitionPair.clear();
transitionPair.push_back(wait_match); transitionPair.push_back(match);
transition_weight[transitionPair] = 1.0;

// Transition weights out of match
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = del_open*(1-ins_open)*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = del_open*(1-ins_open)*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-ins_open)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_end
transitionPair.clear();
transitionPair.push_back(wait_end); transitionPair.push_back(end);
transition_weight[transitionPair] = 1.0;

vector<double> equilibrium = rate_matrix.create_prior();
for (unsigned int i=0; i<states.size(); i++)
    if (state_types[i] == "I") 
        emission_weight_matrix[i] = equilibrium;
}




// 3 - mixture-of-affine gaps branch transducer                                                                      
BranchTrans::BranchTrans(double branch_length_in, Alphabet& alphabet_in, Irrev_EM_matrix& rate_matrix,
                         double ins_open_rate, double del_open_rate, double gap_extend_0, double mix_prior_0, double gap_extend_1, double mix_prior_1, double gap_extend_2, double mix_prior_2)

{
branch_length = branch_length_in;
name = "mixture";

vector<sstring> toks = alphabet_in.tokens();
for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
   alphabet.push_back(string(a->c_str()));

alphabet_size = alphabet.size();

double ins_open = 1-exp(-ins_open_rate*branch_length);
double del_open = 1-exp(-del_open_rate*branch_length);

conditional_sub_matrix = rate_matrix.create_conditional_substitution_matrix(branch_length);

// State identifiers
int start = 0;
int end = 1;
int wait_end = 2;
int match = 3;
int wait_match = 4;
int delete_0 = 5;
int wait_delete_0 = 6;
int insert_0 = 7;
int delete_1 = 8;
int wait_delete_1 = 9;
int insert_1 = 10;
int delete_2 = 11;
int wait_delete_2 = 12;
int insert_2 = 13;

// State list
states.push_back(0);
states.push_back(1);
states.push_back(2);
states.push_back(3);
states.push_back(4);
states.push_back(5);
states.push_back(6);
states.push_back(7);
states.push_back(8);
states.push_back(9);
states.push_back(10);
states.push_back(11);
states.push_back(12);
states.push_back(13);

// State names
state_names.push_back("start");
state_names.push_back("end");
state_names.push_back("wait_end");
state_names.push_back("match");
state_names.push_back("wait_match");
state_names.push_back("delete_0");
state_names.push_back("wait_delete_0");
state_names.push_back("insert_0");
state_names.push_back("delete_1");
state_names.push_back("wait_delete_1");
state_names.push_back("insert_1");
state_names.push_back("delete_2");
state_names.push_back("wait_delete_2");
state_names.push_back("insert_2");

// State types
state_types.push_back("S");
state_types.push_back("E");
state_types.push_back("W");
state_types.push_back("M");
state_types.push_back("W");
state_types.push_back("D");
state_types.push_back("W");
state_types.push_back("I");
state_types.push_back("D");
state_types.push_back("W");
state_types.push_back("I");
state_types.push_back("D");
state_types.push_back("W");
state_types.push_back("I");

// State transitions
vector<state> out;

// Transitions out of insert_2
out.clear();
out.push_back(insert_2);
out.push_back(wait_delete_0);
out.push_back(wait_delete_1);
out.push_back(wait_delete_2);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[insert_2] = out;

// Transitions out of delete_1
out.clear();
out.push_back(insert_0);
out.push_back(insert_1);
out.push_back(insert_2);
out.push_back(wait_match);
out.push_back(wait_delete_1);
out.push_back(wait_end);
outgoing[delete_1] = out;

// Transitions out of insert_0
out.clear();
out.push_back(insert_0);
out.push_back(wait_delete_0);
out.push_back(wait_delete_1);
out.push_back(wait_delete_2);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[insert_0] = out;

// Transitions out of insert_1
out.clear();
out.push_back(insert_1);
out.push_back(wait_delete_0);
out.push_back(wait_delete_1);
out.push_back(wait_delete_2);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[insert_1] = out;

// Transitions out of delete_2
out.clear();
out.push_back(insert_0);
out.push_back(insert_1);
out.push_back(insert_2);
out.push_back(wait_match);
out.push_back(wait_delete_2);
out.push_back(wait_end);
outgoing[delete_2] = out;

// Transitions out of wait_delete_2
out.clear();
out.push_back(delete_2);
outgoing[wait_delete_2] = out;

// Transitions out of start
out.clear();
out.push_back(insert_0);
out.push_back(wait_delete_0);
out.push_back(insert_1);
out.push_back(wait_delete_1);
out.push_back(insert_2);
out.push_back(wait_delete_2);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[start] = out;

// Transitions out of wait_delete_1
out.clear();
out.push_back(delete_1);
outgoing[wait_delete_1] = out;

// Transitions out of wait_delete_0
out.clear();
out.push_back(delete_0);
outgoing[wait_delete_0] = out;

// Transitions out of delete_0
out.clear();
out.push_back(insert_0);
out.push_back(insert_1);
out.push_back(insert_2);
out.push_back(wait_match);
out.push_back(wait_delete_0);
out.push_back(wait_end);
outgoing[delete_0] = out;

// Transitions out of wait_match
out.clear();
out.push_back(match);
outgoing[wait_match] = out;

// Transitions out of match
out.clear();
out.push_back(wait_delete_0);
out.push_back(insert_0);
out.push_back(wait_delete_1);
out.push_back(insert_1);
out.push_back(wait_delete_2);
out.push_back(insert_2);
out.push_back(wait_match);
out.push_back(wait_end);
outgoing[match] = out;

// Transitions out of wait_end
out.clear();
out.push_back(end);
outgoing[wait_end] = out;

//Transition weights
vector<state> transitionPair;

// Transition weights out of insert_2
transitionPair.clear();
transitionPair.push_back(insert_2); transitionPair.push_back(insert_2);
transition_weight[transitionPair] = gap_extend_2;
transitionPair.clear();
transitionPair.push_back(insert_2); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-gap_extend_2)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(insert_2); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-gap_extend_2)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(insert_2); transitionPair.push_back(wait_delete_2);
transition_weight[transitionPair] = (1-gap_extend_2)*del_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(insert_2); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_2)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(insert_2); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of delete_1
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = (1-gap_extend_1)*ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = (1-gap_extend_1)*ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(insert_2);
transition_weight[transitionPair] = (1-gap_extend_1)*ins_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_1)*(1-ins_open);
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = gap_extend_1;
transitionPair.clear();
transitionPair.push_back(delete_1); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of insert_0
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = gap_extend_0;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-gap_extend_0)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-gap_extend_0)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_delete_2);
transition_weight[transitionPair] = (1-gap_extend_0)*del_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_0)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(insert_0); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of insert_1
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = gap_extend_1;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-gap_extend_1)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-gap_extend_1)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_delete_2);
transition_weight[transitionPair] = (1-gap_extend_1)*del_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_1)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(insert_1); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of delete_2
transitionPair.clear();
transitionPair.push_back(delete_2); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = (1-gap_extend_2)*ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(delete_2); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = (1-gap_extend_2)*ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(delete_2); transitionPair.push_back(insert_2);
transition_weight[transitionPair] = (1-gap_extend_2)*ins_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(delete_2); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_2)*(1-ins_open);
transitionPair.clear();
transitionPair.push_back(delete_2); transitionPair.push_back(wait_delete_2);
transition_weight[transitionPair] = gap_extend_2;
transitionPair.clear();
transitionPair.push_back(delete_2); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_delete_2
transitionPair.clear();
transitionPair.push_back(wait_delete_2); transitionPair.push_back(delete_2);
transition_weight[transitionPair] = 1.0;

// Transition weights out of start
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = (1-ins_open)*del_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = (1-ins_open)*del_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(insert_2);
transition_weight[transitionPair] = ins_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_delete_2);
transition_weight[transitionPair] = (1-ins_open)*del_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-ins_open)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(start); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_delete_1
transitionPair.clear();
transitionPair.push_back(wait_delete_1); transitionPair.push_back(delete_1);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_delete_0
transitionPair.clear();
transitionPair.push_back(wait_delete_0); transitionPair.push_back(delete_0);
transition_weight[transitionPair] = 1.0;

// Transition weights out of delete_0
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = (1-gap_extend_0)*ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = (1-gap_extend_0)*ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(insert_2);
transition_weight[transitionPair] = (1-gap_extend_0)*ins_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-gap_extend_0)*(1-ins_open);
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = gap_extend_0;
transitionPair.clear();
transitionPair.push_back(delete_0); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_match
transitionPair.clear();
transitionPair.push_back(wait_match); transitionPair.push_back(match);
transition_weight[transitionPair] = 1.0;

// Transition weights out of match
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_delete_0);
transition_weight[transitionPair] = del_open*(1-ins_open)*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(insert_0);
transition_weight[transitionPair] = ins_open*mix_prior_0;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_delete_1);
transition_weight[transitionPair] = del_open*(1-ins_open)*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(insert_1);
transition_weight[transitionPair] = ins_open*mix_prior_1;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_delete_2);
transition_weight[transitionPair] = del_open*(1-ins_open)*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(insert_2);
transition_weight[transitionPair] = ins_open*mix_prior_2;
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_match);
transition_weight[transitionPair] = (1-ins_open)*(1-del_open);
transitionPair.clear();
transitionPair.push_back(match); transitionPair.push_back(wait_end);
transition_weight[transitionPair] = 1.0;

// Transition weights out of wait_end
transitionPair.clear();
transitionPair.push_back(wait_end); transitionPair.push_back(end);
transition_weight[transitionPair] = 1.0;

vector<double> equilibrium = rate_matrix.create_prior();
for (unsigned int i=0; i<states.size(); i++)
    if (state_types[i] == "I") 
        emission_weight_matrix[i] = equilibrium;
}
