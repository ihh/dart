#include<iostream>
#include "utils.h"
#include <ctype.h>
#include<set>
#include<map>
#include<vector>
#include<string>
#include "exactMatch.h"

using namespace std;


ExactMatch::ExactMatch(string &sequence, node treeNode_in, Alphabet& alphabet_in)
{
  treeNode = treeNode_in; 

  alphabet = string(alphabet_in.nondegenerate_chars());
  
  alphabet_size = alphabet.size();
  float aSize = alphabet_size;

  // NB sequence is assumed in lower-case
  num_delete_states = sequence.length();
  start_state = -1; 
  pre_end_state = num_delete_states;
  end_state = num_delete_states + 1;

  // The states start, pre-end, and end are of 0-type
  state_type_phylogeny[start_state][treeNode] = "";
  state_type_phylogeny[pre_end_state][treeNode] = "";
  state_type_phylogeny[end_state][treeNode] = "";      


  // Initialize binary absorb matrix, adjacency matrix
  vector<int> to; 
  vector<int> from;
  pair<int, int> pairIdx;
  for (int i=0; i<sequence.size(); i++)
	{
	  sequence[i] = tolower(sequence[i]);
	  to.clear();from.clear();
	  //	  to.push_back(i+1); outgoing.push_back(to); // no need to store outgoing transitions (right?)
	  from.push_back(i-1);incoming[i] = from;

	  // This state is of type non-wait
	  state_type_phylogeny[i][treeNode] = stringAt(sequence, i);
	  
	  // This state accounts for position i of the sequence, this is the max and min
	  pairIdx.first=i; pairIdx.second=i; 
	  leaf_coords[i] = pairIdx; 
	  
	  // eventually this ought to be able to handle degenerate characters...for now any non-alphabet symbol
	  // assigned equal weight across the alphabet in the absorb map. 
	  vector<double> delta; 	  
	  if (index(stringAt(sequence,i), alphabet) != -1)
		{
		  for (int charIndex = 0; charIndex < alphabet_size; charIndex++)
			{
			  if (sequence.at(i) == alphabet[charIndex] ) absorb[i].push_back(1);//delta.push_back(1);
			  else absorb[i].push_back(0);//delta.push_back(0);
			}
		  //absorb[i]  = delta;
		}
	  else
		{
		  std::cerr<<"Warning: the character "<<sequence[i]<<" is not in alphabet.  Setting weight for all characters equal\n";
		  for (string::const_iterator omega = alphabet.begin(); omega != alphabet.end(); omega++)
			{
			  delta.push_back(1/aSize);
			}
		  absorb[i]  = delta;
		}
	}

  // Add to/from entries for the WAIT state and END state
  to.clear(); from.clear();
  //  to.push_back(num_delete_states+1); outgoing.push_back(to);
  from.push_back(num_delete_states-1);incoming[pre_end_state] = from;

  to.clear(); from.clear();
  //  to.push_back(num_delete_states+2); outgoing.push_back(to);
  from.push_back(pre_end_state);incoming[end_state] = from;   
}

void ExactMatch::test(void)
{
  for (state s=-1;s<num_delete_states+2; s++)
	{
	  if (s == start_state) std::cerr<<"Start state has no incoming transitions or absorptions\n";
	  else if(s>-1 && s<num_delete_states) 
		{
		  std::cerr<<s<<" is a delete state\n";
		  std::cerr<<"\t Incoming transitions: "; displayVector(incoming.at(s));
		  std::cerr<<"\t Absorption profile: "; displayVector(absorb.at(s));
		  std::cerr<<"\t sum of absorptions: " << sum(absorb.at(s)) << endl;
		  
		  
		}
	  else
		{
		  std::cerr<<s<<" is an ending state\n";
		  std::cerr<<"\t Incoming transitions: "; displayVector(incoming.at(s));
		}
	}
}


// double ExactMatch::absorb_weight(int m, string seqChar)
// {
//   if (m >= num_delete_states) 
// 	{
// 	  std::cout<<"Warning: calling absorption on non-delete state in exact-match transducer . Returning 1\n"; 
// 	  return(1);
// 	}
//   else if(index(seqChar, alphabet) == -1){
// 	std::cout<<"Alphabet character "<<seqChar<<" not found in alphabet!\nExiting...\n";
// 	exit(1);
//   }
//   else
// 	return absorb.at(m).at(index(seqChar, alphabet));
// 	}

// double ExactMatch::absorb_weight(int m, int seqCharIdx)
// {
//   if (m >= num_delete_states) 
// 	{
// 	  std::cout<<"Warning: calling absorption on non-delete state in exact-match transducer . Returning 1\n"; return(1);
// 	}
//   if (seqCharIdx >= alphabet.size()) 
// 	{
// 	  std::cout<<"Character index not within alphabet size\n. Exiting...\n"; 
// 	  exit(1);
// 	}
  
//   else return absorb.at(m).at(seqCharIdx);
// 	}

// double ExactMatch::transition_weight(int m, int mPrime)
// {
//   if (m > num_delete_states+1 || mPrime > num_delete_states+1) 
// 	{
// 	  std::cout<<"Invalid transition called in exact-match transducer. One of "<<m<<" "<<mPrime<<" is out of the allowed range. \nExiting... \n";
// 	  exit(1);
// 	}
//   if (mPrime == m+1) return 1;
//   else return(0);
// }

// vector<state> ExactMatch::get_outgoing(state s)
// {
//   return(outgoing.at(s));
// }

// vector<state> ExactMatch::get_incoming(state s)
// {
//   return(outgoing.at(s));
// }
  
