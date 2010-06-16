#include<iostream>
#include <ctype.h>
#include<set>
#include<map>
#include<vector>
#include<string>

#include "protpal/exactMatch.h"
#include "protpal/utils.h"
#include "util/sstring.h"

using namespace std;


ExactMatch::ExactMatch(string &sequence, node treeNode_in, Alphabet& alphabet_in)
{
  treeNode = treeNode_in; 
  vector<string> single_alphabet;  
  vector<sstring> toks = alphabet_in.tokens(); 
  for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
	alphabet.push_back(string(a->c_str()));
  
  bool hidden = false; 
  if (alphabet[0].size() > 1)
	{
	  std::cerr<<"Using hidden-state alphabet\n";
	  hidden = true; 
	  for (vector<string>::iterator a=alphabet.begin(); a!=alphabet.end(); a++)
		single_alphabet.push_back(stringAt(*a,0));
	}
  else
	single_alphabet = alphabet; // ugh, fix this soon, using DART's alphabet machinery

  alphabet_size = alphabet.size();
  float aSize = alphabet_size;

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
	  // We CAN now handle hidden alphabets, albiet in a primitive sort of way - if the alphabet is seen to have
	  // more than 1 character per token, the 2nd character is assumed to be the 'label' and the first is what matches to 
	  // the sequence
	  vector<double> delta; 
  	  if (index(stringAt(sequence,i), single_alphabet) != -1)
		{
		  for (int charIndex = 0; charIndex < alphabet_size; charIndex++)
			{
			  if (hidden)
				{
				  if (sequence[i] == alphabet[charIndex][0]) absorb[i].push_back(1.0);
				  else absorb[i].push_back(0.0);
				}

			  else
				{
				  if (stringAt(sequence,i) == alphabet[charIndex] ) absorb[i].push_back(1.0);
				  else absorb[i].push_back(0.0);
				}
			}
		}
	  else
		{
		  std::cerr<<"Warning: the character "<<sequence[i]<<" is not in alphabet.  Setting weight for all characters equal\n";
		  for (int charIndex = 0; charIndex!=alphabet_size; charIndex++)
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
		  std::cerr<<"\t Incoming transitions: "; displayVector(incoming[s]);
		  std::cerr<<"\t Absorption profile: "; displayVector(absorb[s]);
		  std::cerr<<"\t sum of absorptions: " << sum(absorb[s]) << endl;
		  
		  
		}
	  else
		{
		  std::cerr<<s<<" is an ending state\n";
		  std::cerr<<"\t Incoming transitions: "; displayVector(incoming[s]);
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
  
