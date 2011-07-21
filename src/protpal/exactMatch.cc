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


ExactMatch::ExactMatch(string &sequence, node treeNode_in, Alphabet& alphabet_in, bool codon_model)
{
  treeNode = treeNode_in; 
  alphabet_size = alphabet_in.size();
  // for alignment envelope 
  MyMap<node, int> coordMap; 
  if ( codon_model )
    {
      if ( sequence.size()%3 )
	{
	  cerr<<"\n\tWarning: adding buffer to end of sequence when using codon model: " << rep(3-sequence.size()%3, "N" ) << endl; 
	  cerr<<"\t(If you are using a codon model, you may want to ensure that all your sequences are of length a multiple of 3)\n"; 
	  sequence += rep(3-sequence.size()%3, "N" );
	}
    }

  /*
  cout<<"Building Alpahbet\n"; 
  Alphabet a("test", 3); 
  a.init_chars("ACG");
  //  a.init_chars("C");
  //  a.init_chars("G");

  cout<<"Alpahbet init'd\n";   
  vector<sstring> toks = a.tokens(); 
  cout<<"Itering thru Alpahbet\n"; 
  for (vector<sstring>::iterator t=toks.begin(); t!=toks.end(); t++)
    cout<<*t<<"\t"; 
  */

  vector<string> seqVector; 
  if ( codon_model)
    treeNode = treeNode; 
    // Nothing
  else
    for (unsigned int i=0; i<sequence.size(); ++i)
      seqVector.push_back(stringAt(sequence, i)); 
  num_delete_states = seqVector.size();

  // Initialize lots of bookkeeping/special state stuff. 
  start_state = -1; 
  coordMap[treeNode] = -1;
  leaf_seq_coords[start_state] = coordMap; 
  pre_end_state = num_delete_states;
  end_state = num_delete_states + 1;

  // The states start, pre-end, and end are of 0-type
  state_type_phylogeny[start_state][treeNode] = "";
  state_type_phylogeny[pre_end_state][treeNode] = "";
  state_type_phylogeny[end_state][treeNode] = "";      

  // Initialize psuedo-binary absorb matrix, adjacency matrix
  vector<int> to; 
  vector<int> from;
  pair<int, int> pairIdx;

  // Initialize the fairly-trivial state structure, along with other necessary variables. 
  // To change
  for ( int i=0; i<num_delete_states; ++i)
	{
	  //	  sequence[i] = tolower(sequence[i]);
	  to.clear();from.clear();
	  from.push_back(i-1);incoming[i] = from;
	  // This state is of type non-wait
	  state_type_phylogeny[i][treeNode] = seqVector[i];
	  
	  // This state accounts for position i of the sequence, this is the max and min
	  pairIdx.first=i; pairIdx.second=i; 
	  leaf_coords[i] = pairIdx; 
	  // for alignment envelope 
	  coordMap[treeNode] = i;
	  leaf_seq_coords[i] = coordMap; 
	}  
  // "Better" way of initializing Felsenstein absorb matrix. OW 7/20/2011
  // This would be  changed for codons...
  Weight_profile seq_profile = alphabet_in.new_seq2weight(sequence);
  if ( num_delete_states != int(seq_profile.size()))
    {
      cerr<< "\nError: number of delete states does not match number of sequence profiles!\n";
      exit(1); 
    }
  for (int i=0; i<num_delete_states; ++i)
    {
      for (int charIndex = 0; charIndex < alphabet_size; charIndex++)
	if ( seq_profile[i].count(charIndex) )
	  absorb[i].push_back( seq_profile[i][charIndex] );
	else
	  //Small but nonzero weight for other chars
	  absorb[i].push_back(1e-10);
    }

  /*  Eventually this might handle codon alphabets...
  for (int i=0; i<num_delete_states; ++i)
    {
      Symbol_weight_map swm = alphabet_in.char2weight(sequence[i]);//seqVector[i]); 
      cerr<<"Char " << seqVector[i].c_str() << " has weights: \n"; 
      for (Symbol_weight_map::iterator wIter = swm.begin(); wIter != swm.end(); wIter++)
	cerr<<"\t" << alphabet_in.int2char(wIter->first) << "\t" << wIter->second << endl; 
    }
  exit(0); 
  */
	
      


  // for alignment envelope 
  coordMap[treeNode] = num_delete_states-1;
  leaf_seq_coords[pre_end_state] = coordMap; 
  
  // Add to/from entries for the WAIT state and END state
  to.clear(); from.clear();
  //  to.push_back(num_delete_states+1); outgoing.push_back(to);
  from.push_back(num_delete_states-1);incoming[pre_end_state] = from;

  to.clear(); from.clear();
  //  to.push_back(num_delete_states+2); outgoing.push_back(to);
  from.push_back(pre_end_state);incoming[end_state] = from;   
  //  std::cerr<<"Created exact match trans from " << sequence << " having " << num_delete_states << " states \n";
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
  
