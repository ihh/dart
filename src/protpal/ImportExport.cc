#include<iostream>
#include<math.h>
#include "protpal/ImportExport.h"
#include "util/sexpr.h"


/// Writing functions 
void AbsorbingTransducer::write_profile(ostream& out, state_path& viterbi_path)
{
  // Writes a profile 
  out.setf(ios::fixed, ios::floatfield);
  out.setf(ios::showpoint);
  out.precision(2);
  out.precision(10); 
  // write transducer to file in a compact format
  pair<state, state> statePair; 

  // Initialize the profile def
  out << "(profile\n";
  add_tag_value_pair(out, "node", treeNode); 
  out << "\n"; 
  // Add the special states - their names are indices,  their types are non-delete
  out << ";; Start/pre-end/end state data\n"; 
  add_basic_state(out, "start", start_state); 
  add_basic_state(out, "wait", pre_end_state); 
  add_basic_state(out, "end", end_state); 
  // Loop over delete states and add them
  out << ";; Delete states' data\n\n"; 
  for (int stateIdx = 0; stateIdx < num_delete_states; stateIdx++)
    // add_delete states handles everything having to do with an individual state, like:
    //   - absorb distribution
    //   - sequence co-ordinates of leaf sequences (for alignment envelope)
    add_delete_state(out, stateIdx); 

  // Loop over transitions and add them
  out << ";; Transition weight data\n\n"; 
  for (MyMap<state, vector<state> >::iterator e = incoming.begin(); e != incoming.end(); e++)
    for (vector<state>::iterator ePrime = e->second.begin(); ePrime != e->second.end(); ePrime++)
      add_transition(out, *ePrime, e->first, get_transition_weight( *ePrime, e->first)); 

  // If viterbi path is non-null, add it after converting each M_id to an integer
  if (viterbi_path.size())
    {
      out << "(viterbi_path "; 
      for (state_path::reverse_iterator vit=viterbi_path.rbegin(); vit!=viterbi_path.rend(); vit++)
	if (mid2int.count(vit->toVector()))
	  out << mid2int[vit->toVector()] << " "; 
      out <<")\n";
    }
  out << ");; end profile for node " << treeNode << "\n";
}

void AbsorbingTransducer::add_tag_value_pair(ostream& out, string tag, string value, bool newline)
{
  out << "(" << tag << " " << value << ")"; 
  if (newline)
    out << "\n"; 
  else
    out << " "; 
}

void AbsorbingTransducer::add_tag_value_pair(ostream& out, string tag, bfloat value, bool newline)
{
  out << "(" << tag << " " << value << ")"; 
  if (newline)
    out << "\n"; 
  else
    out << " "; 

}

void AbsorbingTransducer::add_tag_value_pair(ostream& out, string tag, int value, bool newline)
{
  out << "(" << tag << " " << value << ")"; 
  if (newline)
    out << "\n"; 
  else
    out << " "; 

}

void AbsorbingTransducer::add_tag_value_pair(ostream& out, string tag, double value, bool newline)
{
  out << "(" << tag << " " << value << ")"; 
  if (newline)
    out << "\n"; 
  else
    out << " "; 
}


void AbsorbingTransducer::add_basic_state(ostream& out, string type, int stateIndex)
{
    out << "(state\n"; 
    add_tag_value_pair(out, "name", stateIndex); 
    add_tag_value_pair(out, "type", type); 
    out << ");; end state " << name <<"\n\n";   
}

void AbsorbingTransducer::add_delete_state(ostream& out, int stateIndex)
{
  bool verboseAbsorb = false; 
  // Initialize
  out << "(state \n";   
  add_tag_value_pair(out, "type", "delete"); 
  add_tag_value_pair(out, "name", stateIndex); 
  // add_tag_value_pair(out, "postprob", stored_sampled_profile->post_prob(state2mid[stateIndex])); 

  // Absorption weights
  if (verboseAbsorb)
    {
      for (int charIndex = 0; charIndex < alphabet_size; charIndex++)
	{
	  out << "(absorb "; 
	  add_tag_value_pair(out, "character", alphabet[charIndex], false); 
	  add_tag_value_pair(out, "weight", double(log(get_absorb_weight(stateIndex, charIndex))), false); 
	}
      out << ")\n";
      out << ";; end absorption weights\n"; 
    }
  else
    {
      out << "(absorb "; 
    for (int charIndex = 0; charIndex < alphabet_size; charIndex++)
      out << double(log(get_absorb_weight(stateIndex, charIndex))) << " "; 
    out << ");; end absorb weight info\n\n"; 
    }

  // Leaf indices
  out << "(leaf_indices (" << leaf_coords[stateIndex].first << " " << leaf_coords[stateIndex].second; 
  out << "))\n"; 

  // Leaf sequence coordinates
  out << "(leaf_seq_coords \n"; 
  for (MyMap<node, int>::iterator nodeIter = leaf_seq_coords[stateIndex].begin(); nodeIter != leaf_seq_coords[stateIndex].end(); nodeIter++)
    {
      out << "( "; 
      add_tag_value_pair(out, "node", nodeIter->first, false); 
      add_tag_value_pair(out, "index", nodeIter->second, false); 
      out << ")\n"; 
    }
  out << ") ;; end leaf seq coords\n";

  // Terminate
  out << ");;end delete state " << stateIndex << "\n\n";
}

void AbsorbingTransducer::add_transition(ostream& out, int fromState, int toState, bfloat weight)
{
  out << "(transition "; 
  add_tag_value_pair(out, "from", fromState, false); 
  add_tag_value_pair(out, "to", toState, false); 
  //  add_tag_value_pair(out, "weight", log(weight), false); 
  out << "(weight " << double(log(weight)) <<")"; 
//   cerr<< "Added transition " << weight << " , stored as " << double(log(weight)) << endl; 
//   cerr<< "Added transition " << fromState << " " << toState << " " << weight << " , stored as bfloat " << log(weight) << endl; 
  out << ")\n"; 
}




/// Reading functions 

void AbsorbingTransducer::read_profile(const char* profile_filename, bool logging)
{
  SExpr S, P;
  vector<SExpr*> sVect;

  SExpr_file fh(profile_filename);
  fh.parse_text();
  S.init(fh.text.begin(), fh.text.end());  
  P = S.find_or_die("profile");

  // get tree node that this profile is living at
  treeNode = int(P.find_or_die("node").value().get_atom().to_int());
  if (logging)  
    cerr<<"Treenode deduced as: " << treeNode <<endl; 
  // Find all states in the profile, loop through and parse them in 
  sVect = P.find_all("state"); 
  if (logging)  
    cerr<<"Imported profile has " << sVect.size() << " states" << endl; 
  num_delete_states = sVect.size()-3; // start, pre-end, end are not delete states
  for (vector<SExpr*>::iterator sIter = sVect.begin(); sIter!=sVect.end(); ++sIter)
    {
      if (logging)
	cerr<<"Parsing state...\n"; 
      parse_state(*sIter); // pointer-to-pointer, I think. 
    }

  sVect = P.find_all("transition"); 
  for (vector<SExpr*>::iterator sIter = sVect.begin(); sIter!=sVect.end(); ++sIter)
    {
      if (logging)
	cerr<<"Parsing transition...\n"; 
      parse_transition(*sIter); // pointer-to-pointer, I think. 
    }
  if (logging)
    cerr<<"Parsed the whole profile, done.\n"; 
}


void AbsorbingTransducer::parse_state(SExpr* s)
{
  SExpr S; 
  bfloat weight; 
  try{
    sstring stateType = s->find_or_die("type").value().get_atom();
    int stateIndex = s->find_or_die("name").value().get_atom().to_int();

    if (stateType != "delete")
      {
	if (stateType == "start" )
	  start_state = stateIndex; 
	else if (stateType == "end")
	  end_state = stateIndex; 
	else if (stateType == "wait")
	  pre_end_state = stateIndex; 
      }
    else
      {
	delete_states.push_back(stateIndex); 
	// parse a delete state - in order of hairyness
	// Leaf coords
	vector<sstring> leafCoords = s->find_or_die("leaf_indices").value().atoms_to_strings();
	if (leafCoords.size()  != 2)
	  cerr<<"Error: leafCoords is not length 2!\n";
	else
	  {
	    leaf_coords[stateIndex].first = leafCoords[0].to_int();
	    leaf_coords[stateIndex].second = leafCoords[1].to_int(); 
	  }
	
	// Leaf seq coords
	S = s->find_or_die("leaf_seq_coords");
	int node, seqIndex; 
	for (list<SExpr>::iterator sIter = S.child.begin(); sIter != S.child.end(); ++sIter)
	  {
	    if (sIter == S.child.begin()) 
	      continue;
	    node = sIter->find_or_die("node").value().get_atom().to_int(); 
	    seqIndex = sIter->find_or_die("index").value().get_atom().to_int(); 
	    // cerr<<"In state " << stateIndex << ", node " << node << " had index " << seqIndex << endl; 
	    leaf_seq_coords[stateIndex][node] = seqIndex; 
	  }
	
	// Absorption weight
	S = s->find_or_die("absorb");
	vector<sstring> absorbValues = S.atoms_to_strings();
	if (int(absorbValues.size()) != alphabet_size+1)
	  {
	    cerr<<"ERROR: Imported absorb weights are not the same length as alphabet!\n"; 
	    cerr << "Absorb size:" << absorbValues.size() << " alphabet size: " << alphabet_size+1 << endl;
	    exit(1);
	  }
	  
	for (vector<sstring>::iterator absVal = absorbValues.begin()+1; absVal != absorbValues.end(); absVal++)
	  {
	    weight = bfloat_doubleexp(absVal->to_double());
	    // cerr<< "State " << stateIndex << " absorbs weight, character: " << weight << endl; 
	    absorption_weight[stateIndex].push_back(weight);
	  }
      }
  }
  catch (const Dart_exception& e)
    {
      cerr << "ERROR in parsing state: " << e.what();
    }
}

void AbsorbingTransducer::parse_transition(SExpr* s)
{
  int to, from;
  pair<int, int> transitionPair; 
  double weight;
  bfloat bweight;
  to = s->find_or_die("to").value().get_atom().to_int();  
  from = s->find_or_die("from").value().get_atom().to_int();  
  weight = s->find_or_die("weight").value().get_atom().to_double();  
//   if (isinf(log(bfloat_doubleexp(weight))))
//     {
//       cerr << "Problematic transition weight: " << weight << endl;
//       THROWEXPR("ERROR in reading profile from file: exp(Log(Transition)) is Inf.");
//     }
  bweight = bfloat_doubleexp(weight);// exp(weight);  OW changed 10-27-11 after underflows
  if (!bfloat_is_nonzero(bweight))
    {
      cerr << "Problematic transition weight: " << weight << endl;
      bfloat b = bfloat_doubleexp(weight);
      cerr << b << " bfloat" << endl;
      THROWEXPR("ERROR in reading profile from file: Transition is zero.");
    }
  //  cerr<< "Parsed the transition: from " << from << " to "  << to << " with weight " << " " << bweight << endl; 

  // Set transition weight
  transitionPair.first = from; 
  transitionPair.second = to; 
  transition_weight[transitionPair] = bweight; 

  // Add to incoming map
  incoming[to].push_back(from); 
  
}
void AbsorbingTransducer::verify(bool logging)
{
  //Not yet implemented
  cerr<< "Verification not yet implemented, but important!\n"; 
}
AbsorbingTransducer::AbsorbingTransducer(const char* profile_filename, vector<string> alphabet_in, PHYLIP_tree& tree)
{
  // Initialize an absorbing transducer from a profile stored on disk 

  /*
    Note to self, data that needs to be set:
x    treenode
x    subtreenodes (can use tree)
x    alphabet (from outside)
x    alph.size (above)
x    num_del_states
x    incoming 
x    start_state
x    pre_end_state
x    end_state
x    absorption_weight
    state_type_phylo (optional) (not req'd for reconstruction, but is req'd for alignment...)
x    leaf_coords
x    leaf_seq_coords
x    transition_weight
  */
  alphabet = alphabet_in;
  alphabet_size = alphabet.size();

  incoming.clear(); 
  absorption_weight.clear();
  transition_weight.clear();

  // Read in the profile from a file, using DART's SExpr parser.
  read_profile(profile_filename); 

  // Fill in some extra info about the transducer's position in the tree
  subtreeNodes.clear(); 
  node n;
  for_nodes_pre (tree, treeNode, tree.parent[treeNode], bi)
    {
      const Phylogeny::Branch& b = *bi;
      n = b.second;
      subtreeNodes.push_back(n);
    }  
}


bool AbsorbingTransducer::test_equality(AbsorbingTransducer& other, bool logging, bool drillDown)
{
  // Verify equality of two absorbing transducers
  /*
    x    treenode
    x    subtreenodes (can use tree)
    x    alphabet (from outside)
    x    alph.size (above)
    x    num_del_states
    x    incoming 
    x    start_state
    x    pre_end_state
    x    end_state
    x    absorption_weight
    state_type_phylo (optional) (not req'd for reconstruction, but is req'd for alignment...)
    x    leaf_coords
    x    leaf_seq_coords
    x    transition_weight
  */
  if (logging and drillDown)
    cerr << "Drilling down into transducer differences...\n"; 
  bool same = true; 
  if (transition_weight != other.transition_weight)
    {
      //      if (logging)
      //	cerr<< "The two transducers differ in their  transitions' weights \n"; 
      //      same =  false; 
      if (drillDown)
	{
	  for (map< pair<state, state>, bfloat>::iterator transIter=transition_weight.begin(); transIter!=transition_weight.end(); ++transIter)
	    {
	      if (!other.transition_weight.count(transIter->first))
		cerr<<"Current profile has key that other doesn't: " << transIter->first; 
	      
	      else 
		{
		  if (other.transition_weight[transIter->first] != transIter->second)
		    {
		      cerr<<"Transition from " << transIter->first.first << " to " << transIter->first.second << " in other, current: ";
		      cerr << other.transition_weight[transIter->first]  << " " << transIter->second << endl; 
		    }
		}
	    }
	}
    }
  //  if (absorption_weight != other.absorption_weight)
  //    {
  //      if (logging)
  //	cerr<< "The two transducers differ in their absorption weights \n"; 
  //      same =  false; 
  //    }

  if (start_state != other.start_state)
    {
      if (logging)
	cerr<< "The two transducers differ in their start state \n"; 
      same =  false; 
    }


  if (treeNode != other.treeNode)
    {
      if (logging)
	cerr<< "The two transducers differ in their tree node\n"; 
      same =  false; 
    }

  if (num_delete_states != other.num_delete_states)
    {
      if (logging)
	cerr<< "The two transducers differ in their delete state counts\n"; 
      same =  false; 
    }

  if (incoming != other.incoming)
    {
      if (logging)
	cerr<< "The two transducers differ in their incoming transitions \n"; 
      same =  false; 
    }
  return same; 
}


