#include <queue>
#include <stack>
#include "protpal/CompositePath.h"
#include "tree/phylogeny.h"


using namespace std; 



// ****** CompositePath code ******
//Constructor - two M_ids, assumed to be living at the root.  
CompositePath::CompositePath(M_id m, M_id mPrime, PHYLIP_tree tree_in, double postProb_in, map<node, AbsorbingTransducer>& profiles_in, QTransducer Q_in, bool logging_in)
{
  logging = logging_in; 
  if (logging)
    std::cerr<<"Creating composite path object...\n";

  Q = Q_in; 
  profiles = &profiles_in; 
  postProb = postProb_in;
  tree = tree_in; 
  if (logging)
    std::cerr<<"Creating first composite state object...\n";
  CompositeState s(tree, profiles, &Q_in);

  if (logging)
    std::cerr<<"Filling state tree for first composite state object...\n";
  s.get_profile_states(m, tree.root); 
  path.push_back(s);
  std::cerr<<"Displaying...\n";
  s.display(std::cerr, Q); 

  s.get_profile_states(mPrime, tree.root); 
  path.push_back(s); 
  s.display(std::cerr, Q); 
}

void CompositePath::get_counts(ExpCount& masterCount)
{
  // Add all of the transition counts from this composite path to a master count map
  ExpCount counts; 
  ExpCount::iterator countIter; 

  int j; 

  explode(); 
  for (unsigned int i=0; i<path.size(); i++)
    {
      j = i+1; 
      counts  = count(i,j); 
      for (countIter = counts.begin(); countIter != counts.end(); countIter++)
	masterCount[countIter->first] = (countIter->second)*postProb; 
    }
}


ExpCount CompositePath::count(int i, int j)
{
  ExpCount counts; 
  Node n; 
  pair<state, state> transitionPair; 
    
  // For n in tree.nodes: 
  //   if a real change is detected from state i -> state j at node n:
  //      add component state transition i->j at node n to ExpCount
  // return ExpCount
  for (n = 0; n<tree.nodes(); n++)
    {
      if (component_changed(i, j, n))
	{
	  transitionPair.first = path[i].component_states[n]; 
	  transitionPair.second = path[j].component_states[n]; 
	  counts[transitionPair] += 1.0;	  
	}
    }
  return counts; 
}


bool CompositePath::component_changed(int i, int j, node n)
{
  return true; 
}

void CompositePath::explode(void)
{
  // Expand a composite path which initially consists only of two states.  
  // For each node in depth-first preorder, with the 'right' child coming first,
  // We proceed left-to right in the growing state path.  
  node treeNode = tree.root; 
  vector<CompositeState> newPath; 
  unsigned int pathIdx; 
  stack<node> nodeStack;
  map<int, vector<CompositeState> > nulls; 
  map<int, vector<CompositeState> >::iterator nullsIter; 

  for_rooted_children(tree, treeNode, child)
    nodeStack.push(*child);  // the 1st child is on top, as we need.  It will be popped off first

  while (!nodeStack.empty())
    {
      treeNode = nodeStack.top(); 
      nodeStack.pop();

      for_rooted_children(tree, treeNode, child)
	nodeStack.push(*child);  
      
      // For each state i in the path, expand state i, i+1 at node treeNode, 
      // keeping the vector of null states for later. 
      nulls.clear(); 
      for (pathIdx = path.size(); pathIdx < path.size()-1; pathIdx++)
	nulls[pathIdx] = expand(pathIdx, treeNode); 
      
      // Now insert each of these vectors of null states into the main path
      // This is done by making a new path from scratch, and plunking on the relevant elements,
      // and eventually exchanging the main path with it.
      newPath.clear();
      for (pathIdx = path.size(); pathIdx < path.size(); pathIdx++)
	{
	  newPath.push_back(path[pathIdx]);
	  for ( pathIter = nulls[pathIdx].begin(); pathIter != nulls[pathIdx].end() ; pathIter++)
	    newPath.push_back(*pathIter); 
	}
      path = newPath; 
    }
}
  

string CompositeState::state_type(node n)
{
  bool debug = false; 
  if (n == tree->root)
    {
      std::cerr<<"Error: state_type of root called - root has no incoming branch\n";
      exit(1); 
    }
  else if (profile_states.count(tree->parent[n]) < 1)
    {
      std::cerr<<"Error: node " << tree->node_name[n] << "'s  parent has no profile \n"; 
      std::cerr<<tree->node_name[tree->parent[n]] << "has no profile \n";
      exit(1); 
    }
  M_id m = profile_states[tree->parent[n]];
  if (debug)
    std::cerr<<"Use "<< tree->node_name[tree->parent[n]]; 

  if (m.q_state == placeholder_wait)
    return "W"; 
  vector<node> children; 
  for_rooted_children(*tree, tree->parent[n], child)
    {
      children.push_back(*child); 
    }
  if (n == children[0])  // left child
    {
      if (debug)
	std::cerr<< " L "; 
      return stringAt(Q->get_state_type(m.q_state),2);
    }
  else if (n == children[1])  // right child
    {
      if (debug)
	std::cerr<< " R "; 
      return stringAt(Q->get_state_type(m.q_state),3);
    }
  else 
    {
      std::cerr<<"Child appears to be an orphan, in state_type lookup...\n";
      exit(1);
    }
}

state CompositeState::node_state(node n)
{
  M_id m = profile_states[tree->parent[n]];
  
  if (m.q_state == placeholder_wait)
    return -1; 
  vector<node> children; 
  for_rooted_children(*tree, tree->parent[n], child)
    children.push_back(*child); 
  if (n == children[0]) // left child
    return Q->get_components(m.q_state)[2];
  else if (n == children[1]) // right child
    return Q->get_components(m.q_state)[3];
}


bool CompositeState::are_synced_types(string parentType, string childType)
{
  if (parentType == "S")
    return childType == "S"; 
  if (parentType == "D")
    return childType == "W"; 
  if (parentType == "E")
    return childType == "E"; 

  if (parentType == "W")
    return childType == "W"; 

  if (parentType == "M")
    return childType == "M" || childType == "D"; 

  if (parentType == "I")
    return childType == "M" || childType == "D"; 
  else
    std::cerr<<"Unknown type: " << parentType <<endl; 
}

bool CompositeState::are_synced_nodes(node parent, node child)
{
  string childType=state_type(child), parentType=state_type(parent); 
  return are_synced_types(parentType, childType); 
}
bool CompositeState::synced_below(node n)
{
  node treeNode; 
  for_nodes_pre (*tree, n, tree->parent[n], bi)
    {
      const Phylogeny::Branch& b = *bi;
      treeNode = b.second;
      for_rooted_children(*tree, treeNode, child)
	{
	  if (!are_synced_nodes(treeNode, *child))
	    {
	      std::cerr<<" *** " << endl; 
	      std::cerr<<"Not synced: " << tree->node_name[treeNode] << " " << tree->node_name[*child]<< endl; 
	      std::cerr<<"Not synced: " << state_type(treeNode) << " " << state_type(*child)<< endl; 
	      std::cerr<<" *** " << endl; 
	      return false; 
	    }
	}
    }
  return true; 
}
	    

node CompositeState::active_node(void)
{
  node treeNode, treeNode2; 
  bool synced; 
  vector<node> children; 
  for_nodes_pre (*tree, tree->root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      treeNode = b.second;
      synced = true; 
      
      if (treeNode == tree->root)
	{
	  string Rtype = stringAt(Q->get_state_type(profile_states[tree->root].q_state),0);
	  for_rooted_children(*tree, tree->root, child)
	    {
	      children.push_back(*child); 
	      if (!are_synced_types(Rtype, state_type(*child)) )
		{
		  synced = false; 
		  std::cerr<<"Root not synced with " << tree->node_name[*child] <<endl; 
		}
	    }
	  synced = synced && synced_below(children[0]) && synced_below(children[1])  && Rtype == "I";
	  if (!synced)
	    {
	      if (!synced_below(children[0]))
		std::cerr<<"Not synced below " << tree->node_name[children[0]] <<endl; 
	      else  if (!synced_below(children[1]))
		std::cerr<<"Not synced below " << tree->node_name[children[1]] <<endl; 
	      else
		std::cerr<<"Not synced for unknown reason! " << endl;
	    }
	      
	      
	    
	}
      else
	synced = synced && synced_below(treeNode) && state_type(treeNode) == "I";

      if (synced)
	return treeNode; 
    }
  std::cerr<<"Warning: No active node found\n";
  return 0; 
}


vector<CompositeState> CompositePath::expand(int i, node n)
{
  bool logging = true; 
  // Assume there is a path of CompositeStates, each with a map from tree nodes to M_ids (profile_states).  
  // Taking the i'th and i+1'th of these, we'd like to expand the state sequence at node n, if necessary.
  //  This amounts to getting the summed null states at node n, between the i and i+1'th composite state,
  // and creating a new composite state for each of them, and returning a vector of these states. 
  vector<CompositeState> out; 
  vector<M_id> nulls; 
  vector<M_id>::iterator mIter; 
  M_id m = path[i].profile_states[n];
  M_id mPrime = path[i+1].profile_states[n]; 
  pair<int, int > transitionPair; 
  transitionPair.first = (*profiles)[n].mid2int[m.toVector()]; 
  transitionPair.second = (*profiles)[n].mid2int[mPrime.toVector()]; 

  // if there are no summed null states (e.g. subtree insertions), we're done!
  if ( (*profiles)[n].summed_nulls.count(transitionPair) <1 )
    return out; 
  else 
    nulls = (*profiles)[n].summed_nulls[transitionPair];

  for (mIter = nulls.begin(); mIter != nulls.end(); mIter++)
    {
      CompositeState c(path[i], n, profiles, &Q); // create a new state identical to path[i] except at/below node n
      // fill out the states below it, starting with the null state on top
      c.get_profile_states(*mIter, n);
      out.push_back(c); 
    }
  return out; 
}

// ****** CompositeState code ******
// Constructors
CompositeState::CompositeState(PHYLIP_tree& tree_in, map<node, AbsorbingTransducer>* profiles_in, QTransducer* Q_in)
{
  profiles = profiles_in; 
  tree = &tree_in; 
  prof_states_known = false; 
  component_states_known = false; 
  placeholder_wait = -100; 
  Q = Q_in; 
}

CompositeState::CompositeState(CompositeState& template_state, node n, map<node, AbsorbingTransducer>* profiles_in, QTransducer* Q_in)
{
  profiles = profiles_in; 
  Node treeNode; 
  Phylogeny::Node_vector children; 
  Phylogeny::Node_vector::iterator childIter; 
  tree  = template_state.tree; 
  prof_states_known = false; 
  component_states_known = false; 
  placeholder_wait =-100; 
  Q = Q_in; 

  // fill profile_states until node n
  queue<Node> nodeQueue; 
  nodeQueue.push(tree->root); 

  while (!nodeQueue.empty())
    {
      treeNode = nodeQueue.front(); 
      nodeQueue.pop();
      if (treeNode != n)
	{
	  profile_states[treeNode] = template_state.profile_states[treeNode]; 
	  children = tree->children(treeNode, tree->parent[treeNode]); 
	  for (childIter = children.begin(); childIter != children.end(); childIter++)
	    nodeQueue.push(*childIter); 
	}
    }
}


void CompositeState::get_profile_states(M_id start, node n, bool logging)
{
  // Determine the states of profiles on the tree.  This fills the profile_states map with a state for each node
  // This state is technically related to the transducer on the incoming branch to that node, but for now we are 
  // using the same transducer on each branch

  // Here for now...clunky I know
  int profile_start = -1;
  int profile_wait = 0;
  int profile_delete = 1;
  int profile_pre_end = 2;
  int profile_end = 3;  

  node treeNode, treeNode2; 
  M_id m, mWait; 
  mWait.left_type = profile_wait; 
  mWait.right_type = profile_wait; 
  mWait.q_state = placeholder_wait; // special placeholder Q wait state...
  mWait.left_state = placeholder_wait; // special placeholder wait state...
  mWait.right_state = placeholder_wait; // special placeholder wait state...

  vector<node> children; 
  stack<node> nodeQueue;
  if (logging)
    std::cerr<<"Beginning profile-state filling..\n"; 
  profile_states[n] = start; 
  
  for_rooted_children(*tree, n, child)
    children.push_back(*child); 
  profile_states[children[0]] = (*profiles)[children[0]].state2mid[start.left_state];
  profile_states[children[1]] = (*profiles)[children[1]].state2mid[start.right_state];

  nodeQueue.push(tree->root); 
  while (!nodeQueue.empty())
    {
      treeNode = nodeQueue.top(); 
      nodeQueue.pop();
      if (logging)
	std::cerr<<"Examining node: "<< tree->node_name[treeNode] << endl; 
      if (tree->is_leaf(treeNode))
	{
	  if (logging)
	    std::cerr<<"Node "<< tree->node_name[treeNode] << " is a leaf, nothing to be done. \n"; 
	  continue; 
	}
      else if (profile_states.count(treeNode) < 1)
	{std::cerr<<"Error - no M_id in current node: " << tree->node_name[treeNode]<< endl; exit(1);}
      else m = profile_states[treeNode]; 

      if ( tree->is_internal(treeNode))
	{
	  children.clear(); 
	  for_rooted_children(*tree, treeNode, child)
	    children.push_back(*child); 

	  if (m.left_type != profile_wait)
	    {
	      profile_states[children[0]] = (*profiles)[children[0]].state2mid[m.left_state];
	      nodeQueue.push(children[0]); 
	    }
	  else
	    {
	      if (logging)
		std::cerr<<"Node "<< tree->node_name[children[0]] << " determined to be of type wait \n";
	      if (logging)
		std::cerr<<"Filling wait states for node..."; 
	      for_nodes_pre (*tree, children[0], tree->parent[children[0]], bi)
		{
		  const Phylogeny::Branch& b = *bi;
		  treeNode2 = b.second;
		  profile_states[treeNode2] = mWait; 
		  if (logging)
		    std::cerr<<" " << tree->node_name[treeNode2] << " "; 
		}
	      std::cerr<<"\n";
	    }

	  if (m.right_type != profile_wait)
	    {
	      profile_states[children[1]] = (*profiles)[children[1]].state2mid[m.right_state];
	      nodeQueue.push(children[1]); 
	    }
	  else
	    {
	      if (logging)
		std::cerr<<"Node "<< tree->node_name[children[1]] << " determined to be of type wait \n";
	      if (logging)
		std::cerr<<"Filling wait states for node..."; 
	      for_nodes_pre (*tree, children[1], tree->parent[children[1]], bi)
		{
		  const Phylogeny::Branch& b = *bi;
		  treeNode2 = b.second;
		  profile_states[treeNode2] = mWait; 
		  if (logging)
		    std::cerr<<" " << tree->node_name[treeNode2] << " "; 
		}
	      std::cerr<<"\n";
	    }
	}

    }
}

      
void CompositeState::display(ostream& out, QTransducer& Q )
{
  vector<node> children; 
  map<int, string> int2state;
  int2state[-1]="start";
  int2state[0]="wait";
  int2state[1]="delete";
  int2state[2]="pre-end";
  int2state[3]="end";

  state q_state; 
  string components; 
  out << "\n"; 

  out << "The active node was: " << tree->node_name[active_node()] <<endl; 
  string Rtype = stringAt(Q.get_state_type(profile_states[tree->root].q_state),0);
  out << "Root has state type: " << Rtype << endl; 
  for (node n=0; n<tree->nodes(); n++)
    {
      if (profile_states.count(n)<1 || tree->is_leaf(n))
	continue;
      else
	{
	  children.clear(); 
	  for_rooted_children(*tree, n, child)
	    children.push_back(*child); 
	  q_state = profile_states[n].q_state; 
	  if (q_state == placeholder_wait)
	    components = "WWWW";
	  else
	    components = Q.get_state_type(q_state); 

	  out << "Branch incoming to node "<< tree->node_name[children[0]] << " has state type: " 
	      << components[2] << " check: " << state_type(children[0]) << endl; 
	  out<< "\tThe actual state: "<< Q.B_r.get_state_name(node_state(children[0])) <<endl; 
	  out << "Branch incoming to node "<< tree->node_name[children[1]] << " has state type: " 
	      << components[3] << " check: " << state_type(children[1]) <<endl; 
	  out<< "\tThe actual state: "<< Q.B_r.get_state_name(node_state(children[1])) <<endl; 
	}
    }
  out << "\n"; 
}

void CompositeState::get_component_states(map<node, BranchTrans> transducers)
{
//   map<node, state>::iterator profStateIter; 
//   vector<node> children;
//   node treeNode; 
//   M_id profileState; 
//   state profileStateIdx, Qstate, Lstate, Rstate; 
//   for (profStateIter = profile_states.begin(); profStateIter != profile_states.end(); profStateIter++)
//     {
//       treeNode = profStateIter->first; 
//       profileStateIdx = profStateIter->second; 
//       profileState = (*profiles)[treeNode].state2mid[profileStateIdx];
//       Lstate = profileState.left_state; 
//       Rstate = profileState.right_state; 
      
//       // Same caveat as in set_node
//       for_rooted_children(tree, treeNode, child)
// 	children.push_back(*child);
//       component_states[children[0]] = Lstate; 
//       component_states[children[1]] = Rstate; 
//     }
}


      
      

  

//   for_nodes_pre (*tree, n, tree->parent[n], bi)
//     {
//       const Phylogeny::Branch& b = *bi;
//       treeNode = b.second;
//       if (logging)
// 	std::cerr<<"Considering node: " << tree->node_name[treeNode] <<endl; 
//       if (tree->is_leaf(treeNode))
// 	continue; 
//       else if (profile_states.count(treeNode) < 1)
// 	{std::cerr<<"Error - no M_id in current node: " << tree->node_name[treeNode]<< endl; exit(1);}
//       else m = profile_states[treeNode]; 
//       if ( tree->is_internal(treeNode))
// 	{
// 	  children.clear(); 
// 	  for_rooted_children(*tree, n, child)
// 	    children.push_back(*child); 
// 	  profile_states[children[0]] = (*profiles)[children[0]].state2mid[m.left_state];
// 	  profile_states[children[1]] = (*profiles)[children[1]].state2mid[m.right_state];
// 	}
