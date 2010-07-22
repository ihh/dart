#include <queue>
#include <stack>
#include "protpal/CompositePath.h"
#include "tree/phylogeny.h"


using namespace std; 


int main(void)
{
  const sstring tree_string = "((chimp:.01, ape:.01)CA:.01, human:.1)root;";
  istringstream tree_input (tree_string);
  PHYLIP_tree tree;
  tree.read (tree_input);
  
  CompositeState testState(tree); 
  std::cout<<"Constructed composite state successfully!\n";
}

// ****** CompositePath code ******
//Constructor - two M_ids, assumed to be living at the root.  
CompositePath::CompositePath(M_id m, M_id mPrime, PHYLIP_tree tree_in)
{
  tree = tree_in; 
  CompositeState s(tree);
  s.get_profile_states(m, tree.root); 
  path.push_back(s);
  
  s.get_profile_states(mPrime, tree.root); 
  path.push_back(s); 
}

void CompositePath::explode(void)
{
  // Expand a composite path which initially consists only of two states.  
  // For each node in depth-first preorder, with the 'right' child coming first,
  // We proceed left-to right in the growing state path.  
  node treeNode = tree.root, expandNode;
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
  


vector<CompositeState> CompositePath::expand(int i, node n)
{
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
      CompositeState c(path[i], n); // create a new state identical to path[i] except at/below node n
      // fill out the states below it, starting with the null state on top
      c.get_profile_states(*mIter, n);
      out.push_back(c); 
    }
  return out; 
}

// ****** CompositeState code ******
// Constructors
CompositeState::CompositeState(PHYLIP_tree& tree_in)
{
  tree = &tree_in; 
  prof_states_known = false; 
  component_states_known = false; 
}

CompositeState::CompositeState(CompositeState& template_state, node n )
{
  Node treeNode; 
  Phylogeny::Node_vector children; 
  Phylogeny::Node_vector::iterator childIter; 
  tree  = template_state.tree; 
  prof_states_known = false; 
  component_states_known = false; 

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


void CompositeState::get_profile_states(M_id start, node n)
{
  // Determine the states of profiles on the tree.  This fills the profile_states map with a state for each node
  // This state is technically related to the transducer on the incoming branch to that node, but for now we are 
  // using the same transducer on each branch

  node treeNode; 
  M_id m; 
  vector<node> children; 

  profile_states[n] = start; 
  for_rooted_children(*tree, n, child)
    children.push_back(*child); 
  profile_states[children[0]] = (*profiles)[children[0]].state2mid[start.left_state];
  profile_states[children[1]] = (*profiles)[children[1]].state2mid[start.right_state];

  for_nodes_post (*tree, n, tree->parent[n], bi)
    {
      const Phylogeny::Branch& b = *bi;
      treeNode = b.second;
      if (profile_states.count(treeNode) < 1)
	{ std::cerr<<"Error - no M_id in current node: " << treeNode<< endl; exit(1);}
      else m = profile_states[treeNode]; 
      if ( tree->is_internal(treeNode))
	{
	  children.clear(); 
	  for_rooted_children(*tree, n, child)
	    children.push_back(*child); 
	  profile_states[children[0]] = (*profiles)[children[0]].state2mid[m.left_state];
	  profile_states[children[1]] = (*profiles)[children[1]].state2mid[m.right_state];
	}
    }
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


      
      

  
