#include <queue>
#include <vector>
#include <list>
#include <stack>

#include "protpal/utils.h"
#include "protpal/algebras.h"
#include "protpal/profile.h"
#include "protpal/AlignmentSampler.h"

#include "tree/phylogeny.h"
#include "seq/alignment.h"
#include "util/sstring.h"

using namespace std; 
AlignmentSampler::AlignmentSampler(void)
{
  // Placeholder constructor
}
AlignmentSampler::AlignmentSampler(state_path path_in, node node_in, map<node, Profile>* profiles_in, map<node, AbsorbingTransducer>* AbsProfiles_in, PHYLIP_tree* tree_in)
{

  path = path_in; 
  treeNode = node_in; 
  profiles = profiles_in; 
  AbsProfiles = AbsProfiles_in; 
  tree = tree_in; 
  for (unsigned int i=0; i<path.size(); i++)
    columns.push_back( (*profiles)[treeNode].state_type_phylogeny[path[i].toVector()] );
  profile_start = -1;
  profile_delete = 1;
  profile_pre_end = 2;
}

// sample all subtree alignments, filling columns along the way
void  AlignmentSampler::sample_all(bool viterbi, bool logging_in)
{
  logging = logging_in; 
  have_restored_nulls = true; 
  M_id root_state, source; 
  vector<node> children; 
  string side; 
  colMap subColumns; 
  bool breakout = false; ; 
  Profile* profile = &(*profiles)[treeNode];
  //get children of root
  for_rooted_children(*tree, tree->root, child)
    children.push_back(*child); 
  if (logging)
    std::cerr<<"Beginning recursive sampling at root...\n"; 
  for (unsigned int colIdx = 0; colIdx< columns.size(); colIdx++)
    {
      if (breakout)
	break; 
      root_state = path[colIdx];
      if (root_state.left_type == profile_pre_end && root_state.right_type == profile_pre_end)
	{
	  if (logging)
	    std::cerr<<"Pre-end state reached, ending sampling after this\n"; 
	  breakout = true; 
	}
      if (logging)
	{
	  std::cerr<<"sampling at root from position " << colIdx << "...Root level state:\n";
	  root_state.display(profile->Q); 
	}

      // right side
      if (tree->is_leaf(children[1]))
	{
	  if (logging)
	    std::cerr<<"Not sampling at node " << tree->node_name[children[1]] << ", it is a leaf\n"; 
	}
      else if (root_state.right_type == profile_delete || root_state.right_type == profile_pre_end )
	{
	  side = "right"; 
	  source = find_last_del_state(path, colIdx, side); 
	  subColumns = sample_expanded_path(children[1], source.right_state, root_state.right_state); 
	  add_cols(colIdx, subColumns); 
	}
      
      if (tree->is_leaf(children[0]))
	{
	  if (logging)
	    std::cerr<<"Not sampling at node " << tree->node_name[children[0]] << ", it is a leaf\n"; 
	}
      else if (root_state.left_type == profile_delete || root_state.left_type == profile_pre_end )
	{
	  // left side
	  side = "left"; 
	  source = find_last_del_state(path, colIdx, side); 
	  subColumns = sample_expanded_path(children[0], source.left_state, root_state.left_state); 
	  add_cols(colIdx, subColumns);

	}
    }
}
// Display the alignment in stockholm form
void AlignmentSampler::display(ostream& out, string format)
{
  if (columns.size() == 0)
    return; 
  bool logging = false; 
  map<node, string> alignment; 
  map< node, string>::iterator nodeIter; 
  if (logging)
    {
      std::cerr<<"Displaying " << columns.size() << " alignment columns\n"; 
      std::cerr<<"Column 0 has size: " << columns[0].size() <<endl; 
    }

  if (logging)
    std::cerr<<"Building alignment...\n"; 
  for (colIter = columns.begin(); colIter != columns.end(); colIter++)
    {
      std::cerr<<"Current column has size: " << colIter->size() <<endl; 
      for ( nodeIter = colIter->begin(); nodeIter != colIter->end(); nodeIter++)
	{
	  if (logging)
	    std::cerr<<"Added sequence to alignment from node: " << nodeIter->first << endl; 
	  alignment[nodeIter->first] += nodeIter->second;
	}
    }

  if (logging)
    std::cerr<<"Displaying alignment...\n"; 

  // this is a late add-on, to make alignments a bit more readable...                                                                                                   
  int maxNameLength = 0;
  for (nodeIter = alignment.begin(); nodeIter!= alignment.end(); nodeIter++)
    {
      int nameSize = tree->node_name[nodeIter->first].size();
      maxNameLength = max(maxNameLength, nameSize );
    }

  for (nodeIter = alignment.begin(); nodeIter!= alignment.end(); nodeIter++)
    out << tree->node_name[nodeIter->first] + rep(maxNameLength-tree->node_name[nodeIter->first].size()+4, " ") + nodeIter->second + "\n";
}

// given a subtree alignment, verify that it is flush, then pad remaining nodes (as determined by the tree) with the appropriate gaps
void AlignmentSampler::pad(void)
{
  for (colIter = columns.begin(); colIter != columns.end(); colIter++)
    {
      for (node n = 0; n< tree->nodes(); n++)
	{
	  if (colIter->count(n) < 1 && colIter->size() )
	    (*colIter)[n] = "-"; 
	}
    }
}
// Add all the columns from another alignment to the current alignment, inserting them before the given position. 
void AlignmentSampler::add_cols(unsigned int pos, colMap cols2add)
{
  colMap newColumns; 
  for (unsigned int i=0; i<columns.size(); i++)
    {
      if ( i != pos)
	newColumns.push_back( columns[i] ); 
      else
	{
	  for ( colIter = cols2add.begin(); colIter != cols2add.end(); colIter++)
	    newColumns.push_back(*colIter); 
	  newColumns.push_back(columns[i]);
	}
    }
  columns = newColumns; 
}

M_id AlignmentSampler::find_last_del_state(const state_path& pi, int pathIdx, string side)
{
  pathIdx--; //we know it can't be the first one, since otherwise we wouldn't have called this function. 
  if (side=="left")
    while ( pi[pathIdx].left_type !=  profile_delete && pi[pathIdx].left_type !=  profile_start)
      pathIdx--;
  else if (side=="right")
    while ( pi[pathIdx].left_type !=  profile_delete && pi[pathIdx].left_type !=  profile_start)
      pathIdx--;
  else
    {
      std::cerr<<"unrecognized side: " << side << endl; 
      exit(1); 
    }
  return pi[pathIdx]; 
}


// sample a path through the transducer rooted at node, starting with the state at path[i][node], and ending in path[j][node]
colMap AlignmentSampler::sample_expanded_path(node n, int start, int end, bool viterbi)
{
  AbsorbingTransducer* abs = &((*AbsProfiles)[n]); 
  M_id mstart = abs->state2mid[start], mend = abs->state2mid[end]; 
  return sample_expanded_path(n, mstart, mend, viterbi); 
}

colMap AlignmentSampler::sample_expanded_path(node n, M_id start, M_id end, bool viterbi)
{
  colMap columns, subColumns; 
  DP_container forward; 
  queue<M_id> stateQueue; 
  M_id state; 
  vector<M_id>::iterator mIter; 
  vector<M_id> visited; 
  bfloat emissionWeight; 
  vector<bfloat> weights; 
  pair< vector<int>, vector<int> > transitionPair;  
  bool upstream, deep_logging = false; 
  map< vector<int>, bool> upstream_of_end; 

  // First get pointers to the relevant data structures - incoming/outgoing maps, transition weights, etc. 
  Profile* profile = &((*profiles)[n]); // a little awkward...
  map< vector<int> , vector<M_id> >* incoming = &(profile->incoming); 
  map< vector<int> , vector<M_id> >* outgoing = &(profile->outgoing);
  map< pair<vector<int>, vector<int> >, bfloat >* transition_weight = &(profile->transition_weight);

  // Possibly there's no null states between start and end.  In this case, skip DP and exit this function
//   transitionPair.first = start.toVector();   transitionPair.second = end.toVector(); 
//   if (!profile->between.count(transitionPair))
//     return columns; 

  
  int startIdx = 0, endIdx = profile->backward_states.size(); 
  while ( profile->backward_states[startIdx] != start ) 
    startIdx++; 
  while ( profile->backward_states[endIdx] != end ) 
    endIdx--;   

  if (logging)
    {
      std::cerr<<"Beginning sampling at node " << tree->node_name[n] << " states " << startIdx << " and " << endIdx <<endl; 
      std::cerr<<"Start:\n"; 
      start.display(profile->Q); 
      std::cerr<<"End:\n"; 
      end.display(profile->Q); 
    }
      


  // flag states as being upstream of end or not
  // We will compute DP values for those upstream of 'end' and not for others
  for (int i = profile->backward_states.size(); i>=0; i--)
    {
      state = profile->backward_states[i]; 
      upstream = false; 
      for (mIter = (*outgoing)[state.toVector()].begin(); mIter != (*outgoing)[state.toVector()].end(); mIter++)
	{
	  if ( *mIter == end || upstream_of_end[mIter->toVector()])
	    {
	      upstream = true; 
	      break; 
	    }
	}
      upstream_of_end[state.toVector() ] = upstream; 
    }
  if (logging)
    std::cerr<<"Beginning DP recursion...\n"; 
  // DP recursion 
  forward.set(start, 1.0); 
  for (unsigned int i = startIdx+1; i<profile->backward_states.size(); i++)
    {
      state = profile->backward_states[i]; 
      if (state != end) 
	if (!upstream_of_end[state.toVector()] || profile->is_external(state))
	  continue; 
      if (logging)
	{
	  //std::cerr<< "Position in backward states: " << i << endl; 
	  //	  std::cerr<< "new state \n"; 
	  //	  state.display(profile->Q); 
	}
      
      emissionWeight = profile->compute_emission_weight(state); 
      forward.set( state, 0.0 ); 
      for (mIter = (*incoming)[state.toVector()].begin(); mIter != (*incoming)[state.toVector()].end(); mIter++)
	{
	  transitionPair.first =  mIter->toVector(); 
	  transitionPair.second = state.toVector(); 
	  forward.add( state, forward.get(*mIter) * (*transition_weight)[transitionPair] );
	}
      forward.multiply(state, emissionWeight); 
      if (state == end && logging)
	std::cerr<< "End state reached, DP value: " << forward.get(state) << endl; 
    }
  // End dynamic programming
  if (logging)
    std::cerr<<"Finished DP.\n";   
  // Sample a path via traceback through this DP matrix.  Begin in 'end' and (hopefully!) end in 'start'
  // vector<M_id> visited; M_id state, vector<bfloat> weights; 

  list<M_id> pi; 
  int sampleIdx; 
  pi.push_front(end); 
  state = end; 
  if (logging)
    std::cerr<<"Beginning stochastic traceback.\n";   
  while ( state != start)
    {
      if(deep_logging)
	std::cerr<<"Sampling new state\n"; 
      if (deep_logging)
	state.display(profile->Q); 
      if(deep_logging)
	std::cerr<<"Getting incoming states...\n"; 
      if (! incoming->count(state.toVector() ) )
	{
	  std::cerr<<"Error: state has no incoming transitions!\n";
	  state.display(profile->Q); 
	  exit(1);
	}
      visited = profile->incoming[state.toVector()]; 
      weights.clear(); 
      if(deep_logging)
	std::cerr<<"Getting incoming weights...\n"; 
      for (mIter = visited.begin(); mIter != visited.end(); mIter++)
	{
	  transitionPair.first = mIter->toVector(); 
	  transitionPair.second = state.toVector(); 
	  if (! transition_weight->count(transitionPair) )
	    {
	      std::cerr<<"Error: transition pair was not found in profile's transition matrix.\n"; 
	      exit(1);
	    }
	  weights.push_back( forward.get(*mIter) * (*transition_weight)[transitionPair] ) ; 
	}
      if(deep_logging)
	std::cerr<<"Sampling from weights...\n"; 
      if (viterbi) 
	sampleIdx = maxIndex(weights);
      else 
	sampleIdx = sample(weights);
      if(deep_logging)
	std::cerr<<"Pushing to list\n"; 
      pi.push_front(visited[sampleIdx]); 
      state = visited[sampleIdx]; 
    }
  if(logging)
    std::cerr<<"Final path has length: " << pi.size() << endl;  


  // Now the tree-recursive bit...
  // Between every two states in this newly-sampled path, those subtrees need to be sampled. 
  M_id source; 
  vector<node> children; 
  if(logging)
    std::cerr<< "Beginning to sample in children: "; 
  for_rooted_children(*tree, n, child)
    {
      if(logging)
	std::cerr<< tree->node_name[*child] << " "; 
      children.push_back(*child); 
    }
    if(logging)
      std::cerr<< "\n"; 

  int pathIdx = -1; 
  string side; 
  for (list<M_id>::iterator piIter=pi.begin(); piIter!=pi.end(); piIter++)
    {
      pathIdx++; 
      if (piIter == pi.begin()) continue; 
      // right side
      if (tree->is_leaf(children[1]))
	{
	  if (logging)
	    std::cerr<<"Not sampling at node " << tree->node_name[children[1]] << ", it is a leaf\n"; 
	}
      else if (piIter->right_type == profile_delete || piIter->right_type == profile_pre_end )
	    {
	      side = "right"; 
	      source = find_last_del_state(l2v(pi), pathIdx, side); 
	      subColumns = sample_expanded_path(children[1], source.right_state, piIter->right_state ); 
	      for (colIter = subColumns.begin(); colIter != subColumns.end(); colIter++)
		columns.push_back(padColumn(*colIter, tree)); 
	    }

      if (tree->is_leaf(children[0]))
	{
	  if (logging)
	    std::cerr<<"Not sampling at node " << tree->node_name[children[0]] << ", it is a leaf\n"; 
	}
      else if (piIter->left_type == profile_delete || piIter->left_type == profile_pre_end )
	{
	  // left side
	  side = "left"; 
	  source = find_last_del_state(l2v(pi), pathIdx, side); 
	  subColumns = sample_expanded_path(children[0], source.left_state, piIter->left_state); 
	  for (colIter = subColumns.begin(); colIter != subColumns.end(); colIter++)
	    columns.push_back(padColumn(*colIter, tree)); 
	}
      
      if (pathIdx!=pi.size()-1 && pathIdx != 0)
	{
	  columns.push_back( padColumn(profile->state_type_phylogeny[piIter->toVector()] , tree) ); 
	  if (!profile->state_type_phylogeny.count(piIter->toVector()) )
	    std::cerr<<"STP didn't contain M_id!\n"; 
	  int stpSize = profile->state_type_phylogeny[piIter->toVector()].size(); 
	  std::cerr<<"Size of STP at node " << tree->node_name[n] << " " << stpSize <<endl;
	  std::cerr<<"State:  \n"; 
	  std::cerr<< piIter->q_state << " " << piIter->left_state << " " << piIter->left_type << " " << piIter->right_state << " " <<piIter->right_type << endl; 
	}
      

    }
  if (logging)
    {  
      std::cerr<<"\n\n\n"; 
      std::cerr<<"Alignment sampled at node " << tree->node_name[n] << " " << columns.size()<< " columns :\n"; 
      AlignmentSampler a; 
      a.columns = columns; 
      a.tree = tree; 
      a.display(std::cerr); 
      std::cerr<<"\n\n\n"; 
    }
  return columns; 
}


// *****Functions related to tabulating indel information.  *****
IndelCounter::IndelCounter(Stockholm& stk, PHYLIP_tree* tree_in)
{
  tree = tree_in; 
  Phonebook::iterator seq;
  sstring sequence;
  node n; 
  unsigned int size=0; 
  verySmall = .01; 
  for (seq = stk.row_index.begin(); seq!=stk.row_index.end(); seq++)
    {
      n = index(seq->first, tree->node_name); 
      if ( tree->is_leaf(n) )
	sequence =  stk.get_row_as_string(seq->second);
      else
	sequence = stk.gr_annot[seq->first]["ancrec_CYK_MAP"];
      rows[n] = string(sequence.c_str()); 
    }

  // make sure is flush
  for (rowIter = rows.begin(); rowIter!=rows.end(); rowIter++)
    {
      if (rowIter == rows.begin())
	size = rowIter->second.size(); 
      else if ( rowIter->second.size() != size )
	{
	  std::cerr<<"Error: input alignment appears not to be flush.  Indel counting not possible! Exiting. \n"; 
	  exit(1); 
	}
    }
  L = size; 
}  

IndelCounter::IndelCounter(map<string,string>& sequences, PHYLIP_tree* tree_in)
{
  tree = tree_in; 
  map<string,string>::iterator seq;
  sstring sequence;
  node n; 
  unsigned int size=0; 
  for (seq = sequences.begin(); seq!=sequences.end(); seq++)
    {
      n = index(seq->first, tree->node_name); 
      rows[n] = seq->second;
    }

  // make sure is flush
  for (rowIter = rows.begin(); rowIter!=rows.end(); rowIter++)
    {
      if (rowIter == rows.begin())
	size = rowIter->second.size(); 
      else if ( rowIter->second.size() != size )
	{
	  std::cerr<<"Error: input alignment appears not to be flush.  Indel counting not possible! Exiting. \n"; 
	  exit(1); 
	}
    }
  L = size; 
}  



void IndelCounter::display_indel_info(ostream& out, bool per_branch)
{
  out<<"\n\n### Displaying indel information ###\n";  
  out << "#branch insertion_rate \t insertion_extend \t deletion_rate \t deletion_extend\n"; 
  if (per_branch)
    {
      double insert_rt, delete_rt, delete_ext, insert_ext; 
      node treeNode; 
      for_nodes_pre (*tree, tree->root, -1, bi)
	{
	  const Phylogeny::Branch& b = *bi;
	  treeNode = b.second; 
	  if (treeNode == tree->root)
	    continue; 
	  insert_rt = insert_rate(treeNode);
	  delete_rt = delete_rate(treeNode);
	  
	  insert_ext = insert_extend(treeNode);
	  delete_ext = delete_extend(treeNode);
	
	  out << tree->node_name[treeNode] << "\t" << insert_rt << "\t" << insert_ext << "\t" << delete_rt << "\t" << delete_ext << endl; 
// 	  out << "Insertion rate (insertions per position per unit time) above branch " << tree->node_name[treeNode] << ": " << insert_rt <<endl; 
// 	  if (insert_rt>0.0)
// 	    out << "Insertion extension probability (1/(mean insertion length)) above branch " << tree->node_name[treeNode] << ": " << insert_ext <<endl; 
	  
// 	  out << "Deletion rate (deletions per position per unit time): above branch " << tree->node_name[treeNode] << ": " << delete_rt <<endl; 
// 	  if (delete_rt>0.0)
// 	    out << "Deletion extension probability (1/(mean deletion length)) above branch " << tree->node_name[treeNode] << ": " << delete_ext <<endl; 
	}
    }
  out << "\n### Indel statistics averaged across tree ###\n"; 


  out << "Insertion rate (insertions per position per unit time): " << avg_insert_rate <<endl; 
  out << "Insertion extension probability (1/(mean insertion length)): " << avg_insert_ext <<endl; 
  
  out << "Deletion rate (deletions per position per unit time): " << avg_delete_rate <<endl; 
  out << "Deletion extension probability (1/(mean deletion length)): " << avg_delete_ext <<endl; 
  
  out << "Inserted sequence: ";
  all_insertions(out); 
  out <<"\n"; 
  out << "Deleted sequence: ";
  all_deletions(out); 

  out <<"\n\n\n"; 
}

void IndelCounter::gather_indel_info(bool logging)
{
  // Traverse pairwise alignments, gathering info on each of them - how many matches, inserts, dels, and their lengths
  // For now (?) represent a pairwise alignment with a pair o' strings.  meh
  string parent, child, gap = "-"; 
  node parent_node, child_node;
  unsigned int pair_length, match_count, pos; 
  bool prev_was_ins, prev_was_del, prev_was_match;

  for_nodes_pre (*tree, tree->root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      parent_node = b.first; 
      child_node = b.second; 
      if (child_node == tree->root)
	continue; 
      get_pairwise_alignment(parent_node, child_node, parent, child);
      pair_length = parent.size(); 
      if (logging)
	std::cerr<<"\n\nAbout to start tabulating indels from this pairwise alignment:\n" << 
	  tree->node_name[parent_node] << "  " << parent << "\n" << 
	  tree->node_name[child_node] << "  " << child << "\n"; 


      // start counting indels, matches, etc
      prev_was_ins = prev_was_del = prev_was_match = false; 
      match_count = 0; 
      for (pos=0; pos< pair_length; pos++)
	{
	  if (logging)
	    std::cerr<< "Investigating character pair: " << stringAt(parent, pos) << " " << stringAt(child, pos) << ": ";
	  // match column
	  if (stringAt(parent, pos) != gap && stringAt(child, pos) != gap)
	    {
	      if (logging)
		std::cerr<< " match\n"; 
	      prev_was_match = true; prev_was_ins = prev_was_del = false; 
	      match_count++; 
	    }
	  // delete - parent has a character, but child doesn't
	  else if ( stringAt(parent, pos) != gap && stringAt(child, pos) == gap)
	    {

	      if ( prev_was_del )
		{
		  deletions[child_node][ deletions[child_node].size()-1 ] += stringAt(parent, pos); 
		  if (logging)
		    std::cerr<< " delete-extend\n"; 
		}
	      else
		{
		  deletions[child_node].push_back( stringAt(parent, pos));
		  if (logging)
		    std::cerr<< " delete-open\n"; 
		}
	      prev_was_del = true; prev_was_ins = prev_was_match = false; 
	    }
	  // insert - parent has no character, but child does
	  else if ( stringAt(parent, pos) == gap && stringAt(child, pos) != gap)
	    {
	      if ( prev_was_ins )
		{
		  if (logging)
		    std::cerr<< " insert-extend\n"; 
		  insertions[child_node][ insertions[child_node].size()-1 ] += stringAt(child, pos); 
		}
	      else
		{
		  if (logging)
		    std::cerr<< " insert-open\n"; 
		  insertions[child_node].push_back( stringAt(child, pos));
		}
	      prev_was_ins = true; prev_was_del = prev_was_match = false; 
	    }
	}
      matches[child_node] = match_count; 
      if (logging)
	{
	  std::cerr<<"Done with current alignment. Basic statistics:\n"; 
	  std::cerr<<"Branch length : " << tree->branch_length(parent_node, child_node) << endl; 
	  std::cerr<<"Number of matches: " << match_count << endl; 
	  std::cerr<<"Number of insertions: " << insertions[child_node].size() << endl; 
	  std::cerr<<"Number of deletions: " << deletions[child_node].size() << endl; 
	  std::cerr<<"\n\n"; 
	}

	  
    }
  average_indel_counts(); 
}

void IndelCounter::all_insertions(ostream& out)
{
  node treeNode; 
  for_nodes_pre (*tree, tree->root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      treeNode = b.second; 
      if (treeNode == tree->root)
	continue; 
      for (vector<string>::iterator insIter=insertions[treeNode].begin(); insIter != insertions[treeNode].end(); insIter++)
	out << " " << *insIter << " "; 
    }
}

void IndelCounter::all_deletions(ostream& out)
{
  node treeNode; 
  for_nodes_pre (*tree, tree->root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      treeNode = b.second; 
      if (treeNode == tree->root)
	continue; 
      for (vector<string>::iterator delIter=deletions[treeNode].begin(); delIter != deletions[treeNode].end(); delIter++)
	out << " " << *delIter << " "; 
    }
}

void IndelCounter::average_indel_counts(void)
{
  double insert_rt, delete_rt, delete_ext, insert_ext; 
  double insert_rate_tot = 0.0,  delete_rate_tot = 0.0,  delete_extend_tot = 0.0,  insert_extend_tot = 0.0; 
  node treeNode; 
  double b = double(tree->branches()); 
  int insert_count = 0, delete_count = 0; // actually the count of nodes which have nonzero in/del rates
  for_nodes_pre (*tree, tree->root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      treeNode = b.second; 
      if (treeNode == tree->root)
	continue; 
      insert_rt = insert_rate(treeNode);
      delete_rt = delete_rate(treeNode);

      insert_ext = insert_extend(treeNode);
      delete_ext = delete_extend(treeNode);

      insert_rate_tot +=  insert_rt;
      delete_rate_tot +=  delete_rt;
      
      if (insert_rt > 0.0)
	{
	  insert_extend_tot +=  insert_ext;
	  insert_count++; 
	}
      if (delete_rt > 0.0)
	{
	  delete_extend_tot +=  delete_ext;
	  delete_count++; 
	}
    }
  avg_insert_rate = insert_rate_tot/b; 
  avg_insert_ext = insert_extend_tot/insert_count; 
  
  avg_delete_rate = delete_rate_tot/b; 
  avg_delete_ext = delete_extend_tot/delete_count; 
}

double IndelCounter::insert_extend(node n)
{
  double totalSize = 0; 
  for (vector<string>::iterator inserts = insertions[n].begin(); inserts != insertions[n].end(); inserts++)
    totalSize += inserts->size(); 
  return 1.0 - (double(insertions[n].size()) / totalSize);
}

double IndelCounter::insert_rate(node n)
{
  return insertions[n].size() / ((matches[n] + deletions[n].size() + 1) * max(tree->branch_length(n, tree->parent[n]), verySmall));
}

double IndelCounter::delete_extend(node n)
{
  double totalSize = 0; 
  for (vector<string>::iterator deletes = deletions[n].begin(); deletes != deletions[n].end(); deletes++)
    totalSize += deletes->size(); 
  return 1.0 - (double(deletions[n].size()) / totalSize);
}
  
double IndelCounter::delete_rate(node n)
{
  return deletions[n].size() / ((matches[n] + deletions[n].size() + 1) * max(tree->branch_length(n, tree->parent[n]), verySmall));
}



void IndelCounter::get_pairwise_alignment(node parent, node child, string& parentString, string& childString)
{
  parentString = "";   childString = ""; 
  string gap = "-"; // hmm..maybe should set this somewhere else. 
  for (unsigned int colPos=0; colPos < L ; colPos++)
    {
      if (stringAt(rows[parent],colPos) == gap && stringAt(rows[child],colPos) == gap)
	continue; 
      else
	{
	  parentString += stringAt(rows[parent], colPos); 
	  childString += stringAt(rows[child], colPos); 
	}
    }
}
	


// Lightweight DP container, basic methods for access,adding, etc. 
DP_container::DP_container(void)
{
  // Placeholder constructor
}

void DP_container::set(M_id m, bfloat value)
{
  if (data.count(m.toVector()) > 0)
    {
      std::cerr<<"Error: Attempting to 'set' an already-existing DP cell\n";
      exit(1); 
    }
  else
    data[m.toVector()] = value; 
}
	

void DP_container::add(M_id m, bfloat value)
{
  if (data.count(m.toVector()) < 1)
    {
      std::cerr<<"Error: Attempting to add to  a nonexistentf DP cell\n";
      exit(1); 
    }
  else
    data[m.toVector()] += value; 
}

void DP_container::multiply(M_id m, bfloat value)
{
  if (data.count(m.toVector()) < 1)
    {
      std::cerr<<"Error: Attempting to add to  a nonexistentf DP cell\n";
      exit(1); 
    }
  else
    data[m.toVector()] *= value; 
}


bfloat DP_container::get(M_id m)
{
  if (data.count(m.toVector()) < 1)
    return 0.0; 
  else
    return data[m.toVector()]; 
}

bool all_parents_visited(vector<M_id> parents, vector<M_id> visited)
{
  for (vector<M_id>::iterator parent=parents.begin(); parent!=parents.end(); parent++)
    if ( index(*parent, visited) == -1)
      return false; 
  return true; 
}

map<node, string>  padColumn(map<node, string>& column, PHYLIP_tree* tree )
{
  for (node n = 0; n< tree->nodes(); n++)
    {
      if (column.count(n) < 1 && column.size() )
	column[n] = "-"; 
    }
  if (!column.size())
    std::cerr<<"Warning: returning empty column\n"; 
  return column; 
}
