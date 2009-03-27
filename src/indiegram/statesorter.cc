#include <stack>
#include "indiegram/tripletscfg.h"
#include "indiegram/statesorter.h"

// copy constructor
SCFG_state_sorter::SCFG_state_sorter (const SCFG_state_sorter& base)
  : SCFG_state_typing (base),
    state_graph (base.state_graph)
{ }
 
SCFG_state_sorter::SCFG_state_sorter (const Triplet_SCFG& scfg)
  : SCFG_state_typing (scfg)
{

  // initialize state_graph
  state_graph.resize (num_states());

  // and fill it in
  for (int s = 0; s < num_states(); ++s)
    for (int d = 0; d < num_states(); ++d)
    {
      if (scfg.transition_scores.transition (s, d) > -InfinityScore)
	state_graph[s].insert (make_pair (d, scfg.transition_scores.transition (s, d)));
    }

}


void SCFG_state_sorter::zero_outgoing_bifurcation_transitions()
{
  for (int s = 0; s < num_states(); ++s)
    if (is_bifurc_type (state_type[s]))
      for (int d = 0; d < num_states(); ++d)
	(state_graph[s])[d] = -InfinityScore;
}

void SCFG_state_sorter::add_fake_bifurcation_transitions()
{

  (*this).zero_outgoing_bifurcation_transitions();
  for (int s = 0; s < num_states(); ++s)
    if (is_bifurc_type (state_type[s]))
      {
	state_graph[s].insert (make_pair (bifurc[s].l, (Score) 0));
	state_graph[s].insert (make_pair (bifurc[s].r, (Score) 0));
      }

}


// assignment operator
SCFG_state_sorter& SCFG_state_sorter::operator= (const SCFG_state_sorter& base)
{
  assign_state_typing (base);
  state_graph = base.state_graph;

  return *this;
}


vector<int> SCFG_state_sorter::emit_states() const
{
  vector<int> result;
  for (int s = 0; s < num_states(); ++s)
    if (is_emit_type (state_type[s])) result.push_back (s);
  return result;
}

vector<int> SCFG_state_sorter::nonemit_states_unsorted() const
{
  vector<int> result;
  for (int s = 0; s < num_states(); ++s)
    if (!is_emit_type (state_type[s])) result.push_back (s);
  return result;
}

// to do: 
// - check that indiegram still passes testtripletdp and runs ok on tiny.stk with --testgrammer
// - see if indiegram can run on tiny.stk with tkfst and whether testtkfst works

// first: have it just sort Bifurc states last
// if the add_fake_bifurcation_transitions hack doesn't work,
// then just manually remove them and stuff them all in at the end
// hacky but effective!
vector<int> SCFG_state_sorter::nonemit_states_sorted() const
{
  // see note for Pair_CFG_scores::nonemit_states()
  ((SCFG_state_sorter*) this) -> add_fake_bifurcation_transitions();  // a hacky cast for a hacky technique
  const vector<int> unsorted = nonemit_states_unsorted();
  return topological_sort (unsorted, (*this).state_graph);
}


vector<int> SCFG_state_sorter::topological_sort (const vector<int> states, const State_graph graph) const
{

  // the topological sort algorithm looks like this:
  // while not done
  //  if the graph is empty, then we're done (exit the loop)
  //  pick a node with no predecessors
  //  if no such node exists , then the graph is cyclic (exit the loop)
  //  output that node (number that node)
  //  delete that node from the graph
  // end while

  // find sources and sinks
  vector<int> sources (num_states(), (int) 0);
  vector<vector<int> > sinks (num_states());
  for_const_contents (vector<int>, states, src)
    for (map<int, Score>::const_iterator entry = graph[*src].begin(); entry != graph[*src].end(); ++entry)
      {
	int dest = entry->first;
	++sources[dest];
	sinks[*src].push_back (dest);
      }

  // initialize results and unsorted queue
  vector<int> sorted;  // sorted nodes
  set<int> unsorted (states.begin(), states.end());  // remaining nodes

  // main loop
  while (unsorted.size())
    {
      // check for nodes with no predecessors
      bool found_node = false;
      int node;
      for_const_contents (set<int>, unsorted, n)
	if (sources[*n] == 0)
	  {
	    node = *n;
	    found_node = true;
	    break;
	  }

      // bail if no node found
      if (!found_node)
	break;

      // store the node
      sorted.push_back (node);

      // delete node from graph
      unsorted.erase (node);
      for_const_contents (vector<int>, sinks[node], dest)
	--sources[*dest];
    }

  // check for null cycles
  if (unsorted.size())
    {
      CLOGERR << "Warning: found null cycle(s) during topological sort. Some parses may be missed!\n";
      // put all remaining nodes onto the sorted list
      for_const_contents (set<int>, unsorted, n)
	sorted.push_back (*n);
    }

  // print log message
  if (CTAGGING(0,TOPSORT))
    CL << "Topological sort order of states: (" << sorted << ")\n";

  // return
  return sorted;
}


void SCFG_state_sorter::show_state_graph (State_graph graph, ostream& o) const
{

  for (int s = 0; s < (int) graph.size(); ++s)
    {
      o << s << " -> ";
      for (map<int, Score>::const_iterator entry = graph[s].begin(); entry != graph[s].end(); ++entry)
	{
	  int dest = entry->first;
	  o << dest << ", ";	  
	}
      o << "\n";
    }

}
