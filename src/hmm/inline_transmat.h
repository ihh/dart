#include "util/vector_output.h"
#include "util/map_keys.h"

// templated method code for "transmat.h"

template <class T, class T_matrix>
void Transition_matrix<T,T_matrix>::print_dotfile (ostream& out) const
{
  out << "digraph " << dotfile_graph_name() << " {\n";
  print_dotfile_nodes (out);
  print_dotfile_edges (out);
  out << "}\n";
}

template <class T, class T_matrix>
void Transition_matrix<T,T_matrix>::print_dotfile_nodes (ostream& out) const
{
  out << " node [shape=plaintext]\n";
  print_dotfile_node (out, Start);
  print_dotfile_node (out, End);
  for (int s = 0; s < tm_states(); ++s)
    print_dotfile_node (out, s);
}

template <class T, class T_matrix>
void Transition_matrix<T,T_matrix>::print_dotfile_edges (ostream& out) const
{
  for (int src = Start; src < tm_states(); ++src)
    for (int dest = End; dest < tm_states(); ++dest)
      if (dest != Start)
	if (is_non_null (transition (src, dest)))
	  print_dotfile_edge (out, src, dest);
}

template <class T, class T_matrix>
sstring Transition_matrix<T,T_matrix>::dotfile_graph_name() const
{
  return sstring ("transmat");
}

template <class T, class T_matrix>
sstring Transition_matrix<T,T_matrix>::dotfile_node_id (int state) const
{
  sstring s;
  if (state == Start)
    s << "Start";
  else if (state == End)
    s << "End";
  else
    s << state + 1;
  return s;
}

template <class T, class T_matrix>
sstring Transition_matrix<T,T_matrix>::dotfile_node_label (int state) const
{
  return dotfile_node_id(state);
}

template <class T, class T_matrix>
sstring Transition_matrix<T,T_matrix>::dotfile_edge_weight (int src, int dest) const
{
  sstring s;
  show_element (transition (src, dest), s);
  return s;
}

template <class T, class T_matrix>
void Transition_matrix<T,T_matrix>::print_dotfile_node (ostream& out, int state) const
{
  out << " " << dotfile_node_id(state) << " [";
  typedef map<sstring,sstring> AttrMap;
  AttrMap attr = dotfile_node_attrs(state);
  int n = 0;
  for_const_contents (AttrMap, attr, a)
    out << (n++ ? "," : "") << a->first << '=' << a->second;
  out << "];\n";
}

template <class T, class T_matrix>
map<sstring,sstring> Transition_matrix<T,T_matrix>::dotfile_node_attrs (int state) const
{
  map<sstring,sstring> attr;
  attr[sstring("label")] = dotfile_node_label(state);
  return attr;
}

template <class T, class T_matrix>
void Transition_matrix<T,T_matrix>::print_dotfile_edge (ostream& out, int src, int dest) const
{
  out << " " << dotfile_node_id(src) << " -> " << dotfile_node_id(dest) << " [";
  typedef map<sstring,sstring> AttrMap;
  AttrMap attr = dotfile_edge_attrs(src,dest);
  int n = 0;
  for_const_contents (AttrMap, attr, a)
    out << (n++ ? "," : "") << a->first << '=' << a->second;
  out << "];\n";
}

template <class T, class T_matrix>
map<sstring,sstring> Transition_matrix<T,T_matrix>::dotfile_edge_attrs (int src, int dest) const
{
  map<sstring,sstring> attr;
  sstring ew;
  ew << '"' << dotfile_edge_weight(src,dest) << '"';
  attr[sstring("label")] = ew;
  return attr;
}

template <class T, class T_matrix>
void Transition_matrix<T,T_matrix>::show_transitions (ostream& o) const
{
  int old_prec = o.precision(3);
  save_flags (o);
  right_align (o);

  const int w = element_width();
  const int s = 5;   // width of state names in transition matrix

  o << "Transition matrix " << element_descriptor() << ":\n";
  o.width(s);
  o << "State" << " ";
  for (int i = 0; i < tm_states(); i++) { o.width(w); o << i; o << " "; }
  o.width(w);
  o << "End" << "\n";
  o.width(s);
  o << "Start" << " ";
  for (int i = 0; i < tm_states(); i++) { o.width(w); show_element(start[i],o); o << " "; }
  o.width(w);
  show_element(start_to_end(),o);
  o << "\n";
  for (int i = 0; i < tm_states(); i++)
    {
      o.width(s);
      o << i << " ";
      for (int j = 0; j < tm_states(); j++) { o.width(w); show_element(transition(i,j),o); o << " "; }
      o.width(w);
      show_element(end[i],o);
      o << "\n";
    }
  restore_flags (o);
  o.precision (old_prec);
}

template <class TScoreMatrix>
vector<vector<int> > Transition_methods::incoming_states (const TScoreMatrix& ts)
{
  vector<vector<int> > result (ts.tm_states(), vector<int>());
  for (int j = 0; j < ts.tm_states(); ++j)
    for (int i = 0; i < ts.tm_states(); ++i)
      if (ts.transition(i,j) > -InfinityScore)
	result[j].push_back(i);
  return result;
}

template <class TScoreMatrix>
vector<vector<int> > Transition_methods::selected_outgoing_states (const TScoreMatrix& ts, const vector<int>& selection)
{
  vector<vector<int> > result (ts.tm_states(), vector<int>());
  for (int i = 0; i < ts.tm_states(); ++i)
    for_const_contents (vector<int>, selection, j)
      if (ts.transition(i,*j) > -InfinityScore)
	result[i].push_back(*j);
  return result;
}

template <class TScoreMatrix>
vector<vector<int> > Transition_methods::selected_incoming_states (const TScoreMatrix& ts, const vector<int>& selection)
{
  vector<vector<int> > result (ts.tm_states(), vector<int>());
  for (int j = 0; j < ts.tm_states(); ++j)
    for_const_contents (vector<int>, selection, i)
      if (ts.transition(*i,j) > -InfinityScore)
	result[j].push_back(*i);
  return result;
}

// functor to compare second elements of two pair<int,int>'s
struct CompareSecondIntInPair : binary_function <pair<int,int>, pair<int,int>, bool>
{ result_type operator() (first_argument_type a, second_argument_type b) const { return a.second < b.second; } };

// topological sort
template <class TScoreMatrix>
vector<int> Transition_methods::topological_sort (const TScoreMatrix& ts, const vector<int>& null_states, bool suppress_cycle_warning)
{
  // this code was written to do topological sorts of null states
  // so the word "null" crops up a lot in the variable namespace
  // and i can't be bothered to change it

  // the topological sort algorithm looks like this:
  // while not done
  //  if the graph is empty, then we're done (exit the loop)
  //  pick a node with no predecessors
  //  if no such node exists , then the graph is cyclic (exit the loop)
  //  output that node (number that node)
  //  delete that node from the graph
  // end while

  // find sources and sinks
  vector<int> sources (ts.states(), (int) 0);
  vector<vector<int> > sinks (ts.states());
  for_const_contents (vector<int>, null_states, src)
    for_const_contents (vector<int>, null_states, dest)
    if (*src != *dest && ts.transition (*src, *dest) > -InfinityScore)
      {
	++sources[*dest];
	sinks[*src].push_back (*dest);
      }

  // initialize results and unsorted queue
  vector<int> sorted;  // sorted nodes
  set<int> unsorted (null_states.begin(), null_states.end());  // remaining nodes

  // main loop
  while (unsorted.size())
    {
      // check for nodes with no unsorted predecessors
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
      if (CTAGGING(5,TOPSORT TOPSORT_CYCLE))
	{
	  CL << "Remaining transitions in topological sort graph:\n";
	  for_const_contents (set<int>, unsorted, src)
	    {
	      CL << *src << " ->";
	      for_const_contents (vector<int>, sinks[*src], dest)
		if (unsorted.find(*dest) != unsorted.end())
		  CL << ' ' << *dest;
	      CL << '\n';
	    }
	}
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
