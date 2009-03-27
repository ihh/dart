#ifndef STATESORTER_INCLUDED
#define STATESORTER_INCLUDED

#include "hmm/transmat.h"
#include "indiegram/tripletscfg.h"

typedef vector<map<int, Score> > State_graph;


/// Sorting methods for SCFG states.
/*
 * More robust than those provided by Transition_methods (transmat.h).
 */
struct SCFG_state_sorter : SCFG_state_typing, Grammar_state_enum
{

 public:

  /// Get all emit states.
  vector<int> emit_states() const;

  /// Get all non-emit states.
  vector<int> nonemit_states_unsorted() const;

  /// Get all non-emit states, sorted topologically.
  vector<int> nonemit_states_sorted() const;  // sorted topologically; throws exception if null cycle detected

  /// Perform a topological sort of states according to graph.
  /*
   * Throws an exception if a cycle detected.
   */
  vector<int> topological_sort (const vector<int> states, State_graph graph) const;

  /// Show state graph.
  void show_state_graph (State_graph graph, ostream& o) const;

  // constructors
  SCFG_state_sorter (const SCFG_state_sorter& base); /// copy constructor
  SCFG_state_sorter (const Triplet_SCFG& scfg); /// initialize from a Triplet_SCFG

  SCFG_state_sorter& operator= (const SCFG_state_sorter& base); /// assignment operator

 private:

  void zero_outgoing_bifurcation_transitions();

  void add_fake_bifurcation_transitions();

  /// State graph.
  /*
   * Transition matrix + bifurcations and their children.
   * This representation is memory-inefficient (it takes num_states^2 space)
   * but is fast.
   * state_graph[src_state] is a vector of destination states for src_state.
   */
  State_graph state_graph;

};

#endif /* STATESORTER_INCLUDED */
