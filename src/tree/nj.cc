#include "tree/nj.h"
#include "util/vector_output.h"

#define NJ_NODE_NAME_PREFIX "node_"

NJ_tree::NJ_tree() : PHYLIP_tree() { }

void NJ_tree::build (const Node_name_vec& nn, const array2d<double>& distance_matrix)
{
  // check that there are more than 2 nodes
  if (nn.size() < 2)
    THROWEXPR ("Fewer than 2 nodes; can't make a binary tree");
  // clear the existing tree
  clear();
  node_name.clear();
  parent.clear();
  // copy distance matrix into a conveniently resizable data structure
  vector<vector<double> > dist = distance_matrix.all_columns();
  // estimate tree by neighbor-joining
  // algorithm follows description in Durbin et al, pp170-171
  // first, initialise the list of active nodes
  set<Node> active_nodes;
  for (Node n = 0; n < (int) nn.size(); ++n)
    {
      const Node k = new_node();
      if (k != n)
	THROWEXPR ("new_node() returned unexpected node index " << k << " (expected " << n << ")");
      active_nodes.insert (n);
      node_name.push_back (nn[k]);
      parent.push_back (-1);
    }
  // main loop
  int n_nodes = node_name.size();
  vector<double> avg_dist;
  while (1)
    {
      // get number of active nodes
      const int n_active_nodes = active_nodes.size();
      // loop exit test
      if (n_active_nodes == 2) break;
      if (n_active_nodes < 2)
	THROWEXPR ("Fewer than 2 nodes left -- should never get here");
      // calculate average distances from each node
      avg_dist = vector<double> (n_nodes, (double) 0);
      for_const_contents (set<Node>, active_nodes, node_i_p)
	{
	  double a_i = 0;
	  for_const_contents (set<Node>, active_nodes, node_j_p)
	    if (*node_j_p != *node_i_p)
	      a_i += dist [*node_i_p] [*node_j_p];
	  avg_dist [*node_i_p] = a_i / (double) (n_active_nodes - 2);
	  CTAG(2,NJ) << "Distance correction for node " << *node_i_p << " is " << avg_dist [*node_i_p] << "\n";
	}
      // find minimal compensated distance (with avg distances subtracted off)
      bool first_pair = TRUE;
      double min_dist = 0;
      Node min_i, min_j;
      set<Node>::const_iterator node_i_p, node_j_p, active_nodes_end_p;
      active_nodes_end_p = active_nodes.end();
      for (node_i_p = active_nodes.begin(); node_i_p != active_nodes_end_p; ++node_i_p)
	for (node_j_p = node_i_p, ++node_j_p; node_j_p != active_nodes_end_p; ++node_j_p)
	  {
	    const double compensated_dist = dist [*node_i_p] [*node_j_p] - avg_dist [*node_i_p] - avg_dist [*node_j_p];
	    if (CTAGGING(2,NJ))
	      CL << "Compensated distance from node " << *node_i_p << " to node " << *node_j_p << " is " << compensated_dist << "\n";
	    if (first_pair || compensated_dist < min_dist)
	      {
		min_i = *node_i_p;
		min_j = *node_j_p;
		min_dist = compensated_dist;
		first_pair = FALSE;
	      }
	  }
      // nodes min_i and min_j are neighbors -- join them with new index k
      // first, calculate new distances as per NJ algorithm
      dist.push_back (vector<double> (n_nodes + 1));
      const int k = n_nodes++;
      dist[k][k] = 0;
      const double d_ij = dist[min_i][min_j];
      for (Node m = 0; m < k; ++m)
	dist[m].push_back (dist[k][m] = dist[min_i][m] + dist[min_j][m] - d_ij);
      double d_ik = 0.5 * (d_ij + avg_dist[min_i] - avg_dist[min_j]);
      double d_jk = d_ij - d_ik;
      CTAG(1,NJ) << "Before Kuhner-Felsenstein:\ni=" << min_i << ", j=" << min_j << ", k=" << k << ", d_ij=" << d_ij << ", d_ik=" << d_ik << ", d_jk=" << d_jk << "\nDistances from k to other nodes: " << dist[k] << "\n";
      // apply Kuhner-Felsenstein correction to prevent negative branch lengths
      if (d_ik < 0)
	{
	  d_jk -= d_ik;
	  d_ik = 0;
	}
      if (d_jk < 0)
	{
	  d_ik -= d_jk;
	  d_jk = 0;
	}
      dist[min_i][k] = dist[k][min_i] = d_ik;
      dist[min_j][k] = dist[k][min_j] = d_jk;
      // now update data structures
      const Node my_k = add_node (min_i, d_ik);
      if (my_k != k)
	THROWEXPR ("add_node() returned unexpected node index " << my_k << " (expected " << k << ")");
      CTAG(3,NJ) << "Joining nodes " << min_i << " (" << node_name[min_i] << ") and " << min_j << " (" << node_name[min_j] << ") to common ancestor " << k << " (branch lengths: " << k << "->" << min_i << " = " << d_ik << ", " << k << "->" << min_j << " = " << d_jk << ")\n";
      add_branch (min_j, k, d_jk);
      parent[min_i] = parent[min_j] = k;
      active_nodes.erase (min_i);
      active_nodes.erase (min_j);
      active_nodes.insert (k);
      sstring k_name;
      k_name << NJ_NODE_NAME_PREFIX << k;
      node_name.push_back (k_name);
      parent.push_back (-1);
    }
  // make the root node
  set<Node>::iterator iter = active_nodes.begin();
  const Node i = *iter;
  const Node j = *++iter;
  const double d = max (dist[i][j], 0.);  // don't correct the last node, just keep branch length non-negative
  const Node k = add_node (i, d/2);
  if (k != n_nodes)
    THROWEXPR ("Tree returned unexpected Node index -- CAREFUL; this will probably cause code to break");
  add_branch (j, k, d/2); 
  parent[i] = parent[j] = k;
  root = k;
  node_name.push_back (DART_ROOT_IDENTIFIER);
  parent.push_back (-1);
  CTAG(3,NJ) << "Making root node: joining nodes " << i << " (" << node_name[i] << ") and " <<  j<< " (" << node_name[j] << ") to common ancestor " << k << " (" << node_name[k] << ") (branch length " << d/2 << " to each child)\n";
  // print tree
  if (CTAGGING(4,NJ)) { CL << "Neighbor-joining tree: "; write(CL); }
}
