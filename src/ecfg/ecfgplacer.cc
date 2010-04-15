#include "ecfg/ecfgplacer.h"

ECFG_placer::ECFG_placer (ECFG_scores& ecfg, Stockholm& stock, Tree_alignment& tree_align, double prior_param) :
  ecfg (ecfg),
  stock (stock),
  tree_align (tree_align),
  tree (tree_align.tree),
  prior_param (prior_param)
{
  for (int row = 0; row < tree_align.align.rows(); ++row)
    if (tree_align.row2node[row] < 0)
      unattached_rows.push_back (row);
}

void ECFG_placer::populate_counts()
{
  // reset attach_counts
  attach_counts.clear();
  ECFG_branch_state_counts empty_counts;
  for (int c = 0; c < (int) ecfg.matrix_set.chain.size(); ++c)
    {
      const int states = ecfg.matrix_set.total_states(c);
      empty_counts.push_back (Branch_state_counts (states, states, 0.));
    }
  for_const_contents (vector<int>, unattached_rows, row)
    for_rooted_nodes_post (tree, b)
    {
      const Phylogeny::Node node = (*b).second;
      attach_counts[*row][node] = empty_counts;
    }

  // get CYK matrix & trace
  const int seqlen = tree_align.align.columns();
  const int max_subseq_len = -1;   // unlimit subsequence lengths
  ECFG_auto_envelope env (seqlen, ecfg, max_subseq_len);  // create fold envelope

  Aligned_score_profile asp;
  asp.init (tree_align.align, stock.np, ecfg.alphabet);

  const bool try_fast_prune = false;
  ECFG_CYK_matrix cyk_mx (ecfg, stock, asp, env, try_fast_prune);   // create CYK matrix
  cyk_mx.fill();
  ECFG_cell_score_map cyk_trace = cyk_mx.traceback();  // get traceback

  // create dummy counts
  ECFG_counts dummy_counts (ecfg);

  // loop over nonterminals in parse tree, using reverse order for easy comparison with methods that loop over columns in left-to-right order
  for_const_reverse_contents (ECFG_cell_score_map, cyk_trace, ecsm_ptr)
    {
      const Subseq_coords& coords = ecsm_ptr->first.first;
      const int ecfg_state = ecsm_ptr->first.second;

      // get state info
      const ECFG_state_info& info = ecfg.state_info[ecfg_state];
      if (info.emit_size())
	{
	  // get chain
	  const int chain_idx = info.matrix;
	  const ECFG_chain& chain = ecfg.matrix_set.chain[chain_idx];
	  const int chain_states = ecfg.matrix_set.total_states(chain_idx);

	  // call fill_down
	  cyk_mx.fill_down (dummy_counts, cyk_mx.env.find_subseq_idx (coords.start, coords.len), ecfg_state, 1.);

	  // get posterior state probabilities for each node
	  EM_matrix_base::Column_matrix& colmat = cyk_mx.colmat[ecfg_state];

	  typedef map<Phylogeny::Node,vector<Prob> > Node_state_prob_map;
	  Node_state_prob_map node_state_probs;
	  for_rooted_nodes_post (tree, b)
	    {
	      const Phylogeny::Node n = (*b).second;
	      vector<Prob> state_probs = chain.matrix->create_prior();   // by default (e.g. if gapped), assume prior distribution
	      if (!colmat.gapped[n])
		{
		  for_const_contents (vector<int>, colmat.allowed[n], i)
		    {
		      const Prob node_pp = colmat.node_post_prob (n, *i, tree, *chain.matrix);
		      state_probs[*i] = node_pp;
		    }
		}
	      node_state_probs[n] = state_probs;
	    }

	  // for each unattached row, get state probabilities, take outer product with posterior state probs for each tree node, and accumulate
	  vector<Loge> row_loglike (chain_states);   // this is a dummy vector, so don't re-initialize it inside the loop
	  int gapped = 0;  // flag, set to true by init_row if row is gapped
	  for_const_contents (vector<int>, unattached_rows, row)
	    {
	      vector<Prob> row_prob (chain_states, 0.);
	      info.init_row (ecfg.alphabet, chain.classes, asp, *row, coords, row_loglike, &row_prob, gapped, false);
	      if (!gapped)
		for_const_contents (Node_state_prob_map, node_state_probs, nsp)
		  {
		    const Phylogeny::Node node = nsp->first;
		    const vector<Prob>& node_prob = nsp->second;
		    Branch_state_counts& counts = attach_counts[*row][node][chain_idx];
		    for (int i = 0; i < chain_states; ++i)
		      for (int j = 0; j < chain_states; ++j)
			counts(i,j) += node_prob[i] * row_prob[j];
		  }
	    }
	}
    }
}

ECFG_placer::Attachment_map ECFG_placer::best_attachments (double resolution, double tmax, double tmin)
{
  // populate attach_counts
  populate_counts();

  // set up return data structure
  Attachment_map attach_map;

  // for each unattached row: create ECFG_bell_funcs and find mode for each attachment node
  for_const_contents (vector<int>, unattached_rows, row)
    {
      Loge best_loglike = -InfinityLoge;
      Phylogeny::Node best_node = -1;
      double best_t = 0.;

      for_const_contents (ECFG_row_attachment_counts, attach_counts[*row], node_counts)
	{
	  const Phylogeny::Node node = node_counts->first;
	  const ECFG_branch_state_counts& state_counts = node_counts->second;

	  ECFG_bell_funcs bell_funcs (state_counts, ecfg, prior_param, resolution, tmax, tmin);
	  const double t = bell_funcs.bell_max();
	  const Loge loglike = bell_funcs.func(t);

	  CTAG(5,PLACER) << "Row " << tree_align.align.row_name[*row] << " attachment-node " << tree_align.tree.node_specifier(node) << " branch-length " << t << " log-likelihood " << loglike << "\n";

	  if (loglike > best_loglike || best_node < 0)
	    {
	      best_node = node;
	      best_loglike = loglike;
	      best_t = t;
	    }
	}

      CTAG(6,PLACER) << "Row " << tree_align.align.row_name[*row] << " best-attachment-node " << tree_align.tree.node_specifier(best_node) << " branch-length " << best_t << " log-likelihood " << best_loglike << "\n";

      attach_map[*row] = Attachment_branch (best_node, best_t);
    }

  return attach_map;
}
