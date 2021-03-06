#include "ecfg/ecfgplacer.h"

ECFG_placer::ECFG_placer (ECFG_scores& ecfg, Stockholm& stock, Tree_alignment& tree_align, double prior_param) :
  ecfg (ecfg),
  stock (stock),
  tree_align (tree_align),
  unattached_rows (tree_align.unattached_rows()),
  prior_param (prior_param)
{ }

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
    for_rooted_nodes_post (tree_align.tree, b)
    {
      const Phylogeny::Node node = (*b).second;
      attach_counts[*row][node] = empty_counts;
    }

  // get CYK matrix & trace
  const int max_subseq_len = -1;   // unlimit subsequence lengths
  ECFG_auto_envelope env (stock, ecfg, max_subseq_len);  // create fold envelope

  Aligned_score_profile asp;
  asp.init (tree_align.align, stock.np, ecfg.alphabet);

  const bool try_fast_prune = false;
  ECFG_CYK_matrix cyk_mx (ecfg, stock, asp, env, try_fast_prune);   // create CYK matrix
  cyk_mx.fill();
  ECFG_cell_score_map cyk_trace = cyk_mx.traceback();  // get traceback

  // ensure Tree_alignment's tree stays in sync with ECFG DP matrix's tree
  tree_align.set_tree (cyk_mx.tree);

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
	  cyk_mx.fill_down (cyk_mx.env.find_subseq_idx (coords.start, coords.len), ecfg_state);

	  // get posterior state probabilities for each node
	  EM_matrix_base::Column_matrix& colmat = cyk_mx.colmat[ecfg_state];

	  typedef map<Phylogeny::Node,vector<Prob> > Node_state_prob_map;
	  Node_state_prob_map node_state_probs;
	  for_rooted_nodes_post (cyk_mx.tree, b)
	    {
	      const Phylogeny::Node n = (*b).second;
	      vector<Prob> state_probs (chain_states, 0.);
	      if (colmat.gapped[n])
		state_probs = chain.matrix->create_prior();   // by default (if gapped), assume prior distribution
	      else
		for_const_contents (vector<int>, colmat.allowed[n], i)
		  {
		    const Prob node_pp = colmat.node_post_prob (n, *i, cyk_mx.tree, *chain.matrix);
		    state_probs[*i] = node_pp;
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
		    // CTAG(1,PLACER_DEBUG) << "node="<<cyk_mx.tree.node_specifier(node)<<" row="<<tree_align.align.row_name[*row]<<" node_prob=("<<node_prob<<") row_prob=("<<row_prob<<")\n";
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

void ECFG_placer::attach (double resolution, double tmax, double tmin)
{
  Attachment_map attach_map = best_attachments (resolution, tmax, tmin);
  for_const_contents (Attachment_map, attach_map, row_attachment_branch)
    {
      const int row = row_attachment_branch->first;
      const Attachment_branch& attachment_branch = row_attachment_branch->second;

      // we are updating the Tree_alignment's tree in place here, which is a bit sketchy... IH 4/15/2010
      const Phylogeny::Node new_node = tree_align.tree.add_named_node (tree_align.align.row_name[row], attachment_branch.first, attachment_branch.second);

      tree_align.row2node[row] = new_node;
      tree_align.node2row.push_back (row);
    }

  tree_align.tree_changed();
}

bool ECFG_attachable_tree_alignment_database::has_unattached_rows() const
{
  for_const_contents (list<Tree_alignment>, tree_align_list, tree_align)
    if (tree_align->has_unattached_rows())
      return true;
  return false;
}

void ECFG_attachable_tree_alignment_database::attach_rows (ECFG_scores& ecfg, double prior_param,
							   double time_resolution, double time_max, double time_min)
{
  list<Stockholm>::iterator stock_iter = stock_db->align.begin();
  for_contents (list<Tree_alignment>, tree_align_list, tree_align)
    {
      if (tree_align->has_unattached_rows())
	{
	  ECFG_placer placer (ecfg, *stock_iter, *tree_align, prior_param);
	  placer.attach (time_resolution, time_max, time_min);
	}
      ++stock_iter;
    }
}
