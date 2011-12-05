#include "handel/alitrans.h"
#include "handel/movement.h"

#define SUBST_PRIOR_SAMPLES 100 /* for subst model */

// prefix character for auto-creating PVar's
#define PVAR_PREFIX_CHAR ':'

void Transducer_alignment::sample_subtree (const vector<int>& nodes_to_sample,
					   const vector<PHYLIP_tree>& available_subtrees,
					   double kT,
					   bool viterbi)
{
  // use Redelings-Suchard proposal scheme for this move?
  const bool redsuch = !viterbi && use_Redelings_Suchard && available_subtrees.size() <= 2;
  Loge hastings_ratio_ll = 0.;

  // do a quick sanity check on the input data
  if (available_subtrees.size() == 0)
    THROWEXPR ("In sample_subtree: urk -- need at least one tree here, man");
  for_const_contents (vector<PHYLIP_tree>, available_subtrees, subtree)
    if (subtree->nodes() != (int) nodes_to_sample.size())
      THROWEXPR ("In sample_subtree: whaaaa..?? Subtree has the wrong number of nodes");

  // make indices for nodes in nodes_to_sample set
  vector<int> sample_index (tree.nodes(), -1);
  for (int i = 0; i < (int) nodes_to_sample.size(); ++i)
    sample_index[nodes_to_sample[i]] = i;

  // loop over available subtrees, getting a Forward likelihood for each one & sampling a Forward path
  vector<PHYLIP_tree> move_tree;
  vector<Loge> move_ll;
  vector<Alignment_path> move_path;
  for (int n_subtree = 0; n_subtree < (int) available_subtrees.size(); ++n_subtree)
    {
      // which subtree is this? (for later use by Redelings-Suchard Metropolis-Hastings sampler)
      const bool is_old_subtree = (n_subtree == 0);  // first tree is "old" tree
      const bool is_new_subtree = (n_subtree == ((int) available_subtrees.size()) - 1);  // second tree is "new" tree

      // first, build a PHYLIP_tree by combining the new tree with unsampled branches from the old tree
      // this new tree will have the same node numbering as the old tree
      const PHYLIP_tree& subtree = available_subtrees[n_subtree];
      PHYLIP_tree new_tree;
      vector<int> constrained_phylip_branch, topologically_identical_phylip_branch;
      for (int n = 0; n < tree.nodes(); ++n)
	new_tree.add_node (-1);
      for (int n = 0; n < tree.nodes(); ++n)
	{
	  const int p = tree.parent[n];
	  const bool n_is_root = p < 0;
	  if (!n_is_root)
	    {
	      const bool sample_p = n_is_root ? false : (sample_index[p] >= 0);
	      const bool sample_n = sample_index[n] >= 0;
	      int new_p = p;
	      double new_branch_length = tree.branch_length (p, n);
	      if (sample_p && sample_n)
		{
		  const int new_p_sample_index = subtree.parent[sample_index[n]];
		  if (new_p_sample_index >= 0)
		    {
		      new_p = nodes_to_sample[new_p_sample_index];
		      new_branch_length = subtree.branch_length (new_p_sample_index, sample_index[n]);
		    }
		}

	      new_tree.add_branch (new_p, n, new_branch_length);

	      if (!(sample_p && sample_n))
		constrained_phylip_branch.push_back (n);
	      else if (new_p == p)
		topologically_identical_phylip_branch.push_back (n);
	    }
	}
      new_tree.node_name = tree.node_name;
      new_tree.root = tree.root;
      new_tree.rebuild_parents();

      // print a wee log message
      if (CTAGGING(6,HANDEL))
	{
	  CL << "Sampling a " << nodes_to_sample.size() << "-node transducer";
	  if (available_subtrees.size() > 1)
	    CL << " (tree #" << n_subtree+1 << " of " << available_subtrees.size() << " for this MCMC move)";
	  if (nodes_to_sample.size())
	    for (int n = 0; n < (int) nodes_to_sample.size(); ++n)
	      CL << (n==0 ? ".\nHere are the nodes I'm sampling: " : ", ") << tree.node_specifier (nodes_to_sample[n]);
	  CL << ".\nHere's the tree I'm using:\n";
	  new_tree.write (CL, 0);
	  CL << '\n';
	}

      if (CTAGGING(2,HANDEL_DEBUG))
	{
	  CL << "Subtree (as passed to sample_subtree):\n";
	  subtree.write (CL, 0);
	  CL << '\n';
	}

      // next, build the representation of the composition that Handel_movement uses.
      vector<int> etree2phylip;
      Alignment_path::Decomposition decomp;
      Handel_movement move = prepare_movement (new_tree, constrained_phylip_branch, topologically_identical_phylip_branch, etree2phylip, decomp);

      if (redsuch)
	{
	  move.evaluate_redelings_suchard_inverse_move = is_old_subtree;
	  move.propose_redelings_suchard_move = is_new_subtree;
	  // currently, centroid banding is not implemented for node-flipping moves, where the old & new alignments have different underlying trees
	  // the book-keeping gets a little more complicated in this case:
	  // 1. we would need to handle the new subtree before the old subtree (the opposite of the way we currently do it)
	  // 2. we would need to pass the old & new centroids to the Handel_movement object, rather than having Handel_movement compute them internally
	  if (use_centroid_band && available_subtrees.size() == 1)
	    {
	      move.use_centroid = true;
	      move.centroid_band_width = centroid_band_width;
	    }
	}
      else if (viterbi)
	move.viterbi = true;
      else
	move.nforward = 1;  // fill Forward matrix & sample exactly one Forward alignment

      // run the Handel_movement
      composition_recorder.save (move);  // executes the Handel_movement & dumps to file or stderr
      dotfile_recorder.save (move.ehmm_scores);

      // check for out-of-band old_path
      // (this is just for absolute bands; we check for adaptive(centroid) band overflows below. this should all be in one place. urgghhhhh)
      if (move.old_path_is_outside_hmmoc_band)
	{
	  CLOGERR << "WARNING: current alignment path is outside HMMoC adapter's banding constraint; "
		  << (hmmoc_opts.strict_banded_reversibility ? "skipping this move" : "proceeding anyway")
		  << '\n';
	  if (hmmoc_opts.strict_banded_reversibility)
	    return;
	}

      // finally, build the new Alignment_path and store it
      Alignment_path new_path;
      const vector<int> free_nodes (move.composition.clique[move.composition.free_clique].begin(),
				    move.composition.clique[move.composition.free_clique].end());
      if (redsuch)
	{
	  // get pairwise branch alignment paths
	  if (is_new_subtree)
	    {
	      for (int free_n = 0; free_n < (int) free_nodes.size(); ++free_n)
		{
		  const int free_p = move.peeled_composition.etree.parent[free_n];
		  if (free_p >= 0)
		    {
		      const int etree_n = free_nodes[free_n];
		      const int etree_p = free_nodes[free_p];

		      const int phylip_n = etree2phylip[etree_n];
		      const int phylip_p = etree2phylip[etree_p];

		      const int row_n = node2row[phylip_n];
		      const int row_p = node2row[phylip_p];

		      const Pair_transducer_funcs& branch_trans = move.composition.branch_trans[etree_n];
		      const vector<int>& branch_trace = move.redsuch_path[etree_n];

		      Pairwise_path branch_path;
		      vector<bool> branch_path_column (2);
		      int max_gap_len = 0, current_ins_len = 0, current_del_len = 0;
		      for (int pos = 0; pos < (int) branch_trace.size(); ++pos)
			{
			  const int s = branch_trace[pos];
			  if (s >= 0)
			    {
			      const int t = branch_trans.state_type[s];
			      for (int n = 0; n < 2; ++n)
				branch_path_column[n] = Transducer_state_type_enum::type_node_emit (t, n - 1);   // euch: "n-1" is because this method is really designed for composite transducers, not pair transducers (and composite transducers have an "extra" node at n=0, representing the subroot)
			      // check for gap length
			      const bool is_ins = !branch_path_column[0] && branch_path_column[1];
			      const bool is_del = branch_path_column[0] && !branch_path_column[1];
			      if (is_ins)
				++current_ins_len;
			      else if (is_del)
				++current_del_len;
			      else
				{
				  max_gap_len = max (max (max_gap_len, current_ins_len), current_del_len);
				  current_ins_len = current_del_len = 0;
				}
			      // append column to branch_path
			      branch_path.append_column (branch_path_column);
			    }
			  max_gap_len = max (max (max_gap_len, current_ins_len), current_del_len);
			}

		      // check if max_gap_len is larger than adaptive band
		      if (max_gap_len > centroid_band_width)
			{
			  CLOGERR << "WARNING: largest gap in current alignment path is greater than HMMoC adapter's centroid band width; "
				  << (hmmoc_opts.strict_banded_reversibility ? "skipping this move" : "proceeding anyway")
				  << '\n';
			  if (hmmoc_opts.strict_banded_reversibility)
			    return;
			}

		      // store
		      decomp[Alignment_path::Row_pair (row_p, row_n)] = branch_path;
		    }
		}

	      // log
	      if (CTAGGING(1,HANDEL_REDSUCH))
		for_const_contents (Alignment_path::Decomposition, decomp, dc)
		  {
		    const int row_p = dc->first.first;
		    const int row_n = dc->first.second;

		    const int phylip_p = row2node[row_p];
		    const int phylip_n = row2node[row_n];

		    vector<sstring> row_name;
		    row_name.push_back (new_tree.node_name[phylip_p]);
		    row_name.push_back (new_tree.node_name[phylip_n]);

		    CL << "Alignment " << row_name[0] << "(" << row_p << ") -> " << row_name[1] << "(" << row_n << "):\n";
		    dc->second.show (CL, row_name);
		  }

	      // compose new path
	      new_path.compose_and_log (decomp, false, true);  // don't squash columns; do use transducer gap ordering

	      // store the move
	      move_tree.push_back (new_tree);
	      move_path.push_back (new_path);
	      move_ll.push_back (0.);  // dummy, just to keep move_ll the same size as move_tree & move_path
	    }

	  // accumulate Hastings ratio terms...
	  for_const_contents (vector<Loge>, move.hastings_term, ht)
	    NatsPMulAcc (hastings_ratio_ll, *ht);
	}
      else
	{
	  // get EHMM state path
	  const vector<int>& trace = viterbi ? move.vit_trace : move.fwd_trace.front();
	  if (trace.size())
	    {
	      if (move.ehmm_scores.nodes() - 1 != (int) free_nodes.size())
		THROWEXPR ("free clique fuckup, dude");

	      // get subtree multiple alignment path
	      Alignment_path ehmm_path (move.ehmm_scores.nodes() - 1, trace.size());
	      for (int pos = 0; pos < (int) trace.size(); ++pos)
		{
		  const int s = trace[pos];
		  if (s >= 0)
		    {
		      const int t = move.ehmm_scores.state_type[s];
		      for (int n = 0; n < move.ehmm_scores.nodes() - 1; ++n)
			ehmm_path[n][pos] = Transducer_state_type_enum::type_node_emit (t, n);
		    }
		}

	      // get pairwise branch alignment paths
	      for (int free_n = 0; free_n < (int) free_nodes.size(); ++free_n)
		{
		  const int free_p = move.peeled_composition.etree.parent[free_n];
		  if (free_p >= 0)
		    {
		      const int etree_n = free_nodes[free_n];
		      const int etree_p = free_nodes[free_p];

		      const int phylip_n = etree2phylip[etree_n];
		      const int phylip_p = etree2phylip[etree_p];

		      const int row_n = node2row[phylip_n];
		      const int row_p = node2row[phylip_p];

		      const Pairwise_path branch_path (ehmm_path, free_p, free_n, true);
		      decomp[Alignment_path::Row_pair (row_p, row_n)] = branch_path;
		    }
		}

	      // compose new path
	      new_path.compose_and_log (decomp, false, true);  // don't squash columns; do use transducer gap ordering

	      // store move
	      move_tree.push_back (new_tree);
	      move_path.push_back (new_path);
	      move_ll.push_back (viterbi ? move.vit_ll : move.fwd_ll);

	      // log
	      CTAG(6,HANDEL) << (viterbi ? "Viterbi" : "Forward") << " score is " << Nats2Bits(move_ll.back()) << " bits\n";
	    }
	}
    }

  // randomly choose a move from the available possibilities
  // for Gibbs sampling, this means selecting a subtree & associated path; for Metropolis-Hastings sampling, this means accepting/rejecting the move
  int n_move;
  bool accept_move;
  if (redsuch)
    {
      // Metropolis-Hastings: accept move with probability min(1,exp(hastings_ratio_ll))
      const Prob hastings_ratio = min (1., Nats2Prob (hastings_ratio_ll));

      n_move = 0;
      accept_move = Rnd::decide (hastings_ratio);

      if (move_path.size() != 1)
	THROWEXPR ("Serious clusterfuck, dude");

      // log
      CTAG(6,HANDEL) << "Hastings ratio = " << hastings_ratio << "; move " << (accept_move ? "accepted" : "rejected") << "\n";
    }
  else
    {
      // Gibbs: sample a subtree
      Loge total_ll = -InfinityLoge;
      for_const_contents (vector<Loge>, move_ll, ll)
	NatsPSumAcc (total_ll, *ll);
      Prob p = Rnd::prob();
      for (n_move = 0; n_move < (int) move_ll.size(); ++n_move)
	{
	  p -= Nats2Prob (NatsPMul (move_ll[n_move], -total_ll));
	  if (p <= 0)
	    break;
	}

      // always accept Gibbs moves
      accept_move = true;
    }

  // rebuild
  if (accept_move && n_move < (int) move_tree.size())
    {
      tree = move_tree[n_move];
      align.path = move_path[n_move];
      tree_changed();
      align_changed();
    }
}

Handel_movement Transducer_alignment::prepare_movement (const PHYLIP_tree& new_tree, const vector<int>& constrained_phylip_branch, const vector<int>& topologically_identical_phylip_branch, vector<int>& etree2phylip, Alignment_path::Decomposition& decomp)
{
  // first, make an alphabet for Transducer_SExpr_file
  const vector<sstring>& trans_alph = alphabet().tokens();

  // next, build the representation of the composition that Handel_movement uses.
  // re-order the nodes to ETree node numbering scheme
  ETree etree (0);
  vector<int> phylip2etree (new_tree.nodes(), -1);
  etree2phylip = vector<int> (new_tree.nodes(), -1);
  for_rooted_nodes_pre (new_tree, phylip_b)
    {
      const int
	phylip_p = (*phylip_b).first,
	phylip_n = (*phylip_b).second;
      const int etree_p = phylip_p < 0 ? -1 : phylip2etree[phylip_p];
      const int etree_n = etree.parent.size();
      etree.parent.push_back (etree_p);
      phylip2etree[phylip_n] = etree_n;
      etree2phylip[etree_n] = phylip_n;
    }

  // set up the PScores & Pair_transducer_funcs
  PScores pscores;
  vector<Pair_transducer_funcs> pair_trans (etree.nodes(), Pair_transducer_funcs (0));  // indexed by ETree node index
  const Pair_transducer_scores ppts = prior_pair_trans_sc();
  sstring root_pvar_prefix;
  root_pvar_prefix << safe_tape_name (new_tree, 0) << PVAR_PREFIX_CHAR;
  pair_trans[0] = Pair_transducer_funcs (ppts, pscores, trans_alph, root_pvar_prefix.c_str());

  for (int etree_n = 1; etree_n < etree.nodes(); ++etree_n)
    {
      const int etree_p = etree.parent[etree_n];
      const int phylip_p = etree2phylip[etree_p];
      const int phylip_n = etree2phylip[etree_n];

      const sstring tape_name = safe_tape_name (new_tree, phylip_n);
      CTAG(1,HANDEL_MOVEMENT) << "Preparing transducer for node " << tape_name << '\n';

      const Pair_transducer_scores bpts = branch_pair_trans_sc (new_tree.branch_length (phylip_p, phylip_n));
      sstring pvar_prefix;
      pvar_prefix << tape_name << PVAR_PREFIX_CHAR;
      pair_trans[etree_n] = Pair_transducer_funcs (bpts, pscores, trans_alph, pvar_prefix.c_str());
    }

  // find Insert node for root branch transducer
  int root_ins_state = -1;
  for (int s = 0; s < pair_trans[0].states(); ++s)
    if (pair_trans[0].state_type[s] == TransducerInsertType)
      {
	root_ins_state = s;
	break;
      }

  // make the NodePathMap and the Decomp
  Transducer_SExpr_file::NodePathMap node_path, old_node_path;
  decomp.clear();
  for_const_contents (vector<int>, constrained_phylip_branch, phylip_n)
    {
      const int old_phylip_p = tree.parent[*phylip_n];
      const int row_p = node2row[old_phylip_p];
      const int row_n = node2row[*phylip_n];

      if (row_p < 0)
	THROWEXPR ("Couldn't find " << tree.node_name[old_phylip_p] << " in alignment, but it's in the tree");

      if (row_n < 0)
	THROWEXPR ("Couldn't find " << tree.node_name[*phylip_n] << " in alignment, but it's in the tree");

      const int etree_n = phylip2etree[*phylip_n];
      const Pairwise_path branch_path (align.path, row_p, row_n, true);
      decomp[Alignment_path::Row_pair (row_p, row_n)] = branch_path;

      // make mapping of state type -> index for this branch
      map<int,int> state_index;
      for (int s = 0; s < pair_trans[etree_n].states(); ++s)
	state_index[pair_trans[etree_n].state_type[s]] = s;

      // convert Pairwise_path into state path
      vector<int> state_path;
      int parent_len = 0;
      for (int col = 0; col < branch_path.columns(); ++col)
	{
	  const bool anc = branch_path(0,col), des = branch_path(1,col);
	  const int state_type =
	    anc
	    ? (des ? TransducerMatchType : TransducerDeleteType)
	    : (des ? TransducerInsertType : TransducerWaitType);
	  state_path.push_back (state_index[state_type]);
	  if (anc)
	    ++parent_len;
	}
      node_path[etree_n] = state_path;

      // set root branch implicitly
      if (old_phylip_p == tree.root)
	node_path[0] = vector<int> (parent_len, root_ins_state);
    }

  // if all branches are either constrained or topologically identical to the old tree, then make the NodePathMap for the old path
  // CODE SMELL: this code is more-or-less identical to the previous block (for setting up node_path), and the two blocks should be consolidated
  if ((int) (constrained_phylip_branch.size() + topologically_identical_phylip_branch.size()) == tree.branches())
    for_const_contents (vector<int>, topologically_identical_phylip_branch, phylip_n)
      {
	const int phylip_p = tree.parent[*phylip_n];
	const int row_p = node2row[phylip_p];
	const int row_n = node2row[*phylip_n];
	const int etree_n = phylip2etree[*phylip_n];
	const Pairwise_path branch_path (align.path, row_p, row_n, true);

	// make mapping of state type -> index for this branch
	map<int,int> state_index;
	for (int s = 0; s < pair_trans[etree_n].states(); ++s)
	  state_index[pair_trans[etree_n].state_type[s]] = s;

	// convert Pairwise_path into state path
	vector<int> state_path;
	int parent_len = 0;
	for (int col = 0; col < branch_path.columns(); ++col)
	  {
	    const bool anc = branch_path(0,col), des = branch_path(1,col);
	    const int state_type =
	      anc
	      ? (des ? TransducerMatchType : TransducerDeleteType)
	      : (des ? TransducerInsertType : TransducerWaitType);
	    state_path.push_back (state_index[state_type]);
	    if (anc)
	      ++parent_len;
	  }
	old_node_path[etree_n] = state_path;

	// set root branch implicitly
	if (phylip_p == tree.root)
	  old_node_path[0] = vector<int> (parent_len, root_ins_state);
      }

  // make the NodeProfileMap
  Transducer_SExpr_file::NodeProfileMap node_prof;
  for_iterator (Phylogeny::Node_const_iter, phylip_n, tree.leaves_begin(), tree.leaves_end())
    {
      const int row = node2row[*phylip_n];
      const Score_profile* orig_prof_sc = align.prof[row];
      if (orig_prof_sc)
	node_prof[phylip2etree[*phylip_n]] = Score_profile (*orig_prof_sc);
    }

  // add the banding coefficients
  map<int,double> band_coeff;
  if (use_banding_coefficient)
    {
      const double gr = gap_rate(), gs = mean_gap_size();
      for (int etree_n = 1; etree_n < etree.nodes(); ++etree_n)
	{
	  const int etree_p = etree.parent[etree_n];
	  const int phylip_p = etree2phylip[etree_p];
	  const int phylip_n = etree2phylip[etree_n];
	  const double branch_length = new_tree.branch_length (phylip_p, phylip_n);
	  band_coeff[etree_n] = banding_coefficient * (1. - exp (-gr * branch_length)) * gs;
	}
    }

  // initialize the Handel_movement
  Handel_movement move;
  move.hmmoc_opts = hmmoc_opts;

  move.composition = Transducer_SExpr_file (trans_alph, pscores, etree, pair_trans, node_prof, node_path, old_node_path);
  move.composition.band_coeff = band_coeff;

  move.quiet = false;  // show composition
  move.composite = false;  // don't show composite transducer
  move.acyclic = false;  // don't show acylic transducer

  move.viterbi = false;  // don't do Viterbi
  move.nforward = -1;  // don't do Forward
  move.optacc = false;  // don't do posterior decoding
  move.want_expected_counts = false;  // don't do Backward/peeling
  move.propose_redelings_suchard_move = false;  // don't do Redelings-Suchard proposal
  move.evaluate_redelings_suchard_inverse_move = false;  // don't do Redelings-Suchard inverse proposal

  // set the tape names
  map<int,sstring> new_tape_name;
  vector<sstring> dummy_vec;
  for (int etree_n = 0; etree_n < etree.nodes(); ++etree_n)
    {
      const sstring& nn = new_tree.node_name[etree2phylip[etree_n]];
      if (nn.size())
	new_tape_name[etree_n] = nn;
    }

  // have Transducer_SExpr_file set the branch names & transducer names automatically
  move.composition.rebuild_tree_names (new_tape_name, dummy_vec, dummy_vec);

  // return
  return move;
}

void Transducer_alignment::dump (ostream &out)
{
  // constrain all branches
  vector<int> constrained_phylip_branch, topologically_identical_phylip_branch;
  for_rooted_branches_pre (tree, b)
    constrained_phylip_branch.push_back ((*b).second);

  // build the representation of the composition that Handel_movement uses.
  vector<int> dummy_etree2phylip;
  Alignment_path::Decomposition dummy_decomp;
  Handel_movement move = prepare_movement (tree, constrained_phylip_branch, topologically_identical_phylip_branch, dummy_etree2phylip, dummy_decomp);

  // dump the composition and return
  move.dump_composition (out);
  return;
}

Score Transducer_alignment::conditioned_branch_path_score (const Node_pair& branch) const
{
  Pair_transducer_scores trans = ((Transducer_alignment&)*this).branch_pair_trans_sc (tree.branch_length (branch));  // cast away const
  Pairwise_path path (align.path, node2row[branch.first], node2row[branch.second], true);
  if (sort_indels)
    path.sort_indels();
  return trans.pairwise_path_score (path);
}

void Transducer_alignment::propose_sample_node (Node node, double kT)
{
  propose_sample_or_optimise_node (node, kT, false);
}

bool Transducer_alignment::propose_optimise_node (Node node)
{
  return propose_sample_or_optimise_node (node, 1., true);
}

bool Transducer_alignment::propose_sample_or_optimise_node (Node node, double kT, bool optimise)
{
  const Alignment_path old_path = align.path;

  vector<int> nodes_to_sample;
  PHYLIP_tree subtree;
  // is this node the root?
  if (node == tree.root)
    {
      nodes_to_sample.push_back (node);
      subtree.add_node();
      subtree.node_name.push_back (tree.node_name[node]);

      for_rooted_children (tree, node, c)
	{
	  nodes_to_sample.push_back (*c);
	  subtree.add_node (0, tree.branch_length (node, *c));
	  subtree.node_name.push_back (tree.node_name[*c]);
	}
    }
  else
    {
      nodes_to_sample.push_back (tree.parent[node]);
      subtree.add_node();
      subtree.node_name.push_back (tree.node_name[tree.parent[node]]);

      nodes_to_sample.push_back (node);
      subtree.add_node (0, tree.branch_length (tree.parent[node], node));
      subtree.node_name.push_back (tree.node_name[node]);

      for_rooted_children (tree, node, c)
	{
	  nodes_to_sample.push_back (*c);
	  subtree.add_node (1, tree.branch_length (node, *c));
	  subtree.node_name.push_back (tree.node_name[*c]);
	}
    }

  subtree.root = 0;
  subtree.rebuild_parents();
  const vector<PHYLIP_tree> available_subtrees (1, subtree);

  sample_subtree (nodes_to_sample, available_subtrees, kT, optimise);

  if (CTAGGING(1,HANDEL_PATH))
    {
      CL << "Old path:\n";
      old_path.show (CL, align.row_name);
      CL << "New path:\n";
      align.path.show (CL, align.row_name);
      CL << "I think these are " << (align.path != old_path ? "different" : "the same") << ".\n";
    }
  return align.path != old_path;
}

void Transducer_alignment::propose_sample_branch (const Undirected_pair& branch, double kT)
{
  propose_sample_or_optimise_branch (branch, kT, false);
}

bool Transducer_alignment::propose_optimise_branch (const Undirected_pair& branch)
{
  return propose_sample_or_optimise_branch (branch, 1., true);
}

bool Transducer_alignment::propose_sample_or_optimise_branch (const Undirected_pair& branch, double kT, bool optimise)
{
  const Alignment_path old_path = align.path;

  vector<int> nodes_to_sample (2);
  nodes_to_sample[0] = branch.first;
  nodes_to_sample[1] = branch.second;

  PHYLIP_tree subtree;
  subtree.add_node();
  subtree.node_name.push_back (tree.node_name[branch.first]);

  subtree.add_node (0, tree.branch_length (branch));
  subtree.node_name.push_back (tree.node_name[branch.second]);

  subtree.root = 0;
  subtree.rebuild_parents();
  const vector<PHYLIP_tree> available_subtrees (1, subtree);

  sample_subtree (nodes_to_sample, available_subtrees, kT, optimise);

  if (CTAGGING(1,HANDEL_PATH))
    {
      CL << "Old path:\n";
      old_path.show (CL, align.row_name);
      CL << "New path:\n";
      align.path.show (CL, align.row_name);
      CL << "I think these are " << (align.path != old_path ? "different" : "the same") << ".\n";
    }
  return align.path != old_path;
}

bool Transducer_alignment::propose_sample_branch_swap (Node aunt, Node nephew, Node grumpa, Node dad, double kT)
{
  vector<int> nodes_to_sample;
  PHYLIP_tree old_tree, new_tree;

  // find Ella! (sibling of nephew)
  Node ella = -1;
  for_rooted_children (tree, dad, c)
    if (*c != nephew)
      ella = *c;
  if (ella < 0)
    THROWEXPR ("Couldn't find sibling node");

  // Old tree:          New tree:
  //
  //      0:grumpa         0:grumpa
  //        /  |             /  |
  //    1:dad  aunt:2    1:dad  nephew:2
  //      / |              / |
  // 3:ella nephew:4  3:ella aunt:4

  nodes_to_sample.push_back (grumpa);
  old_tree.add_node();
  new_tree.add_node();
  old_tree.node_name.push_back (tree.node_name[grumpa]);
  new_tree.node_name.push_back (tree.node_name[grumpa]);

  nodes_to_sample.push_back (dad);
  old_tree.add_node (0, tree.branch_length (grumpa, dad));
  new_tree.add_node (0, tree.branch_length (grumpa, dad));
  old_tree.node_name.push_back (tree.node_name[dad]);
  new_tree.node_name.push_back (tree.node_name[dad]);

  nodes_to_sample.push_back (aunt);
  old_tree.add_node (0, tree.branch_length (grumpa, aunt));
  new_tree.add_node (1, tree.branch_length (grumpa, aunt));
  old_tree.node_name.push_back (tree.node_name[aunt]);
  new_tree.node_name.push_back (tree.node_name[aunt]);

  nodes_to_sample.push_back (ella);
  old_tree.add_node (1, tree.branch_length (dad, ella));
  new_tree.add_node (1, tree.branch_length (dad, ella));
  old_tree.node_name.push_back (tree.node_name[ella]);
  new_tree.node_name.push_back (tree.node_name[ella]);

  nodes_to_sample.push_back (nephew);
  old_tree.add_node (1, tree.branch_length (dad, nephew));
  new_tree.add_node (0, tree.branch_length (dad, nephew));
  old_tree.node_name.push_back (tree.node_name[nephew]);
  new_tree.node_name.push_back (tree.node_name[nephew]);

  old_tree.root = new_tree.root = 0;
  old_tree.rebuild_parents();
  new_tree.rebuild_parents();

  vector<PHYLIP_tree> available_subtrees;
  available_subtrees.push_back (old_tree);  // old tree is tree#0
  available_subtrees.push_back (new_tree);  // new tree is tree#1

  sample_subtree (nodes_to_sample, available_subtrees, kT);

  return true;
}

void Transducer_alignment::propose_optimise_branch_length (const Undirected_pair& branch,
							   double tmax, double tres)
{
  THROWEXPR ("Transducer_alignment::propose_optimise_branch_length unimplemented");
}

void Transducer_alignment::align_and_infer_parent (const Score_profile& xseq,
						   const Score_profile& yseq,
						   Node ancestor_node,
						   Alignment_path& axy_path)
{
  THROWEXPR ("Transducer_alignment::align_and_infer_parent unimplemented");
}

sstring Transducer_alignment::safe_tape_name (const PHYLIP_tree& tree, int node)
{
  set<sstring> names;
  sstring name;
  for (int n = 0; n <= node; ++n)
    {
      name.clear();
      if (tree.node_name[n].size())
	for_const_contents (sstring, tree.node_name[n], c)
	  if (*c != '(' && *c != ')' && *c != ';' && !isspace(*c))
	    name << *c;
      else
	name = Transducer_SExpr_file::auto_tape_name (n);
      while (names.find (name) != names.end())
	name << '+';
      names.insert (name);
    }
  return name;
}

void Transducer_alignment_with_subst_model::assign_subst_model (const Transducer_alignment_with_subst_model& ta)
{
  if (ems) delete ems;
  ems = new ECFG_matrix_set (*ta.ems);
  subst_pscores = ta.subst_pscores;
  subst_pcounts = ta.subst_pcounts;
  subst_mutable_pgroups = ta.subst_mutable_pgroups;
  eval_funcs();
}

Substitution_matrix_factory& Transducer_alignment_with_subst_model::submat_factory() const {
  Transducer_alignment_with_subst_model* mutable_this = (Transducer_alignment_with_subst_model*) this;
  mutable_this->eval_funcs();
  return (Substitution_matrix_factory&) (mutable_this->subst_model());
}

EM_matrix_base& Transducer_alignment_with_subst_model::subst_model()
{
  if (ems->chain.size() != 1)
    THROWEXPR("In Transducer_alignment_with_subst_model::subst_model, ECFG_matrix_set is uninitialized");
  return *(ems->chain[0].matrix);
}

sstring Transducer_alignment_with_subst_model::subst_parameter_string() const
{
  return PFunc_builder::mutable_pscores2string (subst_pscores, subst_mutable_pgroups);
}

PScores Transducer_alignment_with_subst_model::propose_subst_params()
{
  PScores proposal (subst_pscores);
  subst_pcounts.randomize (proposal, subst_mutable_pgroups);
  return proposal;
}

void Transducer_alignment_with_subst_model::sample_subst_params()
{
  const int sample_points = SUBST_PRIOR_SAMPLES;

  vector<PScores> x (sample_points + 1);
  vector<Prob> p (sample_points + 1);
  vector<Score> sc (sample_points + 1);

  x[0] = subst_pscores;
  for (int i = 1; i <= sample_points; ++i)
    x[i] = propose_subst_params();

  for (int i = 0; i <= sample_points; ++i)
    {
      subst_pscores = x[i];
      eval_funcs();
      tree_changed();
      sc[i] = alignment_emit_score();

      if (CTAGGING(5,MCMC PARAM_SAMPLE))
	CL << "Subst-param sample " << i << ": " << Score2Bits(sc[i]) << " bits " << subst_parameter_string() << '\n';
    }

  Score min_sc = sc[0];
  for (int i = 1; i <= sample_points; ++i)
    if (sc[i] < min_sc)
      min_sc = sc[i];

  for (int i = 0; i <= sample_points; ++i)
    p[i] = Score2Prob (sc[i] - min_sc);

  const int i = Rnd::choose (p);
  subst_pscores = x[i];
  eval_funcs();
  tree_changed();

  if (CTAGGING(5,MCMC PARAM_SAMPLE))
    CL << "Chose sample " << i << ": " << Score2Bits(sc[i]) << " bits " << subst_parameter_string() << '\n';
}
