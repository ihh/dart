#include "ecfg/ecfgdp.h"
#include "hmm/singletmpl.h"
#include "util/vector_output.h"
#include "ecfg/ecfgsexpr.h"

#if defined(BEAGLE_INCLUDED) && BEAGLE_INCLUDED
#include "libhmsbeagle/beagle.h"
#endif /* BEAGLE_INCLUDED */


#define ECFGDP_REPORT_INTERVAL 100000  /* number of bifurcations between logfile messages during DP */
#define ECFG_MIN_TRACEBACK_ERROR_TOLERANCE  1e-3  /* max permissible error for ECFG traceback */

// heuristic guard against potential underflow errors in Fast_prune:
// reject fast pruning if alignment has over 100 rows.
//
// this was first empirically set to 500, to handle 5S_rRNA, where I first ran into underflow problems
// it was then updated to 100 to be more conservative after I had some problems with smaller alignments
// it's almost certainly not general enough: underflow depends on alphabet size, rate matrix, branch lengths ...
//
// should probably do something a bit more careful
//  -- like trying with fast pruning, then retrying if underflow occurs.
//
// IH, 2/8/2007 (in Melbourne)
inline const bool alignment_safe_from_fast_prune_underflow (const Stockholm& stock)
{
  return stock.rows() < 100;  // prevent Fast_prune underflow
}

ECFG_matrix::ECFG_matrix (const ECFG_scores& ecfg, Stockholm& stock, const ECFG_envelope& env)
  : ecfg (ecfg),
    env (env),
    stock (stock),
    tree (stock),
    final_loglike (-InfinityLoge),
    bifs_done (0),
    last_bifs_done (0),
    total_bifs (env.n_bif()),
    lineage_chain_index (ecfg.matrix_set.chain.size()),
    lineage_matrix (ecfg.matrix_set.chain.size())
{
  // display ECFG
  if (CTAGGING(1,ECFG_SCORES ECFGDP ECFG_SHOW))
    {
      CL << "About to allocate DP matrix for the following ECFG:\n";
      ecfg.show (CL);
    }

  // allocate DP matrix
  if (CTAGGING(3,ALLOC))
    CL << "Allocating " << sizeof(Loge)*env.subseqs()*ecfg.states() << " bytes for DP matrix ("
       << env.subseqs() << " subseqs * " << ecfg.states() << " states)\n";
  cell.resize (env.subseqs(), ecfg.states(), -InfinityLoge);  // default cell score is -infinity

  // prepare ECFG
  ((ECFG_scores&) ecfg).eval_funcs();  // cast away const (lazy hack, to ensure PFunc's are eval'd before any DP is done)

  // initialise ECFG summary data
  emit_states = ecfg.emit_states();
  // nonemit_states are sorted; see Transition_methods::topological_sort()
  nonemit_states = ecfg.nonemit_states();
  // inside algorithm fills null states in reverse topological order
  inside_fill_states = emit_states;
  inside_fill_states.insert (inside_fill_states.end(), nonemit_states.rbegin(), nonemit_states.rend());
  // outside algorithm fills null states in topological order
  outside_fill_states = emit_states;
  outside_fill_states.insert (outside_fill_states.end(), nonemit_states.begin(), nonemit_states.end());
  // incoming & outgoing states
  incoming = Transition_methods::selected_incoming_states (ecfg, outside_fill_states);
  outgoing = Transition_methods::selected_outgoing_states (ecfg, inside_fill_states);
  // incoming & outgoing transition log-likelihoods
  incoming_ll = vector<vector<Loge> > (ecfg.states(), vector<Loge>());
  outgoing_ll = vector<vector<Loge> > (ecfg.states(), vector<Loge>());
  start_ll = vector<Loge> (ecfg.states());
  end_ll = vector<Loge> (ecfg.states());
  for (int s = 0; s < ecfg.states(); ++s)
    {
      for (int i = 0; i < (int) incoming[s].size(); ++i)
	incoming_ll[s].push_back (trans_ll (incoming[s][i], s));
      for (int i = 0; i < (int) outgoing[s].size(); ++i)
	outgoing_ll[s].push_back (trans_ll (s, outgoing[s][i]));
      start_ll[s] = trans_ll (Start, s);
      end_ll[s] = trans_ll (s, End);
    }
  // bifurcations
  left_bifurc = ecfg.left_bifurc();
  right_bifurc = ecfg.right_bifurc();
  // log
  if (CTAGGING(-1,ECFG_SORT_STATES))
    {
      CL << "ECFG emit states: " << emit_states << "\n";
      CL << "ECFG nonemit states: " << nonemit_states << "\n";
      CL << "ECFG inside fill states: " << inside_fill_states << "\n";
      CL << "ECFG outside fill states: " << outside_fill_states << "\n";
    }

  // initialize hybrid chain stuff
  for (int c = 0; c < (int) ecfg.matrix_set.chain.size(); ++c)
    {
      const ECFG_chain& chain = ecfg.matrix_set.chain[c];
      if (chain.matrix)
	{
	  lineage_chain_index[c] = vector<int> (tree.nodes(), (int) c);
	  lineage_matrix[c] = vector<const EM_matrix_base*> (tree.nodes(), chain.matrix);
	}
      else
	{
	  const sstring& default_gs_val = chain.gs_values[0];
	  const int default_chain_index = ((map<sstring,int>&) chain.gs_tag_value_chain_index)[default_gs_val];  // cast away const
	  lineage_chain_index[c] = vector<int> (tree.nodes(), (int) default_chain_index);
	  lineage_matrix[c] = vector<const EM_matrix_base*> (tree.nodes(), ecfg.matrix_set.chain[default_chain_index].matrix);
	  for (int n = 0; n < tree.nodes(); ++n)
	    if (tree.node_name[n].size())
	      {
		sstring gs_val;

		// if the gs_tag is an implicit one, deduce gs_val automagically
		if (chain.gs_tag == ECFG_IMPLICIT_GS_TAG_NODENAME)
		  gs_val = tree.node_name[n];
		else if (chain.gs_tag[0] == ECFG_IMPLICIT_GS_TAG_EQUALS || chain.gs_tag[0] == ECFG_IMPLICIT_GS_TAG_ANCESTOR)
		  {
		    const sstring gs_tag_node_name (chain.gs_tag.begin() + 1, chain.gs_tag.end());
		    bool match;
		    if (chain.gs_tag[0] == ECFG_IMPLICIT_GS_TAG_EQUALS)
		      match = tree.node_name[n] == gs_tag_node_name;
		    else  // ECFG_IMPLICIT_GS_TAG_ANCESTOR
		      {
			match = false;
			for (int p = n; p >= 0 && !match; p = tree.parent[p])
			  if (tree.node_name[p] == gs_tag_node_name)
			    match = true;
		      }
		    gs_val = match ? ECFG_IMPLICIT_GS_VALUE_TRUE : ECFG_IMPLICIT_GS_VALUE_FALSE;
		  }
		else
		  {
		    Stockholm::Annotation gs_annot = stock.gs_annot[tree.node_name[n]];
		    if (gs_annot.find (chain.gs_tag) != gs_annot.end())
		      gs_val = gs_annot[chain.gs_tag];
		  }

		// select the appropriate chain using gs_val
		if (chain.gs_tag_value_chain_index.find (gs_val) != chain.gs_tag_value_chain_index.end())
		  {
		    const int chain_index = ((map<sstring,int>&) chain.gs_tag_value_chain_index)[gs_val];  // cast away const
		    lineage_chain_index[c][n] = chain_index;
		    lineage_matrix[c][n] = ecfg.matrix_set.chain[chain_index].matrix;
		  }
		else
		  CLOGERR << "Warning: no component chain with "
			  << Stockholm_sequence_annotation
			  << " label '" << gs_val
			  << "' in hybrid chain (" << chain.state
			  << "); using default label '" << default_gs_val << "' at node " << tree.node_name[n] << "\n";

	      }
	}
      if (CTAGGING(3,HYBRID_CHAIN))
	CL << "lineage_chain_index[" << c << "] = (" << lineage_chain_index[c] << ")\n";
    }
}

void ECFG_matrix::show (ostream& out) const
{
  save_flags (out);
  out << "States:";
  for (int state = 0; state < ecfg.states(); ++state)
    out << ' ' << ecfg.state_info[state].name << '(' << state << ')';
  out << '\n';
  out << "Subseq       State:";
  for (int state = 0; state < ecfg.states(); ++state)
    {
      out.width(11);
      right_align(out);
      out << ecfg.desc (state);
    }
  out << "\n";
  for (int subseq = 0; subseq < cell.xsize(); ++subseq)
    {
      out << '#';
      sstring subseq_desc;
      subseq_desc << subseq << " (" << env.subseq[subseq].start << '+' << env.subseq[subseq].len << "): ";
      out.width(18);  // = 5+2+4+1+4+3
      left_align(out);
      out << subseq_desc;
      right_align(out);
      for (int state = 0; state < cell.ysize(); ++state)
	{
	  out << ' ';
	  out.width (10);
	  ShowScore (Nats2Score(cell(subseq,state)), out);
	}
      out << " (et)\n";
      show_emit (subseq, out);
    }
  out << " (e) = emit score, (et) = emit score + transition scores\n";
  out << "Final score: " << Nats2Score(final_loglike) << '\n';
  restore_flags (out);
}

void ECFG_matrix::show_progress (const char* algorithm_name, int subseq_idx, bool outside)
{
  bifs_done += env.n_bif (subseq_idx);
  const int total_subseqs = env.subseqs();
  const int subseqs_done = outside ? (total_subseqs - subseq_idx) : (subseq_idx + 1);
  if (bifs_done - last_bifs_done > ECFGDP_REPORT_INTERVAL / tree.nodes()
      || subseqs_done == total_subseqs)
    {
      CTAG(5,ECFGDP) << algorithm_name << ": finished " << subseqs_done << " of " << total_subseqs << " subseqs ("
		     << ((int)(1000.*(double)subseqs_done/(double)total_subseqs))/10.
		     << "%) and " << bifs_done << " of " << total_bifs << " bifurcations ("
		     << ((int)(1000.*bifs_done/total_bifs))/10. << "%)\n";
      last_bifs_done = bifs_done;
    }
}

void ECFG_EM_matrix::use_precomputed (Emit_loglike_matrix& el)
{
  swap (emit_loglike, el);
  fill_up_flag = false;
}

void ECFG_EM_matrix::show_emit (int subseq, ostream& out) const
{
  out.width(19);
  left_align(out);
  out << "";
  right_align(out);
  for (int state = 0; state < ecfg.states(); ++state)
    {
      const int sz = ecfg.state_info[state].emit_size();
      out << ' ';
      out.width (10);
      if (sz && env.subseq[subseq].len >= sz)
	ShowScore (Nats2Score (emit_loglike (subseq, state)), out);
      else
	out << "";
    }
  out << "  (e)\n";
}

ECFG_EM_matrix::ECFG_EM_matrix (const ECFG_scores& ecfg, Stockholm& stock,
				const Aligned_score_profile& asp, const ECFG_envelope& env,
				bool use_fast_prune)
  : ECFG_matrix (ecfg, stock, env),
    asp (asp),
    ll_row (ecfg.states()),
    prob_row (ecfg.states()),
    use_fast_prune (use_fast_prune),
    fast_prune (ecfg.states()),
    fill_up_flag (true)

{
  // debug: switch off fast_prune via logging
  if (CTAGGING(-99,NO_FAST_PRUNE))
    {
      this->use_fast_prune = use_fast_prune = false;
      CL << "Switched off fast_prune\n";
    }

  // allocate emit_loglike
  if (CTAGGING(3,ALLOC))
    CL << "Allocating " << sizeof(Loge)*env.subseqs()*ecfg.states() << " bytes for emit scores ("
       << env.subseqs() << " subseqs * " << ecfg.states() << " states)\n";
  emit_loglike.resize (env.subseqs(), ecfg.states(), 0.);  // default emit_loglike is zero

  // create colmat, ll_row, prob_row
  Column_matrix tmp_colmat;
  colmat = vector<Column_matrix> (ecfg.states(), tmp_colmat);

  // allocate colmat, ll_row, prob_row
  if (CTAGGING(-1,ECFG_EM_MATRIX))
    asp.show (CL);
  for (int s = 0; s < ecfg.states(); ++s)
    if (ecfg.state_info[s].emit_size())
      {
	const int chain_idx = ecfg.state_info[s].matrix;
	const int chain_states = ecfg.matrix_set.total_states (chain_idx);
	colmat[s].alloc (tree.nodes(),
			 chain_states,
			 false); // don't allocate class labels
	if (use_fast_prune)
	  fast_prune[s].prepare (tree, lineage_matrix[chain_idx], colmat[s]);
	ll_row[s] = vector<Loge> (chain_states);
	prob_row[s] = vector<Prob> (chain_states);
      }

  // initialise feature lookup
  set<sstring> feature = ecfg.gc_feature_set();
  align_annot = vector<const char*> ((int) feature.size(), (const char*) 0);
  state_annot = vector<vector<map<sstring,Loge> > > (ecfg.states(), vector<map<sstring,Loge> > ((int) feature.size()));
  int n_feature = 0;
  for_const_contents (set<sstring>, feature, f)  
    {
      CTAG(1,ECFG_EM_MATRIX) << "State annotations for '#=GC " << *f << "':";
      for (int s = 0; s < ecfg.states(); ++s)
	{
	  const ECFG_state_info& info = ecfg.state_info[s];
	  const ECFG_state_info::ECFG_state_annotation::const_iterator a = info.annot.find (*f);
	  if (a != info.annot.end())
	    {
	      for_const_contents (ECFG_state_info::String_prob_dist, a->second, sp) {
		const Score sc = sp->second.eval_sc (ecfg.pscores);
		state_annot[s][n_feature][sp->first] = Score2Nats (sc);
		CL << " '" << sp->first << "' : " << Score2Prob(sc) << " (" << ecfg.state_info[s].name << ")   ";
	      }
	    }
	}
      CL << '\n';
      Stockholm::Annotation::const_iterator row_annot = stock.gc_annot.find (*f);
      if (row_annot != stock.gc_annot.end())
	{
	  align_annot[n_feature] = (row_annot->second.c_str());
	  CL << "Alignment annotation:\n" << align_annot[n_feature] << '\n';
	}
      ++n_feature;
    }
}

void ECFG_EM_matrix::show_coords (ostream& out, int subseq_idx, int state_idx) const
{
  const Subseq_coords& subseq = env.subseq[subseq_idx];
  out << "Subsequence #" << subseq_idx
      << " (coords " << subseq.start << "+" << subseq.len
      << "), state #" << state_idx << " (" << ecfg.desc(state_idx) << ")\n";
  ecfg.state_info[state_idx].show (out);
}

void ECFG_EM_matrix::get_partials (int state, Partial_map& with_context, Partial_map& with_wildcards)
{
  const ECFG_state_info& info = ecfg.state_info[state];
  const ECFG_chain& chain = ecfg.matrix_set.chain[info.matrix];
  const int chain_states = ecfg.matrix_set.total_states (info.matrix);
  const int env_subseqs = env.subseqs();

  // get the Column_matrix, and create a dummy Fast_prune object
  Column_matrix& cm = colmat[state];

  Fast_prune fp;
  fp.prepare (tree, lineage_matrix[info.matrix], cm);

  // init the Partial_map's
  with_context.clear();
  with_wildcards.clear();

  const vector<double> empty_partial_vec (env_subseqs * chain_states, (double) 0.);
  vector<Phylogeny::Node> nodes_with_sequence;
  for_rooted_nodes_pre (tree, b)
    {
      const Phylogeny::Node n = (*b).second;
      if (tree.node2row[n] >= 0)
	{
	  nodes_with_sequence.push_back(n);
	  with_context[n] = empty_partial_vec;
	  if (info.has_context())
	    with_wildcards[n] = empty_partial_vec;
	}
    }

  // make a vector of iterators for quick access
  typedef vector<double>::iterator vec_dbl_iter;
  vector<vec_dbl_iter>
    with_context_vec (tree.nodes()),
    with_wildcards_vec (tree.nodes());
  for_const_contents (vector<Phylogeny::Node>, nodes_with_sequence, n)
    {
      with_context_vec[*n] = with_context[*n].begin();
      if (info.has_context())
	with_wildcards_vec[*n] = with_wildcards[*n].begin();
    }

  // loop over subseqs
  for (int subseq_idx = 0; subseq_idx < env.subseqs(); ++subseq_idx)
    {
      const Subseq_coords& subseq = env.subseq[subseq_idx];
      // call initialise() and copy the results into the partials
      // wildcards
      if (info.has_context())
	{
	  info.initialise (ecfg.alphabet, chain.classes, asp, subseq, stock, tree, cm, &fp, true);  // context columns replaced by wildcard symbols

	  for_const_contents (vector<Phylogeny::Node>, nodes_with_sequence, n)
	    {
	      vec_dbl_iter& n_iter = with_wildcards_vec[*n];
	      vec_dbl_iter F_iter = fp.F[*n].begin();
	      for (int s = 0; s < chain_states; ++s)
		*(n_iter++) = *(F_iter++);
	    }
	}

      // context
      info.initialise (ecfg.alphabet, chain.classes, asp, subseq, stock, tree, cm, &fp, false);  // actual context columns used

      for_const_contents (vector<Phylogeny::Node>, nodes_with_sequence, n)
	{
	  vec_dbl_iter& n_iter = with_context_vec[*n];
	  vec_dbl_iter F_iter = fp.F[*n].begin();
	  for (int s = 0; s < chain_states; ++s)
	    *(n_iter++) = *(F_iter++);
	}
    }
}

void ECFG_EM_matrix::get_branch_transition_matrices_and_root_prior (int state, Transition_matrix_map& branch_transmat, vector<Prob>& root_prior)
{
  const ECFG_state_info& info = ecfg.state_info[state];
  const vector<const EM_matrix_base*>& lineage_mx = lineage_matrix[info.matrix];
  const int chain_states = ecfg.matrix_set.total_states (info.matrix);

  branch_transmat.clear();
  for_rooted_branches_pre (tree, b)
    {
      const array2d<Prob> bt = ((EM_matrix_base*) lineage_mx[(*b).second])->create_conditional_substitution_matrix ((*b).length);  // cast away const
      vector<double> bt_vec (chain_states * chain_states);
      for (int i = 0; i < chain_states; ++i)
	for (int j = 0; j < chain_states; ++j)
	  bt_vec[i*chain_states + j] = bt(i,j);
      branch_transmat[(*b).second].swap (bt_vec);
    }
  root_prior = ((EM_matrix_base*) lineage_mx[tree.root])->create_prior();  // cast away const
}

void ECFG_EM_matrix::use_precomputed_phyloemit (Emit_loglike_matrix& phyloemit)
{
  // calls use_precomputed() then iterates over all cells calling calc_annot_emit_ll() & adding to emit_loglike
  use_precomputed (phyloemit);
  Loge annot_emit_ll;
  for (int subseq_idx = 0; subseq_idx < env.subseqs(); ++subseq_idx)
    for_const_contents (vector<int>, inside_fill_states, s)
      {
	Loge& ell = emit_loglike (subseq_idx, *s);
	if (calc_annot_emit_ll (subseq_idx, *s, annot_emit_ll))
	  NatsPMulAcc (ell, annot_emit_ll);
	else
	  ell = -InfinityLoge;
      }
}

void ECFG_EM_matrix::compute_phylo_likelihoods_with_beagle()
{
#if defined(BEAGLE_INCLUDED) && BEAGLE_INCLUDED

  // set up the outputs & inputs
  const int subseqs = env.subseqs();
  Emit_loglike_matrix emit_loglike_matrix (subseqs, ecfg.states());

  Partial_map with_context, with_wildcards;
  Transition_matrix_map branch_transmat;
  vector<Prob> root_prior;

  // dummy variables for the Beagle-esque aspects of the model & algorithm (pattern weights, categories, etc.)
  const int cumulativeScalingIndex = BEAGLE_OP_NONE;
  const int categoryWeightsIndex = 0;
  const int stateFrequencyIndex = 0;

  const vector<double> patternWeights (subseqs, 1.);

  double dummy_rate[1], dummy_weight[1];
  dummy_rate[0] = dummy_weight[0] = 1.;

  double dummy_loglike_sum;

  // Beagle likes the leaf node indices to be lower than the internal node indices, so reorder the tree...
  vector<int> tree2beagle (tree.nodes());
  int beagle_node = 0;
  for (Phylogeny::Node_const_iter n = tree.leaves_begin(); n != tree.leaves_end(); ++n)
    tree2beagle[*n] = beagle_node++;
  for (Phylogeny::Node_const_iter n = tree.internals_begin(); n != tree.internals_end(); ++n)
    tree2beagle[*n] = beagle_node++;

  const int rootIndex = tree2beagle[tree.root];

  // get Beagle resource list
  BeagleResourceList* rList;
  rList = beagleGetResourceList();

  // loop over states
  for (int state = 0; state < ecfg.states(); ++state)
    if (ecfg.state_info[state].total_size() > 0)
      {
	// get partials (i.e. sequences), transition matrices & prior
	CTAG(5,BEAGLE) << "Preparing partials for Beagle\n";
	get_partials (state, with_context, with_wildcards);

	CTAG(5,BEAGLE) << "Preparing transition matrices for Beagle\n";
	get_branch_transition_matrices_and_root_prior (state, branch_transmat, root_prior);

	const int ctmc_states = (int) root_prior.size();

	// create Beagle instance
	CTAG(5,BEAGLE) << "Requesting Beagle instance for nonterminal " << ecfg.state_info[state].name << " (" << ctmc_states << " states)\n";
	BeagleInstanceDetails instDetails;
	int instance = beagleCreateInstance(tree.leaves(),           /**< Number of tip data elements (input) */
					    tree.nodes(),	           /**< Number of partials buffers to create (input) */
					    0,		           /**< Number of compact state representation buffers to create (input) */
					    ctmc_states,             /**< Number of states in the continuous-time Markov chain (input) */
					    subseqs,                 /**< Number of site patterns to be handled by the instance (input) */
					    1,		           /**< Number of rate matrix eigen-decomposition buffers to allocate (input) -- this must be at least 1, even though we don't use eigen-decompositions, because Beagle needs somewhere to put the state frequencies */
					    tree.nodes(),	           /**< Number of rate matrix buffers (input) */
					    1,                       /**< Number of rate categories (input) */
					    0,                       /**< Number of scaling buffers */
					    NULL,			   /**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
					    0,			   /**< Length of resourceList list (input) */
					    BEAGLE_FLAG_PROCESSOR_GPU,            	   /**< Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input) */
					    0,                       /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
					    &instDetails);
	if (instance < 0)
	  THROWEXPR ("Failed to obtain BEAGLE instance");

	// describe the instance
	int rNumber = instDetails.resourceNumber;
	if (CTAGGING(5,BEAGLE))
	  {
	    sstring beagle_log;
	    beagle_log << "Obtained Beagle instance for resource #" << rNumber;
	    beagle_log << "\n\tName : " << rList->list[rNumber].name;
	    beagle_log << "\n\tDesc : " << rList->list[rNumber].description;
	    beagle_log << "\n\tImpl : " << "GET INFO";
	    beagle_log << "\n\tFlags:";

	    if (instDetails.flags & BEAGLE_FLAG_PROCESSOR_CPU) beagle_log << " PROCESSOR_CPU";
	    if (instDetails.flags & BEAGLE_FLAG_PROCESSOR_GPU) beagle_log << " PROCESSOR_GPU";
	    if (instDetails.flags & BEAGLE_FLAG_PROCESSOR_FPGA) beagle_log << " PROCESSOR_FPGA";
	    if (instDetails.flags & BEAGLE_FLAG_PROCESSOR_CELL) beagle_log << " PROCESSOR_CELL";
	    if (instDetails.flags & BEAGLE_FLAG_PRECISION_DOUBLE) beagle_log << " PRECISION_DOUBLE";
	    if (instDetails.flags & BEAGLE_FLAG_PRECISION_SINGLE) beagle_log << " PRECISION_SINGLE";
	    if (instDetails.flags & BEAGLE_FLAG_COMPUTATION_ASYNCH) beagle_log << " COMPUTATION_ASYNCH";
	    if (instDetails.flags & BEAGLE_FLAG_COMPUTATION_SYNCH)  beagle_log << " COMPUTATION_SYNCH";
	    if (instDetails.flags & BEAGLE_FLAG_EIGEN_REAL)beagle_log << " EIGEN_RAW";
	    if (instDetails.flags & BEAGLE_FLAG_EIGEN_COMPLEX)beagle_log << " EIGEN_COMPLEX";
	    if (instDetails.flags & BEAGLE_FLAG_SCALING_MANUAL)beagle_log << " SCALING_MANUAL";
	    if (instDetails.flags & BEAGLE_FLAG_SCALING_AUTO)beagle_log << " SCALING_AUTO";
	    if (instDetails.flags & BEAGLE_FLAG_SCALING_ALWAYS)beagle_log << " SCALING_ALWAYS";
	    if (instDetails.flags & BEAGLE_FLAG_SCALERS_RAW)beagle_log << " SCALERS_RAW";
	    if (instDetails.flags & BEAGLE_FLAG_SCALERS_LOG)beagle_log << " SCALERS_LOG";
	    if (instDetails.flags & BEAGLE_FLAG_VECTOR_NONE)    beagle_log << " VECTOR_NONE";
	    if (instDetails.flags & BEAGLE_FLAG_VECTOR_SSE)    beagle_log << " VECTOR_SSE";
	    if (instDetails.flags & BEAGLE_FLAG_THREADING_NONE)    beagle_log << " THREADING_NONE";
	    if (instDetails.flags & BEAGLE_FLAG_THREADING_OPENMP)    beagle_log << " THREADING_OPENMP";

	    CL << beagle_log << "\n";
	  }

	// set dummy rate categories
	beagleSetCategoryRates(instance, &dummy_rate[0]);
	beagleSetCategoryWeights(instance, 0, &dummy_weight[0]);

	// set dummy pattern weights
	beagleSetPatternWeights(instance, &patternWeights[0]);

	// allocate space for transitionMatrices
	double* transitionMatrix = (double*) malloc (ctmc_states * ctmc_states * sizeof(double));

	// set transitionMatrices one at a time
	CTAG(5,BEAGLE) << "Sending transition matrices to Beagle instance\n";
	for_rooted_branches_pre (tree, b)
	  {
	    const Phylogeny::Node child = (*b).second;
	    const Transition_matrix_map::const_iterator tmx_iter = branch_transmat.find(child);
	    if (tmx_iter == branch_transmat.end())
	      THROWEXPR("Can't find node " << child);
	    const vector<double>& tmx = tmx_iter->second;
	    beagleSetTransitionMatrix(instance,
				      tree2beagle[child],
				      &tmx[0]);
	  }

	// free space for transitionMatrices
	free (transitionMatrix);

	// set prior for root
	CTAG(5,BEAGLE) << "Sending state frequencies to Beagle instance\n";
	beagleSetStateFrequencies(instance, stateFrequencyIndex, &root_prior[0]);

	// set the sequences for each tip using the partial likelihood arrays with context
	CTAG(5,BEAGLE) << "Sending partials to Beagle instance\n";
	for (Phylogeny::Node_const_iter n = tree.leaves_begin(); n != tree.leaves_end(); ++n)
	  beagleSetTipPartials(instance, tree2beagle[*n], &with_context[*n][0]);

	// create a list of partial likelihood update operations
	// the order is [dest, destScaling, source1, matrix1, source2, matrix2]
	vector<int> ops (BEAGLE_OP_COUNT * tree.internals());
	int k = 0;
	for_rooted_nodes_post (tree, b)
	  {
	    Phylogeny::Node n = (*b).second;
	    if (tree.is_internal(n))
	      {
		Phylogeny::Node_vector kids = tree.children (n, tree.parent[n]);
		if (kids.size() != 2)
		  THROWEXPR ("Beagle can only handle binary trees");
		ops[k] = tree2beagle[n];  // dest
		ops[k+1] = BEAGLE_OP_NONE;  // destinationScaleWrite
		ops[k+2] = BEAGLE_OP_NONE;  // destinationScaleRead
		ops[k+3] = tree2beagle[kids[0]];  // source1
		ops[k+4] = tree2beagle[kids[0]];  // matrix1
		ops[k+5] = tree2beagle[kids[1]];  // source2
		ops[k+6] = tree2beagle[kids[1]];  // matrix2
		k += BEAGLE_OP_COUNT;
	      }
	  }

	// update the partials
	CTAG(5,BEAGLE) << "Scheduling pruning at Beagle instance\n";
	beagleUpdatePartials(instance,          // instance
			     &ops[0],           // operations
			     tree.internals(),  // operationCount
			     BEAGLE_OP_NONE);   // cumulative scaling index

	// calculate the site likelihoods at the root node
	CTAG(5,BEAGLE) << "Scheduling root node calculations at Beagle instance\n";
	beagleCalculateRootLogLikelihoods(instance,                // instance
					  &rootIndex,              // bufferIndices
					  &categoryWeightsIndex,   // weights
					  &stateFrequencyIndex,    // stateFrequencies
					  &cumulativeScalingIndex, // cumulative scaling index
					  1,                       // count
					  &dummy_loglike_sum);     // sum over subseqs

	CTAG(5,BEAGLE) << "Requesting likelihoods from Beagle instance\n";
	vector<double> subseqLogLike (subseqs);
	beagleGetSiteLogLikelihoods(instance,
				    &subseqLogLike[0]);

	// copy subseq log-likelihoods into emit_loglike_matrix
	for (int n = 0; n < subseqs; ++n)
	  emit_loglike_matrix (n, state) = subseqLogLike[n];

	// repeat for with_wildcards; subtract
	if (with_wildcards.size())
	  {
	    // set the sequences for each tip using the partial likelihood arrays with wildcards
	    CTAG(5,BEAGLE) << "Sending partials to Beagle instance (with wildcards)\n";
	    for (Phylogeny::Node_const_iter n = tree.leaves_begin(); n != tree.leaves_end(); ++n)
	      beagleSetTipPartials(instance, tree2beagle[*n], &with_wildcards[*n][0]);

	    // update the partials
	    CTAG(5,BEAGLE) << "Scheduling pruning at Beagle instance (with wildcards)\n";
	    beagleUpdatePartials(instance,          // instance
				 &ops[0],           // operations
				 tree.internals(),  // operationCount
				 BEAGLE_OP_NONE);   // cumulative scaling index

	    // calculate the site likelihoods at the root node
	    CTAG(5,BEAGLE) << "Scheduling root node calculations at Beagle instance (with wildcards)\n";
	    beagleCalculateRootLogLikelihoods(instance,                // instance
					      &rootIndex,              // bufferIndices
					      &categoryWeightsIndex,   // weights
					      &stateFrequencyIndex,    // stateFrequencies
					      &cumulativeScalingIndex, // cumulative scaling index
					      1,                       // count
					      &dummy_loglike_sum);     // sum over subseqs

	    CTAG(5,BEAGLE) << "Requesting likelihoods from Beagle instance (with wildcards)\n";
	    beagleGetSiteLogLikelihoods(instance,
					&subseqLogLike[0]);

	    // subtract subseq log-likelihoods from emit_loglike_matrix
	    for (int n = 0; n < subseqs; ++n)
	      {
		Loge& ell = emit_loglike_matrix (n, state);
		NatsPMulAcc (ell, -subseqLogLike[n]);
	      }
	  }

	// let Beagle clean up
	CTAG(5,BEAGLE) << "Finalizing Beagle instance\n";
	beagleFinalizeInstance(instance);

	CTAG(5,BEAGLE) << "Finalized Beagle instance\n";
      }

  use_precomputed_phyloemit (emit_loglike_matrix);

#else /* BEAGLE_INCLUDED */
  CLOGERR << "An attempt was made to use the Beagle library to compute phylogenetic likelihoods.\nHowever, this program was not compiled with Beagle support.\nTry re-running 'configure' and ensuring that the Beagle library is detected.\n";
#endif /* BEAGLE_INCLUDED */
}

void ECFG_EM_matrix::reconstruct_MAP (Stockholm& stock, const ECFG_cell_score_map& annot, const char* ancrec_tag_cstr, bool annotate_MAP, bool annotate_postprobs, Prob min_reported_postprob)
{
  // #=GR tag
  const sstring ancrec_tag (ancrec_tag_cstr);
  sstring ancrec_tag_pp;
  ancrec_tag_pp << ancrec_tag << '_' << Stockholm_posterior_probability_tag;

  // prepare a list of ancestral nodes to reconstruct, and associated info
  vector<int> nodes_to_build;
  vector<sstring> names_of_nodes_to_build;
  map<int,Alignment_path::Row> row_path;
  map<int,sstring> seq, ppchars;

  // iterate over phylogeny
  for_rooted_nodes_post (tree, b)
    {
      const int node = (*b).second;

      // get node name
      const sstring node_name = tree.node_specifier (node);

      // is row in alignment?
      const bool row_in_alignment = stock.row_index.find (node_name) != stock.row_index.end();

      // get row_path
      if (row_in_alignment)
	{
	  const int row_index = tree.node2row[node];
	  row_path[node] = stock.path.row (row_index);
	}
      else
	row_path[node] = vector<bool> (stock.columns(), true);

      // create ancestral reconstruction string
      seq[node] = sstring (stock.columns(), Alignment::gap_char());
      ppchars[node] = sstring (stock.columns(), Alignment::gap_char());

      // flag this node for building
      nodes_to_build.push_back (node);
      names_of_nodes_to_build.push_back (node_name);
    }

  // print log message
  CTAG(7,ANCREC) << "Reconstructing the following nodes: (" << names_of_nodes_to_build << ")\n";

  // container for postprob annotations
  typedef map<Phylogeny::Node,list<sstring> > Postprob_annot_map;
  Postprob_annot_map annot_map;

  // loop over nonterminals in parse tree
  for_const_contents (ECFG_cell_score_map, annot, ecsm_ptr)
    {
      const Subseq_coords& coords = ecsm_ptr->first.first;
      const int ecfg_state = ecsm_ptr->first.second;

      // get state info
      const ECFG_state_info& info = ecfg.state_info[ecfg_state];
      if (info.emit_size())
	{
	  CTAG(6,ANCREC) << "Reconstructing nonterminal " << info.name << " at subsequence (" << coords.start << ".." << coords.end() << ")\n";

	  // get chain
	  const ECFG_chain& chain = ecfg.matrix_set.chain[info.matrix];

	  // call fill_down
	  fill_down (env.find_subseq_idx (coords.start, coords.len), ecfg_state);

	  // loop over nodes we want to reconstruct
	  vector<int> sym_idx;
	  for_const_contents (vector<int>, nodes_to_build, n)
	    {
	      // find max a posteriori chain state at this node (and record postprobs of all states)
	      typedef map<sstring,Prob> Tuple_prob_map;
	      Tuple_prob_map tuple_postprob;

	      int ml_chain_state = -1;
	      Prob ml_chain_state_prob = 0.;
	      const EM_matrix_base& matrix (chain.matrix ? *chain.matrix : *lineage_matrix[info.matrix][*n]);
	      for (int chain_state = 0; chain_state < matrix.m(); ++chain_state)
		{
		  const Prob chain_state_prob = colmat[ecfg_state].node_post_prob (*n, chain_state, tree, matrix);
		  if (ml_chain_state < 0 || chain_state_prob > ml_chain_state_prob)
		    {
		      ml_chain_state = chain_state;
		      ml_chain_state_prob = chain_state_prob;
		    }

		  // accumulate postprob of this emit tuple
		  if (chain_state_prob >= min_reported_postprob)
		    {
		      sstring tuple;
		      bool found_nongap = false;
		      chain.get_symbol_indices (chain_state, sym_idx);
		      for (int emit_pos = 0; emit_pos < info.l_emit + info.r_emit; ++emit_pos)
			{
			  const int col = info.column_index (coords, emit_pos);
			  if (col >= 0 && col < (int) row_path[*n].size() && row_path[*n][col])
			    {
			      const int emit_sym = sym_idx[emit_pos];
			      const char emit_char = ecfg.alphabet.int2char (emit_sym);
			      tuple << emit_char;
			      found_nongap = true;
			    }
			  else
			    tuple << Alignment::gap_char();
			}
		      if (chain.class_labels.size())
			tuple << chain.class_labels[sym_idx.back()];
		      if (found_nongap)
			tuple_postprob[tuple] += chain_state_prob;
		    }
		}

	      // populate ancestral sequence for this node
	      vector<int> rebuilt_cols;
	      vector<sstring> emit_state;
	      vector<char> emit_chars;
	      const int ppdec = (int) (ml_chain_state_prob * 10 + .5);
	      const char ppchar =
		ppdec < 0
		? (char) '0'
		: (ppdec > 9
		   ? (char) '+'
		   : (char) ('0' + ppdec));
	      for (int emit_pos = 0; emit_pos < info.l_emit + info.r_emit; ++emit_pos)
		{
		  const int col = info.column_index (coords, emit_pos);
		  const int pos = emit_pos + info.l_context;
		  if (col >= 0 && col < (int) row_path[*n].size() && row_path[*n][col])
		    {
		      const int emit_sym = (ml_chain_state / info.mul[emit_pos]) % ecfg.alphabet.size();
		      const char emit_char = ecfg.alphabet.int2char (emit_sym);
		      emit_chars.push_back (emit_char);
		      seq[*n][col] = emit_char;
		      ppchars[*n][col] = ppchar;
		    }
		  rebuilt_cols.push_back (col + 1);
		  emit_state.push_back (chain.state[pos]);
		}

	      // add posterior probability info to Stockholm
	      const sstring node_name = tree.node_specifier (*n);
	      if (annotate_postprobs && !tuple_postprob.empty())
		{
		  sstring gs_val;
		  gs_val << "Nonterm " << info.name << " columns (";
		  for (int i = 0; i < (int) rebuilt_cols.size(); ++i)
		    gs_val << (i > 0 ? "," : "") << rebuilt_cols[i];
		  gs_val << ")";
		  for_const_contents (Tuple_prob_map, tuple_postprob, tp)
		    gs_val << " P(" << tp->first << ")=" << min(1.,tp->second);
		  gs_val << '\n';
		  annot_map[*n].push_front (gs_val);
		}

	      // print log message
	      CTAG(3,ANCREC) << "Reconstructed sequence " << tree.node_specifier(*n) << " columns (" << rebuilt_cols << ") pseudoterminals (" << emit_state << ") chars (" << emit_chars << ") with post.prob. " << ml_chain_state_prob << '\n';
	    }
	}
    }

  // update Stockholm alignment
  for_const_contents (vector<int>, nodes_to_build, n)
    {
      // get node name
      const sstring node_name = tree.node_specifier (*n);

      // is row in alignment?
      const bool row_in_alignment = stock.row_index.find (node_name) != stock.row_index.end();

      // if no row in alignment, add it, so that the #=GR annot will show up
      if (!row_in_alignment)
	{
	  const int new_row_index = stock.add_row (row_path[*n], node_name);
	  stock.update_row_index();
	  tree.node2row[*n] = new_row_index;
	}

      // set #=GS annotations
      if (annotate_postprobs && annot_map.find(*n) != annot_map.end())
	{
	  sstring gs_val = stock.get_gs_annot (node_name, ancrec_tag_pp);
	  for_const_contents (list<sstring>, annot_map[*n], s)
	    gs_val << *s;
	  stock.set_gs_annot (node_name, ancrec_tag_pp, gs_val);
	}

      // set #=GR annotations
      if (annotate_MAP)
	{
	  stock.set_gr_annot (node_name, ancrec_tag, seq[*n]);
	  stock.set_gr_annot (node_name, ancrec_tag_pp, ppchars[*n]);
	}
    }
}

ECFG_inside_matrix::ECFG_inside_matrix (const ECFG_scores& ecfg, Stockholm& stock,
					const Aligned_score_profile& asp, const ECFG_envelope& env,
					bool use_fast_prune)
  : ECFG_EM_matrix (ecfg, stock, asp, env, use_fast_prune)
{ }

void ECFG_inside_matrix::fill()
{
  CTAG(4,INSIDE_MATRIX) << "ECFG_inside_matrix::fill: envelope has " << env.subseqs() << " subsequences.\n";
  for (int source_subseq_idx = 0; source_subseq_idx < env.subseqs(); ++source_subseq_idx)
    {
      env.get_bif_in (env.subseq[source_subseq_idx], bif_in);
      for_const_contents (vector<int>, inside_fill_states, s)
	fill_sum_state (source_subseq_idx, *s);
      show_progress ("Inside", source_subseq_idx, false);
    }

  // do start transitions
  const int final_subseq_idx = env.subseqs() - 1;
  for_const_contents (vector<int>, inside_fill_states, d)
    NatsPSumAcc (final_loglike, NatsPMul (start_ll[*d], cell (final_subseq_idx, *d)));

  // warn if score is negative & infinite
  if (final_loglike <= -InfinityLoge)
    CLOGERR << "Warning: inside algorithm yielded log-likelihood score of -infinity; check your grammar!\n";

  // print matrix to log for debugging
  if (CTAGGING(-1,INSIDE_MATRIX))
    {
      CL << "Inside matrix:\n";
      show(CL);
    }
}

Loge ECFG_inside_matrix::state_inside_ll (int state, const Subseq_coords& subseq) const
{
  const int idx = env.find_subseq_idx (subseq.start, subseq.len);
  if (idx < 0)
    return -InfinityLoge;
  return cell(idx,state);
}

ECFG_outside_matrix::ECFG_outside_matrix (ECFG_inside_matrix& inside, ECFG_counts* counts)
  : ECFG_matrix (inside.ecfg, inside.stock, inside.env),
    inside (inside),
    counts (counts),
    want_substitution_counts (true)
{
  if (counts)
    lineage_stats = inside.get_lineage_stats(*counts);
}

vector<vector<Update_statistics*> > ECFG_EM_matrix::get_lineage_stats (ECFG_counts& counts)
{
  vector<vector<Update_statistics*> > lineage_stats (ecfg.matrix_set.chain.size(), vector<Update_statistics*> (tree.nodes()));
  for (int c = 0; c < (int) ecfg.matrix_set.chain.size(); ++c)
    for (int n = 0; n < tree.nodes(); ++n)
      lineage_stats[c][n] = &counts.stats[lineage_chain_index[c][n]];
  return lineage_stats;
}

void ECFG_outside_matrix::fill()
{
  const Loge inside_final_ll = inside.final_loglike;
  vector<Subseq::Bifurc_out_l> bif_outl;
  vector<Subseq::Bifurc_out_r> bif_outr;
  const char* aa;
  sstring aa_str;
  bool aa_has_wildcards;
  for (int dest_subseq_idx = env.subseqs() - 1; dest_subseq_idx >= 0; --dest_subseq_idx)
    {
      const Subseq_coords& subseq = env.subseq[dest_subseq_idx];
      for_const_contents (vector<int>, outside_fill_states, d)
	{
	  Loge ll = -InfinityLoge;

	  // get info for dest state
	  const ECFG_state_info& info = ecfg.state_info[*d];

	  // find source subseq, check if it's in envelope
	  const int source_subseq_idx = env.find_subseq_idx (subseq.start - info.l_emit, subseq.len + info.emit_size());
	  if (source_subseq_idx >= 0)
	    {
	      // get source subseq
	      const Subseq_coords& src_subseq = env.subseq[source_subseq_idx];
	      const Subseq_coords& dest_emit_coords (src_subseq);  // create alias to help illuminate confusing coords conventions

	      // check source subseq is in minlen/maxlen range for this state; if not, move on to next state
	      if (info.out_of_range (dest_emit_coords))
		continue;

	      // check state info for number of emitted symbols; if source subseq is too short, move on to next state
	      if (info.emit_size() > dest_emit_coords.len)
		continue;

	      // get inside score
	      const Loge inside_minus_final_ll = NatsPMul (inside.cell (source_subseq_idx, *d), -inside_final_ll);

	      // do start transitions
	      if (source_subseq_idx == env.subseqs() - 1)
		{
		  ll = inside.start_ll[*d];
		  if (counts)
		    inside.add_trans_counts (Start, *d, Nats2Prob (NatsPMul (ll, inside_minus_final_ll)), *counts);
		}

	      // loop over incoming transitions
	      const vector<int>& incoming_d = incoming[*d];
	      const vector<Loge>& incoming_ll_d = incoming_ll[*d];
	      for (int incoming_idx = 0; incoming_idx < (int) incoming_d.size(); ++incoming_idx)
		{
		  const int s = incoming_d[incoming_idx];
		  // get emit coords for source state, & check these are in envelope
		  // (this is necessary because emission is outside subsequence for outside DP)
		  const ECFG_state_info& src_info = ecfg.state_info[s];
		  const int src_emit_coords_idx = env.find_subseq_idx (src_subseq.start - src_info.l_emit,
								       src_subseq.len + src_info.emit_size());
		  if (src_emit_coords_idx >= 0)
		    {
		      // calculate transition score; add counts
		      const Loge incoming_ll_s_d = NatsPMul (cell (source_subseq_idx, s), incoming_ll_d[incoming_idx]);
		      if (counts)
			inside.add_trans_counts (s, *d, Nats2Prob (NatsPMul (incoming_ll_s_d, inside_minus_final_ll)), *counts);
		      NatsPSumAcc (ll, incoming_ll_s_d);
		    }
		}

	      // get bifurcation subseq indices
	      // (doing this in the inner loop over outside_fill_states is criminally inefficient... should probably optimize... if gprof confirms this!)
	      env.get_bif_outl (dest_emit_coords, bif_outl);
	      env.get_bif_outr (dest_emit_coords, bif_outr);

	      // loop over incoming bifurcations to left
	      for_const_contents (vector<Subseq::Bifurc_out_l>, bif_outl, b)
		for_const_contents (vector<int>, left_bifurc[*d], s)
		NatsPSumAcc (ll, NatsPMul (cell (b->out, *s), inside.cell (b->l, ecfg.state_info[*s].ldest)));

	      // loop over incoming bifurcations to right
	      for_const_contents (vector<Subseq::Bifurc_out_r>, bif_outr, b)
		for_const_contents (vector<int>, right_bifurc[*d], s)
		NatsPSumAcc (ll, NatsPMul (cell (b->out, *s), inside.cell (b->r, ecfg.state_info[*s].rdest)));

	      // add emit score, *after* calculating post_prob
	      const Prob post_prob = Nats2Prob (NatsPMul (ll, inside_minus_final_ll));
	      NatsPMulAcc (ll, inside.emit_loglike (source_subseq_idx, *d));

	      // call fill_down to accumulate EM_matrix counts
	      if (counts && want_substitution_counts)
		inside.fill_down (source_subseq_idx, *d, counts, &lineage_stats, post_prob);

	      // accumulate annotation counts
	      if (counts)
		for (int a = 0; a < (int) inside.align_annot.size(); ++a)
		  if ((aa = inside.align_annot[a]) != 0)
		    {
		      const ECFG_EM_matrix::String_loglike_dist& sa = inside.state_annot[*d][a];
		      if (sa.size() > 0) {
			get_annot_string (info, aa, src_subseq, aa_str, aa_has_wildcards);
			if (aa_has_wildcards) {
			  Loge sa_ll = -InfinityLoge;
			  for_const_contents (ECFG_EM_matrix::String_loglike_dist, sa, sl) {
			    if (wildcard_match (aa_str, sl->first))
			      NatsPSumAcc (sa_ll, sl->second);
			  }
			  for_const_contents (ECFG_EM_matrix::String_loglike_dist, sa, sl) {
			    if (wildcard_match (aa_str, sl->first))
			      counts->state_annot_count[*d][a][sl->first] += post_prob * Nats2Prob (NatsPMul (sl->second, -sa_ll));
			  }
			} else
			  counts->state_annot_count[*d][a][aa_str] += post_prob;
		      }
		    }

	      // count end transitions
	      if (counts && subseq.len == 0)
		{
		  const Loge end_ll = inside.end_ll[*d];
		  inside.add_trans_counts (*d, End, Nats2Prob (NatsPMul3 (ll, end_ll, -inside_final_ll)), *counts);
		}
	    }
	  // store
	  cell (dest_subseq_idx, *d) = ll;
	}
      show_progress ("Outside", dest_subseq_idx, true);
    }

  // print matrix to log for debugging
  if (CTAGGING(-1,OUTSIDE_MATRIX))
    {
      CL << "Outside matrix:\n";
      show(CL);
    }
}

ECFG_inside_outside_matrix::ECFG_inside_outside_matrix (const ECFG_scores& ecfg, Stockholm& stock, const Aligned_score_profile& asp, const ECFG_envelope& env, ECFG_counts* counts)
  : inside (ecfg, stock, asp, env),
    outside (inside, counts)
{ }

void ECFG_inside_outside_matrix::fill()
{
  inside.fill();
  outside.fill();
}

Loge ECFG_inside_outside_matrix::post_transition_ll (int src_state, int dest_state, int subseq_idx) const
{
  const Subseq_coords& dest_subseq = inside.env.subseq[subseq_idx];
  Subseq_coords src_subseq = dest_subseq;

  Loge inout_loglike = -InfinityLoge;
  bool coords_ok = false;
  if (src_state == Start)
    {
      if (subseq_idx == inside.env.subseqs() - 1)
	{
	  coords_ok = true;
	  // get inside-outside log-likelihood
	  inout_loglike = inside.cell (subseq_idx, dest_state);
	}
    }
  else
    {
      // Find actual emit coords of src_state, and check it's in envelope
      const ECFG_state_info& src_state_info = inside.ecfg.state_info[src_state];
      src_subseq.start -= src_state_info.l_emit;
      src_subseq.len += src_state_info.emit_size();

      if (inside.env.find_subseq_idx (src_subseq.start, src_subseq.len) >= 0)
	{
	  coords_ok = true;
	  // get inside-outside log-likelihood
	  inout_loglike = outside.cell (subseq_idx, src_state);
	  if (dest_state != End)
	    NatsPMulAcc (inout_loglike, inside.cell (subseq_idx, dest_state));
	}
    }

  return coords_ok
    ? NatsPMul3 (inout_loglike,
		 inside.trans_ll (src_state, dest_state),
		 -inside.final_loglike)
    : -InfinityLoge;
}

Loge ECFG_inside_outside_matrix::post_state_ll (int dest_state, int subseq_idx) const
{
  Loge state_ll = -InfinityLoge;

  // sum over transitions
  for_const_contents (vector<int>, inside.incoming[dest_state], s)
    NatsPSumAcc (state_ll, post_transition_ll (*s, dest_state, subseq_idx));

  // debug
  if (Nats2Prob(state_ll) > 2)
    {
      for_const_contents (vector<int>, inside.incoming[dest_state], s)
	{
	  Prob tp = Nats2Prob(post_transition_ll (*s, dest_state, subseq_idx));
	  if (tp > 2)
	    {
	      CLOGERR<<"s="<<*s<<" dest="<<dest_state<<" subseq_idx="<<subseq_idx<<" post_trans="<<tp<<"\n";

	      Score out_sc = Nats2Score(outside.cell (subseq_idx, *s));
	      Score in_sc = Nats2Score(inside.cell (subseq_idx, dest_state));
	      Score trans_sc = Nats2Score(inside.trans_ll (*s, dest_state));
	      Score final_sc = Nats2Score(inside.final_loglike);
	      CLOGERR<<"out="<<out_sc<<" in="<<in_sc<<" trans="<<trans_sc<<" final="<<final_sc<<"\n";
	    }
	}
    }

  if (subseq_idx == inside.env.subseqs() - 1)
    NatsPSumAcc (state_ll, post_transition_ll (Start, dest_state, subseq_idx));

  // sum over bifurcations
  const Subseq_coords& subseq = inside.env.subseq[subseq_idx];
  Loge bif_ll = -InfinityLoge;

  vector<Subseq::Bifurc_out_l> bif_outl;
  inside.env.get_bif_outl (subseq, bif_outl);
  for_const_contents (vector<Subseq::Bifurc_out_l>, bif_outl, b)
    for_const_contents (vector<int>, inside.left_bifurc[dest_state], s)
    NatsPSumAcc (bif_ll, NatsPMul (outside.cell (b->out, *s), inside.cell (b->l, inside.ecfg.state_info[*s].ldest)));

  vector<Subseq::Bifurc_out_r> bif_outr;
  inside.env.get_bif_outr (subseq, bif_outr);
  for_const_contents (vector<Subseq::Bifurc_out_r>, bif_outr, b)
    for_const_contents (vector<int>, inside.right_bifurc[dest_state], s)
    NatsPSumAcc (bif_ll, NatsPMul (outside.cell (b->out, *s), inside.cell (b->r, inside.ecfg.state_info[*s].rdest)));

  NatsPSumAcc (state_ll, NatsPMul3 (bif_ll, inside.cell (subseq_idx, dest_state), -inside.final_loglike));

  // return
  return state_ll;
}

Loge ECFG_inside_outside_matrix::post_transition_ll (int src_state, int dest_state, const Subseq_coords& subseq) const
{
  return post_transition_ll (src_state, dest_state, inside.env.find_subseq_idx (subseq.start, subseq.len));
}

Loge ECFG_inside_outside_matrix::post_state_ll (int dest_state, const Subseq_coords& subseq) const
{
  return post_state_ll (dest_state, inside.env.find_subseq_idx (subseq.start, subseq.len));
}

void ECFG_inside_outside_matrix::annotate (Stockholm& stock, GFF_list& gff_list, const sstring& seqname, const ECFG_cell_score_map& annot, const sstring& logpostprob_tag) const
{
  sstring ppbycol (stock.columns(), ECFG_annotation_wildcard);
  for_const_contents (ECFG_cell_score_map, annot, ss)
    {
      const Subseq_coords& subseq = ss->first.first;
      const int state = ss->first.second;
      const ECFG_state_info& info = inside.ecfg.state_info[state];

      const Loge state_ll = post_state_ll (state, subseq);
      const Prob state_pp = Nats2Prob (state_ll);

      annotate_post_prob (gff_list, seqname, state, subseq);

      const char pp_char = '0' + minmax ((int) (10.*state_pp), 0, 9);
      if (info.emit_size())  // emit state: mark up postprob annotation columns with a postprob summary char (0-9)
	{
	  for (int pos = 0; pos < info.l_emit; ++pos)
	    ppbycol[subseq.start + pos] = pp_char;
	  for (int pos = 0; pos < info.r_emit; ++pos)
	    ppbycol[subseq.end() - info.r_emit + pos] = pp_char;
	}
    }
  stock.set_gc_annot (logpostprob_tag, ppbycol);
}

void ECFG_inside_outside_matrix::annotate_all_post_state_ll (GFF_list& gff_list, const sstring& seqname, const ECFG_cell_score_map& annot, const sstring& annot_tag, double min_postprob) const
{
  for (int subseq_idx = 0; subseq_idx < inside.env.subseqs(); ++subseq_idx)
    {
      const Subseq_coords& subseq = inside.env.subseq[subseq_idx];
      if (subseq.len > 0)
	for_const_contents (vector<int>, inside.inside_fill_states, s)
	  {
	    ECFG_cell_score_map::const_iterator annot_state_iter = annot.find (ECFG_subseq_state (subseq, *s));
	    const bool annot_found = annot_state_iter != annot.end();
	    const Prob pp = Nats2Prob (post_state_ll (*s, subseq_idx));
	    if (pp >= min_postprob || annot_found)
	      {
		map<sstring,sstring> extra;
		if (annot_found)
		  extra[annot_tag] = sstring("1");
		annotate_post_prob (gff_list, seqname, *s, subseq, &extra);
	      }
	  }
    }
}

void ECFG_inside_outside_matrix::annotate_post_prob (GFF_list& gff_list, const sstring& seqname, int state, const Subseq_coords& subseq, const map<sstring,sstring>* extra_annotations) const
{
  const int subseq_idx = inside.env.find_subseq_idx (subseq.start, subseq.len);
  const Log2 inout_lg = Nats2Bits (post_state_ll (state, subseq_idx));
  const Log2 inside_lg = Nats2Bits (inside.cell (subseq_idx, state));

  GFF gff;
  gff.seqname = seqname.size() ? seqname : inside.stock.get_name();
  gff.source = inside.ecfg.name;
  gff.feature = inside.ecfg.state_info[state].name;
  gff.start = subseq.start + 1;
  gff.end = subseq.end();
  gff.score = inout_lg;

  const sstring unique_id = gff_list.create_unique_id();
  gff.set_value (GFF_ID_tag, unique_id.c_str());

  sstring inout_tag, inout_score_str;
  inout_score_str << inout_lg;
  gff.set_value (ECFG_GFF_LogPostProb_tag, inout_score_str.c_str());

  sstring inside_tag, inside_score_str;
  inside_score_str << inside_lg;
  gff.set_value (ECFG_GFF_LogInsideProb_tag, inside_score_str.c_str());

  if (extra_annotations)
    {
      typedef map<sstring,sstring> AnnotMap;
      for_const_contents (AnnotMap, *extra_annotations, tv)
	gff.set_value (tv->first.c_str(), tv->second.c_str());
    }

  gff_list.push_back (gff);
}

void ECFG_inside_outside_matrix::annotate_hidden_classes (Stockholm& stock, const ECFG_cell_score_map& annot)
{
  set<sstring> feature_names;
  for_const_contents (vector<ECFG_chain>, inside.ecfg.matrix_set.chain, chain)
    if (chain->classes > 1)
      feature_names.insert (chain->class_row);

  const int align_cols = stock.columns();
  const sstring default_annot (align_cols, ECFG_annotation_wildcard);
  for_rooted_nodes_post (inside.tree, n)
    {
      const int align_row = inside.tree.node2row[(*n).second];
      if (align_row >= 0)
	{
	  const sstring& row_name = stock.row_name[align_row];
	  for_const_contents (set<sstring>, feature_names, feature_name)
	    stock.gr_annot[row_name][*feature_name] = default_annot;
	}
    }

  for_const_contents (ECFG_cell_score_map, annot, ss)
    {
      const Subseq_coords& subseq = ss->first.first;
      const int state = ss->first.second;
      const ECFG_state_info& info = inside.ecfg.state_info[state];

      if (info.emit_size())
	{
	  const ECFG_chain& chain = inside.ecfg.matrix_set.chain[info.matrix];
	  if (chain.class_labels.size())
	    {
	      const int subseq_idx = inside.env.find_subseq_idx (subseq.start, subseq.len);
	      inside.fill_up (subseq_idx, state);
	      inside.fill_down (subseq_idx, state);

	      Column_matrix& colmat (inside.colmat[state]);
	      const int observed_states = inside.ecfg.matrix_set.observed_states (info.matrix);
	      for_rooted_nodes_post (inside.tree, b)
		{
		  const Phylogeny::Node node = (*b).second;
		  const int align_row = inside.tree.node2row[node];
		  if (align_row >= 0)
		    {
		      // get column indices
		      vector<int> cols;
		      for (int pos = 0; pos < info.l_emit; ++pos)
			cols.push_back (subseq.start + pos);
		      for (int pos = 0; pos < info.r_emit; ++pos)
			cols.push_back (subseq.end() - info.r_emit + pos);
		      // check for gaps
		      int n_gaps = 0;
		      for_const_contents (vector<int>, cols, col)
			if (!stock.path (align_row, *col))
			  ++n_gaps;
		      if (!colmat.gapped[node] && n_gaps < info.emit_size())
			{
			  const sstring& row_name = stock.row_name[align_row];
			  sstring& gr = stock.gr_annot[row_name][chain.class_row];

			  const EM_matrix_base& matrix (chain.matrix ? *chain.matrix : *inside.lineage_matrix[info.matrix][node]);
			  vector<double> class_post_prob (chain.classes, 0.);
			  for (int c = 0; c < chain.classes; ++c)
			    for (int s = c * observed_states; s < (c+1) * observed_states; ++s)
			      class_post_prob[c] += colmat.node_post_prob (node, s, inside.tree, matrix);

			  const int ml_class = max_element (class_post_prob.begin(), class_post_prob.end()) - class_post_prob.begin();
			  const sstring& class_label = chain.class_labels[ml_class];
			  for (int pos = 0; pos < info.emit_size(); ++pos)
			    gr[cols[pos]] = class_label[pos];
			}
		    }
		}
	    }
	}
    }
}

ECFG_CYK_matrix::ECFG_CYK_matrix (const ECFG_scores& ecfg, Stockholm& stock,
				  const Aligned_score_profile& asp, const ECFG_envelope& env, bool try_fast_prune)
  : ECFG_EM_matrix (ecfg, stock, asp, env, try_fast_prune && alignment_safe_from_fast_prune_underflow (stock)),  // guard against underflow
    final_state (UndefinedState)
{ }

void ECFG_CYK_matrix::fill()
{
  // handle case of no subsequences
  if (env.subseqs() == 0)
	final_state = End;

  // main DP loop
  for (int source_subseq_idx = 0; source_subseq_idx < env.subseqs(); ++source_subseq_idx)
    {
      env.get_bif_in (env.subseq[source_subseq_idx], bif_in);
      for_const_contents (vector<int>, inside_fill_states, s)
	if (ecfg.state_info[*s].sum_state)
	  fill_sum_state (source_subseq_idx, *s);  // sum state: fill by summing probs (like Inside)
	else
	  fill_max_state (source_subseq_idx, *s);  // not a sum state: fill by taking max's
      show_progress ("CYK", source_subseq_idx, false);
    }

  // do start transitions
  const int final_subseq_idx = env.subseqs() - 1;
  for_const_contents (vector<int>, inside_fill_states, d)
  {
    const Loge f_ll =  NatsPMul (start_ll[*d], cell (final_subseq_idx, *d));
    if (f_ll > final_loglike)
      {
	final_loglike = f_ll;
	final_state = *d;
      }
  }

  // warn if score is negative & infinite
  if (final_loglike <= -InfinityLoge)
    CLOGERR << "Warning: CYK algorithm yielded log-likelihood score of -infinity; check your grammar!\n";

  // print matrix to log for debugging
  if (CTAGGING(-1,CYK_MATRIX))
    {
      CL << "CYK matrix:\n";
      show(CL);
    }
}

ECFG_cell_score_map ECFG_CYK_matrix::traceback()
{
  // sanity check on CYK matrix
  if (final_state == UndefinedState)
    THROWEXPR ("Traceback start co-ords undefined");

  if (final_loglike <= -InfinityLoge)
    THROWEXPR ("Final score is -infinity; traceback likely to break");

  // init traceback vars
  int subseq_idx = env.subseqs() - 1;
  int traceback_state = final_state;
  stack<int> subseq_idx_stack;
  stack<int> state_stack;
  ECFG_cell_score_map trace;

  // print log message
  CTAG(4,CYK_TRACEBACK) << "Starting CYK traceback at subseq #" << subseq_idx << ", state #" << traceback_state << "\n";

  // start traceback
  while (true)
  {
    // remember coords for later log message
    const int last_state = traceback_state;
    const int last_subseq_idx = subseq_idx;

    // terminate this branch of the parse tree if we've reached the End state
    if (traceback_state == End)
    {
      if (state_stack.empty()) break;

      subseq_idx = subseq_idx_stack.top();
      traceback_state = state_stack.top();

      subseq_idx_stack.pop();
      state_stack.pop();

      continue;
    }

    // get subseq & state info
    const Subseq_coords& subseq = env.subseq[subseq_idx];
    const ECFG_state_info& info = ecfg.state_info[traceback_state];

    // print log message
    if (CTAGGING(3,CYK_TRACEBACK))
      CL << "CYK traceback: at subsequence (" << subseq.start << ',' << subseq.end() << ") state " << info.name << '\n';

    // set up traceback score vars
    const Loge traceback_ll = cell (subseq_idx, traceback_state);
    Loge test_ll, best_ll = -InfinityLoge;

    // annotate
    trace[ECFG_subseq_state (subseq, traceback_state)] = traceback_ll;

    // check for bifurcation
    if (info.bifurc)
      {
	const int l = info.ldest;
	const int r = info.rdest;
	int best_l_idx = -1, best_r_idx = -1;

	env.get_bif_in (subseq, bif_in);
	for_const_contents (vector<Subseq::Bifurc_in>, bif_in, b)
	  if ((test_ll = NatsPMul (cell (b->l, l), cell (b->r, r))) > best_ll)
	    {
	      best_l_idx = b->l;
	      best_r_idx = b->r;
	      best_ll = test_ll;
	    }

	if (best_l_idx < 0)
	  THROWEXPR ("ECFG traceback failed: couldn't find bifurcation at subseq #" << subseq_idx << ", state #" << traceback_state);

	subseq_idx_stack.push (best_r_idx);
	subseq_idx = best_l_idx;
	state_stack.push (r);
	traceback_state = l;

      }
    else  // not a bifurcation
      {
	// get the skinny on this cell
	const int dest_subseq_idx = env.find_subseq_idx (subseq.start + info.l_emit, subseq.len - info.emit_size());
	if (dest_subseq_idx >= 0)
	  {
	    // get dest subseq
	    const Subseq_coords& dest = env.subseq[dest_subseq_idx];

	    // retrieve emit score calculated in fill()
	    const Loge emit_ll = emit_loglike (subseq_idx, traceback_state);

	    // do end transitions
	    int next_state = UndefinedState;
	    if (dest.len == 0 && (test_ll = NatsPMul (end_ll[traceback_state], emit_ll)) > best_ll)
	      {
		next_state = End;
		best_ll = test_ll;
	      }

	    // loop over outgoing transitions
	    for_const_contents (vector<int>, outgoing[traceback_state], d)
	      if ((test_ll = NatsPMul3 (cell (dest_subseq_idx, *d),
					 trans_ll (traceback_state, *d),
					 emit_ll)) > best_ll)
		  {
		    next_state = *d;
		    best_ll = test_ll;
		  }

	    if (next_state == UndefinedState)
	      THROWEXPR ("ECFG traceback failed: couldn't find dest state at subseq #" << subseq_idx << ", state #" << traceback_state);

	    traceback_state = next_state;
	    subseq_idx = dest_subseq_idx;
	  }
	else
	  THROWEXPR ("ECFG traceback failed: couldn't find dest subseq at subseq #" << subseq_idx << ", state #" << traceback_state);
      }		// end if (info.bifurc)

    if (!info.sum_state && abs (best_ll - traceback_ll) > ECFG_MIN_TRACEBACK_ERROR_TOLERANCE)
      CLOGERR << "Warning: ECFG traceback log-likelihood at subseq #" << last_subseq_idx << ", state #" << last_state
	      << " is " << best_ll << "; expected " << traceback_ll << "\n";
  }			// end while (1)

  return trace;
}

void ECFG_trainer::do_dp (bool do_outside)
{
  loglike = 0.;
  if (do_outside)
    counts.clear (pseud_init, pseud_mutate, pseud_wait);  // reset counts to pseudocounts
  for (int n = 0; n < (int) stock_db.size(); ++n)
    {
      Stockholm* stock = stock_db[n];

      CTAG(4,ECFG_EM) << "Starting dynamic programming for alignment #" << n+1
		      << " (" << stock->rows() << " rows, " << stock->columns() << " columns)\n";

      const Aligned_score_profile& asp = asp_vec[n];
      ECFG_auto_envelope env (*stock, ecfg, max_subseq_len);
      ECFG_inside_matrix inside (ecfg, *stock, asp, env, false);  // don't ever use fast pruning during training (not even post-iteration; likelihoods could turn out differently)
      inside.fill();	// calls fill_up() for each allowed subsequence, state pair

      CTAG(4,ECFG_EM) << "Alignment #" << n+1 << " inside log-likelihood is "
		      << Nats2Bits (inside.final_loglike) << " bits\n";

      if (inside.final_loglike <= -InfinityLoge)
	{
	  CLOGERR << "Skipping zero-probability alignment (this is not good. Check your grammar)\n";
	  continue;
	}
      NatsPMulAcc (loglike, inside.final_loglike);
      if (do_outside)
	{
	  ECFG_outside_matrix outside (inside, &counts);
	  outside.want_substitution_counts = ecfg.update_rates;
	  outside.fill();	// calls fill_up() and fill_down()
	}
    }

  // warn if not all matrices filled down (e.g., short alignment)
  if (do_outside && ecfg.update_rates)
    for (int i = 0; i < (int) counts.stats.size(); ++i)
      if (ecfg.matrix_set.chain[i].matrix && !counts.filled_down[i])
	CLOGERR << "Warning: insufficient data to optimize rate matrix #" << i+1 << "\n";
}

void ECFG_trainer::iterate_quick_EM (int forgive)
{
  // bail if nothing to update
  const bool ecfg_parametric = ecfg.has_parametric_functions();
  if (!ecfg.update_rates && !ecfg.update_rules && !ecfg_parametric)
    {
      CLOGERR << "Warning: tried to run EM without updating rates, rules _or_ params. This has no effect";
      return;
    }

  // if no parameter pseudocounts have been specified, borrow the global ones for chains.
  // this prevents some unpleasant side-effects of training without priors
  // (basically, rates can become huge, diagonalization can fail, log-likelihoods can even go positive...)
  // IH, 2/19/2007
  if (ecfg_parametric && ecfg.pcounts.groups() == 0)
    {
      CLOGERR << "Warning: no pseudocounts for grammar parameters. Using chain pseudocounts.\n";
      ecfg.pcounts = PCounts (ecfg.pscores);
      for (int g = 0; g < ecfg.pscores.groups(); ++g)
	{
	  const int s = ecfg.pscores.group_size(g);
	  if (s == 1)
	    {
	      ecfg.pcounts[PVar(g,0)] = pseud_mutate;
	      ecfg.pcounts.wait[g] = pseud_wait;
	    }
	  else
	    for (int v = 0; v < s; ++v)
	      ecfg.pcounts[PVar(g,v)] = pseud_init;
	}
    }

  // set up EM vars
  ECFG_matrix_set matrix_set = (ECFG_matrix_set &)ecfg.matrix_set;	// cast away const
  best_loglike = 0;
  int dec = 0;

  // the parameters for an ECFG_scores are currently split between several places...
  // ...this is clearly terrible design, and should be intelligently updated
  // (probably in a generic EM algorithm)
  Concrete_sparse_transition_scores
    best_transmat (ecfg.states(), -InfinityScore),
    old_transmat (ecfg.states(), -InfinityScore);  // the ECFG's transition matrix (null rule probabilities)
  vector <EM_matrix_params> best_params, old_params;  // the initial distributions & rate matrices for the Markov chains
  PScores best_pscores, old_pscores;  // the parameter values, for ECFGs in which the transitions/chains are parametric

  bool reached_max_iter = false;
  for (int iter = 0; ; ++iter)
  {
    // check for max iterations
    if (em_max_iter >= 0 && iter >= em_max_iter)
    {
      CTAG(7,ECFG_EM ECFG_EM_PROGRESS) << "EM hit " << em_max_iter << " iterations; stopping\n";
      reached_max_iter = true;
      break;
    }

    // save old scores
    old_transmat.assign_transition_matrix (ecfg);
    old_params.clear();
    for_contents (vector<ECFG_chain>, ecfg.matrix_set.chain, em_mat)
      if (em_mat->matrix)
	old_params.push_back ((EM_matrix_params&) *em_mat->matrix);
    old_pscores = ecfg.pscores;

    // get new counts
    do_dp (true);

    // write grammar with counts to training log
    if (training_log)
      {
	ECFG_scores ecfg_copy (ecfg);
	sstring iter_str, loglike_str, best_loglike_str;
	iter_str << EG_EM_UPDATE << " " << iter;
	loglike_str << EG_EM_BITS << " " << -Nats2Bits(loglike);
	best_loglike_str << EG_EM_BEST_BITS << " " << -Nats2Bits(loglike);
	ecfg_copy.meta.insert(ecfg_copy.meta.begin(),SExpr(best_loglike_str.begin(),best_loglike_str.end()));
	ecfg_copy.meta.insert(ecfg_copy.meta.begin(),SExpr(loglike_str.begin(),loglike_str.end()));
	ecfg_copy.meta.insert(ecfg_copy.meta.begin(),SExpr(iter_str.begin(),iter_str.end()));
	ECFG_builder::ecfg2stream (*training_log, ecfg_copy.alphabet, ecfg_copy, iter>0 ? &counts : 0);
	training_log->flush();
      }

    // update ECFG
    counts.update_ecfg (ecfg);

    if (CTAGGING(5,ECFG_EM ECFG_COUNTS))
      counts.show (ecfg, CL);

    // dump ECFG to logfile if appropriate
    if (CTAGGING(1,ECFG_EM ECFG_EM_SEXPR))
      {
	CL << "Grammar after EM iteration #" << iter+1 << ":\n";
	ECFG_builder::ecfg2stream (CL, ecfg.alphabet, ecfg);
      }

    // check for increase in log-likelihood
    const Loge prev_best = best_loglike;
    if (iter == 0 || loglike > best_loglike)
    {
      best_loglike = loglike;
      best_transmat = old_transmat;
      best_params = old_params;
      best_pscores = old_pscores;

      save_grammar();
    }

    // print log message
    CTAG(6,ECFG_EM ECFG_EM_PROGRESS) << "EM iteration #" << iter+1
				     << ": log-likelihood = " << Nats2Bits(loglike) << " bits\n";

    // check for decrease in log-likelihood
    if (iter > 0)
    {
      const double inc = (loglike - prev_best) / (abs(prev_best) < TINY ? 1. : abs(prev_best));  // IH, 4/20/2005
      CTAG(3,ECFG_EM ECFG_EM_PROGRESS) << "(previous best = " << Nats2Bits(prev_best)
				       << " bits; fractional improvement = " << inc << ")\n";
      if (inc < em_min_inc)
      {
        if (loglike < prev_best)
          CTAG(7,ECFG_EM ECFG_EM_PROGRESS) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_best)
					   << " to " << Nats2Bits(loglike) << " bits during EM\n";
        if (++dec > forgive)
        {
          CTAG(7,ECFG_EM ECFG_EM_PROGRESS) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
          break;
        }
      }
      else
        dec = 0;
    }	// end if (iter > 0)
  }	// end for (int iter = 0; ; ++iter)

  // since we only update the best log-likelihood the round *after* the new params are calculated,
  // we need to check for an increase one last time.
  // However, if we've maxed out our total number of iterations, then we DON'T do this last DP run,
  // because the user has waited long enough (this also keeps runtime down for the case when em_max_iter==0)
  if (!reached_max_iter)
    {
      CTAG(6,ECFG_EM ECFG_EM_PROGRESS) << "Checking post-iteration log-likelihood\n";
      do_dp (false);
      CTAG(6,ECFG_EM ECFG_EM_PROGRESS) << "Post-iteration log-likelihood = " << Nats2Bits(loglike) << " bits\n";
      // if log-likelihood hasn't improved (and we did at least one EM iteration), restore the best parameters
      if (em_max_iter != 0 && loglike < best_loglike)
	{
	  ecfg.assign_transition_matrix (best_transmat);
	  for (int i = 0, j = 0; i < (int) ecfg.matrix_set.chain.size(); ++i)
	    if (ecfg.matrix_set.chain[i].matrix)
	      {
		ecfg.matrix_set.chain[i].matrix->assign_matrix_params (best_params[j++]);
		ecfg.matrix_set.chain[i].matrix->update();
	      }
	  ecfg.pscores = best_pscores;

	  save_grammar();  // hook to save intermediate grammar to file
	}
      else  // log-likelihood did improve (or we didn't do any EM iterations, so have nothing to improve on)
	best_loglike = loglike;
    }

  // and return
  CTAG(6,ECFG_EM) << "Best log-likelihood: " << Nats2Bits(best_loglike) << "\n";
}
