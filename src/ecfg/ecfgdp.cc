#include "ecfg/ecfgdp.h"
#include "hmm/singletmpl.h"
#include "util/vector_output.h"
#include "ecfg/ecfgsexpr.h"

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
			  << "); using default label '" << default_gs_val << "'\n";

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

    use_fast_prune (use_fast_prune),
    fast_prune (ecfg.states()),
    fill_up_flag (true),

    src_gap_profile (tree.nodes()),
    dest_gap_profile (tree.nodes()),
    link_extend_sc (ecfg.states()),
    link_end_sc (ecfg.states()),
    ins_event_sc (ecfg.states(), tree.nodes()),
    del_event_sc (ecfg.states(), tree.nodes()),
    match_event_sc (ecfg.states(), tree.nodes())

{
  // debug: switch off fast_prune via logging
  if (CTAGGING(-99,NO_FAST_PRUNE))
    {
      this->use_fast_prune = use_fast_prune = false;
      CL << "Switched off fast_prune\n";
    }

  // create colmat vector of Column_matrix objects
  Column_matrix tmp_colmat;
  colmat = vector<Column_matrix> (ecfg.states(), tmp_colmat);

  // allocate emit_loglike
  if (CTAGGING(3,ALLOC))
    CL << "Allocating " << sizeof(Loge)*env.subseqs()*ecfg.states() << " bytes for emit scores ("
       << env.subseqs() << " subseqs * " << ecfg.states() << " states)\n";
  emit_loglike.resize (env.subseqs(), ecfg.states(), 0.);  // default emit_loglike is zero

  // allocate Column_matrix's
  if (CTAGGING(-1,ECFG_EM_MATRIX))
    asp.show (CL);
  for (int s = 0; s < ecfg.states(); ++s)
    if (ecfg.state_info[s].emit_size())
      {
	const int chain_idx = ecfg.state_info[s].matrix;
	colmat[s].alloc (tree.nodes(),
			 ecfg.matrix_set.total_states (chain_idx),
			 false); // don't allocate class labels
	if (use_fast_prune)
	  fast_prune[s].prepare (tree, lineage_matrix[chain_idx], colmat[s]);
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

  // initialise precomputed scores for indel events
  total_branch_len = 0.;
  for_rooted_branches_pre (tree, b)
    total_branch_len += (*b).length;
  for (int state = 0; state < ecfg.states(); ++state)
    {
      const ECFG_state_info& info = ecfg.state_info[state];

      link_extend_sc[state] = Prob2Score (info.link_extend);
      link_end_sc[state] = Prob2Score (info.link_end);

      const Loge no_ins_ll = -total_branch_len * info.ins_rate;
      const Score ins_sc = Prob2Score (1. - Nats2Prob (no_ins_ll));

      ins_event_sc (state, tree.root) = Nats2Score (no_ins_ll);

      for_rooted_branches_pre (tree, b)
	{
	  const int node = (*b).second;
	  const double length = (*b).length;

	  const Score node_sc = Prob2Score (length / total_branch_len);
	  const Loge no_del_ll = -length * info.del_rate;

	  ins_event_sc (state, node) = ScorePMul (ins_sc, node_sc);
	  del_event_sc (state, node) = Prob2Score (1. - Nats2Prob (no_del_ll));
	  match_event_sc (state, node) = Nats2Score (no_del_ll);
	}
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

// helpers for fill_up
// tests if a==b, allowing both a & b to contain wildcards
bool wildcard_match (const sstring& a, const sstring& b) {
  for (int pos = 0; pos < (int) a.size(); ++pos)
    if (a[pos] != b[pos] && a[pos] != ECFG_annotation_wildcard && b[pos] != ECFG_annotation_wildcard)
      return false;
  return true;
}

// get annotation string
void get_annot_string (const ECFG_state_info& info, const char* aa, const Subseq_coords& subseq, sstring& aa_str, bool& has_wildcards) {
  aa_str.clear();
  has_wildcards = false;

  for (int pos = 0; pos < info.l_emit; ++pos) {
    const char aa_char = aa[pos + subseq.start];
    aa_str << aa_char;
    if (aa_char == ECFG_annotation_wildcard)
      has_wildcards = true;
  }

  for (int pos = 0; pos < info.r_emit; ++pos) {
    const char aa_char = aa[subseq.end() - info.r_emit + pos];
    aa_str << aa_char;
    if (aa_char == ECFG_annotation_wildcard)
      has_wildcards = true;
  }
}

bool ECFG_EM_matrix::calc_annot_emit_ll (int subseq_idx, int state_idx, Loge& annot_emit_ll)
{
  annot_emit_ll = 0;

  const ECFG_state_info& info = ecfg.state_info[state_idx];
  const Subseq_coords& subseq = env.subseq[subseq_idx];

  // if the alignment annotation is incompatible with this emit state, assign probability zero and bail
  const char *aa;
  sstring aa_str;
  bool aa_has_wildcards;
  for (int a = 0; a < (int) align_annot.size(); ++a)
    if ((aa = align_annot[a]) != 0)
      {
	const String_loglike_dist& sa = state_annot[state_idx][a];
	if (sa.size() > 0) {
	  get_annot_string (info, aa, subseq, aa_str, aa_has_wildcards);

	  if (aa_has_wildcards) {
	    Loge sa_ll = -InfinityLoge;
	    for_const_contents (String_loglike_dist, sa, sl) {
	      if (wildcard_match (aa_str, sl->first))
		NatsPSumAcc (sa_ll, sl->second);
	    }
	    NatsPMulAcc (annot_emit_ll, sa_ll);

	  } else {
	    String_loglike_dist::const_iterator sa_iter = sa.find (aa_str);
	    if (sa_iter == sa.end())
	      {
		annot_emit_ll = -InfinityLoge;
		return false;
	      }
	    NatsPMulAcc (annot_emit_ll, sa_iter->second);
	  }
	}
      }

  return annot_emit_ll > -InfinityLoge;
}

// fill_up
bool ECFG_EM_matrix::fill_up (int subseq_idx, int state_idx, bool condition_on_context)
{
  // make refs to key variables
  Loge& emit_ll = emit_loglike (subseq_idx, state_idx);
  const ECFG_state_info& info = ecfg.state_info[state_idx];
  const Subseq_coords& subseq = env.subseq[subseq_idx];

  // if this is a bifurc or null state, assign probability one and bail
  if (info.bifurc || info.emit_size() == 0)  // bifurc or null?
    {
      emit_ll = 0.;  // assign probability one, so DP doesn't incur a penalty
      return false;  // return false, because this is not an emit state
    }

  // if the alignment annotation is incompatible with this emit state, assign probability zero and bail
  Loge annot_emit_ll;
  if (!calc_annot_emit_ll (subseq_idx, state_idx, annot_emit_ll))
    {
      emit_ll = -InfinityLoge;
      return false;
    }

  // print log message describing the subsequence co-ords & the state info
  if (CTAGGING(0,FILL_UP))
  {
    CL << "fill_up:\n";
    show_coords (CL, subseq_idx, state_idx);
  }

  // initialise the Column_matrix and call fill_up to do Felsenstein pruning
  Column_matrix& cm = colmat[state_idx];
  Fast_prune& fp = fast_prune[state_idx];

  // first, compute Felsenstein likelihood with emitted columns replaced with wildcards, in order to condition on the context
  Loge marginal_ll = 0.;
  const ECFG_chain& chain = ecfg.matrix_set.chain[info.matrix];
  const bool use_context = condition_on_context && (info.l_context > 0 || info.r_context > 0);
  if (use_context)
    {
      if (info.initialise (ecfg.alphabet, chain.classes, asp, subseq, stock, tree, cm, use_fast_prune ? &fp : (Fast_prune*) 0, true))
	{
	  if (use_fast_prune)
	    marginal_ll = Prob2Nats (fp.prune());
	  else
	    {
	      cm.fill_up (lineage_matrix[info.matrix], tree, subseq_idx);
	      marginal_ll = cm.total_log_likelihood();
	    }
	  if (CTAGGING(3,FILL_UP))
	    CL << "Marginal log likelihood (context): " << Nats2Bits(marginal_ll) << " bits\n";
	}
      else
	CTAG(3,FILL_UP) << "[context: initialise returned false]\n";
    }

  // now use the actual emitted columns
  if (info.initialise (ecfg.alphabet, chain.classes, asp, subseq, stock, tree, cm, use_fast_prune ? &fp : (Fast_prune*) 0, false))
    {
      Loge joint_ll = 0;
      if (use_fast_prune)
	joint_ll = Prob2Nats (fp.prune());
      else
	{
	  cm.fill_up (lineage_matrix[info.matrix], tree, subseq_idx);
	  joint_ll = cm.total_log_likelihood();
	}
      emit_ll = NatsPMul3 (annot_emit_ll, joint_ll, -marginal_ll);  // add in annot_emit_ll
      if (CTAGGING(3,FILL_UP))
	{
	  CL << "(Subseq " << subseq.start << "+" << subseq.len << ", state " << state_idx << ") ";
	  if (use_context)
	    {
	      CL << "Joint log likelihood (context,emit): " << Nats2Bits(joint_ll) << " bits\n";
	      CL << "Conditional log likelihood (emit|context): " << Nats2Bits(emit_ll) << " bits\n";
	    }
	  else
	    CL << "Emit log likelihood: " << Nats2Bits(emit_ll) << " bits\n";
	}
      return true;  // success
    }
  else
    CTAG(3,FILL_UP) << "[initialise returned false]\n";

  // initialise failed; assign probability zero
  emit_ll = -InfinityLoge;
  return false;
}

void ECFG_EM_matrix::fill_down (ECFG_counts& counts, int subseq_idx, int state_idx, double weight)
{
  // set up lineage_stats
  vector<vector<Update_statistics*> > lineage_stats (ecfg.matrix_set.chain.size(), vector<Update_statistics*> (tree.nodes()));
  for (int c = 0; c < (int) ecfg.matrix_set.chain.size(); ++c)
    for (int n = 0; n < tree.nodes(); ++n)
      lineage_stats[c][n] = &counts.stats[lineage_chain_index[c][n]];
  // check state info
  const ECFG_state_info& info = ecfg.state_info[state_idx];
  if (info.bifurc || info.total_size() == 0)
    return;
  // call fill_up
  if (fill_up (subseq_idx, state_idx))   // only proceed if fill_up succeeded
    {
      // print fill_down log message
      if (CTAGGING(0,FILL_DOWN))
	CL << "fill_down: weight=" << weight << "\n";
      // do pruning
      Column_matrix& cm = colmat[state_idx];
      // NB lineage-dependent models: the following two lines need to set up vector<EM_matrix_base*> and vector<Update_statistics*> and set *all* relevant filled_down[] flags
      cm.fill_down (lineage_matrix[info.matrix], tree, lineage_stats[info.matrix], subseq_idx, weight);
      for (int n = 0; n < tree.nodes(); ++n)
	counts.filled_down[lineage_chain_index[info.matrix][n]] = true;		// pk for EM - remember that this matrix has been filled down
    }
}

void ECFG_EM_matrix::fill_sum_state (int source_subseq_idx, int source_state_idx)
{
  const Subseq_coords& subseq = env.subseq[source_subseq_idx];
  const ECFG_state_info& info = ecfg.state_info[source_state_idx];
  Loge ll = -InfinityLoge;

  // if subseq is out of range for this state, probability stays at zero
  if (!info.out_of_range (subseq))
    {
      // bifurcation or "regular" state?
      if (info.bifurc)
	{
	  const int l = info.ldest;
	  const int r = info.rdest;

	  for_const_contents (vector<Subseq::Bifurc_in>, bif_in, b)
	    NatsPSumAcc (ll, NatsPMul (cell (b->l, l), cell (b->r, r)));
	}
      else
	{
	  const int dest_subseq_idx = env.find_subseq_idx (subseq.start + info.l_emit, subseq.len - info.emit_size());
	  if (dest_subseq_idx >= 0)
	    {
	      const Subseq_coords& dest = env.subseq[dest_subseq_idx];

	      if (dest.len == 0)  // do end transitions
		NatsPSumAcc (ll, effective_trans (source_state_idx, End, subseq, dest));
	      // loop over outgoing transitions
	      for_const_contents (vector<int>, outgoing[source_state_idx], d)
		NatsPSumAcc (ll, NatsPMul (cell (dest_subseq_idx, *d), effective_trans (source_state_idx, *d, subseq, dest)));

	      // add emit score (even if fill_up fails, in which case emit score might be -inf)
	      if (fill_up_flag)
		fill_up (source_subseq_idx, source_state_idx);
	      NatsPMulAcc (ll, emit_loglike (source_subseq_idx, source_state_idx));
	    }
	}
    }

  // store
  cell (source_subseq_idx, source_state_idx) = ll;
}

void ECFG_EM_matrix::fill_max_state (int source_subseq_idx, int source_state_idx)
{
  const Subseq_coords& subseq = env.subseq[source_subseq_idx];
  const ECFG_state_info& info = ecfg.state_info[source_state_idx];
  Loge ll = -InfinityLoge;

  // if subseq is out of range for this state, probability stays at zero
  if (!info.out_of_range (subseq))
    {
      // bifurcation or "regular" state?
      if (info.bifurc)
	{
	  const int l = info.ldest;
	  const int r = info.rdest;

	  for_const_contents (vector<Subseq::Bifurc_in>, bif_in, b)
	    ll = max (ll, NatsPMul (cell (b->l, l), cell (b->r, r)));
	}
      else
	{
	  const int dest_subseq_idx = env.find_subseq_idx (subseq.start + info.l_emit, subseq.len - info.emit_size());
	  if (dest_subseq_idx >= 0)
	    {
	      const Subseq_coords& dest = env.subseq[dest_subseq_idx];

	      if (dest.len == 0)  // do end transitions
		ll = max (ll, effective_trans (source_state_idx, End, subseq, dest));
	      // loop over outgoing transitions
	      for_const_contents (vector<int>, outgoing[source_state_idx], d)
		ll = max (ll, NatsPMul (cell (dest_subseq_idx, *d), effective_trans (source_state_idx, *d, subseq, dest)));

	      // add emit score (even if fill_up fails, in which case emit score might be -inf)
	      if (fill_up_flag)
		fill_up (source_subseq_idx, source_state_idx);
	      NatsPMulAcc (ll, emit_loglike (source_subseq_idx, source_state_idx));
	    }
	}
    }

  // store
  cell (source_subseq_idx, source_state_idx) = ll;
}

void ECFG_EM_matrix::reconstruct_MAP (Stockholm& stock, const ECFG_cell_score_map& annot, const char* ancrec_tag_cstr, bool annotate_postprobs)
{
  // create dummy counts
  ECFG_counts dummy_counts (ecfg);

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
	  fill_down (dummy_counts, env.find_subseq_idx (coords.start, coords.len), ecfg_state, 1.);

	  // loop over nodes we want to reconstruct
	  for_const_contents (vector<int>, nodes_to_build, n)
	    {
	      // find max a posteriori chain state at this node
	      int ml_chain_state = -1;
	      Prob ml_chain_state_prob = 0.;
	      const EM_matrix_base& matrix (chain.matrix ? *chain.matrix : *lineage_matrix[info.matrix][*n]);
	      for (int chain_state = 0; chain_state < chain.matrix->m(); ++chain_state)
		{
		  const Prob chain_state_prob = colmat[ecfg_state].node_post_prob (*n, chain_state, tree, matrix);
		  if (ml_chain_state < 0 || chain_state_prob > ml_chain_state_prob)
		    {
		      ml_chain_state = chain_state;
		      ml_chain_state_prob = chain_state_prob;
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
		  if (row_path[*n][col])
		    {
		      const int emit_sym = (ml_chain_state / info.mul[emit_pos]) % ecfg.alphabet.size();
		      const char emit_char = ecfg.alphabet.int2char (emit_sym);
		      seq[*n][col] = emit_char;
		      ppchars[*n][col] = ppchar;
		      emit_chars.push_back (emit_char);
		    }
		  rebuilt_cols.push_back (col + 1);
		  emit_state.push_back (chain.state[pos]);
		}

	      // add posterior probability info to Stockholm
	      if (annotate_postprobs)
		{
		  const sstring node_name = tree.node_specifier (*n);
		  sstring gs_val = stock.get_gs_annot (node_name, ancrec_tag_pp);
		  gs_val << "State " << info.name
			 << " columns (" << rebuilt_cols
			 << ") chars (" << emit_chars
			 << ") postprob " << ml_chain_state_prob
			 << '\n';
		  stock.set_gs_annot (node_name, ancrec_tag_pp, gs_val);
		}

	      // print log message
	      CTAG(5,ANCREC) << "Reconstructed sequence " << tree.node_specifier(*n) << " columns (" << rebuilt_cols << ") pseudoterminals (" << emit_state << ") chars (" << emit_chars << ") with post.prob. " << ml_chain_state_prob << '\n';
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

      // set #=GR annotations
      stock.set_gr_annot (node_name, ancrec_tag, seq[*n]);
      stock.set_gr_annot (node_name, ancrec_tag_pp, ppchars[*n]);
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
  const Subseq_coords& final_subseq = env.subseq[final_subseq_idx];
  for_const_contents (vector<int>, inside_fill_states, d)
    NatsPSumAcc (final_loglike, NatsPMul (effective_trans (Start, *d, final_subseq, final_subseq), cell (final_subseq_idx, *d)));

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
{ }

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
		  ll = inside.effective_trans (Start, *d, dest_emit_coords, dest_emit_coords);
		  if (counts)
		    inside.add_trans_counts (Start, *d, dest_emit_coords, dest_emit_coords,
					     Nats2Prob (NatsPMul (ll, inside_minus_final_ll)), *counts);
		}

	      // loop over incoming transitions
	      for_const_contents (vector<int>, incoming[*d], s)
		{
		  // get emit coords for source state, & check these are in envelope
		  // (this is necessary because emission is outside subsequence for outside DP)
		  const ECFG_state_info& src_info = ecfg.state_info[*s];
		  const int src_emit_coords_idx = env.find_subseq_idx (src_subseq.start - src_info.l_emit,
								       src_subseq.len + src_info.emit_size());
		  if (src_emit_coords_idx >= 0)
		    {
		      // calculate transition score; add counts
		      const Subseq_coords& src_emit_coords = env.subseq[src_emit_coords_idx];
		      const Loge incoming_ll = NatsPMul (cell (source_subseq_idx, *s),
							 inside.effective_trans (*s, *d, src_emit_coords, dest_emit_coords));
		      if (counts)
			inside.add_trans_counts (*s, *d, src_emit_coords, dest_emit_coords,
						 Nats2Prob (NatsPMul (incoming_ll, inside_minus_final_ll)), *counts);
		      NatsPSumAcc (ll, incoming_ll);
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
		inside.fill_down (*counts, source_subseq_idx, *d, post_prob);

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
		  const Loge end_ll = inside.effective_trans (*d, End, dest_emit_coords, subseq);
		  inside.add_trans_counts (*d, End, dest_emit_coords, subseq,
					   Nats2Prob (NatsPMul3 (ll, end_ll, -inside_final_ll)), *counts);
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
		 inside.effective_trans (src_state, dest_state, src_subseq, dest_subseq),
		 -inside.final_loglike)
    : -InfinityLoge;
}

Loge ECFG_inside_outside_matrix::post_state_ll (int dest_state, int subseq_idx) const
{
  Loge state_ll = -InfinityLoge;

  // sum over transitions
  for_const_contents (vector<int>, inside.incoming[dest_state], s)
    NatsPSumAcc (state_ll, post_transition_ll (*s, dest_state, subseq_idx));

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

void ECFG_inside_outside_matrix::annotate_all_post_state_ll (GFF_list& gff_list, const sstring& seqname, const ECFG_cell_score_map& annot, const sstring& annot_tag) const
{
  for (int subseq_idx = 0; subseq_idx < inside.env.subseqs(); ++subseq_idx)
    {
      const Subseq_coords& subseq = inside.env.subseq[subseq_idx];
      if (subseq.len > 0)
	for_const_contents (vector<int>, inside.inside_fill_states, s)
	  {
	    ECFG_cell_score_map::const_iterator annot_state_iter = annot.find (ECFG_subseq_state (subseq, *s));
	    const bool annot_found = annot_state_iter != annot.end();

	    annotate_post_prob (gff_list, seqname, *s, subseq);
	    if (annot_found)
	      gff_list.back().set_value (annot_tag.c_str(), "1");
	  }
    }
}

void ECFG_inside_outside_matrix::annotate_post_prob (GFF_list& gff_list, const sstring& seqname, int state, const Subseq_coords& subseq) const
{
  const int subseq_idx = inside.env.find_subseq_idx (subseq.start, subseq.len);

  GFF gff;
  gff.seqname = seqname.size() ? seqname : inside.stock.get_name();
  gff.source = inside.ecfg.name;
  gff.feature = inside.ecfg.state_info[state].name;
  gff.start = subseq.start + 1;
  gff.end = subseq.end();
  gff.score = Nats2Bits (inside.cell (subseq_idx, state));

  sstring score_str;
  score_str << Nats2Bits (post_state_ll (state, subseq_idx));
  gff.set_value (ECFG_GFF_LogPostProb_tag, score_str.c_str());

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

  ECFG_counts dummy_counts (inside.ecfg);
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
	      dummy_counts.stats[info.matrix].states = 0;  // signal that we don't want substitution counts
	      inside.fill_down (dummy_counts, subseq_idx, state, 1.);

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
  const Subseq_coords& final_subseq = env.subseq[final_subseq_idx];
  for_const_contents (vector<int>, inside_fill_states, d)
  {
    const Loge f_ll =  NatsPMul (effective_trans (Start, *d, final_subseq, final_subseq), cell (final_subseq_idx, *d));
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
	    if (dest.len == 0 && (test_ll = NatsPMul (effective_trans (traceback_state, End, subseq, dest), emit_ll)) > best_ll)
	      {
		next_state = End;
		best_ll = test_ll;
	      }

	    // loop over outgoing transitions
	    for_const_contents (vector<int>, outgoing[traceback_state], d)
	      if ((test_ll = NatsPMul3 (cell (dest_subseq_idx, *d),
					 effective_trans (traceback_state, *d, subseq, dest),
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
      ECFG_envelope env (stock->columns(), max_subseq_len);
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

  // the parameters for an ECFG_scores are currently split between four places...
  // ...this is clearly terrible design, and should be intelligently updated
  // (probably in a generic EM algorithm)
  Concrete_sparse_transition_scores
    best_transmat (ecfg.states(), -InfinityScore),
    old_transmat (ecfg.states(), -InfinityScore);  // the ECFG's transition matrix (null rule probabilities)
  vector <EM_matrix_params> best_params, old_params;  // the initial distributions & rate matrices for the Markov chains
  vector <ECFG_state_info> best_info, old_info;  // the gap models of the ECFG's emit states
  PScores best_pscores, old_pscores;  // the parameter values, for ECFGs in which the transitions/chains/indel models are parametric

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
    old_info = ecfg.state_info;
    old_pscores = ecfg.pscores;

    // get new counts
    do_dp (true);

    // write grammar with counts to training log
    if (training_log) {
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
      best_info = old_info;
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
	  ecfg.state_info = best_info;
	  ecfg.pscores = best_pscores;

	  save_grammar();  // hook to save intermediate grammar to file
	}
      else  // log-likelihood did improve (or we didn't do any EM iterations, so have nothing to improve on)
	best_loglike = loglike;
    }

  // and return
  CTAG(6,ECFG_EM) << "Best log-likelihood: " << Nats2Bits(best_loglike) << "\n";
}
