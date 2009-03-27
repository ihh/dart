#include "handel/movement.h"
#include "handel/hmmoc_adapter.h"

void Handel_movement::read_composition (const vector<sstring>& filenames)
{
  SExpr_file sexpr_file (filenames);
  composition = Transducer_SExpr_file (sexpr_file.sexpr);
}

void Handel_movement::dump_composition (ostream& out)
{
  // set up some flags
  fill_forward = nforward >= 0;
  fill_backward = want_expected_counts || optacc;
  redelings_suchard = propose_redelings_suchard_move || evaluate_redelings_suchard_inverse_move;
  const bool do_multidimensional_alignment = viterbi || fill_forward || fill_backward;
  do_alignment = do_multidimensional_alignment || redelings_suchard;

  // set up Redelings-Suchard proposal distribution, if none was specified
  if (redelings_suchard && composition.proposal_branches.size() == 0)
    composition.autopropose();

  // if old_path is specified for all unconstrained branches, then compute pre-resampling score
  bool got_old_loglike = false;
  vector<sstring> pathless_branch_names;
  if (composition.free_clique >= 0)
    {
      for (int n = 0; n < composition.etree.nodes(); ++n)
	if (composition.path.find (n) == composition.path.end()
	    && composition.old_path.find (n) == composition.old_path.end())
	  pathless_branch_names.push_back (composition.branch_name[n]);

      if (pathless_branch_names.size() == 0)
	{
	  Transducer_SExpr_file::NodePathMap path (composition.path);
	  path.insert (composition.old_path.begin(), composition.old_path.end());

	  old_loglike = composition.get_path_loglike (path);
	  got_old_loglike = true;
	}
    }

  // if simulation specified, then simulate paths & sequences
  if (simulate)
    composition.simulate();

  // if alignment file specified, dump alignment
  if (stockfile.size())
    {
      ofstream stockfile_stream (stockfile.c_str());
      Stockade stockade = composition.stockade();
      stockade.align.write_Stockholm (stockfile_stream);
    }

  // if fully constrained, get the score
  fully_constrained = composition.free_clique < 0;
  vector<Pair_transducer_scores> branch_sc;
  if (fully_constrained)
    {
      branch_sc = composition.make_branch_sc();
      composition.peel_clique (0, branch_sc, peeling_loglike, 0);
    }

  // peel back constrained branches
  peeling_required = composition.path.size() > 0 && !fully_constrained;
  peeled_composition = 
    peeling_required
    ? composition.peel_constrained (peeling_loglike, normalize_peeled_profiles)
    : composition;

  // if the old path was fully specified, AND we are using the HMMoC adapter (currently the only implementation that does banding),
  // then check to see whether the old path is within the banding constraint; and, if it's not, then flag this so that the move can be rejected for MCMC (for detailed balance reasons)
  if (got_old_loglike && (use_hmmoc_for_forward || use_hmmoc_for_viterbi))
    {
      Alignment_path::Decomposition old_decomp;
      for (int n = 1; n < peeled_composition.etree.nodes(); ++n)
	old_decomp[Alignment_path::Row_pair (peeled_composition.etree.parent[n], n)] = peeled_composition.branch_trans[n].state2align_path (peeled_composition.old_path[n], true);

      Alignment_path old_align_path;
      old_align_path.compose (old_decomp, false, true);

      long band_diameter = 0;
      if (peeled_composition.has_all_banding_coefficients())
	band_diameter = peeled_composition.get_band_diameter();

      if (!HMMoC_adapter_options::path_is_in_band (old_align_path, band_diameter))  // if band_diameter==0, will use upper-limit banding constraints
	old_path_is_outside_hmmoc_band = true;
    }

  // output the input (minimal filter behavior)
  // note that we do this *after* peeling constrained branches,
  // which may introduce a new singleton transducer
  if (!quiet)
    {
      composition.show_defs (out);
      composition.show_tree (out, 0, show_constrained_composite_paths);
      out.flush();
    }

  // show the log-likelihood (if the path is constrained on all branches)
  if (fully_constrained)
    {
      out << '\n';
      dump_raw_loglike (out, TSEXPR_FULL_SCORE, peeling_loglike);
      out.flush();
    }

  // show the old log-likelihood (if the old path is specified on all unconstrained branches)
  if (got_old_loglike)
    {
      out << '\n';
      dump_raw_loglike (out, TSEXPR_OLD_SCORE, old_loglike);
      out.flush();
    }

  // show the peeled composition
  const bool show_peeled_composition = (peeling_required || fully_constrained) && !quiet;
  if (show_peeled_composition)
    {
      // show transducer definitions if they were previously suppressed
      // commented out by IH, 3/3/2008, because it's confusing to print transducers when "--quiet" was explicitly requested to turn them off.
      // however, I leave the commented-out code here as an illustration that
      // peeled_composition is self-contained: its show_defs method prints
      // all the branch transducer definitions required to do an independent
      // phylocomposer run.

      // if (quiet)
      //  {
      //   peeled_composition.show_defs (out);
      //   out.flush();
      //  }

      // show peeled composition
      out << "\n(" << TSEXPR_PEELED;
      peeled_composition.show_tree (out, 1, false);
      out << ' ';
      dump_score (out, TSEXPR_PSCORE, (Score) 0);  // dump_score will add peeling_loglike to the supplied score of zero
      out << ")\n";
      out.flush();
    }

  // show reconstructed ancestral tapes
  if (composition.nodes_to_reconstruct.size())
    {
      if (!fully_constrained)
	THROWEXPR ("Can't reconstruct ancestral nodes unless state path is specified on every branch");
      for_const_contents (set<int>, composition.nodes_to_reconstruct, n)
	{
	  // since fully_constrained==true, branch_sc is already initialized
	  Loge dummy_loglike = 0.;
	  const Score_profile prof_sc = composition.peel_clique (0, branch_sc, dummy_loglike, *n);
	  out << '(' << TSEXPR_RECONS
	      << "\n (" << TSEXPR_NAME << ' ' << composition.tape_name[*n]
	      << ")\n ";
	  composition.show_prof (prof_sc, out, 1);
	  out << ")\n";
	}
    }

  // are we doing alignment? if so, we'll need sequences & emit labels
  if (do_alignment && peeled_composition.prof.empty())
    THROWEXPR ("Alignment was requested, but no sequence data supplied");
  if (do_alignment && fully_constrained)
    THROWEXPR ("Alignment was requested, but state paths were specified on all branches, so alignment is fully constrained");

  // see if we're using the composite transduer as an HMM later
  use_composite_as_HMM = acyclic || do_multidimensional_alignment || (hmmoc_opts.hmmoc_filename_prefix.size() > 0);

  // do we need to do anything else?
  const bool ehmm_funcs_defined = composite || use_composite_as_HMM || dotfile.size();
  if (!ehmm_funcs_defined && !redelings_suchard)
    return;

  // use HMMoC?
  use_hmmoc_for_forward = hmmoc_opts.try_to_use_hmmoc_adapter && fill_forward && !fill_backward;
  use_hmmoc_for_viterbi = hmmoc_opts.try_to_use_hmmoc_adapter && viterbi;

  // build the EHMM_transducer_funcs
  if (ehmm_funcs_defined)
    ehmm_funcs = peeled_composition.build_composite();

  // output the EHMM_transducer_funcs, if specified
  if (composite)
    {
      peeled_composition.show_composite_name_format (out);
      peeled_composition.show_composite (ehmm_funcs, out);
      out.flush();
    }

  // alignment requires emit labels...
  if (do_multidimensional_alignment && !ehmm_funcs.has_emit_labels())
    {
      CLOGERR << "In the following transducer:\n";
      peeled_composition.show_composite (ehmm_funcs, CLOGERR);
      THROWEXPR ("Alignment was requested, but no transducer states have emit labels. How do you want to score matches?");
    }

  // if we're using the composite transduer as an HMM later, check that it is jointly normalized
  if (use_composite_as_HMM && !ehmm_funcs.is_joint())
    THROWEXPR ("Oops -- composite transducer isn't jointly normalized -- bailing out!\n\n"
	       << "To do alignment (or null state elimination), the composite transducer\n"
	       << "must be jointly normalized (like an HMM). In other words, the input sequence\n"
	       << "must be forced to have zero length. To enforce this, you need to place\n"
	       << "a transducer on the root branch that has no match or delete states.");

  // save to dotfile, if specified
  if (dotfile.size())
    {
      ofstream dotfile_stream (dotfile.c_str());
      if (!dotfile_stream)
	THROWEXPR ("Couldn't open dotfile");
      ehmm_funcs.print_dotfile (dotfile_stream);
    }

  // do we go any further, using the composite transducer as an HMM for alignment?
  if (use_composite_as_HMM)
    {
      // build the EHMM_transducer_scores
      ehmm_scores = ehmm_funcs.eval_sc (peeled_composition.pscores);

      // get the list of observed sequences & "null" states that don't emit to observed sequences
      const vector<ENode> obs = peeled_composition.observed_nodes();
      vector<int> emit_states;
      ehmm_scores.get_emit_and_null_states (obs, emit_states, states_to_eliminate);

      // figure out which states are inaccessible
      // Note that this has already been attempted at the EHMM_transducer_funcs stage, but we do it again here in case any states are newly recognizable as being inaccessible based on the numerical label values
      // (if any transition labels evaluate to zero, some states might be inaccessible even though they were not recognized as such by EHMM_transducer_funcs).

      // FIXME: for parity with EHMM_transducer_funcs, this code should probably be moved to a separate method EHMM_transducer_scores::inaccessible_states
      // (NB beware of counting null states twice... once for being null, and once for being inaccessible... this may be a danger if a standalone EHMM_transducer_scores::inaccessible_states method is created.)

      // first look for paths from Start
      vector<int> start_accessible (ehmm_scores.states(), (int) 0);
      set<int> visited;
      stack<int> src_states;
      src_states.push (Grammar_state_enum::Start);
      while (src_states.size())
	{
	  const int src = src_states.top();
	  src_states.pop();

	  if (visited.find(src) == visited.end())
	    {
	      visited.insert (src);
	      for (int dest = 0; dest < ehmm_scores.states(); ++dest)
		if (ehmm_scores.transition (src, dest) > -InfinityScore)
		  {
		    start_accessible[dest] = 1;
		    src_states.push (dest);
		  }
	    }
	}

      // now look for paths from End
      vector<int> end_accessible (ehmm_scores.states(), (int) 0);
      visited.clear();
      stack<int> dest_states;
      dest_states.push (Grammar_state_enum::End);
      while (dest_states.size())
	{
	  const int dest = dest_states.top();
	  dest_states.pop();

	  if (visited.find(dest) == visited.end())
	    {
	      visited.insert (dest);
	      for (int src = 0; src < ehmm_scores.states(); ++src)
		if (ehmm_scores.transition (src, dest) > -InfinityScore)
		  {
		    end_accessible[src] = 1;
		    dest_states.push (src);
		  }
	    }
	}

      // add any inaccessible states to the elimination list
      for_const_contents (vector<int>, emit_states, s)    // if this code is moved to EHMM_transducer_scores, this should be a loop over ALL states, not just emit states
	if (!start_accessible[*s] || !end_accessible[*s])
	  states_to_eliminate.push_back (*s);

      // build the Eliminated_EHMM_transducer_scores
      elim_ehmm_scores = Eliminated_EHMM_transducer_scores (ehmm_scores, states_to_eliminate);

      // output the Eliminated_EHMM_transducer_scores, if requested
      if (acyclic)
	{
	  out << '\n';
	  elim_ehmm_scores.show_sexpr (out, peeled_composition.composite_name, &peeled_composition.pscores, false);
	  out.flush();
	}

      // dump the Eliminated_EHMM_transducer_scores to a HMMoC file, if specified
      if (hmmoc_opts.hmmoc_filename_prefix.size())
	{
	  HMMoC_adapter hmmoc_adapter (*this);
	  hmmoc_adapter.dump_hmmoc_model_to_file (hmmoc_opts.hmmoc_filename_prefix.c_str(),
						  hmmoc_opts.hmmoc_filename_prefix.c_str(),
						  hmmoc_opts.hmmoc_filename_prefix.c_str());
	}

      // sample Viterbi alignment?
      seq_vec = peeled_composition.make_seq_vec();
      if (viterbi)
	{
	  // use HMMoC?
	  bool do_viterbi_yourself = true;
	  if (use_hmmoc_for_viterbi)
	    {
	      Score hmmoc_vit_sc;
	      vector<int> hmmoc_vit_path;

	      HMMoC_adapter hmmoc_adapter (*this);
	      const bool hmmoc_succeeded = hmmoc_adapter.exec_Viterbi (hmmoc_vit_sc, hmmoc_vit_path);

	      if (hmmoc_succeeded)
		{
		  // output Viterbi likelihood & store
		  out << '\n';
		  dump_score (out, TSEXPR_VSCORE, hmmoc_vit_sc);
		  vit_ll = unpeeled_ll (hmmoc_vit_sc);

		  // output traceback
		  // we sample eliminated states here, but do not do the same after using the "native" Viterbi implementation
		  // this is because the hmmoc adapter always uses the acyclic (eliminated) EHMM, whereas our native Viterbi uses the composite (uneliminated)
		  dump_and_store_viterbi_path (out, elim_ehmm_scores.sample_eliminated (hmmoc_vit_path), vit_ll);

		  // no need to do it ourselves
		  do_viterbi_yourself = false;
		}
	    }

	  if (do_viterbi_yourself)
	    {
	      // populate Viterbi matrix
	      vit.trans_sc = &ehmm_scores;
	      vit.seq = seq_vec;

	      // allocate, fill
	      vit.alloc();
	      vit.fill();

	      // output Viterbi likelihood & store
	      out << '\n';
	      dump_score (out, TSEXPR_VSCORE, vit.end_sc);
	      vit_ll = unpeeled_ll (vit.end_sc);

	      // output Viterbi trace
	      dump_and_store_viterbi_path (out, vit.traceback(), vit_ll);
	    }
	}

      // sample Forward alignment(s)?
      if (fill_forward)
	{
	  // use HMMoC?
	  bool do_forward_yourself = true;
	  if (use_hmmoc_for_forward)
	    {
	      Score hmmoc_fwd_sc;
	      vector<vector<int> > hmmoc_fwd_path;

	      HMMoC_adapter hmmoc_adapter (*this);
	      const bool hmmoc_succeeded = hmmoc_adapter.exec_forward (hmmoc_fwd_sc, hmmoc_fwd_path, nforward);

	      if (hmmoc_succeeded)
		{
		  // output Forward likelihood & store
		  out << '\n';
		  dump_score (out, TSEXPR_FSCORE, hmmoc_fwd_sc);
		  fwd_ll = unpeeled_ll (hmmoc_fwd_sc);

		  // output tracebacks
		  for (int n = 0; n < nforward; ++n)
		    dump_and_store_forward_path (out, hmmoc_fwd_path[n]);

		  // no need to do it ourselves
		  do_forward_yourself = false;
		}
	    }

	  if (do_forward_yourself)
	    {
	      // populate forward matrix
	      fwd.trans_sc = &elim_ehmm_scores;
	      fwd.seq = seq_vec;

	      // allocate, fill
	      fwd.alloc();
	      fwd.fill();

	      // output Forward likelihood & store
	      out << '\n';
	      dump_score (out, TSEXPR_FSCORE, fwd.end_sc);
	      fwd_ll = unpeeled_ll (fwd.end_sc);

	      // sample tracebacks
	      for (int n = 0; n < nforward; ++n)
		{
		  const vector<int> acyclic_fwd_trace = fwd.sample_traceback();
		  dump_and_store_forward_path (out, acyclic_fwd_trace);
		}
	    }
	}

      // fill backward matrix
      if (fill_backward)
	{
	  // populate the forward-backward matrix
	  back.want_counts = want_expected_counts;
	  back.trans_sc = &elim_ehmm_scores;
	  back.seq = seq_vec;

	  // allocate & fill forward & backward DP matrices
	  if (fill_forward)  // forward already filled?
	    back.init_fwd (fwd);
	  back.alloc();
	  back.fill();

	  // show counts?
	  if (want_expected_counts)
	    {
	      if (peeling_required)
		CLOGERR << "Warning: posterior expected label usage counts do not yet include\n"
			<< "transitions and emissions from the constrained part of the tree...\n";
	      pcounts = PCounts (peeled_composition.pscores);
	      back.inc_var_counts (ehmm_funcs, pcounts, peeled_composition.pscores, &elim_ehmm_scores);
	      PFunc_builder::pcounts2stream (out, pcounts, TSEXPR_POSTEXPECT, (PCounts*) 0, false);
	      out.flush();
	    }

	  // do optimal accuracy alignment?
	  if (optacc)
	    {
	      // allocate & fill optimal accuracy DP matrix; get traceback alignment
	      back.init_sumpairs_reward();
	      back.alloc_optacc();
	      back.fill_optacc();
	      optacc_path = back.optacc_traceback();

	      // print optimal accuracy alignment path
	      out << "\n(" << TSEXPR_OPTACC << " (" << TSEXPR_TYPE << " (";
	      for (int col = 0; col < optacc_path.columns(); ++col)
		{
		  vector<sstring> match_node;
		  for (int row = 0; row < optacc_path.rows(); ++row)
			if (optacc_path (row, col))
			  match_node.push_back (peeled_composition.tape_name[back.observed_seqs[row]]);
		  out << (col > 0 ? " (" : "(") << match_node << ')';
		}
	      out << ")))\n";
	      out.flush();
	    }
	}
    }

  // Redelings-Suchard MCMC kernel
  if (redelings_suchard)
    {
      Transducer_SExpr_file& pc (peeled_composition);  // shorthand

      // check validity
      if (pc.path.size() > 1 || (pc.path.size() == 1 && pc.path.find(0) == pc.path.end()))  // this test should always be false, due to earlier code
	THROWEXPR ("Redelings-Suchard kernel only implemented for peeled composition");  // should be unreachable

      if (!pc.proposal_branches.size())
	THROWEXPR ("Redelings-Suchard kernel requested, but no proposal schedule specified!");

      if (evaluate_redelings_suchard_inverse_move && !got_old_loglike)
	{
	  if (composition.free_clique < 0)
	    THROWEXPR ("Redelings-Suchard MCMC kernel requires at least some branch paths unspecified");

	  THROWEXPR ("Redelings-Suchard MCMC kernel requires complete specification of old state, but the following branch paths were unspecified:\n" << pathless_branch_names);
	}

      // loop over proposal compositions
      Transducer_SExpr_file::NodePathMap new_path = pc.path;
      set<int> node_set, branches_seen;
      for_const_contents (Transducer_SExpr_file::RedSuchSchedule, pc.proposal_branches, branches)
	{
	  // figure out which branches are constrained and which are unconstrained; also which nodes are included, and which is the proposal root
	  vector<int> cons_branch, uncons_branch;
	  int root = pc.etree.nodes();
	  for_const_contents (vector<int>, *branches, b)
	    {
	      const int parent = pc.etree.parent[*b];
	      node_set.insert (*b);
	      if (*b < root)
		root = *b;
	      if (parent >= 0)
		{
		  node_set.insert (parent);
		  if (parent < root)
		    root = parent;
		}
	      if (branches_seen.find (*b) == branches_seen.end())
		{
		  uncons_branch.push_back (*b);
		  branches_seen.insert (*b);
		}
	      else
		cons_branch.push_back (*b);
	    }
	  if (root == pc.etree.nodes())
	    THROWEXPR ("Couldn't find root of proposal subtree");

	  // build node index maps between full<-->proposal trees
	  const vector<int> node_prop2full (node_set.begin(), node_set.end());
	  vector<int> node_full2prop;
	  for (int n_full = 0; n_full < pc.etree.nodes(); ++n_full)
	    {
	      const int n_prop = find (node_prop2full.begin(), node_prop2full.end(), n_full) - node_prop2full.begin();
	      node_full2prop.push_back (n_prop == (int) node_prop2full.size() ? -1 : n_prop);
	    }

	  // build ETree
	  ETree prop_etree = pc.etree.subtree (node_prop2full);

	  // populate branch_trans, branch_trans_name, prof, path, branch_name & tape_name for proposal composition
	  // also create path & old_path for "inverse move" composition (other side of Hastings ratio)
	  Transducer_SExpr_file::NodePathMap prop_path, inv_path, inv_old_path;
	  const Transducer_SExpr_file::NodePathMap dummy_path;  // always empty
	  Transducer_SExpr_file::NodeProfileMap prop_prof;
	  vector<sstring> prop_branch_trans_name (prop_etree.nodes());
	  vector<Pair_transducer_funcs> prop_branch_trans (prop_etree.nodes());
	  vector<sstring> prop_branch_name (prop_etree.nodes());
	  map<int,sstring> prop_tape_name;
	  map<int,double> prop_band_coeff;

	  // prop_prof
	  for (int n_prop = 0; n_prop < (int) node_prop2full.size(); ++n_prop)
	    if (pc.prof.find (node_prop2full[n_prop]) != pc.prof.end())
	      prop_prof[n_prop] = pc.prof[node_prop2full[n_prop]];

	  // prop_path: constrained branches
	  for_const_contents (vector<int>, cons_branch, b_full)
	    prop_path[node_full2prop[*b_full]] = new_path[*b_full];

	  // inv_path: constrained branches
	  for_const_contents (vector<int>, cons_branch, b_full)
	    if (pc.old_path.find (*b_full) != pc.old_path.end())
	      inv_path[node_full2prop[*b_full]] = pc.old_path[*b_full];

	  // inv_old_path == prop_path ?  (i.e. is new_path identical to old_path on the constrained branches?)
	  // if so, this is a "virgin proposal", allowing us to re-use some calculations between prop_move and inv_move
	  const bool virgin_proposal = inv_old_path == prop_path;
	  const bool calc_old_path_ll_during_prop_step = virgin_proposal && propose_redelings_suchard_move && evaluate_redelings_suchard_inverse_move;

	  // inv_old_path (or inv_path if virgin proposal): unconstrained branches
	  for_const_contents (vector<int>, uncons_branch, b_full)
	    if (pc.old_path.find (*b_full) != pc.old_path.end())
	      (calc_old_path_ll_during_prop_step ? inv_path : inv_old_path)[node_full2prop[*b_full]] = pc.old_path[*b_full];

	  // prop_band_coeff
	  for_const_contents (vector<int>, uncons_branch, b_full)
	    if (pc.band_coeff.find (*b_full) != pc.band_coeff.end())
	      prop_band_coeff[node_full2prop[*b_full]] = pc.band_coeff[*b_full];

	  // prop_branch_trans
	  for_const_contents (vector<int>, *branches, b_full)
	    {
	      const int b_prop = node_full2prop[*b_full];
	      if (b_prop >= 0)
		{
		  prop_branch_trans[b_prop] = pc.branch_trans[*b_full];
		  prop_branch_trans_name[b_prop] = pc.branch_trans_name[*b_full];
		}
	    }

	  // prop_branch_name, prop_tape_name
	  for (int n_prop = 1; n_prop < prop_etree.nodes(); ++n_prop)
	    {
	      prop_branch_name[n_prop] = pc.branch_name[node_prop2full[n_prop]];
	      prop_tape_name[n_prop] = pc.tape_name[node_prop2full[n_prop]];
	    }

	  // if we have a new root, place singleton transducer on root branch; otherwise, use root transducer from full tree
	  if (root == 0)
	    {
	      prop_branch_trans[0] = pc.branch_trans[0];
	      prop_branch_trans_name[0] = pc.branch_trans_name[0];
	      prop_branch_name[0] = pc.branch_name[0];
	      prop_tape_name[0] = pc.tape_name[0];
	      // copy across path and/or old_path for root branch
	      if (pc.path.find (0) != pc.path.end())
		prop_path[0] = inv_path[0] = pc.path[0];
	      else if (pc.old_path.find (0) != pc.old_path.end())
		inv_old_path[0] = pc.old_path[0];
	    }
	  else
	    {
	      Singleton_transducer_funcs sing_funcs (pc.alphabet_size);
	      prop_branch_trans[0] = sing_funcs;
	      prop_branch_trans_name[0] = sstring (Singleton_transducer_name);
	      // place an all-insert path on the root branch, so it's properly constrained
	      if (prop_prof.find (0) != prop_prof.end())
		prop_path[0] = inv_path[0] = vector<int> (prop_prof[0].size(), sing_funcs.insert_state());
	    }

	  // proposal move
	  if (propose_redelings_suchard_move)
	    {
	      // build composition
	      // for economy, if this was a "virgin proposal", then use this composition to evaluate the score of inv_old_path
	      Transducer_SExpr_file prop_comp (pc.alphabet, pc.pscores, prop_etree, prop_branch_trans, prop_prof, prop_path, calc_old_path_ll_during_prop_step ? inv_old_path : dummy_path);
	      prop_comp.rebuild_tree_names (prop_tape_name, prop_branch_name, prop_branch_trans_name);
	      prop_comp.band_coeff = prop_band_coeff;

	      if (CTAGGING(1,REDSUCH_DEBUG))
		{
		  typedef map<int,sstring> NameMap;
		  for_const_contents (NameMap, prop_tape_name, ptn)
		    CL << "prop_tape_name[" << ptn->first << "] = " << ptn->second << '\n';
		  for_const_contents (NameMap, prop_comp.tape_name, ptn)
		    CL << "prop_comp.tape_name[" << ptn->first << "] = " << ptn->second << '\n';
		}

	      // do the move
	      Handel_movement prop_move;
	      prop_move.composition = prop_comp;
	      prop_move.nforward = 1;  // compute forward score & sample 1 path
	      prop_move.hmmoc_opts = hmmoc_opts;
	      prop_move.use_hmmoc_for_forward = use_hmmoc_for_forward;
	      prop_move.composite = true;  // we will want to see the composite transducer in the event of an error

	      sstring prop_comp_output;
	      try {
		prop_move.dump_composition (prop_comp_output);
	      } catch (const Dart_exception& e) {
		CLOGERR << e.what() << "Proposal composition:\n" << prop_comp_output;
		THROWEXPR ("Proposal error");
	      }
	      if (CTAGGING(1,REDSUCH_PROPOSAL))
		CL << "Proposal composition:\n" << prop_comp_output;

	      // check for out-of-band old_path
	      if (calc_old_path_ll_during_prop_step)
		old_path_is_outside_hmmoc_band = old_path_is_outside_hmmoc_band || prop_move.old_path_is_outside_hmmoc_band;

	      // get the trace & stick it in new_path
	      const vector<int>& prop_trace = prop_move.fwd_trace[0];
	      const vector<vector<int> > prop_branch_path = prop_move.ehmm_scores.branch_paths (prop_trace);

	      vector<int> node_prop2free (node_prop2full.size(), -1);
	      int total_free = 0;
	      for_const_contents (set<int>, prop_move.composition.free_clique_set(), n)
		node_prop2free[*n] = total_free++;

	      for_const_contents (vector<int>, uncons_branch, n_full)
		{
		  const int n_prop = node_full2prop[*n_full];
		  const int n_free = node_prop2free[n_prop];
		  const vector<int>& n_path = prop_branch_path[n_free];
		  if (CTAGGING(1,REDSUCH_PATH))
		    CL << "Copying path for branch " << prop_branch_name[n_prop] << " (" << n_path << ")\n";
		  new_path[*n_full] = n_path;
		}

	      // store Hastings ratio term (denominator)
	      Transducer_SExpr_file& ppc = prop_move.peeled_composition;
	      Transducer_SExpr_file::NodePathMap prop_path_map = ppc.get_path_map (prop_move.ehmm_funcs, prop_trace);
	      const Loge prop_path_ll = NatsPMul (ppc.get_path_loglike (prop_path_map), prop_move.peeling_loglike);
	      hastings_term.push_back (-NatsPMul (prop_path_ll, -prop_move.fwd_ll));

	      // store Hastings ratio term (numerator)
	      if (calc_old_path_ll_during_prop_step)
		hastings_term.push_back (NatsPMul (prop_move.old_loglike, -prop_move.fwd_ll));

	      // display prop_comp
	      out << "\n(" << TSEXPR_PROPOSAL;
	      prop_comp.show_tree (out, 1, false);
	      out << ' ';
	      if (calc_old_path_ll_during_prop_step)
		prop_move.dump_raw_loglike (out, TSEXPR_OLD_SCORE, prop_move.old_loglike);
	      prop_move.dump_raw_loglike (out, TSEXPR_FSCORE, prop_move.fwd_ll);
	      prop_move.dump_peeled_composite_trace (out, TSEXPR_FPATH, prop_trace, prop_path_ll);
	      out << ")\n";
	      out.flush();
	    }

	  // inverse move
	  if (evaluate_redelings_suchard_inverse_move && !calc_old_path_ll_during_prop_step)
	    {
	      // build composition
	      Transducer_SExpr_file inv_comp (pc.alphabet, pc.pscores, prop_etree, prop_branch_trans, prop_prof, inv_path, inv_old_path);
	      inv_comp.rebuild_tree_names (prop_tape_name, prop_branch_name, prop_branch_trans_name);
	      inv_comp.band_coeff = prop_band_coeff;

	      // do the move
	      Handel_movement inv_move;
	      inv_move.composition = inv_comp;
	      inv_move.nforward = 0;  // compute forward score, but don't sample any paths
	      inv_move.hmmoc_opts = hmmoc_opts;
	      inv_move.use_hmmoc_for_forward = use_hmmoc_for_forward;
	      inv_move.composite = true;  // we will want to see the composite transducer in the event of an error

	      sstring inv_comp_output;
	      try {
		inv_move.dump_composition (inv_comp_output);
	      } catch (const Dart_exception& e) {
		CLOGERR << e.what() << "Inverse proposal composition:\n" << inv_comp_output;
		THROWEXPR ("Inverse proposal error");
	      }
	      if (CTAGGING(1,REDSUCH_PROPOSAL))
		CL << "Inverse proposal composition:\n" << inv_comp_output;

	      // check for out-of-band old_path
	      old_path_is_outside_hmmoc_band = old_path_is_outside_hmmoc_band || inv_move.old_path_is_outside_hmmoc_band;

	      // store Hastings ratio term (numerator)
	      hastings_term.push_back (NatsPMul (inv_move.old_loglike, -inv_move.fwd_ll));

	      // display inv_comp
	      out << "\n(" << TSEXPR_INVERSE;
	      inv_comp.show_tree (out, 1, false);
	      out << ' ';
	      inv_move.dump_raw_loglike (out, TSEXPR_OLD_SCORE, inv_move.old_loglike);
	      inv_move.dump_raw_loglike (out, TSEXPR_FSCORE, inv_move.fwd_ll);
	      out << ")\n";
	      out.flush();
	    }
	}

      // store proposed path, print score, etc
      if (propose_redelings_suchard_move)
	{
	  // map node indices from peeled_composition to main composition
	  const set<int>& free_clique_set = composition.free_clique_set();
	  const vector<int> free_clique_vec (free_clique_set.begin(), free_clique_set.end());
	  for (int n = 0; n < pc.etree.nodes(); ++n)
	    {
	      if (new_path.find (n) == new_path.end())
		THROWEXPR ("Failed to propose a path for node " << n);
	      if (CTAGGING(1,REDSUCH_PATH))
		CL << "Copying path from peeled-node #" << n << " (" << pc.tape_name[n]
		   << ") to unpeeled-node #" << free_clique_vec[n] << " (" << composition.tape_name[free_clique_vec[n]]
		   << "): (" << new_path[n] << ")\n";
	      redsuch_path[free_clique_vec[n]] = new_path[n];
	    }

	  // get score of new_path
	  redsuch_ll = NatsPMul (pc.get_path_loglike (new_path), peeling_loglike);

	  // display path (creating composite path only if we have an EHMM with which to look up composite states)
	  if (ehmm_funcs_defined)
	    {
	      // compose new_path; place in redsuch_trace
	      vector<vector<int> > branch_path;
	      for (int n = 0; n < pc.etree.nodes(); ++n)
		branch_path.push_back (redsuch_path[n]);

	      redsuch_trace = ehmm_funcs.composite_path (branch_path);

	      // debug: check that composite_path worked
	      if (CTAGGING(1,REDSUCH_COMPOSITE_PATH))
		{
		  const vector<vector<int> > bp = ehmm_funcs.branch_paths (redsuch_trace);
		  CTAG(1,REDSUCH_COMPOSITE_PATH) << "Comparing original branch paths to decomposed composite path (" << pc.etree.nodes() << " nodes):\n";
		  for (int n = 0; n < pc.etree.nodes(); ++n)
		    {
		      CL << "  redsuch_path[" << n << "]:";
		      for_const_contents (vector<int>, redsuch_path[n], s)
			CL << ' ' << ehmm_funcs.branch_transducer[n].get_state_name(*s);
		      CL << '\n';
		      CL << "           bp[" << n << "]:";
		      for_const_contents (vector<int>, bp[n], s)
			CL << ' ' << ehmm_funcs.branch_transducer[n].get_state_name(*s);
		      CL << '\n';
		      CL << (bp[n] == redsuch_path[n] ? " (identical)\n" : " (different!)\n");
		    }
		}

	      // display redsuch_trace
	      dump_peeled_composite_trace (out, TSEXPR_RSPATH, redsuch_trace, redsuch_ll);
	    }
	  else
	    {
	      if (CTAGGING(1,REDSUCH_PATH))
		for_const_contents (Transducer_SExpr_file::NodePathMap, redsuch_path, np)
		  CL << "redsuch_path[" << np->first << "] = (" << np->second << ")\n";
	      dump_unpeeled_node_path_map (out, TSEXPR_RSPATH, redsuch_path, redsuch_ll);
	    }

	  // store Hastings ratio term (numerator)
	  hastings_term.push_back (redsuch_ll);
	}

      // store Hastings ratio term for old path (denominator)
      if (evaluate_redelings_suchard_inverse_move)
	hastings_term.push_back (-old_loglike);

      // display Hastings ratio terms
      if (hastings_term.size())
	{
	  vector<Log2> hastings_term_bits;
	  Loge log_hastings_ratio = 0.;
	  for_const_contents (vector<Loge>, hastings_term, h_ll)
	    {
	      NatsPMulAcc (log_hastings_ratio, *h_ll);
	      hastings_term_bits.push_back (Nats2Bits (*h_ll));
	    }
	  out << "\n(" << TSEXPR_HTERM << " (" << hastings_term_bits
	      << "))\n";

	  // if we have all the necessary parts, then compute & print Hastings ratio
	  if (propose_redelings_suchard_move && evaluate_redelings_suchard_inverse_move)
	    out << '(' << TSEXPR_HASTINGS << ' ' << -Nats2Bits (min (log_hastings_ratio, 0.)) << ")\n";
	}
    }
}

Loge Handel_movement::unpeeled_ll (Score peeled_sc)
{
  const Loge peeled_ll = Score2Nats (peeled_sc);
  return
    peeling_required || fully_constrained
    ? NatsPMul (peeling_loglike, peeled_ll)
    : peeled_ll;
}

void Handel_movement::dump_score (ostream& out, const char* tag, Score sc)
{
  dump_raw_loglike (out, tag, unpeeled_ll (sc));
}

void Handel_movement::dump_raw_loglike (ostream& out, const char* tag, Loge ll)
{
  out << '(' << tag << ' ' <<
    peeled_composition.score_sexpr (ll)
      << ")\n";
  out.flush();
}

void Handel_movement::dump_peeled_composite_trace (ostream& out, const char* tag, const vector<int>& path, Loge path_ll)
{
  out << "\n(" << tag << '\n';
  peeled_composition.show_path_with_tree (ehmm_funcs, path, path_ll, out, 1);
  out << ")\n";
  out.flush();
}

void Handel_movement::dump_unpeeled_node_path_map (ostream& out, const char* tag, const Transducer_SExpr_file::NodePathMap& unpeeled_path, Loge path_ll)
{
  const set<int>& free_clique_set = composition.free_clique_set();
  const vector<int> free_clique_vec (free_clique_set.begin(), free_clique_set.end());

  Transducer_SExpr_file::NodePathMap peeled_path;
  for (int n_free = 0; n_free < (int) free_clique_vec.size(); ++n_free)
    peeled_path[n_free] = ((Transducer_SExpr_file::NodePathMap&)unpeeled_path)[free_clique_vec[n_free]];  // cast away const

  out << "\n(" << tag << '\n';
  peeled_composition.show_path_with_tree (peeled_path, path_ll, out, 1);
  out << ")\n";
  out.flush();
}

void Handel_movement::dump_and_store_forward_path (ostream& out, const vector<int>& acyclic_path)
{
  const vector<int> cyclic_path = elim_ehmm_scores.sample_eliminated (acyclic_path);

  Transducer_SExpr_file& pc = peeled_composition;
  Transducer_SExpr_file::NodePathMap path_map = pc.get_path_map (ehmm_funcs, cyclic_path);
  const Loge path_ll = NatsPMul (pc.get_path_loglike (path_map), peeling_loglike);

  fwd_trace.push_back (cyclic_path);
  dump_peeled_composite_trace (out, TSEXPR_FPATH, cyclic_path, path_ll);
}

void Handel_movement::dump_and_store_viterbi_path (ostream& out, const vector<int>& composite_path, Loge path_ll)
{
  vit_trace = composite_path;
  dump_peeled_composite_trace (out, TSEXPR_VPATH, vit_trace, path_ll);
}
