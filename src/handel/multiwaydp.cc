#include "handel/multiwaydp.h"

// Transducer_tree_emission
Transducer_tree_emission::Transducer_tree_emission (const vector<const Symbol_score_map*>& node_scores,
						    const vector<int>& branch_trans_state,
						    ENode inserter,
						    const vector<ENode>& deleters,
						    const vector<ENode>& matchers)
  : node_scores (node_scores),
    branch_trans_state (branch_trans_state),
    inserter (inserter),
    deleters (deleters),
    matchers (matchers)
{
  state_name << '(' << branch_trans_state << ')';
}


// Transducer_peeler
void Transducer_peeler::init (int a, const ETree& e, const vector<Pair_transducer_scores>& bt)
{
  alphabet_size = a;
  etree = e;
  branch_transducer = bt;

  if (alphabet_size == 0)
    THROWEXPR ("Zero alphabet size!");

  subtrees = ESubtrees (etree);

  F_matrix = vector<vector<Score> > (etree.nodes(), vector<Score> (alphabet_size));
  E_matrix = vector<vector<Score> > (etree.nodes(), vector<Score> (alphabet_size));
  G_matrix = vector<vector<Score> > (etree.nodes(), vector<Score> (alphabet_size));

  branch_trans_counts = vector<Pair_transducer_counts*> (branch_transducer.size(), (Pair_transducer_counts*) 0);
}

void Transducer_peeler::init (int a, const ETree& e, const vector<Pair_transducer_scores>& bt, vector<Pair_transducer_counts>& btc)
{
  init (a, e, bt);
  for (int node = 0; node < (int) btc.size(); ++node)
    branch_trans_counts[node] = &btc[node];
}

void Transducer_peeler::prune (const Transducer_tree_emission& emit)
{
  if (emit.inserter < 0)
    THROWEXPR ("In state (" << emit.branch_trans_state << "): root branch is in an absorbing state; can't calculate emit score by pruning");

  // Pre-fill F_matrix from node_scores, initializing allowed symbols at all nodes
  for (int n = 0; n < (int) emit.node_scores.size(); ++n)
    if (emit.node_scores[n])
      {
	fill (F_matrix[n].begin(), F_matrix[n].end(), -InfinityScore);
	for_const_contents (Symbol_score_map, *emit.node_scores[n], sym_sc)
	  F_matrix[n][sym_sc->first] = sym_sc->second;
      }
    else
      fill (F_matrix[n].begin(), F_matrix[n].end(), 0);

  // propagate [sym-->] deleters upwards through F_matrix
  for_const_contents (vector<ENode>, emit.deleters, n)
    {
      Pair_transducer_scores& dtrans = branch_transducer[*n];
      const int dstate = emit.branch_trans_state[*n];
      const ENode p = etree.parent[*n];

      for (int sym = 0; sym < alphabet_size; ++sym)
	{
	  // build E_matrix
	  const Score del_sc = dtrans.delete_val (dstate, sym);
	  // store E_matrix, accumulate in F_matrix
	  E_matrix[*n][sym] = del_sc;
	  ScorePMulAcc (F_matrix[p][sym], del_sc);
	}
    }

  // propagate [xsym-->ysym] matchers upwards through F_matrix
  for_const_contents (vector<ENode>, emit.matchers, n)
    {
      Pair_transducer_scores& mtrans = branch_transducer[*n];
      const int mstate = emit.branch_trans_state[*n];
      const ENode p = etree.parent[*n];

      for (int xsym = 0; xsym < alphabet_size; ++xsym)
	{
	  // build E_matrix
	  Score branch_sc = -InfinityScore;
	  for (int ysym = 0; ysym < alphabet_size; ++ysym)
	    ScorePSumAcc (branch_sc, ScorePMul (F_matrix[*n][ysym],
						mtrans.match_val (mstate, xsym, ysym)));

	  // store E_matrix, accumulate in F_matrix
	  E_matrix[*n][xsym] = branch_sc;
	  ScorePMulAcc (F_matrix[p][xsym], branch_sc);
	}
    }

  // terminate recursion at [-->ysym] inserter, summarizing F_matrix at pruning_sc
  Pair_transducer_scores& itrans = branch_transducer[emit.inserter];
  const int istate = emit.branch_trans_state[emit.inserter];

  pruning_sc = -InfinityScore;
  for (int ysym = 0; ysym < alphabet_size; ++ysym)
    ScorePSumAcc (pruning_sc, ScorePMul (F_matrix[emit.inserter][ysym],
					 itrans.insert_val (istate, ysym)));

  // log
  if (CTAGGING(1,TRANSDUCER_F_MATRIX))
    {
      CL << "Transducer pruning matrix for state " << emit.state_name
	 << ": matchers=(" << emit.matchers << ") deleters=(" << emit.deleters << ") inserter=" << emit.inserter << '\n';
      for_const_contents (vector<ENode>, emit.matchers, n)
	{
	  CL << (find (emit.deleters.begin(), emit.deleters.end(), *n) == emit.deleters.end() ? "Match " : "Delete")
	     << " node " << *n << ':';
	  for (int sym = 0; sym < alphabet_size; ++sym)
	    {
	      CL << ' ';
	      ShowScore (F_matrix[*n][sym], CL);
	    }
	  if (emit.node_scores[*n])
	    {
	      CL << "    (node scores:";
	      for_const_contents (Symbol_score_map, *emit.node_scores[*n], sym_sc)
		CL << ' ' << sym_sc->first << "=>" << sym_sc->second;
	      CL << ')';
	    }
	  CL << '\n';
	}
      CL << "Insert node " << emit.inserter << ':';
      for (int sym = 0; sym < alphabet_size; ++sym)
	{
	  CL << ' ';
	  ShowScore (F_matrix[emit.inserter][sym], CL);
	}
      if (emit.node_scores[emit.inserter])
	{
	  CL << "    (node scores:";
	  for_const_contents (Symbol_score_map, *emit.node_scores[emit.inserter], sym_sc)
	    CL << ' ' << sym_sc->first << "=>" << sym_sc->second;
	  CL << ')';
	}
      CL << '\n';
      CL << "Final score: " << pruning_sc << '\n';
    }
}

void Transducer_peeler::peel (const Transducer_tree_emission& emit, double weight)
{
  // --- Begin copied definitions from comments in Transducer_peeler class prototype

  // Let N be an ENode tree index and I,J,K be alphabet symbols.
  // Let T[N] be all observations in subtree rooted at node N.
  // Let Parent[N] be the parent of N.
  // F_matrix[N][I] = Prob2Score P({observations in T[N]} | state at N is I)
  // E_matrix[N][J] = Prob2Score P({observations in T[N]} | state at Parent[N] is J)
  // G_matrix[N][K] = Prob2Score P({observations NOT in T[N]} AND state at N is K)

  // --- End copied definitions

  // Relevant posterior probabilities for computing expected counts:

  // M(N,I,J) = P(state at Parent[N] is I, state at matcher N is J | observations)
  // D(N,J) = P(state at deleter N is J | observations)
  // I(N,J) P(state at inserter N is J | observations)

  // Let E(N,I) = Score2Prob(E_matrix[N][I]), & similarly for F_matrix, G_matrix.
  // Define the G-E product:
  // GE(N,I) = P({observations NOT in T[N]} AND state at Parent[N] is I)
  //         = G(Parent[N],I) * \product_{siblings S of N} E(S,I)

  // Let transMatch(N,I,J), transInsert(N,J) and transDelete(N,I)
  //  be the subscripted transducer match, insert and delete emit labels.

  // Then:
  // M(N,I,J) = GE(N,I) * transMatch(N,I,J) * F(N,J) / P(observations)
  //   D(N,I) = GE(N,I) * transDelete(N,I)           / P(observations)
  //   I(N,J) =           transInsert(N,J)  * F(N,J) / P(observations)

  // Because of the way we encode observation constraints in the initialization of G,
  // I(N,J) is actually calculated using the following (equivalent) formula:

  //   I(N,J) =                     G(N,J)  * F(N,J) / P(observations)

  // Now the actual code...

  // pre-fill G_matrix with zeroes
  for (int n = 0; n < (int) emit.node_scores.size(); ++n)
    /* used to fill G_matrix with node_scores, like F_matrix; now think this is a buggy mistake
    if (emit.node_scores[n])
      {
	fill (G_matrix[n].begin(), G_matrix[n].end(), -InfinityScore);
	for_const_contents (Symbol_score_map, *emit.node_scores[n], sym_sc)
	  G_matrix[n][sym_sc->first] = sym_sc->second;
      }
    else
    */
    fill (G_matrix[n].begin(), G_matrix[n].end(), 0);

  // initialize G_matrix at inserter [-->ysym]
  Pair_transducer_scores& itrans = branch_transducer[emit.inserter];
  const int istate = emit.branch_trans_state[emit.inserter];

  Pair_transducer_counts* btc = branch_trans_counts[emit.inserter];
  for (int ysym = 0; ysym < alphabet_size; ++ysym)
    {
      ScorePMulAcc (G_matrix[emit.inserter][ysym], itrans.insert_val (istate, ysym));

      // accumulate ysym insert counts
      if (btc)
	btc->insert_val(istate,ysym) += weight * Score2Prob (ScorePMul3 (G_matrix[emit.inserter][ysym],
									 F_matrix[emit.inserter][ysym],
									 -pruning_sc));
    }

  // propagate [psym-->nsym] matchers downwards through G_matrix
  for_const_reverse_contents (vector<ENode>, emit.matchers, n)
    {
      Pair_transducer_scores& mtrans = branch_transducer[*n];
      const int mstate = emit.branch_trans_state[*n];
      const ENode p = etree.parent[*n];
      const ENodeVec& siblings = subtrees.children[p];
      Pair_transducer_counts* btc = branch_trans_counts[*n];

      vector<Score> branch_sc (alphabet_size, -InfinityScore);
      for (int psym = 0; psym < alphabet_size; ++psym)
	{
	  // Take product of G_matrix and E_matrix:
	  // GE[psym] = G[parent][psym] \prod_{siblings} E[sibling][psym]
	  // = parent_sc
	  Score parent_sc = G_matrix[p][psym];
	  for_const_contents (vector<ENode>, siblings, s)
	    if (*s != *n)
	      ScorePMulAcc (parent_sc, E_matrix[*s][psym]);

	  // accumulate GE[psym]*match[psym->nsym] messages in G_matrix[nsym]
	  for (int nsym = 0; nsym < alphabet_size; ++nsym)
	    {
	      const Score GE_match_sc = ScorePMul (parent_sc,
						   mtrans.match_val (mstate, psym, nsym));
	      ScorePSumAcc (branch_sc[nsym], GE_match_sc);

	      // accumulate psym->nsym match counts
	      if (btc)
		btc->match_val(mstate,psym,nsym) += weight * Score2Prob (ScorePMul3 (GE_match_sc,
										     F_matrix[*n][nsym],
										     -pruning_sc));
	    }
	}

      for (int nsym = 0; nsym < alphabet_size; ++nsym)
	ScorePMulAcc (G_matrix[*n][nsym], branch_sc[nsym]);
    }

  // accumulate counts at [psym-->] deleters
  for_const_contents (vector<ENode>, emit.deleters, n)
    {
      Pair_transducer_counts* btc = branch_trans_counts[*n];
      if (btc)
	{
	  Pair_transducer_scores& dtrans = branch_transducer[*n];
	  const int dstate = emit.branch_trans_state[*n];
	  const ENode p = etree.parent[*n];
	  const ENodeVec& siblings = subtrees.children[p];

	  for (int psym = 0; psym < alphabet_size; ++psym)
	    {
	      // Take product of G_matrix and E_matrix:
	      // GE[psym] = G[parent][psym] \prod_{siblings} E[sibling][psym]
	      // = parent_sc
	      Score parent_sc = G_matrix[p][psym];
	      for_const_contents (vector<ENode>, siblings, s)
		if (*s != *n)
		  ScorePMulAcc (parent_sc, E_matrix[*s][psym]);

	      // accumulate psym delete counts
	      btc->delete_val(dstate,psym) += weight * Score2Prob (ScorePMul3 (parent_sc,
									       dtrans.delete_val(dstate,psym),
									       -pruning_sc));
	    }
	}
    }

  // log
  if (CTAGGING(1,TRANSDUCER_G_MATRIX))
    {
      CL << "Transducer peeling matrix for state " << emit.state_name
	 << ": matchers=(" << emit.matchers << ") deleters=(" << emit.deleters << ") inserter=" << emit.inserter << '\n';
      for_const_contents (vector<ENode>, emit.matchers, n)
	{
	  CL << (find (emit.deleters.begin(), emit.deleters.end(), *n) == emit.deleters.end() ? "Match " : "Delete")
	     << " node " << *n << ':';
	  for (int sym = 0; sym < alphabet_size; ++sym)
	    {
	      CL << ' ';
	      ShowScore (G_matrix[*n][sym], CL);
	    }
	  CL << '\n';
	}
      CL << "Insert node " << emit.inserter << ':';
      for (int sym = 0; sym < alphabet_size; ++sym)
	{
	  CL << ' ';
	  ShowScore (G_matrix[emit.inserter][sym], CL);
	}
      CL << '\n';
    }
}

Symbol_score_map Transducer_peeler::post_sc (ENode node, bool normalize) const
{
  Symbol_score_map ssm;
  for (int sym = 0; sym < alphabet_size; ++sym)
    {
      const Score sc = ScorePMul3 (F_matrix[node][sym], G_matrix[node][sym], (normalize ? -pruning_sc : 0));
      if (sc > -InfinityScore)
	ssm[sym] = sc;
    }
  return ssm;
}

// Transducer_DP_base
Transducer_DP_base::Transducer_DP_base()
  : trans_sc (0), want_counts (false), alphabet_size (0), kT (1.), ehmm_transition_counts (0, 0.)
{ }

void Transducer_DP_base::alloc()
{
  if (CTAGGING(1,TRANSDUCER TRANSDUCER_ALLOC))
    {
      CL << "Transducer_DP_base::alloc() called for the following transducer:\n";
      trans_sc->show (CL);
    }

  inserter = vector<ENode> (trans_sc->states());
  deleters = vector<vector<ENode> > (trans_sc->states(), vector<ENode>());
  matchers = vector<vector<ENode> > (trans_sc->states(), vector<ENode>());
  for (int s = 0; s < trans_sc->states(); ++s)
    trans_sc->summarize_emission (s, inserter[s], deleters[s], matchers[s]);

  branch_trans_counts.clear();
  if (want_counts)
    {
      ehmm_transition_counts = Transducer_counts (*trans_sc, 0.);
      // initialize peeler to collect counts
      for (int node = 0; node <= trans_sc->max_node_index(); ++node)
	branch_trans_counts.push_back (Pair_transducer_counts (trans_sc->branch_transducer[node]));
      peeler.init (trans_sc->alphabet_size, trans_sc->etree, trans_sc->branch_transducer, branch_trans_counts);
    }
  else
    peeler.init (trans_sc->alphabet_size, trans_sc->etree, trans_sc->branch_transducer);  // init peeler without counts

  node_scores = vector<const Symbol_score_map*> (seq.size(), (const Symbol_score_map*) 0);

  observed_seqs.clear();
  for (int i = 0; i < (int) seq.size(); ++i)
    if (seq[i])
      observed_seqs.push_back (i);

  incoming = Transition_methods::incoming_states (*trans_sc);

  vector<int> all_states (trans_sc->states());
  for (int i = 0; i < trans_sc->states(); ++i)
    all_states[i] = i;
  outgoing = Transition_methods::selected_outgoing_states (*trans_sc, all_states);
}

Score Transducer_DP_base::calc_cell_emit (int state, const vector<int>& seq_coords, Prob weight)
{
  // get state type; bail with zero score if a wait state
  const State_type emit_type = trans_sc->state_type[state];
  if (emit_type == 0)
    return 0;

  // populate node_scores
  for (int i = 0, n = 0; n < (int) seq.size(); ++n)  // n is index into node_scores[] and seq[]; i is index into seq_coords[] and observed_seq[]
    {
      if (type_node_emit (emit_type, n))
	if (seq[n])
	  if (seq_coords[i] > 0)
	    node_scores[n] = &((*seq[n])[seq_coords[i] - 1]);
	  else
	    return -InfinityScore;   // bail out with likelihood zero if we're too close to the edge of the matrix
	else  // emission at node n, but no Score_profile
	  node_scores[n] = 0;
      if (seq[n])
	  ++i;
    }

  // delegate pruning to Transducer_peeler
  Transducer_tree_emission tree_emission (node_scores, trans_sc->branch_trans_states[state], inserter[state], deleters[state], matchers[state]);
  tree_emission.state_name = trans_sc->sexpr_state_name (state);
  tree_emission.state_name << " at coords (" << seq_coords << ')';
  peeler.prune (tree_emission);
  if (want_counts)
    peeler.peel (tree_emission, weight);
  const Score pruning_sc = peeler.pruning_sc;

  // do the temperature correction
  return (Score) (((double) pruning_sc) / kT);
}


// Transducer_DP_matrix
double Transducer_DP_matrix::cells() const
{
  double n = states();
  for_const_contents (vector<int>, dim, l)
    n *= *l;
  return n;
}

void Transducer_DP_matrix::alloc()
{
  Transducer_DP_base::alloc();

  dim.clear();
  for_const_contents (vector<int>, observed_seqs, i)
    dim.push_back (seq[*i]->size() + 1);
  seqzero = vector<int> (dim.size(), (int) 0);

  if (CTAGGING(5,TRANSDUCER TRANSDUCER_ALLOC))
    {
      CL << "Allocating " << states();
      for_const_contents (vector<int>, dim, l)
	CL << " * " << *l;
      CL << " = " << cells() << " cells for transducer dynamic programming matrix\n";
    }

  // efficiently allocate cell_sc
  cell_sc.clear();
  cell_sc = vector<Cell_multi_array> (states(), Cell_multi_array());
  for (int s = 0; s < states(); ++s)
    {
      Cell_multi_array tmp_cell_sc (dim, -InfinityScore);
      cell_sc[s].swap (tmp_cell_sc);
    }

  // initialize end_sc
  end_sc = -InfinityScore;

  // set up internal vars
  n_types = 1 << trans_sc->nodes();
  seq_delta = array2d<int> (states(), sequences(), (int) 0);
  for (int s = 0; s < states(); ++s)
    for (int d = 0; d < sequences(); ++d)
      if (type_node_emit (trans_sc->state_type[s], observed_seqs[d]))
	seq_delta(s,d) = 1;

  if (CTAGGING(-3,TRANSDUCER))
    CL << "Transducer_DP_matrix::observed_seqs = (" << observed_seqs
       << ")\nTransducer_DP_matrix::trans_sc->state_type = (" << trans_sc->state_type
       << ")\nTransducer_DP_matrix::seq_delta:\n" << seq_delta;
}

void Transducer_DP_matrix::show (ostream& o)
{
  CL.save_logfile_state();
  sstring out;

  out << "Dimensions: (" << dim << ")\n";
  out << "Cell\tState:Score:Emit  (scores of -infinity are omitted)\n";
  for (Counter seq_coords (seqzero, seqzero, dim);
       seq_coords < seq_coords.end();
       ++seq_coords)
    {
      out << seq_coords << '\t';
      for (int state = 0; state < states(); ++state)
	{
	  const Score sc = cell_sc[state][seq_coords];
	  if (sc > -InfinityScore)
	    out << trans_sc->get_state_name (state)
		<< ':' << sc
		<< ':' << calc_cell_emit (state, seq_coords)
		<< ' ';
	}
      out << '\n';
    }
  out << "End\t" << end_sc << "\n";

  CL.restore_logfile_state();
  o << out;
}

// Transducer_forward_matrix
Transducer_forward_matrix::Transducer_forward_matrix()
{ }

void Transducer_forward_matrix::fill()
{
  CTAG(4,TRANSDUCER TRANSDUCER_FILL) << "Filling " << cells() << " cells in transducer forward dynamic programming matrix\n";

  // do outgoing transitions from start state
  vector<int> seq_coords (seqzero);
  for (int state = 0; state < states(); ++state)
    {
      // check to see if destination cell is outside the matrix
      bool accessible = true;
      for (int d = 0; d < sequences(); ++d)
	if ((seq_coords[d] = seq_delta (state, d)) >= dim[d])
	  {
	    accessible = false;
	    break;
	  }
      // if not outside, then add transition score, and store
      if (accessible)
	cell_sc[state][seq_coords] = trans_sc->start[state];
    }

  // set up sequence coord vectors
  vector<int> src_seq_coords (seqzero);
  // loop over all destination sequence coords
  for (Counter dest_seq_coords (seqzero, seqzero, dim);
       dest_seq_coords < dest_seq_coords.end();
       ++dest_seq_coords)
    {
      // loop over all destination states (they're all emit states, so no sorting is necessary)
      for (int dest_state = 0; dest_state < states(); ++dest_state)
	{
	  // find source cell; check if it's outside the matrix
	  bool accessible = true;
	  for (int d = 0; d < sequences(); ++d)
	    if ((src_seq_coords[d] = dest_seq_coords[d] - seq_delta(dest_state,d)) < 0)
	      {
		accessible = false;
		break;
	      }
	  if (accessible)
	    {
	      // get address of current cell; calculate emit score
	      Score& dest_sc = cell_sc[dest_state][dest_seq_coords];

	      // loop over incoming source states, updating cell score
	      for_const_contents (vector<int>, incoming[dest_state], src_state)
		ScorePSumAcc (dest_sc, ScorePMul (cell_sc[*src_state][src_seq_coords],
						  trans_sc->transition (*src_state, dest_state)));

	      // add emit score
	      const Score emit_sc = calc_cell_emit (dest_state, dest_seq_coords);
	      ScorePMulAcc (dest_sc, emit_sc);
	    }
	}
    }

  // check for all zero-length sequences
  bool all_seqs_empty = true;
  for (int d = 0; d < sequences(); ++d)
    if ((seq_coords[d] = dim[d] - 1) > 0)
      all_seqs_empty = false;

  // do end state
  if (all_seqs_empty)
    ScorePSumAcc (end_sc, trans_sc->start_to_end());
  else
    for (int state = 0; state < states(); ++state)
      ScorePSumAcc (end_sc, ScorePMul (cell_sc[state][seq_coords], trans_sc->end[state]));

  // log
  if (CTAGGING(2,TRANSDUCER TRANSDUCER_DP_MATRIX))
    {
      CL << "Transducer forward DP matrix:\n";
      show (CL);
    }
}

vector<int> Transducer_forward_matrix::sample_traceback()
{
  vector<int> trace;

  if (end_sc <= -InfinityScore)
    CLOGERR << "Forward likelihood is zero; skipping traceback\n";
  else
    {
      vector<int> dest_seq_coords (dim);
      for (int d = 0; d < sequences(); ++d)
	--dest_seq_coords[d];
      int dest_state = Grammar_state_enum::End;

      vector<vector<int> > trace_coords;  // for traceback-out-of-bounds debugging only
      do
	{
	  trace.push_back (dest_state);
	  trace_coords.push_back (dest_seq_coords);

	  vector<int> src_seq_coords (dest_seq_coords);
	  int src_state = dest_state;

	  Score emit_sc = 0;
	  if (dest_state != Grammar_state_enum::End)
	    {
	      emit_sc = calc_cell_emit (dest_state, dest_seq_coords);
	      for (int d = 0; d < sequences(); ++d)
		if ((src_seq_coords[d] -= seq_delta(dest_state,d)) < 0)
		  {
		    CLOGERR << "Transducer_forward_matrix: traceback out of bounds\n";
		    CL << "Transducer:\n";
		    trans_sc->show (CL);
		    show (CL);
		    CL << "Traceback path:\n";
		    for_const_contents (vector<vector<int> >, trace_coords, tc)
		      CL << *tc << "\n";
		    THROWEXPR ("Traceback went out of bounds at " << src_seq_coords);
		  }
	    }

	  bool start_ok = true;
	  for (int d = 0; start_ok && d < sequences(); ++d)
	    if (src_seq_coords[d] > 0)
	      start_ok = false;

	  vector<Score> sc;
	  for (src_state = 0; src_state < states(); ++src_state)
	    sc.push_back (ScorePMul3 (cell_sc[src_state][src_seq_coords],
				      trans_sc->transition (src_state, dest_state),
				      emit_sc));

	  if (start_ok)
	    sc.push_back (ScorePMul (trans_sc->transition (Grammar_state_enum::Start, dest_state),
				     emit_sc));

	  Score total_sc = -InfinityScore;
	  for_const_contents (vector<Score>, sc, s)
	    ScorePSumAcc (total_sc, *s);

	  Prob p = Rnd::prob();
	  for (src_state = 0; src_state < states(); ++src_state)
	    if ((p -= Score2Prob (ScorePMul (sc[src_state], -total_sc))) <= 0)
	      break;

	  if (start_ok && src_state == states())
	    src_state = Grammar_state_enum::Start;

	  dest_seq_coords = src_seq_coords;
	  dest_state = src_state;
	}
      while (dest_state != Grammar_state_enum::Start);
      trace.push_back (dest_state);

      reverse (trace.begin(), trace.end());
    }

  return trace;
}

Transducer_backward_matrix::Transducer_backward_matrix (bool wc)
  : Transducer_DP_matrix(),
    fwd_initialized(false)
{
  want_counts = wc;
}

void Transducer_backward_matrix::init_fwd (Transducer_forward_matrix& existing_fwd)
{
  fwd.end_sc = existing_fwd.end_sc;
  fwd.cell_sc.swap (existing_fwd.cell_sc);

  fwd_initialized = true;
}

void Transducer_backward_matrix::alloc()
{
  if (!fwd_initialized)
    {
      fwd.trans_sc = trans_sc;
      fwd.seq = seq;
      fwd.alloc();
    }

  Transducer_DP_matrix::alloc();

  n_compressed_types = 1 << sequences();
  cell_match_by_type = vector<multi_array<Prob> > (n_compressed_types, multi_array<Prob> (dim, 0.));

  set<State_type> valid_uncompressed_type_set;
  for_const_contents (vector<State_type>, trans_sc->state_type, t)
    valid_uncompressed_type_set.insert (*t);

  compressed_state_type = vector<unsigned int> (n_types);
  for_const_contents (set<State_type>, valid_uncompressed_type_set, t)
    {
      unsigned int c = 0;
      for (int i = 0; i < (int) sequences(); ++i)
	if (type_node_emit (*t, observed_seqs[i]))
	  c |= 1 << i;
      compressed_state_type[*t] = c;

      cell_match_by_type[c] = multi_array<Prob> (dim, 0.);
    }
}

void Transducer_backward_matrix::fill()
{
  // fill the forward matrix
  if (!fwd_initialized)
    fwd.fill();

  if (fwd.end_sc <= -InfinityScore)
    CLOGERR << "Warning: forward score is -infinity. Skipping backward matrix fill.";
  else
    {
      // log
      CTAG(4,TRANSDUCER TRANSDUCER_FILL) << "Filling " << cells() << " cells in transducer backward dynamic programming matrix\n";

      // handle start-to-end transition
      if (dim == vector<int> (dim.size(), (int) 1))
	ehmm_transition_counts.start_to_end() += Score2Prob (ScorePMul (trans_sc->start_to_end(),
									-fwd.end_sc));

      // do incoming transitions from end state
      vector<int> dest_seq_coords (dim), src_seq_coords (dim);
      for_contents (vector<int>, dest_seq_coords, c)
	--*c;
      for (int state = 0; state < states(); ++state)
	{
	  cell_sc[state][dest_seq_coords] = trans_sc->end[state];
	  if (want_counts)
	    {
	      const Prob post_transition_prob = Score2Prob (ScorePMul3 (fwd.cell_sc[state][dest_seq_coords],
									trans_sc->end[state],
									-fwd.end_sc));
	      ehmm_transition_counts.end[state] += post_transition_prob;
	    }
	}

      // loop over all source sequence coords
      for (Counter cell_counter (seqzero, seqzero, dim);
	   cell_counter < cell_counter.end();
	   ++cell_counter)
	{
	  // set up src_seq_coords
	  for (int i = 0; i < (int) src_seq_coords.size(); ++i)
	    src_seq_coords[i] = dim[i] - cell_counter[i] - 1;

	  // loop over all source states (they're all emit states, so no sorting is necessary)
	  for (int src_state = 0; src_state < states(); ++src_state)
	    {
	      // get source cell score
	      Score& src_sc = cell_sc[src_state][src_seq_coords];

	      // subtract final score from forward score
	      const Score fwd_minus_final_sc = ScorePMul (fwd.cell_sc[src_state][src_seq_coords], -fwd.end_sc);

	      // loop over outgoing dest states
	      for_const_contents (vector<int>, outgoing[src_state], dest_state)
		{
		  // find dest cell; check if it's outside the matrix
		  bool accessible = true;
		  for (int d = 0; d < sequences(); ++d)
		    if ((dest_seq_coords[d] = src_seq_coords[d] + seq_delta(*dest_state,d)) >= dim[d])
		      {
			accessible = false;
			break;
		      }
		  if (accessible)
		    {
		      const Score src_dest_sc = ScorePMul (cell_sc[*dest_state][dest_seq_coords],
							   trans_sc->transition (src_state, *dest_state));
		      ScorePSumAcc (src_sc, src_dest_sc);

		      // add post prob for this transition to accumulated counts
		      if (want_counts)
			{
			  const Prob post_transition_prob = Score2Prob (ScorePMul (src_dest_sc, fwd_minus_final_sc));
			  ehmm_transition_counts.transition (src_state, *dest_state) += post_transition_prob;
			}
		    }
		}

	      // get post prob of being in this state
	      const Prob post_state_prob = Score2Prob (ScorePMul (src_sc, fwd_minus_final_sc));

	      // update cell_match_by_type
	      cell_match_by_type[compressed_state_type[trans_sc->state_type[src_state]]][src_seq_coords] += post_state_prob;

	      // get emit score & accumulate counts for emission
	      const Score emit_sc = calc_cell_emit (src_state, src_seq_coords, post_state_prob);

	      // add emit_sc
	      ScorePMulAcc (src_sc, emit_sc);
	    }
	}

      // do counts from start state
      // also calculate final score for consistency check during debugging
      vector<int> seq_coords (seqzero);
      for (int state = 0; state < states(); ++state)
	{
	  // check to see if destination cell is outside the matrix
	  bool accessible = true;
	  for (int d = 0; d < sequences(); ++d)
	    if ((seq_coords[d] = seq_delta (state, d)) >= dim[d])
	      {
		accessible = false;
		break;
	      }
	  // if not outside, then add transition score, and store
	  if (accessible)
	    {
	      const Score start_sc = ScorePMul (cell_sc[state][seq_coords], trans_sc->start[state]);
	      ScorePSumAcc (end_sc, start_sc);
	      if (want_counts)
		{
		  const Prob post_transition_prob = Score2Prob (ScorePMul (start_sc, -fwd.end_sc));
		  ehmm_transition_counts.start[state] += post_transition_prob;
		}
	    }
	}

      // log
      if (CTAGGING(2,TRANSDUCER TRANSDUCER_DP_MATRIX))
	{
	  CL << "Transducer backward DP matrix:\n";
	  show (CL);
	}

      if (CTAGGING(2,TRANSDUCER TRANSDUCER_OPTACC))
	{
	  CL << "Transducer optimal accuracy postprob matrix:\n";

	  CL << "Dimensions: (" << dim << ")\n";
	  CL << "Cell\tType:Postprob\n";
	  for (Counter seq_coords (seqzero, seqzero, dim);
	       seq_coords < seq_coords.end();
	       ++seq_coords)
	    {
	      CL << seq_coords << '\t';
	      for (unsigned int c = 1; c < n_compressed_types; ++c)
		if (cell_match_by_type[c].dim().size())
		  CL << c << ':' << cell_match_by_type[c][seq_coords] << ' ';
	      CL << '\n';
	    }
	}
    }
}

void Transducer_backward_matrix::init_sumpairs_reward()
{
  reward = vector<double> (n_compressed_types, 0.);
  for (unsigned int c = 1; c < n_compressed_types; ++c)
    {
      unsigned int tmp = c;
      int seqs = 0;
      while (tmp)
	{
	  if (tmp & 1)
	    ++seqs;
	  tmp = tmp >> 1;
	}
      reward[c] = seqs * (seqs - 1) / 2;
      if (CTAGGING(2,TRANSDUCER TRANSDUCER_OPTACC))
	CL << "Reward for state type " << c << " (" << seqs << " sequences) is " << reward[c] << "\n";
    }
}

void Transducer_backward_matrix::inc_var_counts (const EHMM_transducer_funcs& ehmm_funcs,
						 PCounts& var_counts,
						 const PScores& var_scores,
						 const Eliminated_EHMM_transducer_scores* elim_sc,
						 Prob weight)
{
  if (!want_counts)
    THROWEXPR ("Get counts. What counts? I didn't collect any counts. You didn't want any counts!");

  if (elim_sc)
    {
      if (trans_sc != (EHMM_transducer_scores*) elim_sc)
	THROWEXPR ("In Transducer_backward_matrix:\n 'elim_ehmm_sc' needs to be the same object as my 'trans_sc', kid.");

      Transducer_counts corrected_ehmm_transition_counts = elim_sc->count_eliminated (ehmm_transition_counts);
      Transducer_methods::inc_transition_counts (corrected_ehmm_transition_counts, ehmm_funcs, var_counts, var_scores, weight);
    }
  else
    Transducer_methods::inc_transition_counts (ehmm_transition_counts, ehmm_funcs, var_counts, var_scores, weight);

  for (int node = 0; node <= trans_sc->max_node_index(); ++node)
    {
      CTAG(3,TRANSDUCER_COUNTS) << "Calling inc_emit_label_counts at tape #" << node << " (" << trans_sc->get_tape_name(node) << ")\n";
      Transducer_methods::inc_emit_label_counts (branch_trans_counts[node], ehmm_funcs.branch_transducer[node], var_counts, var_scores, weight);
    }
}


void Transducer_backward_matrix::alloc_optacc()
{
  optacc = multi_array<double> (dim);
}

void Transducer_backward_matrix::fill_optacc()
{
  if (fwd.end_sc <= -InfinityScore)
    CLOGERR << "Warning: forward score is -infinity. Skipping optimal accuracy matrix fill.";
  else
    {
      CTAG(4,TRANSDUCER TRANSDUCER_FILL) << "Filling " << cells() << " cells in transducer optimal accuracy dynamic programming matrix\n";

      // set up sequence coord vectors
      vector<int> src_seq_coords (seqzero);
      // loop over all destination sequence coords
      for (Counter dest_seq_coords (seqzero, seqzero, dim);
	   dest_seq_coords < dest_seq_coords.end();
	   ++dest_seq_coords)
	{
	  // get address of cell
	  double& cell_score = optacc[dest_seq_coords];

	  // loop over compressed state types (each compressed type is an incoming direction)
	  for (unsigned int c = 1; c < n_compressed_types; ++c)
	    {
	      // find source cell, check it's in bounds
	      bool accessible = true;
	      for (int i = 0; i < (int) sequences(); ++i)
		if ((src_seq_coords[i] = dest_seq_coords[i] - (c & (1<<i) ? 1 : 0)) < 0)
		  {
		    accessible = false;
		    break;
		  }
	      if (accessible)
		cell_score = max (cell_score, optacc[src_seq_coords] + reward[c] * cell_match_by_type[c][dest_seq_coords]);
	    }
	}

      // log
      if (CTAGGING(1,TRANSDUCER_OPTACC))
	show_optacc (CL);
    }
}

Alignment_path Transducer_backward_matrix::optacc_traceback()
{
  Alignment_path path (sequences(), 0);

  if (fwd.end_sc <= -InfinityScore)
    CLOGERR << "Warning: forward score is -infinity. Skipping optimal accuracy traceback.\n";
  else
    {

      vector<int> src_seq_coords (seqzero), dest_seq_coords (dim);
      for (int d = 0; d < sequences(); ++d)
	--dest_seq_coords[d];

      while (1)
	{
	  bool at_zero = true;
	  for (int d = 0; d < sequences(); ++d)
	    if (dest_seq_coords[d])
	      {
		at_zero = false;
		break;
	      }
	  if (at_zero)
	    break;

	  unsigned int best_dir = 0;
	  double best_score = 0;
	  for (unsigned int c = 1; c < n_compressed_types; ++c)
	    {
	      // find source cell, check it's in bounds
	      bool accessible = true;
	      for (int i = 0; i < (int) sequences(); ++i)
		if ((src_seq_coords[i] = dest_seq_coords[i] - (c & (1<<i) ? 1 : 0)) < 0)
		  {
		    accessible = false;
		    break;
		  }
	      if (accessible)
		{
		  const double trial_score = optacc[src_seq_coords] + reward[c] * cell_match_by_type[c][dest_seq_coords];
		  if (trial_score > best_score || best_dir == 0)
		    {
		      best_score = trial_score;
		      best_dir = c;
		    }
		}
	    }

	  if (best_dir == 0)
	    THROWEXPR ("Optimal-accuracy traceback failed at (" << dest_seq_coords << ")");

	  vector<bool> column (sequences(), false);
	  for (int r = 0; r < sequences(); ++r)
	    {
	      const int delta = best_dir & (1 << r) ? 1 : 0;
	      column[r] = delta;
	      dest_seq_coords[r] -= delta;
	    }
	  path.append_column (column);
	}

      for (int r = 0; r < path.rows(); ++r)
	reverse (path.row(r).begin(), path.row(r).end());
    }

  return path;
}

void Transducer_backward_matrix::show_optacc (ostream& out)
{
  CL.save_logfile_state();
  sstring s;

  out << "Optimal accuracy matrix dimensions: (" << dim << ")\n";
  out << "Cell\tScore\n";
  for (Counter seq_coords (seqzero, seqzero, dim);
       seq_coords < seq_coords.end();
       ++seq_coords)
    s << seq_coords << '\t' << optacc[seq_coords] << '\n';

  CL.restore_logfile_state();
  out << s;
}

Transducer_Viterbi_matrix::Transducer_Viterbi_matrix()
{ }

void Transducer_Viterbi_matrix::fill()
{
  CTAG(4,TRANSDUCER TRANSDUCER_FILL) << "Filling " << cells() << " cells in transducer Viterbi dynamic programming matrix\n";

  // do outgoing transitions from start state
  vector<int> seq_coords (seqzero);
  for (int state = 0; state < states(); ++state)
    {
      // check to see if destination cell is outside the matrix
      bool accessible = true;
      for (int d = 0; d < sequences(); ++d)
	if ((seq_coords[d] = seq_delta (state, d)) >= dim[d])
	  {
	    accessible = false;
	    break;
	  }
      // if not outside, then add transition score, and store
      if (accessible)
	cell_sc[state][seq_coords] = trans_sc->start[state];
    }

  // set up sequence coord vectors
  vector<int> src_seq_coords (seqzero);

  // sort destination states
  vector<int> null_states, emit_states;
  trans_sc->get_emit_and_null_states (peeler.subtrees.leaves, emit_states, null_states);
  Transition_methods::topological_sort (*trans_sc, null_states, true);
  vector<int> sorted_states;
  sorted_states.insert (sorted_states.end(), emit_states.begin(), emit_states.end());
  sorted_states.insert (sorted_states.end(), null_states.begin(), null_states.end());

  if (CTAGGING(1,TRANSDUCER TRANSDUCER_VITERBI))
    CL << "Viterbi fill order: " << sorted_states << '\n';

  // loop over all destination sequence coords
  for (Counter dest_seq_coords (seqzero, seqzero, dim);
       dest_seq_coords < dest_seq_coords.end();
       ++dest_seq_coords)
    {
      // loop over all destination states (they're all emit states, so no sorting is necessary)
      for_const_contents (vector<int>, sorted_states, ds_ptr)
	{
	  const int dest_state = *ds_ptr;

	  // find source cell; check if it's outside the matrix
	  bool accessible = true;
	  for (int d = 0; d < sequences(); ++d)
	    if ((src_seq_coords[d] = dest_seq_coords[d] - seq_delta(dest_state,d)) < 0)
	      {
		accessible = false;
		break;
	      }
	  if (accessible)
	    {
	      // get address of current cell; calculate emit score
	      Score& dest_sc = cell_sc[dest_state][dest_seq_coords];

	      // loop over incoming source states, updating cell score
	      for_const_contents (vector<int>, incoming[dest_state], src_state)
		dest_sc = max (dest_sc, ScorePMul (cell_sc[*src_state][src_seq_coords],
						   trans_sc->transition (*src_state, dest_state)));

	      // add emit score
	      const Score emit_sc = calc_cell_emit (dest_state, dest_seq_coords);
	      ScorePMulAcc (dest_sc, emit_sc);
	    }
	}
    }

  // check for all zero-length sequences
  bool all_seqs_empty = true;
  for (int d = 0; d < sequences(); ++d)
    if ((seq_coords[d] = dim[d] - 1) > 0)
      all_seqs_empty = false;

  // do end state
  if (all_seqs_empty)
    end_sc = max (end_sc, trans_sc->start_to_end());
  else
    for (int state = 0; state < states(); ++state)
      end_sc = max (end_sc, ScorePMul (cell_sc[state][seq_coords], trans_sc->end[state]));

  // log
  if (CTAGGING(2,TRANSDUCER TRANSDUCER_DP_MATRIX))
    {
      CL << "Transducer Viterbi DP matrix:\n";
      show (CL);
    }
}

vector<int> Transducer_Viterbi_matrix::traceback()
{
  vector<int> trace;

  if (end_sc <= -InfinityScore)
    CLOGERR << "Viterbi likelihood is zero; skipping trace\n";
  else
    {
      vector<int> dest_seq_coords (dim);
      for (int d = 0; d < sequences(); ++d)
	--dest_seq_coords[d];
      int dest_state = Grammar_state_enum::End;

      vector<vector<int> > trace_coords;  // for traceback-out-of-bounds debugging only
      do
	{
	  trace.push_back (dest_state);
	  trace_coords.push_back (dest_seq_coords);

	  vector<int> src_seq_coords (dest_seq_coords);
	  int src_state = dest_state;

	  Score emit_sc = 0;
	  if (dest_state != Grammar_state_enum::End)
	    {
	      emit_sc = calc_cell_emit (dest_state, dest_seq_coords);
	      for (int d = 0; d < sequences(); ++d)
		if ((src_seq_coords[d] -= seq_delta(dest_state,d)) < 0)
		  {
		    CLOGERR << "Transducer_Viterbi_matrix: traceback out of bounds\n";
		    CL << "Transducer:\n";
		    trans_sc->show (CL);
		    show (CL);
		    CL << "Traceback path:\n";
		    for_const_contents (vector<vector<int> >, trace_coords, tc)
		      CL << *tc << "\n";
		    THROWEXPR ("Traceback went out of bounds at " << src_seq_coords);
		  }
	    }

	  bool start_ok = true;
	  for (int d = 0; start_ok && d < sequences(); ++d)
	    if (src_seq_coords[d] > 0)
	      start_ok = false;

	  vector<Score> sc;
	  for (src_state = 0; src_state < states(); ++src_state)
	    sc.push_back (ScorePMul3 (cell_sc[src_state][src_seq_coords],
				      trans_sc->transition (src_state, dest_state),
				      emit_sc));

	  if (start_ok)
	    sc.push_back (ScorePMul (trans_sc->transition (Grammar_state_enum::Start, dest_state),
				     emit_sc));

	  src_state = max_element (sc.begin(), sc.end()) - sc.begin();

	  if (start_ok && src_state == states())
	    src_state = Grammar_state_enum::Start;

	  dest_seq_coords = src_seq_coords;
	  dest_state = src_state;
	}
      while (dest_state != Grammar_state_enum::Start);
      trace.push_back (dest_state);

      reverse (trace.begin(), trace.end());
    }

  return trace;
}
