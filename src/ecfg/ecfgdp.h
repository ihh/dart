#ifndef ECFGDP_INCLUDED
#define ECFGDP_INCLUDED

#include "ecfg/ecfg.h"
#include "ecfg/ecfgenv.h"
#include "ecfg/fastprune.h"

// Stockholm #=GR tag for ancestral sequence reconstructions:
// maximum a posteriori (MAP) residue states for each sequence,
// conditional on MAP Cocke-Younger-Kasami (CYK) parse tree.
#define CYK_MAP_reconstruction_tag "ancrec_CYK_MAP"

// min postprob for ancrec reporting
#define DEFAULT_MIN_ANCREC_POSTPROB 0.01

// DP matrix base class
// Cell scores are stored in this class, but intermediate emission log-likelihoods (as e.g. computed by libHMSBeagle) are stored in a subclass
struct ECFG_matrix : ECFG_enum, Stream_saver
{
  // data
  const ECFG_scores& ecfg;
  const ECFG_envelope& env;
  Stockholm& stock;
  Stockholm_tree tree;
  array2d<Loge> cell;   // cell(subseq_idx,state)
  Loge final_loglike;
  double bifs_done, last_bifs_done, total_bifs;  // progress tracking

  // ECFG summary data
  vector<int> emit_states, nonemit_states, inside_fill_states, outside_fill_states;
  vector<vector<int> > incoming, outgoing, left_bifurc, right_bifurc;
  vector<vector<Loge> > incoming_ll, outgoing_ll;
  vector<Loge> start_ll, end_ll;

  // data for hybrid chains
  vector<vector<int> > lineage_chain_index;  // lineage_chain_index[chain][node] is a vector tracking assignments of branches to component chain indices for chain #n
  vector<vector<const EM_matrix_base*> > lineage_matrix;  // lineage_matrix[n] is a vector tracking assignments of branches to EM_matrix_base* pointers for chain #n

  // constructor
  ECFG_matrix (const ECFG_scores& ecfg, Stockholm& stock, const ECFG_envelope& env);
  // virtual destructor, to keep gcc happy
  virtual ~ECFG_matrix() { }

  // ECFG_subseq_score_calculator method
  Loge state_ll (int state, const Subseq_coords& subseq) const;

  // transition log-likelihood
  inline Loge trans_ll (int src_state, int dest_state) const;

  // helpers for calculating annotation log-likelihoods
  static inline bool wildcard_match (const sstring& a, const sstring& b);
  static inline void get_annot_string (const ECFG_state_info& info, const char* aa, const Subseq_coords& subseq, sstring& aa_str, bool& has_wildcards);

  // show methods
  void show (ostream& out) const;
  virtual void show_emit (int subseq, ostream& out) const { }
  void show_progress (const char* algorithm_name, int subseq_idx, bool outside);
};

// DP matrix with EM-related data
// This class contains the emit log-likelihood table that could be optimized using libHMSBeagle
struct ECFG_EM_matrix : ECFG_matrix
{
  // alignment info
  const Aligned_score_profile& asp;
  vector<const char*> align_annot;   // accessed as align_annot[featureIndex]
  typedef map<sstring,Loge> String_loglike_dist;
  vector<vector<String_loglike_dist> > state_annot;  // accessed as state_annot[stateIndex][featureIndex][columnAnnotationString]

  // emit log-likelihoods for each (subseq,state) pair.
  // filling this data structure is expensive, but the structure can be re-used between CYK, Inside & Outside algorithms.
  // in theory the algorithm to fill the structure can be parallelized....
  // Filling this array2d could be optimized using libHMSBeagle (NB there is additional code involved too, due to probabilistic annotation tracks and context-dependent emissions)
  typedef array2d<Loge> Emit_loglike_matrix;
  Emit_loglike_matrix emit_loglike;  // accessed as emit_loglike(subseq_idx,state)
  void use_precomputed (Emit_loglike_matrix& emit_loglike);  // swaps in emit_loglike, clears fill_up_flag

  // workspace
  vector<Column_matrix> colmat;   // "scratch" Column_matrix's, indexed by state
  bool use_fast_prune;  // set this to use fast pruning algorithm
  vector<Fast_prune> fast_prune;  // faster implementation of the pruning algorithm; only used if use_fast_prune is TRUE
  bool fill_up_flag;  // clear this flag if emit_loglike() is precomputed
  vector<Subseq::Bifurc_in> bif_in;

  // constructor
  ECFG_EM_matrix (const ECFG_scores& ecfg, Stockholm& stock, const Aligned_score_profile& asp, const ECFG_envelope& env, bool use_fast_prune = false);

  // EM methods
  inline bool calc_annot_emit_ll (int subseq_idx, int state_idx, Loge& annot_emit_ll);  // calculates likelihood due to probabilistic/deterministic annotation tracks, places in annot_emit_ll, returns true if annot_emit_ll > -InfinityLoge
  inline bool fill_up (int subseq_idx, int state_idx, bool condition_on_context = true);  // returns TRUE if fill_up was called on Column_matrix (i.e. kosher emit state)
  inline void fill_down (ECFG_counts& counts, int subseq_idx, int state_idx, double weight);  // updates stats, returns log-likelihood

  // Beagle helper methods
  typedef map<Phylogeny::Node,vector<double> > Partial_map;   // indexing: Partial_map[node][subseq*STATES + state]
  typedef map<Phylogeny::Node,vector<double> > Transition_matrix_map;  // indexing: Transition_matrix_map[childNode][i*STATES + j]
  void get_partials (int state, Partial_map& with_context, Partial_map& with_wildcards);  // both Partial_map's are cleared; with_wildcards is only filled if state has emit context
  void get_branch_transition_matrices_and_root_prior (int state, Transition_matrix_map& branch_transmat, vector<Prob>& root_prior);
  void use_precomputed_phyloemit (Emit_loglike_matrix& phyloemit);  // calls use_precomputed() then iterates over all cells calling calc_annot_emit_ll() & adding to emit_loglike

  // main Beagle method called by users of this class
  void compute_phylo_likelihoods_with_beagle();  // call before fill()

  // methods to fill cells
  // These methods assume that env.get_bif_in() has already been called by the fill routine.
  inline void fill_sum_state (int subseq_idx, int state_idx);
  inline void fill_max_state (int subseq_idx, int state_idx);

  // method for accumulating transition counts
  inline void add_trans_counts (int src_state, int dest_state, double weight, ECFG_counts& counts) const;

  // ancestral state reconstruction
  void reconstruct_MAP (Stockholm& stock, const ECFG_cell_score_map& annot, const char* ancrec_tag = CYK_MAP_reconstruction_tag, bool annotate_map = true, bool annotate_postprobs = false, Prob min_reported_postprob = DEFAULT_MIN_ANCREC_POSTPROB);

  // display
  void show_emit (int subseq, ostream& out) const;
  void show_coords (ostream& out, int subseq_idx, int state_idx) const;
};

// ECFG trainer
struct ECFG_trainer
{
  // data
  ECFG_scores& ecfg;
  const vector<Stockholm*>& stock_db;
  const vector<Aligned_score_profile>& asp_vec;
  ECFG_counts counts;
  Loge loglike;  // log-likelihood of last iteration
  Loge best_loglike;  // best log-likelihood during training
  int max_subseq_len;
  ostream* training_log;

  // EM control
  double em_min_inc;
  int em_max_iter;

  // pseudocounts
  double pseud_init;
  double pseud_mutate;
  double pseud_wait;

  // constructor
  ECFG_trainer (ECFG_scores& ecfg, const vector<Stockholm*>& stock_db, const vector<Aligned_score_profile>& asp_vec,
		int max_subseq_len, ostream* training_log = 0)
    : ecfg (ecfg), stock_db (stock_db), asp_vec (asp_vec), counts (ecfg), max_subseq_len (max_subseq_len), training_log (training_log),
      em_min_inc (.001), em_max_iter (-1),
      pseud_init(0.), pseud_mutate(0.), pseud_wait(0.)
  { }

  // virtual destructor
  virtual ~ECFG_trainer() { }

  // training methods
  void do_dp (bool do_outside);
  void iterate_quick_EM (int forgive);

  // virtual hook for saving grammars during training
  virtual void save_grammar() { }
};

// Inside matrix
// In this matrix, cell(SUBSEQ,STATE) is the sum-of-scores for all parse subtrees rooted in STATE and emitting SUBSEQ
struct ECFG_inside_matrix : ECFG_EM_matrix, ECFG_inside_calculator
{
  // constructor
  ECFG_inside_matrix (const ECFG_scores& ecfg, Stockholm& stock, const Aligned_score_profile& asp, const ECFG_envelope& env, bool use_fast_prune = false);
  // build method
  void fill();
  // virtual ECFG_inside_calculator method
  Loge state_inside_ll (int state, const Subseq_coords& subseq) const;
};

// Outside matrix
// In the Outside matrix, cell(STATE,SUBSEQ) is the sum-of-scores for all parse sub-trees ending in STATE,
// where the subtree from STATE upwards has emitted everything outside SUBSEQ
// (so the emission from STATE is also outside SUBSEQ).
struct ECFG_outside_matrix : ECFG_matrix
{
  // data
  ECFG_inside_matrix& inside;
  ECFG_counts* counts;
  bool want_substitution_counts;  // if false, EM down-fill won't be called
  // constructor
  ECFG_outside_matrix (ECFG_inside_matrix& inside, ECFG_counts* counts);
  // build method
  void fill();
};

// Inside-outside matrix
struct ECFG_inside_outside_matrix : Grammar_state_enum, ECFG_posterior_probability_calculator
{
  // data
  ECFG_inside_matrix inside;
  ECFG_outside_matrix outside;

  // constructor
  ECFG_inside_outside_matrix (const ECFG_scores& ecfg, Stockholm& stock, const Aligned_score_profile& asp, const ECFG_envelope& env, ECFG_counts* counts = 0);

  // fill method
  void fill();

  // Posterior log-probability calculation methods.
  // These give the posterior log-probability that a parse tree contains a subtree whose root is labeled DEST,
  // where the subtree from DEST downwards has emitted SUBSEQ (including the emission from DEST).
  // (In other words, the subseq-indexing convention here follows Inside/CYK, not Outside/KYC.)
  // In the case of post_transition_sc(), it's also required that the parent of the node labeled DEST must be labeled SRC.
  Loge post_transition_ll (int src_state, int dest_state, int subseq_idx) const;
  Loge post_state_ll (int dest_state, int subseq_idx) const;

  Loge post_transition_ll (int src_state, int dest_state, const Subseq_coords& subseq) const;
  Loge post_state_ll (int dest_state, const Subseq_coords& subseq) const;

  // Annotation methods
  void annotate (Stockholm& stock, GFF_list& gff_list, const sstring& seqname, const ECFG_cell_score_map& annot, const sstring& logpostprob_tag) const;
  void annotate_all_post_state_ll (GFF_list& gff_list, const sstring& seqname, const ECFG_cell_score_map& annot, const sstring& annot_tag) const;
  void annotate_post_prob (GFF_list& gff_list, const sstring& seqname, int state, const Subseq_coords& subseq) const;
  void annotate_hidden_classes (Stockholm& stock, const ECFG_cell_score_map& annot);
};

// CYK matrix
// In this matrix, cell(SUBSEQ,STATE) is the max score for any parse subtree rooted in STATE and emitting SUBSEQ
struct ECFG_CYK_matrix : ECFG_EM_matrix
{
  // data
  int final_state;  // state of final cell, for traceback

  // constructor
  ECFG_CYK_matrix (const ECFG_scores& ecfg, Stockholm& stock, const Aligned_score_profile& asp, const ECFG_envelope& env, bool try_fast_prune = true);

  // build method
  void fill();

  // Traceback/annotation method
  ECFG_cell_score_map traceback();  // automatically annotates Stockholm alignment
};



// inline method defs

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
  const bool use_context = condition_on_context && info.has_context();
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
      // do peeling
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
		NatsPSumAcc (ll, end_ll[source_state_idx]);
	      // loop over outgoing transitions
	      const vector<int>& outgoing_s = outgoing[source_state_idx];
	      const vector<Loge>& outgoing_ll_s = outgoing_ll[source_state_idx];
	      for (int outgoing_idx = 0; outgoing_idx < (int) outgoing_s.size(); ++outgoing_idx)
		{
		  const int d = outgoing_s[outgoing_idx];
		  NatsPSumAcc (ll, NatsPMul (cell (dest_subseq_idx, d), outgoing_ll_s[outgoing_idx]));
		}

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
		ll = max (ll, end_ll[source_state_idx]);
	      // loop over outgoing transitions
	      const vector<int>& outgoing_s = outgoing[source_state_idx];
	      const vector<Loge>& outgoing_ll_s = outgoing_ll[source_state_idx];
	      for (int outgoing_idx = 0; outgoing_idx < (int) outgoing_s.size(); ++outgoing_idx)
		{
		  const int d = outgoing_s[outgoing_idx];
		  ll = max (ll, NatsPMul (cell (dest_subseq_idx, d), outgoing_ll_s[outgoing_idx]));
		}

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

Loge ECFG_matrix::trans_ll (int src_state, int dest_state) const
{
  return Score2Nats (ecfg.transition (src_state, dest_state));
}

void ECFG_EM_matrix::add_trans_counts (int src_state, int dest_state,
				       double weight, ECFG_counts& counts) const
{
  counts.transition (src_state, dest_state) += weight;
  return;
}

// helpers for fill_up
// tests if a==b, allowing both a & b to contain wildcards
bool ECFG_matrix::wildcard_match (const sstring& a, const sstring& b) {
  for (int pos = 0; pos < (int) a.size(); ++pos)
    if (a[pos] != b[pos] && a[pos] != ECFG_annotation_wildcard && b[pos] != ECFG_annotation_wildcard)
      return false;
  return true;
}

// get annotation string
void ECFG_matrix::get_annot_string (const ECFG_state_info& info, const char* aa, const Subseq_coords& subseq, sstring& aa_str, bool& has_wildcards) {
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

#endif /* ECFGDP_INCLUDED */
