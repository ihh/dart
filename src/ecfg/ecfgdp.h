#ifndef ECFGDP_INCLUDED
#define ECFGDP_INCLUDED

#include "ecfg/ecfg.h"
#include "ecfg/ecfgenv.h"
#include "ecfg/fastprune.h"

// Stockholm #=GR tag for ancestral sequence reconstructions:
// maximum a posteriori (MAP) residue states for each sequence,
// conditional on MAP Cocke-Younger-Kasami (CYK) parse tree.
#define CYK_MAP_reconstruction_tag "ancrec_CYK_MAP"

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

  // data for hybrid chains
  vector<vector<int> > lineage_chain_index;  // lineage_chain_index[chain][node] is a vector tracking assignments of branches to component chain indices for chain #n
  vector<vector<const EM_matrix_base*> > lineage_matrix;  // lineage_matrix[n] is a vector tracking assignments of branches to EM_matrix_base* pointers for chain #n

  // constructor
  ECFG_matrix (const ECFG_scores& ecfg, Stockholm& stock, const ECFG_envelope& env);
  // virtual destructor, to keep gcc happy
  virtual ~ECFG_matrix() { }

  // ECFG_subseq_score_calculator method
  Loge state_ll (int state, const Subseq_coords& subseq) const;

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

  // indel pseudo-events: precomputed scores & workspace variables
  vector<int> src_gap_profile, dest_gap_profile;  // workspace for storing gap profiles while calculating effective transition scores
  double total_branch_len;  // total branch length of tree
  vector<Score> link_extend_sc, link_end_sc;  // indexed by [state]
  array2d<Score> ins_event_sc, del_event_sc, match_event_sc;  // indexed by (state, node)

  // constructor
  ECFG_EM_matrix (const ECFG_scores& ecfg, Stockholm& stock, const Aligned_score_profile& asp, const ECFG_envelope& env, bool use_fast_prune = false);

  // EM methods
  bool calc_annot_emit_ll (int subseq_idx, int state_idx, Loge& annot_emit_ll);  // calculates likelihood due to probabilistic/deterministic annotation tracks, places in annot_emit_ll, returns true if annot_emit_ll > -InfinityLoge
  bool fill_up (int subseq_idx, int state_idx, bool condition_on_context = true);  // returns TRUE if fill_up was called on Column_matrix (i.e. kosher emit state)
  void fill_down (ECFG_counts& counts, int subseq_idx, int state_idx, double weight);  // updates stats, returns log-likelihood

  // Beagle helper methods
  // TODO: implement these!
  typedef map<Phylogeny::Node,Score_profile> Partial_map;
  typedef map<Phylogeny::Node,array2d<Prob> > Transition_matrix_map;  // indexed by child node
  void get_partials (int state, Partial_map& with_context, Partial_map& with_wildcards) { }  // both Partial_map's are cleared; with_wildcards is only filled if state has emit context
  void get_branch_transition_matrices_and_root_prior (int state, Transition_matrix_map& branch_transmat, vector<double>& root_prior) { }
  void use_precomputed_phyloemit (Emit_loglike_matrix& phyloemit) { }  // calls use_precomputed() then iterates over all cells calling calc_annot_emit_ll() & adding to emit_loglike

  // methods to fill cells
  // These methods assume that env.get_bif_in() has already been called by the fill routine.
  void fill_sum_state (int subseq_idx, int state_idx);
  void fill_max_state (int subseq_idx, int state_idx);

  // methods for calculating rate-gap transitions
  // method for testing if a gap profile is illegal (i.e. doesn't contain exactly one clique)
  inline bool gap_profile_illegal (const vector<int>& gapped, int root) const;

  // effective transition log-likelihood
  inline Loge effective_trans (int src_state, int dest_state,
			       const Subseq_coords& src_subseq, const Subseq_coords& dest_subseq) const;

  // dest_gap_state event score
  inline Score dest_event_sc (int dest_state, const vector<int>& dest_gap_profile, int dest_root) const;

  // method for accumulating transition counts
  inline void add_trans_counts (int src_state, int dest_state, const Subseq_coords& src_subseq, const Subseq_coords& dest_subseq,
				double weight, ECFG_counts& counts) const;

  // method for accumulating indel counts
  inline void add_event_counts (int dest_state, const vector<int>& dest_gap_profile, int dest_root,
				double weight, ECFG_counts& counts) const;

  // ancestral state reconstruction
  void reconstruct_MAP (Stockholm& stock, const ECFG_cell_score_map& annot, const char* ancrec_tag = CYK_MAP_reconstruction_tag, bool annotate_postprobs = false);

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

Loge ECFG_EM_matrix::effective_trans (int src_state, int dest_state,
				      const Subseq_coords& src_subseq, const Subseq_coords& dest_subseq) const
{
  // get state info
  const ECFG_state_info* src_info = src_state >= 0 ? &ecfg.state_info[src_state] : (ECFG_state_info*) 0;
  const ECFG_state_info* dest_info = dest_state >= 0 ? &ecfg.state_info[dest_state] : (ECFG_state_info*) 0;

  // figure out if src & dest states use pseudo-indel modeling
  const bool src_indels = src_info ? src_info->indels : false;
  const bool dest_indels = dest_info ? dest_info->indels : false;

  // get base transition score
  const Score base_trans_sc = ecfg.transition (src_state, dest_state);

  // dispose of quick & easy case: transitions between non-indel states
  if (!src_indels && !dest_indels)  // src_indels == dest_indels == false
    return Score2Nats (base_trans_sc);

  // get src & dest gap profiles
  int src_root, dest_root;
  vector<int>& src_gap = ((ECFG_EM_matrix*)this)->src_gap_profile;  // cast away const
  vector<int>& dest_gap = ((ECFG_EM_matrix*)this)->dest_gap_profile;  // cast away const

  if (src_indels)
    {
      src_info->get_gap_profile (stock, tree, src_subseq, src_gap, src_root);
      if (gap_profile_illegal (src_gap, src_root))
	return -InfinityLoge;
    }

  if (dest_indels)
    {
      dest_info->get_gap_profile (stock, tree, dest_subseq, dest_gap, dest_root);
      if (gap_profile_illegal (dest_gap, dest_root))
	return -InfinityLoge;
    }

  // handle transitions involving pseudo-indel states
  Score trans_sc = -InfinityScore;
  if (src_indels)
    {
      if (dest_indels)  // src_indels == true, dest_indels == true
	{
	  const Score new_link_sc = ScorePMul (link_end_sc[src_state], base_trans_sc);
	  const Score event_sc = dest_event_sc (dest_state, dest_gap_profile, dest_root);
	  trans_sc = ScorePMul (new_link_sc, event_sc);
	  if (src_state == dest_state && src_gap == dest_gap)
	    {
	      const Score same_link_sc = link_extend_sc[src_state];
	      ScorePSumAcc (trans_sc, same_link_sc);
	    }
	}

      else  // src_indels == true, dest_indels == false
	trans_sc = ScorePMul (link_end_sc[src_state], base_trans_sc);
    }

  else  // src_indels == false, dest_indels == true
    trans_sc = ScorePMul (base_trans_sc, dest_event_sc (dest_state, dest_gap, dest_root));

  // return
  return Score2Nats (trans_sc);
}

void ECFG_EM_matrix::add_trans_counts (int src_state, int dest_state,
				       const Subseq_coords& src_subseq, const Subseq_coords& dest_subseq,
				       double weight, ECFG_counts& counts) const
{
  // get state info
  const ECFG_state_info* src_info = src_state >= 0 ? &ecfg.state_info[src_state] : (ECFG_state_info*) 0;
  const ECFG_state_info* dest_info = dest_state >= 0 ? &ecfg.state_info[dest_state] : (ECFG_state_info*) 0;

  // figure out if src & dest states use pseudo-indel modeling
  const bool src_indels = src_info ? src_info->indels : false;
  const bool dest_indels = dest_info ? dest_info->indels : false;

  // get base transition score
  const Score base_trans_sc = ecfg.transition (src_state, dest_state);

  // dispose of quick & easy case: transitions between non-indel states
  if (!src_indels && !dest_indels)  // src_indels == dest_indels == false
    {
      counts.transition (src_state, dest_state) += weight;
      return;
    }

  // get src & dest gap profiles
  int src_root, dest_root;
  vector<int>& src_gap = ((ECFG_EM_matrix*)this)->src_gap_profile;  // cast away const
  vector<int>& dest_gap = ((ECFG_EM_matrix*)this)->dest_gap_profile;  // cast away const

  if (src_indels)
    {
      src_info->get_gap_profile (stock, tree, src_subseq, src_gap, src_root);
      if (gap_profile_illegal (src_gap, src_root))
	return;
    }

  if (dest_indels)
    {
      dest_info->get_gap_profile (stock, tree, dest_subseq, dest_gap, dest_root);
      if (gap_profile_illegal (dest_gap, dest_root))
	return;
    }

  // print log message showing gap profiles
  if (CTAGGING(2,ECFG_TRANS_COUNTS))
    {
      CL << "Counting transitions from state " << src_state << ", subseq " << src_subseq.start << "+" << src_subseq.len;
      if (src_indels)
	{
	  CL << " (";
	  for (int n = 0; n < tree.nodes(); ++n)
	    CL << (src_gap[n] ? '-' : '*');
	  CL << ")";
	}

      CL << " to state " << dest_state << ", subseq " << dest_subseq.start << "+" << dest_subseq.len;
      if (dest_indels)
	{
	  CL << " (";
	  for (int n = 0; n < tree.nodes(); ++n)
	    CL << (dest_gap[n] ? '-' : '*');
	  CL << ")";
	}
      CL << ", weight " << weight << "\n";
    }

  // handle transitions involving pseudo-indel states
  if (src_indels)
    {
      if (dest_indels)  // src_indels == true, dest_indels == true
	{
	  const Score trans_sc = Nats2Score (effective_trans (src_state, dest_state, src_subseq, dest_subseq));
	  const Score new_link_sc = ScorePMul (link_end_sc[src_state], base_trans_sc);
	  const Score event_sc = dest_event_sc (dest_state, dest_gap, dest_root);

	  const double new_link_weight = weight * Score2Prob (ScorePMul3 (new_link_sc, event_sc, -trans_sc));
	  counts.link_end_count[src_state] += new_link_weight;
	  counts.transition (src_state, dest_state) += new_link_weight;
	  add_event_counts (dest_state, dest_gap, dest_root, new_link_weight, counts);

	  if (src_state == dest_state && src_gap == dest_gap)
	    {
	      const Score same_link_sc = link_extend_sc[src_state];
	      const double same_link_weight = weight * Score2Prob (ScorePMul (same_link_sc, -trans_sc));
	      counts.link_extend_count[src_state] += same_link_weight;
	    }
	}

      else  // src_indels == true, dest_indels == false
	{
	  counts.link_end_count[src_state] += weight;
	  counts.transition (src_state, dest_state) += weight;
	}
    }

  else  // src_indels == false, dest_indels == true
    {
      counts.transition (src_state, dest_state) += weight;
      add_event_counts (dest_state, dest_gap, dest_root, weight, counts);
    }
}

Score ECFG_EM_matrix::dest_event_sc (int dest_state, const vector<int>& dest_gap, int dest_root) const
{
  Score sc = ins_event_sc (dest_state, dest_root);
  for_branches_post (tree, dest_root, tree.parent[dest_root], b)
    if (!dest_gap[(*b).first])
      ScorePMulAcc (sc, (dest_gap[(*b).second] ? del_event_sc : match_event_sc) (dest_state, (*b).second));
  return sc;
}

void ECFG_EM_matrix::add_event_counts (int dest_state, const vector<int>& dest_gap, int dest_root,
				       double weight, ECFG_counts& counts) const
{
  counts.ins_wait[dest_state] += weight * total_branch_len;
  if (dest_root != tree.root)
    counts.ins_count[dest_state] += weight;

  for_branches_post (tree, dest_root, tree.parent[dest_root], b)
    if (!dest_gap[(*b).first])
      {
	counts.del_wait[dest_state] += weight * (*b).length;
	if (dest_gap[(*b).second])
	  counts.del_count[dest_state] += weight;
      }
}

bool ECFG_EM_matrix::gap_profile_illegal (const vector<int>& gapped, int root) const
{
  if (root < 0)
    return true;

  for_branches_post (tree, root, tree.parent[root], b)
    if (gapped[(*b).first] && !gapped[(*b).second])
      return true;

  return false;
}

#endif /* ECFGDP_INCLUDED */
