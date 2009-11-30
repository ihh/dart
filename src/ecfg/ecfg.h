#ifndef ECFG_INCLUDED
#define ECFG_INCLUDED

#include "hsm/em_matrix.h"
#include "hmm/transmat.h"
#include "scfg/foldenv.h"
#include "tree/tree_alignment.h"
#include "irrev/pfunc_em_matrix.h"
#include "ecfg/ecfgenv.h"
#include "ecfg/fastprune.h"
#include "seq/gff.h"

// Miscellaneous #define's for input/output
#define ECFG_annotation_wildcard Default_annotation_wildcard_char
#define ECFG_default_name        "ECFG"
#define ECFG_default_nonterminal "S"

// GFF tags
#define ECFG_GFF_default_seqname   "Alignment"
#define ECFG_GFF_LogPostProb_tag   "lgPost"      /* log_2 P(Parse Tree Uses Nonterminal | Alignment)  */
#define ECFG_GFF_LogInsideProb_tag "lgInside"    /* log_2 P(Inside Annotation, Inside Alignment | Parse Tree Rooted At Nonterminal)  */

// ECFG enums & typedefs
struct ECFG_enum : Grammar_state_enum
{
  typedef EM_matrix::Column_matrix Column_matrix;
  typedef EM_matrix::Update_statistics Update_statistics;
  enum Update_policy { Rev = 0, Irrev = 1, Rind = 2, Parametric = 3, Hybrid = 4, TotalPolicies = 5 };
};

// Markov chain for ECFG substitution model
struct ECFG_chain : ECFG_enum
{
  // data
  EM_matrix_base* matrix;  // NB lineage-dependent models: this should be null for lineage-dependent models
  int word_len;
  Update_policy type;
  vector<sstring> state;  // badly named; this is actually the set of chain pseudoterminals

  // hidden classes
  sstring class_row;
  vector<sstring> class_labels;
  int classes;

  // PFunc stuff
  bool is_parametric;
  EM_matrix_funcs* matrix_funcs;

  // data for hybrid chains
  string gs_tag;  // the "row" value in the grammar file
  vector<sstring> gs_values;  // the "label" values in the grammar file
  map<sstring,int> gs_tag_value_chain_index;  // gs_tag_value_chain_index[gs_value] = index of ECFG_chain for branches to nodes labeled with "#=GS gs_tag gs_value"
};

// Set of substitution matrices for an Evolutionary CFG
struct ECFG_matrix_set : ECFG_enum
{
  // general data
  const Alphabet& alphabet;  // singlet alphabet (may have 4 symbols, i.e. "ACGU", or 5, i.e. "ACGU-")
  // matrix descriptors
  vector<ECFG_chain> chain;  // the matrices
  // constructors
  ECFG_matrix_set (const ECFG_matrix_set& ems);
  ECFG_matrix_set (const Alphabet& alphabet);
  // destructor -- deletes all EM_matrix's in chain
  ~ECFG_matrix_set();
  // builder method
  int add_matrix (int wordlen, Update_policy = Rev, int n_classes = 1, double tres = DEFAULT_TIMEPOINT_RES);  // returns index of new matrix

  // accessors
  int total_states (int chain_idx) const;
  int observed_states (int chain_idx) const;
  int observed_states_by_word_len (int word_len) const;
};

// state type info for an emit state of an Evolutionary CFG
struct ECFG_state_info : ECFG_enum
{
  // state name
  sstring name;

  // bifurcation parameters
  bool bifurc;  // flag for whether this is a bifurcation state; if FALSE, state is an emit state
  int ldest, rdest;  // left and right destination states (bifurcation states only)

  // emit parameters
  int matrix;  // index of ECFG_chain (emit states only)
  int l_context, l_emit, r_emit, r_context;  // number of columns in (respectively) L-context, L-emission, R-emission, R-context
  vector<int> mul;  // multipliers for context & emitted symbols
  vector<bool> comp;  // flags for complementing context & emitted symbols
  int min_len, max_len;  // valid subsequence lengths; NB min_len must be >= emit_size()
  bool infix;  // if true, then max_len is set equal to ECFG_main::max_subseq_len
  bool prefix, suffix;  // if true, then suppress warnings if max_len > ECFG_main::max_subseq_len (has no other effect other than to suppress warnings)

  // gap parameters
  bool gaps_ok;  // if true, then gaps are allowed (contingent on wild_gaps); if false, no gaps are allowed at all
  bool wild_gaps;  // if true, then partial gaps are treated as wildcards; if false, then partial gaps aren't allowed
  bool indels;  // if true, then state emits a link whose gap profile is determined by indel pseudo-events
  Prob link_extend, link_end;  // probability of extending/ending an existing link
  double ins_rate, del_rate;  // link insertion/deletion rate (NB not a "true" indel model like e.g. TKF91)

  // PFunc's
  bool has_parametric_indels;
  PFunc link_extend_func, link_end_func, ins_rate_func, del_rate_func;

  // By-column annotation; e.g. "<>" for an emitlr state, "012" for a codon state
  typedef map<sstring,PFunc> String_prob_dist;
  typedef map<sstring,String_prob_dist> ECFG_state_annotation;
  ECFG_state_annotation annot;  // indexed by Stockholm "#=GC" tag.
  bool sum_state;  // true if this is scored as a sum state by CYK

  // GFF annotations
  list<GFF> gff;

  // constructor
  ECFG_state_info();
  ECFG_state_info (int l_emit, int r_emit);
  ECFG_state_info (int l_context, int l_emit, int r_emit, int r_context);

  // accessors
  inline int total_size() const { return l_context + l_emit + r_emit + r_context; }
  inline int emit_size() const { return l_emit + r_emit; }
  inline bool has_context() const { return l_context > 0 || r_context > 0; }

  // method to check if subseq is out of range
  inline bool out_of_range (const Subseq_coords& subseq) const
  { return subseq.len < min_len || (max_len >= 0 && subseq.len > max_len ); }

  // Column_matrix initialisation methods: return FALSE if failed
  inline bool init_row (const Alphabet& alphabet, int hidden_classes,
			const Aligned_score_profile& align, int row, const Subseq_coords& subseq,
			vector<Loge>& colmat_row, vector<Prob>* fastprune_row,
			int& gapped, bool use_wildcards_for_emit_columns) const;
  inline bool initialise (const Alphabet& alphabet, int hidden_classes,
			  const Aligned_score_profile& align, const Subseq_coords& subseq,
			  const Stockholm& stock, const Stockholm_tree& tree,
			  Column_matrix& colmat, Fast_prune* fast_prune,
			  bool use_wildcards_for_emit_columns) const;

  // method to get gap profile for a subseq
  inline void get_gap_profile (const Stockholm& stock, const Stockholm_tree& tree, const Subseq_coords& subseq,
			       vector<int>& gapped, int& root) const;

  // method to test if a node is a wildcard due to gaps, or missing row info
  inline bool node_is_wild (const Stockholm& stock, const Stockholm_tree& tree, const Subseq_coords& subseq, int node) const;

  // method to test if a node should be pruned
  inline bool node_is_redundant (const Stockholm& stock, const Stockholm_tree& tree, const Subseq_coords& subseq,
				 const vector<int>& gapped, const Phylogeny::Branch& b) const;

  // method to get column index for a position in a multi-column emission
  inline int column_index (const Subseq_coords& subseq, int col_index) const;

  // method to test if a column index represents an emit column
  inline bool is_emit_column (int col_index) const;

  // output
  void show (ostream& out) const;
};

// some special subclasses of ECFG_state_info
// ECFG_codon
// ECFG_codon_rev

struct ECFG_null_state_info : ECFG_state_info
{ ECFG_null_state_info(); };

struct ECFG_emitl_state_info : ECFG_state_info
{ ECFG_emitl_state_info (int matrix_idx); };

struct ECFG_emitr_state_info : ECFG_state_info
{ ECFG_emitr_state_info (int matrix_idx); };

struct ECFG_emitlr_state_info : ECFG_state_info
{ ECFG_emitlr_state_info (int matrix_idx, int alphabet_size, const char* tag); };

struct ECFG_bifurc_state_info : ECFG_state_info
{ ECFG_bifurc_state_info (int l, int r); };

template<class T>
struct ECFG : Transition_matrix <T, array2d <T, array2d_sparse_vector<T> > >, ECFG_enum
{
  // data
  const Alphabet& alphabet;

  // constructors
  ECFG (const Alphabet& alphabet);
  ECFG (const Alphabet& alphabet, int states);
  ECFG (const Alphabet& alphabet, int states, T t);
  ECFG (const ECFG<T>& ecfg);

  template<class S>
  ECFG (const ECFG<S>& ecfg, T t);

  // helpers
  int states() const { return this->tm_states(); }  // PK 2/15/05 see http://gborg.postgresql.org/pipermail/libpqxx-general/2004-July/000568.html
  virtual sstring desc (int state) const
  {
    sstring my_desc;
    my_desc << state;
    return my_desc;
  }
};

// Annotation
typedef pair<Subseq_coords,int> ECFG_subseq_state;
typedef map<ECFG_subseq_state,Loge> ECFG_cell_score_map;

// ECFG_funcs
struct ECFG_funcs : ECFG<PFunc>
{
  // constructors
  ECFG_funcs (const Alphabet& alphabet, int states)
    : ECFG<PFunc> (alphabet, states, PFunc())
  { }

  ECFG_funcs (const Alphabet& alphabet)
    : ECFG<PFunc> (alphabet)
  { }

  // display methods
  const char* element_descriptor() const { return "probability functions"; }
  int  element_width() const { return 16; }
  void show_element (const PFunc& element, ostream& o) const { element.show(o); }
};

// abstract interface for retrieving inside probabilities
struct ECFG_inside_calculator
{
  // virtual destructor
  virtual ~ECFG_inside_calculator() { }

  // virtual methods implemented by ECFG_inside_matrix
  virtual Loge state_inside_ll (int state, const Subseq_coords& subseq) const = 0;
};

// abstract interface for calculating state & transition posterior log-probabilities
struct ECFG_posterior_probability_calculator
{
  // virtual destructor
  virtual ~ECFG_posterior_probability_calculator() { }

  // virtual methods implemented by ECFG_inside_outside_matrix
  virtual Loge post_transition_ll (int src_state, int dest_state, int subseq_idx) const = 0;
  virtual Loge post_transition_ll (int src_state, int dest_state, const Subseq_coords& subseq) const = 0;
  virtual Loge post_state_ll (int dest_state, int subseq_idx) const = 0;
  virtual Loge post_state_ll (int dest_state, const Subseq_coords& subseq) const = 0;
};

// ECFG_scores
struct ECFG_scores : ECFG<Score>
{
  // data
  ECFG_matrix_set matrix_set;
  vector<ECFG_state_info> state_info;
  sstring name;
  list<SExpr> meta;  // contains "(meta ...)" expressions, which are preserved but ignored
  list<SExpr> transient_meta;  // contains "(meta ...)" expressions, which are ignored and NOT preserved (they're output, but won't be read back in again)

  // parametric stuff
  bool has_parametric_transitions;  // if true, then trans_funcs are used
  ECFG_funcs trans_funcs;  // transition PFunc's
  PScores pscores;  // parameter space & scores
  PCounts pcounts;  // pseudocounts for pscores
  set<int> mutable_pgroups;  // list of non-constant PGroup indices to be updated during EM

  // update flags
  bool update_rates, update_rules;

  // constructors
  ECFG_scores (const Alphabet& alphabet)
    : ECFG<Score> (alphabet), matrix_set (alphabet), state_info(), name (ECFG_default_name),
      has_parametric_transitions (false), trans_funcs (alphabet),
      update_rates (true), update_rules (true)
  { }

  ECFG_scores (const Alphabet& alphabet, int states)
    : ECFG<Score> (alphabet, states, -InfinityScore), matrix_set (alphabet), state_info (states), name (ECFG_default_name),
      has_parametric_transitions (false), trans_funcs (alphabet, states),
      update_rates (true), update_rules (true)
  {
    init_default_state_names();
  }

  ECFG_scores (const ECFG_scores& ecfg)
    : ECFG<Score> (ecfg), matrix_set (ecfg.matrix_set), state_info (ecfg.state_info), name (ecfg.name),
      has_parametric_transitions (ecfg.has_parametric_transitions), trans_funcs (ecfg.trans_funcs),
      pscores (ecfg.pscores), pcounts (ecfg.pcounts), mutable_pgroups (ecfg.mutable_pgroups),
      update_rates (ecfg.update_rates), update_rules (ecfg.update_rules)
  { }

  template<class T> ECFG_scores (const ECFG<T>& ecfg)
    : ECFG<Score> (ecfg), matrix_set (ecfg.alphabet), state_info (ecfg.states), name (ECFG_default_name),
      has_parametric_transitions (false), trans_funcs (alphabet, ecfg.states),
      update_rates (ecfg.update_rates), update_rules (ecfg.update_rules)
  {
    init_default_state_names();
  }

  void init_default_state_names();

  // methods to test whether this is a left- or a right-regular grammar (if so, a compact envelope can be used)
  bool is_left_regular() const;
  bool is_right_regular() const;

  // method to test if ECFG has *any* parametric functions whatsoever (transitions, chains, gap models...)
  bool has_parametric_functions() const;

  // method to test if ECFG has any chains with hidden classes
  bool has_hidden_classes() const;

  // method to test if there are any GFF annotations
  bool has_GFF() const;

  // method to set max_len for infix states, and issue warnings for non-prefix and non-suffix states with excessive max_len
  void set_infix_len (int max_subseq_len);

  // method to return a pointer to the first single-pseudoterminal chain in this grammar, or null if there is no such chain
  const ECFG_chain* first_single_pseudoterminal_chain() const;

  // helpers
  set<sstring> gc_feature_set() const;
  vector<int> nonemit_states_unsorted() const;
  vector<int> nonemit_states() const;  // sorted topologically; throws exception if null cycle detected
  vector<int> emit_states() const;
  vector<int> bifurc_states() const;
  vector<vector<int> > left_bifurc() const;  // left_bifurc()[rdest] = list of states bifurcating to (*,rdest)
  vector<vector<int> > right_bifurc() const;  // right_bifurc()[ldest] = list of states bifurcating to (ldest,*)

  // bifurcation pseudo-transitions to fake out the topological sort
  void add_fake_bifurcation_transitions();
  void remove_fake_bifurcation_transitions();

  // I/O
  void write (ostream& out) const;
  void read (istream& in);

  // annotation
  void annotate (Stockholm& stock, const ECFG_cell_score_map& annot) const;
  void make_GFF (GFF_list& gff_list,
		 const ECFG_cell_score_map& annot,
		 const char* seqname = ECFG_GFF_default_seqname,
		 ECFG_posterior_probability_calculator* pp_calc = 0,
		 ECFG_inside_calculator* ins_calc = 0) const;

  // evaluate PFunc's for parametric rates & probabilities
  void eval_funcs();

  // display methods
  void show (ostream& out) const;

  const char* element_descriptor() const { return "scores"; }
  int         element_width() const { return 10; }
  void        show_element (const Score& element, ostream& o) const { ShowScore (element, o); }
};

// ECFG_counts: container for accumulating expected counts during E-step of EM
struct ECFG_counts : ECFG<Prob>
{
  // data
  const ECFG_scores* ecfg;
  vector<Update_statistics> stats;  // EM_matrix update statistics
  vector<bool> filled_down;  // flag indicating whether, for each matrix, fill_down() has been called at least once

  // indel update statistics (pseudo-indel states only)
  vector<double> ins_count, del_count, ins_wait, del_wait;
  vector<double> link_extend_count, link_end_count;

  // annotation counts
  typedef map<sstring,Prob> String_counts;
  vector<vector<String_counts> > state_annot_count;  // indexed as state_annot_count[stateIndex][gcFeatureTag][columnAnnotationString]

  // parameter counts
  PCounts var_counts;

  // constructor
  ECFG_counts (const ECFG_scores& ecfg);

  // clear method
  void clear (double pseud_init = 0., double pseud_mutate = 0., double pseud_wait = 0.);

  // methods for doing M-step of EM
  void update_ecfg (ECFG_scores& ecfg);
  void update_ecfg_rates (ECFG_scores& ecfg);
  void update_ecfg_rules (ECFG_scores& ecfg);
  void update_ecfg_params (ECFG_scores& ecfg);

  // transmit counts to PCounts
  void inc_var_counts (PCounts&           var_counts,
		       const ECFG_scores& ecfg,
		       const Prob         weight = 1.0);

  // show method
  void show (const ECFG_scores& ecfg, ostream& out) const;

  // inherited output methods
  const char* element_descriptor() const { return "counts"; }
  int         element_width() const { return 10; }
  void        show_element (const Prob& element, ostream& o) const { o << element; }
};

// automatically initialized ECFG_envelope
struct ECFG_auto_envelope : ECFG_envelope {
  ECFG_auto_envelope (int seqlen, const ECFG_scores& ecfg, int max_subseq_len = -1) {
    init (seqlen, ecfg.is_left_regular() || ecfg.is_right_regular() ? 0 : max_subseq_len);
  }
};

// class to sample a multiple alignment, given a tree
struct ECFG_simulation : Stockade, Grammar_state_enum
{
  // data
  const ECFG_scores& ecfg;
  ECFG_counts counts;

  // constructor
  ECFG_simulation (const ECFG_scores& ecfg, const PHYLIP_tree& tree);

  // method to add the chain counts into the alignment as a #=GF field
  void add_counts_to_stockade();
};

// inline & templated method defs

// ECFG_state_info
bool ECFG_state_info::init_row (const Alphabet& alphabet, int hidden_classes, const Aligned_score_profile& align, int row,
				const Subseq_coords& subseq, vector<Loge>& colmat_row, vector<Prob>* fastprune_row,
				int& gapped, bool use_wildcards_for_emit_columns) const
{
  // get iterators for all positions
  int count = 1;  // total number of combinations
  vector<Symbol_score_iterator> begin, end;
  const int size = total_size();
  const int alph_size = alphabet.size();
  begin.reserve (size);
  end.reserve (size);
  int gaps = 0;
  for (int i = 0; i < size; ++i)
    {
      const bool is_emit = is_emit_column (i);
      const int col = column_index (subseq, i);
      bool wildcard = is_emit ? use_wildcards_for_emit_columns : (col < 0 || col >= align.columns());
      if (!wildcard)
	{
	  const Symbol_score_map* ssm = align (row, col);
	  if (ssm == 0)  // gap
	    {
	      if ((!is_emit) || (gaps_ok && wild_gaps))
		wildcard = true;  // put in a wildcard if this is a context (rather than an emit) column, or if wild_gaps is true
	      else  // if we get here, then is_emit==true, and either gaps_ok==false or wild_gaps==false
		++gaps;
	    }
	  else  // not a gap
	    {
	      begin.push_back (ssm->begin());
	      end.push_back (ssm->end());
	      count *= ssm->size();
	    }
	}
      if (wildcard)
	{
	  begin.push_back (alphabet.wild_ssm.begin());
	  end.push_back (alphabet.wild_ssm.end());
	  count *= alph_size;
	}
    }

  // in the event that wild_gaps==false, decide how to handle gaps
  if (gaps_ok && gaps == emit_size())
    {
      // all gaps, so return OK but leave colmat_row and fastprune_row untouched
      gapped = true;
      return true;
    }
  else if (gaps)
    {
      // partial gaps are not OK
      gapped = true;
      return false;
    }

  // calculate multiplier for hidden class index
  int class_mul = 1;
  for (int i = 0; i < size; ++i)
    class_mul *= alph_size;

  // iterate through all symbols at all positions
  vector<Symbol_score_iterator> iter (begin);
  while (count-- > 0)
    {
      // compute symbol index
      int sym = 0;
      Score sc = 0;
      for (int i = 0; i < size; ++i)
	{
	  int c = (*iter[i]).first;
	  if (comp[i])
	    c = alphabet.complement (c);
	  sym += mul[i] * c;
	  ScorePMulAcc (sc, (*iter[i]).second);
	}

      // write (symbol,prob) pair into Column_matrix and Fast_prune objects for all classes
      const Loge ll = Score2Nats (sc);
      const Prob p = Score2Prob (sc);
      for (int hidden_class = 0; hidden_class < hidden_classes; ++hidden_class)
	{
	  const int i = sym + hidden_class * class_mul;
	  colmat_row[i] = ll;
	  if (fastprune_row)
	    (*fastprune_row)[i] = p;
	}

      // increment iterators
      for (int j = 0; j < size; ++j)
	if (++iter[j] == end[j])
	  iter[j] = begin[j];
	else
	  break;
    }

  // returned OK
  gapped = false;
  return true;
}

bool ECFG_state_info::node_is_wild (const Stockholm& stock, const Stockholm_tree& tree, const Subseq_coords& subseq, int node) const
{
  const int row = tree.node2row[node];
  if (row < 0)
    return true;

  // node is wild iff wild_gaps is true, and all children are gapped
  const int sz = emit_size();
  for (int i = 0; i < sz; ++i)
    if (stock.path (row, column_index (subseq, l_context + i)))
      return false;
  return wild_gaps;
}

bool ECFG_state_info::node_is_redundant (const Stockholm& stock, const Stockholm_tree& tree, const Subseq_coords& subseq,
					 const vector<int>& gapped, const Phylogeny::Branch& b) const
{
  const int p = b.first;
  const int n = b.second;
  bool prune = false;

  if (p < 0 ? true : gapped[p])  // only proceed if at root, or parent is gapped
    if (node_is_wild (stock, tree, subseq, n))  // only proceed if node is a wildcard
      {
	int n_ungapped_children = 0;
	for_children (tree, n, p, c)
	  if (!gapped[*c])  // loop over ungapped children
	    ++n_ungapped_children;
	if (n_ungapped_children == 1)
	  prune = true;
      }

  return prune;
}

int ECFG_state_info::column_index (const Subseq_coords& subseq, int pos) const
{
  pos -= l_context;
  return pos < l_emit ? subseq.start + pos : subseq.end() + pos - (l_emit + r_emit);
}

bool ECFG_state_info::is_emit_column (int pos) const
{
  pos -= l_context;
  return pos >= 0 && pos < l_emit + r_emit;
}

bool ECFG_state_info::initialise (const Alphabet& alphabet, int hidden_classes,
				  const Aligned_score_profile& align, const Subseq_coords& subseq,
				  const Stockholm& stock, const Stockholm_tree& tree,
				  Column_matrix& colmat, Fast_prune* fast_prune,
				  bool use_wildcards_for_emit_columns) const
{
  if (out_of_range (subseq))
    return false;

  colmat.clear();
  if (fast_prune)
    fast_prune->clear();
  for_rooted_nodes_post (tree, b)
    {
      const int p = (*b).first;
      const int n = (*b).second;
      const int row = tree.node2row[n];

      if (row < 0)  // row has no alignment path info; treat as a wildcard iff at least one child ungapped
	{
	  colmat.gapped[n] = true;
	  for_children (tree, n, p, c)
	    if (!colmat.gapped[*c])
	      {
		// wildcard
		colmat.gapped[n] = false;
		fill (colmat.U[n].begin(), colmat.U[n].end(), 0.);
		if (fast_prune)
		  fill (fast_prune->F[n].begin(), fast_prune->F[n].end(), 1.);
		break;
	      }
	}
      else
	{
	  // if a Fast_prune object was supplied, then populate that
	  vector<Prob>* fastprune_row = fast_prune ? &fast_prune->F[n] : (vector<Prob>*) 0;
	    
	  if (!init_row (alphabet, hidden_classes, align, row, subseq, colmat.U[n], fastprune_row, colmat.gapped[n], use_wildcards_for_emit_columns))
	    return false;
	}
    }

  // prune lineages of *'s
  for_rooted_nodes_pre (tree, b)
    {
      const int n = (*b).second;
      if (node_is_redundant (stock, tree, subseq, colmat.gapped, *b))
	{
	  colmat.gapped[n] = true;
	  fill (colmat.U[n].begin(), colmat.U[n].end(), -InfinityLoge);
	  if (fast_prune)
	    fill (fast_prune->F[n].begin(), fast_prune->F[n].end(), 0.);
	}
    }

  // init signposts and return
  colmat.init_signposts (tree);
  return colmat.clique.size() == 1;  // must have exactly one clique
}

void ECFG_state_info::get_gap_profile (const Stockholm& stock, const Stockholm_tree& tree, const Subseq_coords& subseq,
				       vector<int>& gapped, int& root) const
{
  root = -1;
  if (out_of_range (subseq))
    {
      // ensure that out-of-range subseqs have a well-defined gap profile
      for_contents (vector<int>, gapped, g)
	*g = true;
      return;
    }

  const Alignment_path& path = stock.path;
  for_rooted_nodes_post (tree, b)
    {
      const int p = (*b).first;
      const int n = (*b).second;
      const int row = tree.node2row[n];

      if (row < 0)  // row has no alignment path info; treat as a wildcard iff at least one child is ungapped
	{
	  gapped[n] = true;
	  for_children (tree, n, p, c)
	    if (!gapped[*c])
	      {
		// wildcard
		gapped[n] = false;
		break;
	      }
	}
      else
	{
	  gapped[n] = false;
	  for (int i = 0; i < l_emit + r_emit; ++i)
	    if (!path (row, column_index (subseq, l_context + i)))
	      {
		gapped[n] = true;
		break;
	      }
	}
    }

  // prune redundant nodes
  for_rooted_nodes_pre (tree, b)
    if (node_is_redundant (stock, tree, subseq, gapped, *b))
      gapped[(*b).second] = true;

  // find root
  for_rooted_nodes_post (tree, b)
    if (!gapped[(*b).second])
      root = (*b).second;
}

// ECFG
template<class T>
ECFG<T>::ECFG (const Alphabet& alph)
  : Transition_matrix <T, array2d <T, array2d_sparse_vector<T> > > (),
    alphabet (alph)
{ }

template<class T>
ECFG<T>::ECFG (const Alphabet& alph, int states)
  : Transition_matrix <T, array2d <T, array2d_sparse_vector<T> > > (states),
    alphabet (alph)
{ }

template<class T>
ECFG<T>::ECFG (const Alphabet& alph, int states, T t)
  : Transition_matrix <T, array2d <T, array2d_sparse_vector<T> > > (states, t),
    alphabet (alph)
{ }

template<class T>
ECFG<T>::ECFG (const ECFG<T>& ecfg)
  : Transition_matrix <T, array2d <T, array2d_sparse_vector<T> > > (ecfg),
    alphabet (ecfg.alphabet)
{ }

template<class T>
template<class S>
ECFG<T>::ECFG (const ECFG<S>& ecfg, T t)
  : Transition_matrix <T, array2d <T, array2d_sparse_vector<T> > > (ecfg.states(), t),
    alphabet (ecfg.alphabet)
{ }

#endif /* ECFG_INCLUDED */
