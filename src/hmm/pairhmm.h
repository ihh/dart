#ifndef PAIRHMM_INCLUDED
#define PAIRHMM_INCLUDED

#include "seq/alignment.h"
#include "hmm/transmat.h"
#include "util/strsaver.h"

// default names
#define Pair_HMM_default_name         "Pair_HMM"
#define Pair_HMM_default_state_prefix "Pair_HMM_state_"

// Pair_HMM state type enum
struct Pair_HMM_state_type_enum
{
  enum State_type { EmitX = 1, EmitY = 2, EmitXY = 3, Null = 0 };
  static int type_dx (State_type t) { return t & 1; }
  static int type_dy (State_type t) { return (t & 2) >> 1; }
};

// Templated Pair_HMM class is the base class for all pair HMM parameter classes
//
template<class T>
struct Pair_HMM : Transition_matrix<T>, Pair_HMM_state_type_enum
{
  // typedefs
  typedef T entry_type;

  // names
  sstring name;
  vector<sstring> state_name;

  // state type & emit data
  vector<State_type>   state_type;
  vector<vector<T> >   single_emit;
  vector<array2d<T> >  pair_emit;

  // alphabet
  const Alphabet*      alphabet;

  // constructors
  Pair_HMM (int states, const Alphabet* alphabet = 0);
  Pair_HMM (int states, T t, const Alphabet* alphabet = 0);

  // builders
  void init_emit (int state, State_type type, const T& val);

  // accessors
  int states() const { return state_type.size(); }   // excludes start & end states

  // dx/dy helpers
  int dx (int state) const { return type_dx (state_type[state]); }
  int dy (int state) const { return type_dy (state_type[state]); }
  static State_type dxdy_type (int dx, int dy) { return dx==1 ? (dy==1 ? EmitXY : EmitX) : (dy==1 ? EmitY : Null); }

  // show method
  void show (ostream& o) const;

  // method to check if same dimensions as another Pair HMM
  template<class S>
  bool same_dimensions (const Pair_HMM<S>& s) const;
};

// Pair_HMM_scores contains transition and emission scores
//
struct Pair_HMM_scores : Pair_HMM<Score>
{
  Pair_HMM_scores (int states = 0, const Alphabet* alphabet = 0) : Pair_HMM<int> (states, -InfinityScore, alphabet) { }

  vector<int>   find_optimal_state_path         (const Pairwise_path& path) const;      // looks for optimal state path by DP, throws exception if none compatible
  Pairwise_path convert_state_path_to_alignment (const vector<int>& state_path) const;

  Score path_transition_score (const vector<int>& state_path) const;
  Score path_emit_score (const vector<int>& state_path, const Score_profile& xseq, const Score_profile& yseq) const;
  Score path_score (const vector<int>& state_path, const Score_profile& xseq, const Score_profile& yseq) const;

  Score pairwise_path_score (const Score_profile& xseq, const Score_profile& yseq, const Pairwise_path& path);  // returns Forward score

  void scale_all_scores (double beta);     // multiplies all scores by beta

  // Regexps for parsing BLAST matrices
  static Regexp header_line_regexp;
  static Regexp nonempty_line_regexp;

  // method to read a BLAST matrix into the pair_emit array for a state
  void read_BLAST_matrix (istream& in, int state);
  
  // sparseness methods
  vector<int>           null_states_unsorted() const;   // null states
  vector<int>           null_states() const;            // null states, sorted topologically (croaks on null cycles)
  vector<int>           emit_states() const;            // non-null states

  // null state elimination method
  void eliminate_null_states();

  // show methods
  const char* element_descriptor() const { return "scores"; }
  int  element_width() const { return 10; }
  void show_element (const Score& element, ostream& o) const { ShowScore (element, o); }
};

// Pair_HMM_derivatives contains derivatives (w.r.t. some variable) of transition probabilities
//
struct Pair_HMM_derivatives : Pair_HMM<double>
{
  sstring wrt;  // derivative is "with respect to" (wrt) the parameter described in this string (for display only)
  Pair_HMM_derivatives (int states, const Alphabet* alphabet = 0);
  Pair_HMM_derivatives (int states, const char* x, const Alphabet* alphabet = 0);   // x is a string describing the parameter w.r.t. which the derivatives are taken

  const char* element_descriptor() const { return wrt.c_str(); }
  int  element_width() const { return 10; }
  void show_element (const double& element, ostream& o) const { o << element; }
};

// Pair_HMM_mask structure is used to restrict Pair_HMM_scores parameter updates (e.g. during expectation maximization).
//
struct Pair_HMM_mask : Pair_HMM<int>
{
  Pair_HMM_mask (const Pair_HMM_scores& hmm, bool init_flag = 1);
  void set_incoming_transitions (int state);
  void clear_incoming_transitions (int state);
  void set_outgoing_transitions (int state);
  void clear_outgoing_transitions (int state);
  void set_emit (int state);
  void clear_emit (int state);
  void set_all_start();
  void clear_all_start();
  void set_all_end();
  void clear_all_end();
  void set_all_transitions();
  void clear_all_transitions();
  void set_all_emit();
  void clear_all_emit();
  const char* element_descriptor() const { return "update flags"; }
  int  element_width() const { return 10; }
  void show_element (const double& element, ostream& o) const { o << element; }
};

// Pair_HMM_counts structure contains expected transition-usage counts
//
// NB:    Expected usage count for transition T   =   T * dP/dT    where P is the likelihood
//
struct Pair_HMM_counts : Pair_HMM<double>
{
  double  log_likelihood;    //  log(P)

  Pair_HMM_counts (const Pair_HMM_scores& hmm);

  void   add_counts (const Pair_HMM_counts& counts);
  void   subtract_counts (const Pair_HMM_counts& counts);
  void   update_HMM_scores (Pair_HMM_scores& hmm, const Pair_HMM_mask& mask, bool sample = 0, double kT = 1) const;    // if sample==TRUE, probabilities are Dirichlet with pseudocounts scaled by 1/kT

  double dloglike_dx (const Pair_HMM_scores& hmm, const Pair_HMM_derivatives& deriv) const;         // deriv structure contains derivatives dT/dX of transition & emission probabilities w.r.t. X

  Score   add_counts_from_unaligned_sequences (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq);  // returns Forward score
  Score   add_counts_from_state_path (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const vector<int>& state_path);  // returns score of state path

  const char* element_descriptor() const { return "counts"; }
  int  element_width() const { return 10; }
  void show_element (const Prob& element, ostream& o) const { o << element; }

  Score add_counts_from_aligned_sequences (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const Pairwise_path& path);  // returns alignment score
};

// DP classes
//
struct Pair_emit_score_calculator
{
  const Pair_HMM_scores&  hmm;
  const Score_profile&    xseq;
  const Score_profile&    yseq;

  vector<Score_profile*>  xseq_evolved;         // use pointers for speed, because vector is sparse

  Pair_emit_score_calculator (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq);
  ~Pair_emit_score_calculator();

  inline Score calc_single_emit_score (const vector<Score>& single, const Symbol_score_map& ssmap) const;
  inline Score calc_pair_emit_score (const array2d<Score>& pair, const Symbol_score_map& ssmapx, const Symbol_score_map& ssmapy) const;
  inline Score calc_symbol_score_map_inner_product (const Symbol_score_map& ssmapx, const Symbol_score_map& ssmapy) const;
  inline Score calc_emit_score (int state, int xpos, int ypos) const;

  // in the following functions, posterior_sc does *not* include emit score
  // (thus in fact it is misnamed -- it should be posterior_minus_emit_sc or something like that)
  inline void accum_single_emit_score (Score posterior_sc, Pair_HMM_counts& counts, int state, const Symbol_score_map& ssmap);
  inline void accum_pair_emit_score (Score posterior_sc, Pair_HMM_counts& counts, int state, const Symbol_score_map& ssmapx, const Symbol_score_map& ssmapy);
  inline void accum_emit_score (Score posterior_sc, Pair_HMM_counts& counts, int state, int xpos, int ypos);
};

struct Pair_DP_matrix_base : Pair_emit_score_calculator, Stream_saver
{
  const Alphabet*       alphabet;
  
  array2d<int> pair_env;  // pair_env(i,j)=1 => can align subseqs 1..i of x and 1..j of y
  vector<array2d<Score> > cell;

  const int             xsize;      // equals xseq.size() + 1
  const int             ysize;      // equals yseq.size() + 1
  const int             states;     // excludes start & end states
  Score                 final_score;
  
  inline Pair_DP_matrix_base (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq);
  inline Pair_DP_matrix_base (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const array2d<int>& pair_env);

  void show (ostream& o) const;
  void show_sparse (ostream& o) const;  // only prints cells for which pair_env==1
  sstring hmm_dump() const;

  inline char score2char (const Symbol_score_map& m) const;
};


class Pair_Viterbi_DP_matrix : public Pair_DP_matrix_base
{
public:
  Pair_Viterbi_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq);

  vector<int> optimal_state_path() const;  // excludes Start, End

private:
  inline Score calc_score (int dest_st, int xpos, int ypos, Score sc) const;
  void fill();
};


class Pair_forward_DP_matrix : public Pair_DP_matrix_base
{
public:
  Pair_forward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq);
  Pair_forward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const array2d<int>& pair_env);

  vector<int> sample_state_path() const;   // excludes Start, End

private:
  inline Score calc_score (int dest_st, int xpos, int ypos, Score sc) const;
  void fill();
};

class Pair_forward_backward_DP_matrix : public Pair_DP_matrix_base
{
  // after backward recursion, cell [s](x,y) contains probability of emitting
  // final (xseq.size()+dx-x) symbols of xseq & final (yseq.size()+dy-y) symbols of yseq
  // starting from state s (with the state's emission included in the backward subseq),
  // where dx,dy are emit vectors from state s
  //  ... i.e. subseqs x+1-dx..L and y+1-dy..M (if L==xseq.size() and M==yseq.size()) ...
  // NB xsize = xseq.size() + 1, ysize = yseq.size() + 1
public:
  Pair_forward_DP_matrix forward_matrix;
  Pair_HMM_counts counts;

  const vector<array2d<Score> >& forward_cell;  // reference to forward_matrix.cell
  const Score& forward_score;  // reference to forward_matrix.final_score

  Pair_forward_backward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq);
  Pair_forward_backward_DP_matrix (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const array2d<int>& pair_env);

  Score post_state_score (int state, int x, int y) const;  // posterior probability of being in cell [state] after emitting (x) symbols of xseq & (y) symbols of yseq

  void dump_forward_cell (ostream& o) const;  // hack method for displaying preserved_forward_cell

private:

  inline Score back_calc_score (int src_st, int xpos, int ypos);
  void fill();
};

///////////////////////////////////////////////////////////////////////
//                                                                   //
//                template & inline method defs                      //
//                                                                   //
///////////////////////////////////////////////////////////////////////


template<class T>
Pair_HMM<T>::Pair_HMM (int states, const Alphabet* alphabet) :
  Transition_matrix<T> (states),
  name (Pair_HMM_default_name),
  state_name (states, sstring (Pair_HMM_default_state_prefix)),
  state_type (states, (State_type) Null),
  single_emit (states, vector<T>()),
		 pair_emit (states, array2d<T>()),
			      alphabet (alphabet)
{
  for (int s = 0; s < states; ++s)
    state_name[s] << (s + 1);
}

template<class T>
Pair_HMM<T>::Pair_HMM (int states, T t, const Alphabet* alphabet) :
  Transition_matrix<T> (states, t),
  name (Pair_HMM_default_name),
  state_name (states, sstring (Pair_HMM_default_state_prefix)),
  state_type (states, (State_type) Null),
  single_emit (states, vector<T>()),
		 pair_emit (states, array2d<T>()),
			      alphabet (alphabet)
{
  for (int s = 0; s < states; ++s)
    state_name[s] << (s + 1);
}

template<class T>
void Pair_HMM<T>::init_emit (int state, State_type type, const T& val)
{
  if (alphabet == 0) THROWEXPR ("Null alphabet in call to init_emit");
  const int A = alphabet->size();
  state_type[state] = type;
  if (type == EmitXY)
    pair_emit[state] = array2d<T> (A, A, val);
  else if (type == EmitX || type == EmitY)
    single_emit[state] = vector<T> (A, val);
}

template<class T>
void Pair_HMM<T>::show (ostream& o) const
{
  int old_prec = o.precision(3);
  this->save_flags (o);
  this->right_align (o);

  const int w = this->element_width();

  this->show_transitions (o);

  o << "Emission profile " << this->element_descriptor() << ":\n";
  if (alphabet)
    {
      o << "                  ";
      for (int s = 0; s < alphabet->size(); s++) { o << " "; o.width(w); o << alphabet->int2char(s); }
      o << "\n";
    }
  for (int i = 0; i < states(); i++)
    {
      o.width(2);
      o << i;
      switch (state_type[i])
	{
	case Null:   o << " (type null)    "; break;
	case EmitX:  o << " (type x-emit)  "; break;
	case EmitY:  o << " (type y-emit)  "; break;
	case EmitXY: o << " (type xy-emit) "; break;
	default: o << " (type unknown)\n"; break;
	}
      if (state_type[i] == EmitX || state_type[i] == EmitY)
	{
	  o << "(";
	  for (int x = 0; x < (int) single_emit[i].size(); x++)
	    { o.width(w); Pair_HMM<T>::show_element(single_emit[i][x],o); if (x < (int) single_emit[i].size()-1) o << " "; }
	  o << ")\n";
	}
      else if (state_type[i] == EmitXY)
	{
	  for (int y = 0; y < pair_emit[i].ysize(); y++)
	    {
	      if (y > 0) o << "                  ";
	      o << "(";
	      for (int x = 0; x < pair_emit[i].xsize(); x++) { o.width(w); Pair_HMM<T>::show_element(pair_emit[i](x,y),o); if (x < pair_emit[i].xsize()-1) o << " "; }
	      o << ")";
	      if (alphabet) o << " " << alphabet->int2char(y);
	      o << "\n";
	    }
	}
      else
	o << "\n";
    }
  this->restore_flags (o);
  o.precision (old_prec);
}

template<class T>
template<class S>
bool Pair_HMM<T>::same_dimensions (const Pair_HMM<S>& s) const
{
  if (states() != s.states()) return 0;
  for (int i = 0; i < states(); ++i)
    {
      if (state_type[i] != (State_type) s.state_type[i]) return 0;
      switch (state_type[i])
	{
	case EmitXY:
	  if (pair_emit[i].xsize() != s.pair_emit[i].xsize() || pair_emit[i].ysize() != s.pair_emit[i].ysize()) return 0;
	  break;
	case EmitX:
	case EmitY:
	  if (single_emit[i].size() != s.single_emit[i].size()) return 0;
	  break;
	case Null:
	  break;
	default:
	  THROW Standard_exception ("Unknown state type!");
	  break;
	}
    }
  return 1;
}


Score Pair_emit_score_calculator::calc_single_emit_score (const vector<Score>& single, const Symbol_score_map& ssmap) const
{
  Score sc = -InfinityScore;
  for_const_contents (Symbol_score_map, ssmap, ss)
    ScorePSumAcc (sc, ScorePMul (single[(*ss).first], (*ss).second));
  return sc;
}

Score Pair_emit_score_calculator::calc_pair_emit_score (const array2d<Score>& pair, const Symbol_score_map& ssmapx, const Symbol_score_map& ssmapy) const
{
  Score sc = -InfinityScore;
  for_const_contents (Symbol_score_map, ssmapx, ssx)
    {
      for_const_contents (Symbol_score_map, ssmapy, ssy)
	ScorePSumAcc (sc, ScorePMul3 (pair ((*ssx).first, (*ssy).first), (*ssx).second, (*ssy).second));
    }
  return sc;
}

Score Pair_emit_score_calculator::calc_symbol_score_map_inner_product (const Symbol_score_map& ssmapx, const Symbol_score_map& ssmapy) const
{
  Score sc = -InfinityScore;
  for_const_contents (Symbol_score_map, ssmapy, ssy)
    ScorePSumAcc (sc, (*ssy).second + ((Symbol_score_map&) ssmapx) [(*ssy).first]);
  return sc;
}

Score Pair_emit_score_calculator::calc_emit_score (int state, int xpos, int ypos) const
{
  Score sc;
  switch (hmm.state_type[state])
    {
    case Pair_HMM_scores::Null: sc = 0; break;
    case Pair_HMM_scores::EmitX: sc = calc_single_emit_score (hmm.single_emit[state], xseq[xpos-1]); break;
    case Pair_HMM_scores::EmitY: sc = calc_single_emit_score (hmm.single_emit[state], yseq[ypos-1]); break;
    case Pair_HMM_scores::EmitXY: sc = calc_symbol_score_map_inner_product ((*xseq_evolved[state])[xpos-1], yseq[ypos-1]); break;
      // this is what the old 'case' used to be. calc_pair_emit_score() is slower than calc_symbol_score_map_inner_product().
      //
      //	case Pair_HMM_scores::EmitXY: sc = calc_pair_emit_score (hmm.pair_emit[state], xseq[xpos-1], yseq[ypos-1]); break;
    default: sc = -InfinityScore; break;
    }
  return sc;
}

// in the following functions, posterior_sc does *not* include emit score
// (thus in fact it is misnamed -- it should be posterior_minus_emit_sc or something like that)
void Pair_emit_score_calculator::accum_single_emit_score (Score posterior_sc, Pair_HMM_counts& counts, int state, const Symbol_score_map& ssmap)
{
  const vector<Score>&     single         =  hmm.single_emit[state];
  vector<double>&          single_counts  =  counts.single_emit[state];
  for_const_contents (Symbol_score_map, ssmap, ss)
    single_counts[(*ss).first] += Score2Prob (ScorePMul3 ((*ss).second, single[(*ss).first], posterior_sc));
}

void Pair_emit_score_calculator::accum_pair_emit_score (Score posterior_sc, Pair_HMM_counts& counts, int state, const Symbol_score_map& ssmapx, const Symbol_score_map& ssmapy)
{
  const array2d<Score>&    pair          =  hmm.pair_emit[state];
  array2d<double>&         pair_counts   =  counts.pair_emit[state];
  for_const_contents (Symbol_score_map, ssmapx, ssx)
    {
      for_const_contents (Symbol_score_map, ssmapy, ssy)
	pair_counts((*ssx).first, (*ssy).first) += Score2Prob (ScorePMul (ScorePMul ((*ssx).second, (*ssy).second),
									  ScorePMul (pair((*ssx).first,(*ssy).first), posterior_sc)));
    }
}

void Pair_emit_score_calculator::accum_emit_score (Score posterior_sc, Pair_HMM_counts& counts, int state, int xpos, int ypos)
{
  switch (hmm.state_type[state])
    {
    case Pair_HMM_scores::Null:
      break;
    case Pair_HMM_scores::EmitX:
      accum_single_emit_score (posterior_sc, counts, state, xseq[xpos-1]);
      break;
    case Pair_HMM_scores::EmitY:
      accum_single_emit_score (posterior_sc, counts, state, yseq[ypos-1]);
      break;
    case Pair_HMM_scores::EmitXY:
      accum_pair_emit_score (posterior_sc, counts, state, xseq[xpos-1], yseq[ypos-1]);
      break;
    default:
      break;
    }
}

Pair_DP_matrix_base::Pair_DP_matrix_base (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq) :
  Pair_emit_score_calculator (hmm, xseq, yseq),
  alphabet (hmm.alphabet),
  pair_env (xseq.size() + 1, yseq.size() + 1, (int) 1),
  cell (hmm.states(), array2d<int> (xseq.size() + 1, yseq.size() + 1, -InfinityScore)),
  xsize (xseq.size() + 1), ysize (yseq.size() + 1), states (hmm.states()),
  final_score (-InfinityScore)
{
  if (CTAGGING(4,ALLOC PAIRENV)) CL << "Allocated Pair HMM DP matrix:  " << xsize << " * " << ysize << "  cells by  " << states << "  states; complete pair envelope\n";
  if (CTAGGING(-3,PAIRHMMDP_MODEL))
    { CL << "Pair HMM (in DP matrix constructor):\n"; hmm.show (CL); }
}

Pair_DP_matrix_base::Pair_DP_matrix_base (const Pair_HMM_scores& hmm, const Score_profile& xseq, const Score_profile& yseq, const array2d<int>& pair_env) :
  Pair_emit_score_calculator (hmm, xseq, yseq),
  alphabet (hmm.alphabet),
  pair_env (pair_env),
  cell (hmm.states(), array2d<int> (xseq.size() + 1, yseq.size() + 1, -InfinityScore)),
  xsize (xseq.size() + 1), ysize (yseq.size() + 1), states (hmm.states()),
  final_score (-InfinityScore)
{
  if (pair_env.xsize() != xsize || pair_env.ysize() != ysize)
    THROWEXPR ("Pair envelope is " << pair_env.xsize() << "*" << pair_env.ysize() << ", DP matrix is " << xsize << "*" << ysize << "\n");
  if (CTAGGING(4,ALLOC PAIRENV))
    {
      CL << "Allocated Pair HMM DP matrix:  " << xsize << " * " << ysize << "  cells by  " << states << "  states; partial pair envelope\n";
      if (CTAGGING(1,PAIRENV)) CL << "Pair envelope:\n" << pair_env;
    }
  if (CTAGGING(-3,PAIRHMMDP_MODEL))
    { CL << "Pair HMM (in DP matrix constructor):\n"; hmm.show (CL); }
}

char Pair_DP_matrix_base::score2char (const Symbol_score_map& m) const
{
  if (alphabet) return alphabet->score2char (m);
  int i = Score_profile::score_map_consensus(m);
  if (i < 10) return i + '0';
  return i-10 + 'a';
}

Score Pair_Viterbi_DP_matrix::calc_score (int dest_st, int xpos, int ypos, Score sc) const
{
  int prevx = xpos - hmm.dx(dest_st);
  int prevy = ypos - hmm.dy(dest_st);
  for (int src_st = 0; src_st < states; src_st++)
    sc = max (sc, ScorePMul (cell[src_st] (prevx, prevy), hmm.transition (src_st, dest_st)));
  return ScorePMul (sc, calc_emit_score (dest_st, xpos, ypos));
}

Score Pair_forward_DP_matrix::calc_score (int dest_st, int xpos, int ypos, Score sc) const
{
  int prevx = xpos - hmm.dx(dest_st);
  int prevy = ypos - hmm.dy(dest_st);
  for (int src_st = 0; src_st < states; src_st++)
    ScorePSumAcc (sc, cell[src_st] (prevx, prevy) + hmm.transition (src_st, dest_st));
  return ScorePMul (sc, calc_emit_score (dest_st, xpos, ypos));
}

Score Pair_forward_backward_DP_matrix::back_calc_score (int src_st, int xpos, int ypos)
{
  Score back_sc = -InfinityScore;
      
  const Pair_HMM_scores::State_type t = hmm.state_type[src_st];
  if (t == Pair_HMM_scores::Null
      || (t == Pair_HMM_scores::EmitX && xpos > 0)
      || (t == Pair_HMM_scores::EmitY && ypos > 0)
      || (t == Pair_HMM_scores::EmitXY && xpos > 0 && ypos > 0))
    {
      if (xpos == xsize-1 && ypos == ysize-1)
	{
	  back_sc = hmm.end[src_st];
	  counts.end[src_st] = Score2Prob (ScorePMul3 (forward_cell[src_st] (xsize-1,ysize-1), back_sc, -forward_score));
	}
      for (int dest_st = 0; dest_st < states; dest_st++)
	{
	  const int dest_x = xpos + hmm.dx(dest_st);
	  const int dest_y = ypos + hmm.dy(dest_st);
	  if (dest_x < xsize && dest_y < ysize)
	    {
	      const Score dest_sc = cell[dest_st] (dest_x, dest_y);
	      const Score trans_sc = hmm.transition (src_st, dest_st);
	      ScorePSumAcc (back_sc, ScorePMul (dest_sc, trans_sc));
	      counts.transition(src_st, dest_st) += Score2Prob (ScorePMul (ScorePMul3 (forward_cell[src_st] (xpos,ypos), trans_sc, dest_sc), -forward_score));
	    }
	}
      const Score emit_sc = calc_emit_score(src_st,xpos,ypos);
      accum_emit_score (ScorePMul (ScorePMul3 (forward_cell[src_st] (xpos,ypos), -emit_sc, back_sc), -forward_score), counts, src_st, xpos, ypos);   // have to subtract off emit_sc so it doesn't get counted twice
      ScorePMulAcc (back_sc, emit_sc);
    }
  return back_sc;
}

#endif
