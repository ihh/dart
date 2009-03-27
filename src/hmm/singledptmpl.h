#ifndef SINGLE_DP_TMPL_INCLUDED
#define SINGLE_DP_TMPL_INCLUDED

#include <time.h>
#include "hmm/singletmpl.h"

// abstract DP interfaces

struct Single_Viterbi_interface
{
  virtual Score       get_Viterbi_score() const;
  virtual vector<int> optimal_state_path() const;

  virtual ~Single_Viterbi_interface();
};

struct Single_forward_interface
{
  virtual Score       get_forward_score() const;
  virtual vector<int> sample_state_path() const;

  virtual ~Single_forward_interface();
};

struct Single_forward_backward_interface : Stream_saver
{
  virtual Score                    get_forward_score() const;
  virtual const Single_HMM_counts& get_expected_counts() const;
  virtual const vector<Metaprob>&  get_expected_metacounts() const;

  virtual ~Single_forward_backward_interface();

  void show_metacounts (ostream& o) const;
};

// matrix factory interface

struct Single_matrix_factory
{
  virtual Single_Viterbi_interface*          new_Viterbi (const Single_HMM_scores& hmm, const Named_profile& np) = 0;
  virtual Single_forward_interface*          new_forward (const Single_HMM_scores& hmm, const Named_profile& np) = 0;
  virtual Single_forward_backward_interface* new_forward_backward (const Single_HMM_scores& hmm, const Named_profile& np) = 0;
  virtual ~Single_matrix_factory() { }
};

// templates for Viterbi, Forward, Forward-Backward
//
// the template class EmitScoreCalculator must provide the following methods...
//
//     EmitScoreCalculator (const Single_HMM_scores& hmm, const Named_profile& seq);   // constructor
//     Score calc_emit_score (int state, int seqpos) const;
//     Score accum_emit_score (int posterior_sc, Single_HMM_counts& counts, vector<Metaprob>& metacounts, int state, int seqpos);
//     Score path_score (const vector<int>& state_path) const;   // only used by sample_state_path()
//     static const char* calc_type();  // identifies this EmitScoreCalculator
//
// ...and the following data members...
//
//     Named_profile& seq;            // the sequence
//     const Single_HMM_scores& hmm;  // the HMM topology & scores
//     char char_at_pos (int x);      // display method

// base class for all matrices
template <class EmitScoreCalculator>
class Single_DP_matrix_base : public EmitScoreCalculator, protected Stream_saver
{
public:
  array2d<Score>   cell;        // the matrix itself

  const int        size;        // equals seq.size() + 1
  const int        states;      // excludes start & end states
  int              final_score;
  
  // base constructor
  Single_DP_matrix_base (const Single_HMM_scores& hmm, const Named_profile& seq) :
    EmitScoreCalculator (hmm, seq),
    cell (hmm.states(), seq.size() + 1),
    size (seq.size() + 1),
    states (hmm.states()),
    final_score (-InfinityScore)
    {
      this->alphabet = &hmm.alphabet();
      if (CTAGGING(4,ALLOC)) CL << "Allocated DP matrix:  " << size << "  cells by  " << states << "  states\n";
      if (CTAGGING(-3,DP_HMM)) { CL << "HMM passed to DP matrix:\n"; hmm.show(CL); }
    }

  // display methods
  void show (ostream& o) const;
  sstring hmm_dump() const;

protected:
  // method to calculate cell + transition score
  inline Score calc_transition_score (int src_st, int dest_st, int src_column) const
  { return ScorePMul (cell (src_st, src_column), this->hmm.transition (src_st, dest_st)); }
  // timing stuff
  void report_speed (clock_t start_time, const char* matrix_type);
};

// Viterbi matrix
template <class EmitScoreCalculator>
class Single_Viterbi_matrix_template
  : public Single_DP_matrix_base <EmitScoreCalculator>, public Single_Viterbi_interface
{
public:
  Single_Viterbi_matrix_template (const Single_HMM_scores& hmm, const Named_profile& seq);
  
  Score get_Viterbi_score() const;
  vector<int> optimal_state_path() const;

private:
  // generic loop-over-incoming-states calc method for Viterbi matrix
  inline Score calc_incoming (const int dest_st, const int src_column, const vector<int>& incoming) const
  {
    Score sc = -InfinityScore;
    for_const_contents (vector<int>, incoming, src_st)
      sc = max (sc, this->calc_transition_score (*src_st, dest_st, src_column));
    return sc;
  }
  // we have one calc method for emit states & one for null states
  inline Score calc_score (const int dest_st, const int column, const vector<int>& incoming) const
    {
      Score sc = calc_incoming (dest_st, column - 1, incoming);
      ScorePMulAcc (sc, this->calc_emit_score (dest_st, column - 1));
      return sc;
    }
  inline Score calc_null_score (const int dest_st, const int column, const vector<int>& incoming) const
    {
      const Score sc = calc_incoming (dest_st, column, incoming);
      return sc;
    }
};

// forward matrix
template <class EmitScoreCalculator>
class Single_forward_matrix_template : public Single_DP_matrix_base <EmitScoreCalculator>, public Single_forward_interface
{
public:
  Single_forward_matrix_template (const Single_HMM_scores& hmm, const Named_profile& seq);

  Score get_forward_score() const;
  vector<int> sample_state_path() const;

private:
  // generic loop-over-incoming-states calc method for forward matrix
  inline Score calc_incoming (const int dest_st, const int src_column, const vector<int>& incoming) const
  {
    Score sc = -InfinityScore;
    for_const_contents (vector<int>, incoming, src_st)
      ScorePSumAcc (sc, this->calc_transition_score (*src_st, dest_st, src_column));
    return sc;
  }
  // we have one calc method for emit states & one for null states
  inline Score calc_score (const int dest_st, const int column, const vector<int>& incoming) const
    {
      Score sc = calc_incoming (dest_st, column - 1, incoming);
      ScorePMulAcc (sc, this->calc_emit_score (dest_st, column - 1));
      return sc;
    }
  inline Score calc_null_score (const int dest_st, const int column, const vector<int>& incoming) const
    {
      Score sc = calc_incoming (dest_st, column, incoming);
      return sc;
    }
};

// forward-backward matrix
template <class EmitScoreCalculator>
class Single_forward_backward_matrix_template
  : private Single_forward_matrix_template <EmitScoreCalculator>, public Single_forward_backward_interface
{
public:
  Single_HMM_counts        counts;
  vector<Metaprob>         metacounts;
  Score                    get_forward_score() const;
  const Single_HMM_counts& get_expected_counts() const;
  const vector<Metaprob>&  get_expected_metacounts() const;

  Single_forward_backward_matrix_template (const Single_HMM_scores& hmm, const Named_profile& seq);

private:
  const array2d<Score>& forward_cell;

  // backwards version of calc_transition_score
  inline Score back_calc_transition_score (const int src_st, const int dest_st, const int dest_col) const
  { return ScorePMul (this->cell (dest_st, dest_col), this->hmm.transition (src_st, dest_st)); }

  // backwards version of calc_incoming
  inline void back_calc_outgoing (Score& sc, const int src_st, const int src_col, const int dest_col, const vector<int>& outgoing)
  {
    for_const_contents (vector<int>, outgoing, dest_st)
      {
	const Score trans_sc = back_calc_transition_score (src_st, *dest_st, dest_col);
	ScorePSumAcc (sc, trans_sc);
	counts.transition (src_st, *dest_st) += Score2Prob (ScorePMul3 (forward_cell (src_st, src_col), trans_sc, -this->final_score));
      }
  }

  // wrapper for accum_emit_score
  inline void back_calc_emit_score (Score& back_sc, const int src_st, const int src_col)
  {
    const Score posterior_sc = ScorePMul3 (forward_cell (src_st, src_col), back_sc, -this->final_score);
    const Score emit_sc = this->accum_emit_score (posterior_sc, counts, metacounts, src_st, src_col - 1);
    ScorePMulAcc (back_sc, emit_sc);
  }

  // method to update end counts
  inline void accum_end_score (const int src_st)
  {
    counts.end[src_st] += Score2Prob (ScorePMul3 (forward_cell (src_st, this->size-1), this->hmm.end[src_st], -this->final_score));
  }
  
  // method to calculate start score
  inline void accum_start_score (Score& start_sc, const int dest_st)
  {
    const int dest_col = this->hmm.state_type[dest_st] == Single_HMM_scores::Emit ? 1 : 0;
    if (dest_col < this->size)
      {
	const Score state_start_sc = ScorePMul (this->hmm.start[dest_st], this->cell (dest_st, dest_col));
	counts.start[dest_st] += Score2Prob (state_start_sc - this->final_score);
	ScorePSumAcc (start_sc, state_start_sc);
      }
  }
  
  // backwards score calc method for emit states
  inline Score back_calc_score (int src_st, int column, const vector<int>& outgoing_emit, const vector<int>& outgoing_null)
    {
      Score back_sc = -InfinityScore;

      // transitions to null states
      back_calc_outgoing (back_sc, src_st, column, column, outgoing_null);

      // transitions to emit states
      back_calc_outgoing (back_sc, src_st, column, column + 1, outgoing_emit);

      // emit score
      back_calc_emit_score (back_sc, src_st, column);

      return back_sc;
    }

  // score calc method for null states
  inline Score back_calc_null_score (int src_st, int column, const vector<int>& outgoing_emit, const vector<int>& outgoing_null)
    {
      Score back_sc = -InfinityScore;

      // transitions to null states
      back_calc_outgoing (back_sc, src_st, column, column, outgoing_null);

      // transitions to emit states
      back_calc_outgoing (back_sc, src_st, column, column + 1, outgoing_emit);

      return back_sc;
    }

  // special backwards score calc methods for the final column, placed here for symmetry despite only being called at the end boundary

  // backwards score calc methods for final-column emit states
  inline Score back_calc_end_score (int src_st, const vector<int>& outgoing_null)
    {
      const int column = this->size - 1;
      Score back_sc = this->hmm.end[src_st];

      // transitions to null states
      back_calc_outgoing (back_sc, src_st, column, column, outgoing_null);

      // emit score
      back_calc_emit_score (back_sc, src_st, column);

      // update counts
      accum_end_score (src_st);

      return back_sc;
    }

  // backwards score calc method for final-column null states
  inline Score back_calc_null_end_score (int src_st, const vector<int>& outgoing_null)
    {
      const int column = this->size - 1;
      Score back_sc = this->hmm.end[src_st];

      // transitions to null states
      back_calc_outgoing (back_sc, src_st, column, column, outgoing_null);

      // update counts
      accum_end_score (src_st);

      return back_sc;
    }
};


// methods code

// method to display a matrix
template <class EmitScoreCalculator>
void Single_DP_matrix_base <EmitScoreCalculator>::show (ostream& o) const
{
  save_flags (o);
  
  o << "    ";
  o << "        0:*";
  for (int x = 0; x < size - 1; x++)
    {
      o.width(9);
      right_align (o);
      o << x + 1 << ':' << this->char_at_pos(x);
    }
  o << "\n";
  for (int s = 0; s < states; s++)
    {
      left_align (o);
      o.width(4);
      o << s;
      for (int x = 0; x < size; x++)
	{
	  o << " ";
	  right_align (o);
	  o.width(10);
	  ShowScore (cell(s,x), o);
	}
      o << "\n";
    }
  o << "Final score: ";
  ShowScore (final_score, o);
  o << "\n";

  restore_flags (o);
}

// debugging HMM + matrix output
template <class EmitScoreCalculator>
sstring Single_DP_matrix_base <EmitScoreCalculator>::hmm_dump() const
{
  // until sstring::operator<<(int) bug is fixed, we use the following hack
  ostream& dump = CLOGERR;
  //  sstring dump;
  dump << "Single HMM:\n";
  this->hmm.show (dump);
  dump << "Single HMM DP matrix:\n";
  show (dump);
  //  return dump;
  return sstring();
}

// timing method
template <class EmitScoreCalculator>
void Single_DP_matrix_base <EmitScoreCalculator>::report_speed (clock_t start_time, const char* matrix_type)
{
  if (CTAGGING (2, CELLS_PER_SEC))
    {
      clock_t end_time = clock();
      const int n_cells = cell.xsize() * cell.ysize();
      const double n_secs = ((double) (end_time - start_time)) / (double) CLOCKS_PER_SEC;
      const double cells_per_sec = ((double) n_cells) / n_secs;
      CL << matrix_type << " " << this->calc_type() << " fill: " << n_cells << " cells in " << n_secs << " seconds (" << cells_per_sec << " cells/sec)\n";
    }
}

// constructor for a Viterbi matrix
template <class EmitScoreCalculator>
Single_Viterbi_matrix_template <EmitScoreCalculator>::Single_Viterbi_matrix_template
(const Single_HMM_scores& hmm, const Named_profile& seq)
  : Single_DP_matrix_base <EmitScoreCalculator> (hmm, seq)
{
  const vector<int> null_states = hmm.null_states();
  const vector<int> emit_states = hmm.emit_states();
  
  const vector<vector<int> > incoming = hmm.incoming_states();
  
  clock_t start_time = clock();
  
  // initialise zeroth column
  //
  for_const_contents (vector<int>, emit_states, s)
    this->cell(*s,0) = -InfinityScore;
  
  for_const_contents (vector<int>, null_states, s)
    this->cell(*s,0) = max (hmm.start[*s], calc_null_score (*s, 0, incoming[*s]));
  
  // initialise first column
  //
  if (this->size > 1)
    {
      for_const_contents (vector<int>, emit_states, s)
	this->cell(*s,1) = max (ScorePMul (this->hmm.start[*s], this->calc_emit_score (*s, 0)), calc_score (*s, 1, incoming[*s]));
      
      for_const_contents (vector<int>, null_states, s)
	this->cell(*s,1) = calc_null_score (*s, 1, incoming[*s]);
    }
  
  // fill the matrix
  //
  for (int x = 2; x < this->size; x++)
    {
      for_const_contents (vector<int>, emit_states, s)
	this->cell(*s,x) = calc_score (*s, x, incoming[*s]);
      
      for_const_contents (vector<int>, null_states, s)
	this->cell(*s,x) = calc_null_score (*s, x, incoming[*s]);
    }
  
  // get the final score
  //
  for (int s = 0; s < this->states; s++)
    this->final_score = max (this->final_score, ScorePMul (this->cell(s,this->size-1), this->hmm.end[s]));
  if (this->size == 1)
    this->final_score = max (this->final_score, this->hmm.start_to_end());
  
  this->report_speed (start_time, "Viterbi");
  
  if (CTAGGING(-1,DP_MATRIX VITERBI_MATRIX)) { CL << "Viterbi matrix:\n"; this->show(CL); }

  CTAG(4,DP_SCORE VITERBI_SCORE) << "Viterbi log-likelihood is " << Score2Bits(this->final_score) << " bits\n";
}


// accessor for Viterbi score
template <class EmitScoreCalculator>
int Single_Viterbi_matrix_template <EmitScoreCalculator>::get_Viterbi_score() const
{ return this->final_score; }

// traceback method for Viterbi
template <class EmitScoreCalculator>
vector<int> Single_Viterbi_matrix_template <EmitScoreCalculator>::optimal_state_path() const
{
  if (this->final_score == -InfinityScore)
    THROWEXPR ("Final score is -infinity; traceback likely to break");

  vector<int> state_path;
  
  if (this->size > 1 || this->final_score != this->hmm.start_to_end())
    {
      int x = this->size - 1;
      int current_state = -1;
      for (int s = 0; s < this->hmm.states(); s++)
	if (ScorePMul (this->cell(s,x), this->hmm.end[s]) == this->final_score) { current_state = s; break; }
      if (current_state < 0) THROWEXPR (this->hmm_dump() << "Single HMM: traceback failed at first hurdle");

      while (1)
	{
	  state_path.push_back (current_state);

	  Score current_score = this->cell (current_state, x);
	  if (this->hmm.state_type[current_state] == Single_HMM_scores::Emit)
	    {
	      if (--x < 0) THROWEXPR (this->hmm_dump() << "Single HMM: traceback out of bounds at x=" << x << ", state=" << current_state);
	      current_score -= this->calc_emit_score (current_state, x);
	    }
	  if (x == 0 && current_score == this->hmm.start[current_state]) break;
	  
	  int prev_state = -1;
	  for (int s = 0; s < this->hmm.states(); s++)
	    if (current_score == this->calc_transition_score (s, current_state, x)) { prev_state = s; break; }
	  if (prev_state < 0) THROWEXPR (this->hmm_dump() << "Single HMM: traceback failed at x=" << x << ", state=" << current_state);

	  current_state = prev_state;
	}
  
      reverse (state_path.begin(), state_path.end());
    }

  return state_path;
}

// constructor for a forward matrix
template <class EmitScoreCalculator>
Single_forward_matrix_template <EmitScoreCalculator>::Single_forward_matrix_template
(const Single_HMM_scores& hmm, const Named_profile& seq)
  : Single_DP_matrix_base <EmitScoreCalculator> (hmm, seq)
{
  const vector<int> null_states = hmm.null_states();
  const vector<int> emit_states = hmm.emit_states();
  
  const vector<vector<int> > incoming = hmm.incoming_states();
  
  clock_t start_time = clock();
  
  // initialise zeroth column
  //
  for_const_contents (vector<int>, emit_states, s)
    this->cell(*s,0) = -InfinityScore;
  
  for_const_contents (vector<int>, null_states, s)
    this->cell(*s,0) = ScorePSum (hmm.start[*s], calc_null_score (*s, 0, incoming[*s]));
  
  // initialise first column
  //
  if (this->size > 1)
    {
      for_const_contents (vector<int>, emit_states, s)
	this->cell(*s,1) = ScorePSum (ScorePMul (hmm.start[*s], this->calc_emit_score (*s, 0)), calc_score (*s, 1, incoming[*s]));
      
      for_const_contents (vector<int>, null_states, s)
	this->cell(*s,1) = calc_null_score (*s, 1, incoming[*s]);
    }
  
  // fill the matrix
  //
  for (int x = 2; x < this->size; x++)
    {
      for_const_contents (vector<int>, emit_states, s)
	this->cell(*s,x) = calc_score (*s, x, incoming[*s]);
      
      for_const_contents (vector<int>, null_states, s)
	this->cell(*s,x) = calc_null_score (*s, x, incoming[*s]);
    }
  
  // get the final score
  //
  for (int s = 0; s < this->states; s++)
    ScorePSumAcc (this->final_score, ScorePMul (this->cell(s,this->size-1), this->hmm.end[s]));
  if (this->size == 1) ScorePSumAcc (this->final_score, hmm.start_to_end());
  
  this->report_speed (start_time, "Forward");
  
  if (CTAGGING(-1,DP_MATRIX FORWARD_MATRIX)) { CL << "Forward matrix:\n"; this->show(CL); }

  CTAG(4,DP_SCORE FORWARD_SCORE) << "Forward log-likelihood is " << Score2Bits(this->final_score) << " bits\n";
}


// accessor for forward score
template <class EmitScoreCalculator>
Score Single_forward_matrix_template <EmitScoreCalculator>::get_forward_score() const
{ return this->final_score; }

// sampling traceback for forward matrix
template <class EmitScoreCalculator>
vector<int> Single_forward_matrix_template <EmitScoreCalculator>::sample_state_path() const
{
  vector<int> state_path;

  vector<int> sc (this->hmm.states() + 1);
  sc[0] = this->size == 1 ? this->hmm.start_to_end() : -InfinityScore;
  int x = this->size - 1;
  for (int s = 0; s < this->hmm.states(); s++)
    sc[s+1] = ScorePMul (this->cell(s,x), this->hmm.end[s]);
  int current_state = Rnd::choose (Score2ProbVecNorm (sc)) - 1;

  if (current_state > -1)
    {
      while (1)
	{
	  state_path.push_back (current_state);
	  if (this->hmm.state_type[current_state] == Single_HMM_scores::Emit)
	    if (--x < 0) break;

	  for (int s = 0; s < this->hmm.states(); s++)
	    sc[s] = this->calc_transition_score (s, current_state, x);
	  vector<double> weight = Score2ProbVecNorm(sc);
	  current_state = Rnd::choose(weight);
	}
  
      reverse (state_path.begin(), state_path.end());
    }

  if (CTAGGING(5,DP_SCORE))
    CL << "Sampled single-HMM path posterior log-probability = " << Score2Bits (this->path_score (state_path) - this->final_score) << " bits\n";
  return state_path;
}

// constructor for a forward-backward matrix
template <class EmitScoreCalculator>
Single_forward_backward_matrix_template <EmitScoreCalculator>::Single_forward_backward_matrix_template
(const Single_HMM_scores& hmm, const Named_profile& seq) :
  Single_forward_matrix_template <EmitScoreCalculator> (hmm, seq),
  counts (hmm),
  forward_cell (this->cell)
{
  if (hmm.uses_metascores())   // don't initialise the metacounts vector unless we have to. things are slow enough already
    metacounts = vector<Metaprob> (hmm.max_metascore_idx() + 1, Metaprob (seq.size(), (Prob) 0));

  counts.log_likelihood = Score2Nats(this->final_score);

  if (this->final_score > -InfinityScore)  // skip the backward phase if probability is zero
    {
      // we want the null states sorted in *reverse* topological order
      vector<int> _null_states = hmm.null_states();
      reverse (_null_states.begin(), _null_states.end());
      
      // make const references for null_states & emit_states
      const vector<int>& null_states = _null_states;
      const vector<int> emit_states = hmm.emit_states();
      
      // get incoming & outgoing staets
      const vector<vector<int> > outgoing_emit = hmm.selected_outgoing_states (emit_states);
      const vector<vector<int> > outgoing_null = hmm.selected_outgoing_states (null_states);
      
      clock_t start_time = clock();
      
      // do final column
      //
      for_const_contents (vector<int>, null_states, s)
	this->cell (*s, this->size - 1) = back_calc_null_end_score (*s, outgoing_null[*s]);
      
      if (this->size > 1)
	for_const_contents (vector<int>, emit_states, s)
	  this->cell (*s, this->size - 1) = back_calc_end_score (*s, outgoing_null[*s]);
      else
	for_const_contents (vector<int>, emit_states, s)
	  this->cell (*s, this->size - 1) = -InfinityScore;
      
      // fill the matrix
      //
      for (int x = this->size - 2; x >= 1; --x)
	{
	  for_const_contents (vector<int>, null_states, s)
	    this->cell (*s, x) = back_calc_null_score (*s, x, outgoing_emit[*s], outgoing_null[*s]);
	  
	  for_const_contents (vector<int>, emit_states, s)
	    this->cell (*s, x) = back_calc_score (*s, x, outgoing_emit[*s], outgoing_null[*s]);
	}
      
      // do zeroth column
      //
      if (this->size > 1)
	for_const_contents (vector<int>, null_states, s)
	  this->cell (*s, 0) = back_calc_null_score (*s, 0, outgoing_emit[*s], outgoing_null[*s]);
      
      // do start transitions
      //
      int start_sc = -InfinityScore;
      for (int s = 0; s < this->states; s++)
	accum_start_score (start_sc, s);
      
      if (this->size == 1)
	counts.start_to_end() = Score2Prob (ScorePMul (hmm.start_to_end(), -this->final_score));
      
      this->report_speed (start_time, "Backward");
      
      if (CTAGGING(-1,DP_MATRIX BACKWARD_MATRIX)) { CL << "Backward matrix:\n"; this->show(CL); }
      CTAG(3,DP_SCORE BACKWARD_SCORE) << "Backward log-likelihood is " << Score2Bits(start_sc) << " bits\n";
    }
  
  if (CTAGGING(4,DP_COUNTS)) { CL << "Single-HMM counts from forward-backward matrix:\n"; counts.show(CL); show_metacounts(CL); }
}

// accessor for forward score in a forward-backward matrix
template <class EmitScoreCalculator>
Score Single_forward_backward_matrix_template <EmitScoreCalculator>::get_forward_score() const
{ return this->final_score; }

// accessor for HMM counts in a forward-backward matrix
template <class EmitScoreCalculator>
const Single_HMM_counts& Single_forward_backward_matrix_template <EmitScoreCalculator>::get_expected_counts() const
{ return counts; }

// accessor for metacounts in a forward-backward matrix
template <class EmitScoreCalculator>
const vector<Metaprob>& Single_forward_backward_matrix_template <EmitScoreCalculator>::get_expected_metacounts() const
{ return metacounts; }

#endif
