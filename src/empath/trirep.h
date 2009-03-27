#ifndef TRIREP_INCLUDED
#define TRIREP_INCLUDED

#include "empath/trainer.h"

struct Trirep_set
{
  int max_len;
  const Alphabet& alphabet;
  set<Digitized_biosequence> reps;  // sorted by rep length
  bool contains (Digitized_biosequence& dsq);  // checks for cyclic variants
  void add_degenerates();  // adds two-base degenerate chars (assumes DNA)
  Trirep_set (int max_motif_len, int max_palindrome_len, const Alphabet& alphabet = DNA_alphabet);
  friend ostream& operator<< (ostream& o, const Trirep_set& trirep_set);
};

struct Trirep_model : public Trainable
{
  const Trirep_set& trirep_set;
  int min_reps;
  int shape;
  static int n_states (const Trirep_set& trirep_set, int units);
  Boolean_group match;
  Boolean_group ins;
  Boolean_group del;
  struct Repdata
  {
    Weight_profile rep_profile;
    vector<int> start_state;
    PGroup start;
    Boolean_group loop;
    Repdata (const Digitized_biosequence& dsq, const Alphabet& alph);
  };
  vector<Repdata> repdata;
  vector<sstring> repname;
  PGroup which_rep;
  // constructor
  Trirep_model (Trirep_set& trirep_set, PScores& pscore, int min_reps, int shape, double p_mismatch, double p_ins, double p_del);
  // I/O
  void write_repdist (ostream& out) const;
};

#endif  /* TRIREP_INCLUDED */
