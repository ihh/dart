#include <fstream>
#include <deque>
#include "hmm/pairhmm.h"
#include "tree/substitution_matrix_factory.h"

// Wildcard_transducer
// Pair HMM without emissions
// Used to represent a Handel indel model
struct Wildcard_transducer: Pair_HMM_scores
{
  // constructor
  Wildcard_transducer (int states) : Pair_HMM_scores (states, &Dummy_alphabet) { }
};

// Handel_model
// A Substitution_matrix_factory that also represents indel probabilities
// using a Wildcard_transducer
struct Handel_model : Substitution_matrix_factory
{
  // gamma()
  // Returns the geometric probability distribution parameter, gamma.
  //   P(sequence length at equilibrium = L) = (1-gamma) * gamma^L
  virtual Prob gamma() = 0;

  // transducer(t)
  // Returns the conditionally normalised Pair HMM for branch length t (no emissions)
  virtual Wildcard_transducer transducer (double t) = 0;

  // method to make conditional Pair HMM with emissions
  Pair_HMM_scores conditional_pair_HMM (double t);

  // method to make joint Pair HMM with emissions
  Pair_HMM_scores joint_pair_HMM (double t);
};
