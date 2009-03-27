#ifndef PFOLD_INCLUDED
#define PFOLD_INCLUDED

#include "hmm/transmat.h"
#include "ecfg/ecfg.h"

// PFOLD rate matrices
struct PFOLD_loop_matrix : EM_matrix
{
  PFOLD_loop_matrix();
};

struct PFOLD_stem_matrix : EM_matrix
{
  PFOLD_stem_matrix();
};

// PFOLD grammar (with a little tweaking)
// States: S F L B(=Bifurcation) U(=Unpaired)
//
// Start -> S
// S  -> B (0.868534) | L (0.131466)
// F' -> F (0.787640) | B (0.212360)
// L  -> U (0.894603) | F (0.105397)
// B  -> L S
// U' -> End
// U  -> x U'
// F  -> x F' y

struct PFOLD_state_enum : Grammar_state_enum
{
  enum { pfoldNull = UndefinedState - 1 };
  enum PFOLD_state { pfoldS = 0, pfoldF = 1, pfoldL = 2, pfoldB = 3, pfoldU = 4, pfoldTotalStates = 5 };
  // display methods
  static const char* PFOLD_state_string (PFOLD_state state);
};

struct PFOLD_ECFG : ECFG_scores, PFOLD_state_enum
{
  // constructor
  PFOLD_ECFG (const char* tag = "PFOLD");
  // EM_matrix accessors
  EM_matrix_base& stem_matrix() { return *matrix_set.chain[state_info[pfoldF].matrix].matrix; }
  EM_matrix_base& loop_matrix() { return *matrix_set.chain[state_info[pfoldU].matrix].matrix; }
  // display method
  sstring desc (int state) const { return sstring (PFOLD_state_string ((PFOLD_state) state)); }
};

// PFOLD grammar with pseudo-indel states
struct IFOLD_ECFG : PFOLD_ECFG
{
  IFOLD_ECFG();
};

// Null model grammar for any alphabet
struct Null_ECFG : ECFG_scores
{
  // constructor
  Null_ECFG (const Alphabet& alph, bool rev = true);
};

// Null model grammar for RNA
struct Null_RNA_ECFG : Null_ECFG
{
  // constructor
  Null_RNA_ECFG (bool rev = true);
};

// Null model grammar for RNA, with indels
struct Null_indel_RNA_ECFG : Null_RNA_ECFG
{
  // constructor
  Null_indel_RNA_ECFG();
};

// Nearest-neighbor model
struct Nearest_neighbor_ECFG : ECFG_scores
{
  // pair alphabet
  Alphabet pair_alphabet;
  // constructor
  Nearest_neighbor_ECFG();
};

// Codon model
struct Codon_ECFG : ECFG_scores
{
  // triplet alphabet
  Alphabet triplet_alphabet;
  // constructor
  Codon_ECFG (bool both_strands = true);
};

#endif /* PFOLD_INCLUDED */
