#ifndef ENERGY21_INCLUDED
#define ENERGY21_INCLUDED

#include "util/score.h"
#include "util/array2d.h"

// RNA duplex lookup table

struct Duplex_lookup
{
  // duplex symbol <--> index lookup tables
  enum { Symbols = 4, Duplexes = 6 };
  vector<int> duplex_sym5;  // A, C, G, U, G, U
  vector<int> duplex_sym3;  // U, G, C, A, U, G
  array2d<int> duplex_idx;  // (A,U)=0, (C,G)=1 ... (U,G)=5, (A,A)=-1 ...
  Duplex_lookup();
};

// RNA base pair stacking energies
// Parent class -- does nothing -- just a container
struct RNA_energy
{
  // energies
  array2d<array2d<Score> > stacked;
  array2d<array2d<Score> > hairpin;
  array2d<array2d<Score> > interior;
  vector<Score> interior_loop;
  vector<Score> bulge_loop;
  vector<Score> hairpin_loop;
  // constructor
  RNA_energy();
  // helper method to do an affine fit (y = A + B*x) to a segment of a curve
  static void affine_fit (const vector<Score>& loop, int min_len, int max_len, Score& A_ret, Score& B_ret);
};

// RNA base pair stacking energies at 21 degrees Centigrade -- courtesy of Turner lab
// Integer multiplier is -1000
// I guess the energies are in Kcal/mole
struct Energy21 : RNA_energy
{
  Energy21();
};

// tetraloops (-2.1 kcal/mole bonus)
// GAAA GCAA GGAA GUAA GAGA GCGA GUGA
// UUCG UACG UCCG
// UUUA CUUG AUUU

#endif /* ENERGY21_INCLUDED */
