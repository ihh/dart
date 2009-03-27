#ifndef DISTMAT_INCLUDED
#define DISTMAT_INCLUDED

#include "util/array2d.h"
#include "util/math_fn.h"
#include "seq/alignment.h"

// function to return a distance for two integer indices
struct Dist_func : binary_function<int,int,double>
{
  virtual double operator() (int, int);
  virtual ~Dist_func();
};

// abstract factory for making Dist_func's from multiple alignments
struct Dist_func_factory
{
  virtual Dist_func* create_dist_func (const Alignment& align) = 0;
  virtual ~Dist_func_factory() { }
};

// class to return a list of the N indices (from 1 to M) nearest to a given index k
struct Nearest : vector<int>
{
  Nearest (Dist_func& f, int N, int M, int k);
};

// general distance matrix class
struct Dist_matrix : array2d<double>
{
  // constructor
  Dist_matrix (int size, Dist_func& f);
  // method to return candidate similar indices
  vector<int> similar (int idx, int n_sim) const;  // returns indices of n_sim entries most similar to idx
  int most_similar (int idx) const;  // wrapper for similar() with n_sim=1
};

#endif /* DISTMAT_INCLUDED */
