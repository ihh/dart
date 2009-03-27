#ifndef SUB_DISTMAT_INCLUDED
#define SUB_DISTMAT_INCLUDED

#include "seq/distmat.h"
#include "tree/substitution_matrix_factory.h"

// substitution counts
struct Substitution_counts
{
  // the set of all residue pairs & their frequencies
  typedef pair<int,int> Residue_pair;
  typedef map<Residue_pair,double> Residue_pair_count;
  Residue_pair_count res_pair_count;
  // constructor
  Substitution_counts (const Alignment& align, int row1, int row2);
  // show method
  void show (ostream& out) const;
};

// substitution-based time log-likelihood function
struct Subst_log_like : unary_function <double, Loge>
{
  // data
  Substitution_matrix_factory& submat_factory;
  const Substitution_counts& subst_counts;
  // constructor
  Subst_log_like (Substitution_matrix_factory& submat_factory, const Substitution_counts& subst_counts);
  // eval method
  Loge operator() (double t);
};

// time derivative of substitution log-likelihood function
struct Subst_log_like_dt : unary_function <double, double>
{
  // data
  Substitution_matrix_factory& submat_factory;
  const Substitution_counts& subst_counts;
  // constructor
  Subst_log_like_dt (Substitution_matrix_factory& submat_factory, const Substitution_counts& subst_counts);
  // eval method
  double operator() (double t);
};

// concrete Dist_func_factory for substitution models
struct Subst_dist_func_factory : Dist_func_factory
{
  // Dist_func class
  struct Subst_dist_func : Dist_func
  {
    // data
    const Subst_dist_func_factory& factory;
    const Alignment& align;
    // constructor
    Subst_dist_func (const Subst_dist_func_factory& factory, const Alignment& align) : factory (factory), align (align) { }
    // method to evaluate pairwise distance
    double operator() (int i, int j);
  };
  // data
  Substitution_matrix_factory& submat_factory;
  double tres;
  double tmax;
  // constructor
  Subst_dist_func_factory (Substitution_matrix_factory& submat_factory,
			   double tres = .01, double tmax = 10.);
  // factory method
  Subst_dist_func* create_dist_func (const Alignment& align);
};

#endif /* SUB_DISTMAT_INCLUDED */
