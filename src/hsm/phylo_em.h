#ifndef PHYLO_EM_INCLUDED
#define PHYLO_EM_INCLUDED

#include "util/math_fn.h"
#include "util/array2d.h"

// a "clean" implementation of phylo-EM using templates to allow for complex numbers.
// given a transition i-->j, eigenvalues & left/right eigenvectors of a rate matrix, and a time parameter t,
// returns a matrix whose on-diagonal elements correspond to wait times, and off-diagonals are usage counts.
template <class NUM>
struct Phylo_EM
{
  // data
  vector<NUM> mu;  // mu[k] = k'th eigenvalue
  array2d<NUM> U;  // U(i,k) = k'th right eigenvector, i'th entry
  array2d<NUM> Uinv;  // Uinv(k,i) = k'th left eigenvector, i'th entry

  // constructor
  Phylo_EM (vector<NUM> mu, array2d<NUM> U, array2d<NUM> Uinv)
    : mu(mu), U(U), Uinv(Uinv)
  { }

  // size accessor
  inline int states() const { return (int) mu.size(); }

  // method
  array2d<NUM> get_counts (int src, int dest, double T);
};

template<class NUM>
array2d<NUM> Phylo_EM<NUM>::get_counts (int a, int b, double T)
{
  // calculate Q
  array2d<NUM> Q (states(), states(), (NUM) 0);
  for (int i = 0; i < states(); ++i)
    for (int j = 0; j < states(); ++j)
      for (int k = 0; k < states(); ++k)
	Q(i,j) += U(i,k) * mu[k] * Uinv(k,j);

  // log Q
  if (CTAGGING(-1,PHYLO_EM_DEBUG))
    {
      CL << "Q matrix calculated by Phylo_EM::get_counts:\n";
      for (int i = 0; i < states(); ++i)
	{
	  for (int j = 0; j < states(); ++j)
	    CL << Q(i,j) << ' ';
	  CL << '\n';
	}
    }

  // calculate exp(mu[k]*T) for all k
  vector<NUM> exp_mu_T (states());
  for (int k = 0; k < states(); ++k)
    exp_mu_T[k] = exp (mu[k] * T);

  // calculate M_ab
  NUM M_ab = 0;
  for (int k = 0; k < states(); ++k)
    M_ab += U(a,k) * exp_mu_T[k] * Uinv(k,b);

  // calculate J
  array2d<NUM> J (states(), states());
  for (int k = 0; k < states(); ++k)
    for (int l = 0; l < states(); ++l)
      if (k == l || abs (mu[k] - mu[l]) < TINY)
	J(k,l) = T * exp_mu_T[k];
      else
	J(k,l) = (exp_mu_T[k] - exp_mu_T[l]) / (mu[k] - mu[l]);

  // calculate C
  array2d<NUM> C (states(), states());
  for (int i = 0; i < states(); ++i)
    for (int j = 0; j < states(); ++j)
      {
	NUM I_ab_ij = 0;
	for (int k = 0; k < states(); ++k)
	  for (int l = 0; l < states(); ++l)
	    I_ab_ij += U(a,k) * Uinv(k,i) * U(j,l) * Uinv(l,b) * J(k,l);
	C(i,j) = I_ab_ij / M_ab;
	// multiply off-diagonal elements by Q_ij
	if (i != j)
	  C(i,j) *= Q(i,j);
      }

  return C;
}

#endif /* PHYLO_EM_INCLUDED */
