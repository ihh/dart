#ifndef RATE_MATRIX_INCLUDED
#define RATE_MATRIX_INCLUDED

#include <iostream>
#include "util/sstring.h"
#include <vector>
#include "newmat/newmat.h"
#include "util/dexception.h"
#include "tree/substitution_matrix_factory.h"
#include "util/opts_list.h"

class Rate_matrix : public Substitution_matrix_factory
{
 private:
  Simple_alphabet a;
  vector<double>  frequency;    // equilibrium frequencies of residues
  vector<double>  sqrt_freq;    // square roots of residue frequencies
  DiagonalMatrix  D;            // diagonal elements are eigenvalues of symmetrised rate matrix
  Matrix          U;            // column vectors are eigenvectors of symmetrised rate matrix

 public:
  typedef Standard_exception Format_exception;

  array2d<double> create_conditional_substitution_matrix (double time);
  array2d<double> differentiate_conditional_substitution_matrix (double time);
  vector<double>  create_prior() { return frequency; }
  const Alphabet& alphabet() { return a; }

  void read (istream& in);

  Rate_matrix() : a("user",""), D(), U() {}
  Rate_matrix (const Alphabet& a, const vector<double>& freq, const Matrix& evec_cols, const vector<double>& eval);

  static bool rate_matrix_help (Opts_list* ol);
};

#endif
