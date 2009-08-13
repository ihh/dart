#ifndef MATRIX_ADAPTOR_INCLUDED
#define MATRIX_ADAPTOR_INCLUDED

#include <vector>

#include "newmat/newmat.h"
#include "newmat/newmatap.h"

#include "util/multi_array.h"
#include "util/array2d.h"
#include "util/dexception.h"
#include "util/logfile.h"


// Jacobi wrapper

void Jacobi_catch(const SymmetricMatrix&, DiagonalMatrix&, Matrix&);

// NB that array2d<> and multi_array<> use (x,y) addressing, whereas Matrix uses (row,col).
// Converting between these types thus looks like transposition.
//

class Matrix_to_multi_array_adaptor : public multi_array<double>
{
 public:
  Matrix_to_multi_array_adaptor (const Matrix& matrix) : multi_array<double> (vector<int> (matrix.Ncols(), matrix.Nrows()))
    {
      for (int i = 0; i < matrix.Nrows(); i++)
	for (int j = 0; j < matrix.Ncols(); j++)
	  (*this) (j, i) = ((Matrix&) matrix) (i+1, j+1);
    }
};

class multi_array_to_Matrix_adaptor : public Matrix
{
 public:
  multi_array_to_Matrix_adaptor (const multi_array<double>& array) :
    Matrix (array.rank() == 2 ? array.dim()[1] : 1, array.rank() == 2 ? array.dim()[0] : 1)
    {
      if (array.rank() != 2) THROWEXPR ("Dimensionality mismatch");
      for (int i = 0; i < array.dim()[1]; i++)
	for (int j = 0; j < array.dim()[0]; j++)
	  (*this) (i+1, j+1) = ((multi_array<double>&) array) (j, i);
    }
};

class Matrix_to_array2d_adaptor : public array2d<double>
{
 public:
  Matrix_to_array2d_adaptor (const Matrix& matrix) : array2d<double> (matrix.Ncols(), matrix.Nrows())
    {
      for (int i = 0; i < matrix.Nrows(); i++)
	for (int j = 0; j < matrix.Ncols(); j++)
	  (*this) (j, i) = ((Matrix&) matrix) (i+1, j+1);
    }
};

class array2d_to_Matrix_adaptor : public Matrix
{
 public:
  array2d_to_Matrix_adaptor (const array2d<double>& array) :
    Matrix (array.ysize(), array.xsize())
    {
      for (int i = 0; i < array.ysize(); i++)
	for (int j = 0; j < array.xsize(); j++)
	  (*this) (i+1, j+1) = ((array2d<double>&) array) (j, i);
    }
};

class vector_to_RowVector_adaptor : public RowVector
{
 public:
  vector_to_RowVector_adaptor (const vector<double>& v) : RowVector (v.size())
    { for (int i = 0; i < (int) v.size(); i++) (*this) (i+1) = ((vector<double>&) v) [i]; }
};

class vector_to_ColumnVector_adaptor : public ColumnVector
{
 public:
  vector_to_ColumnVector_adaptor (const vector<double>& v) : ColumnVector (v.size())
    { for (int i = 0; i < (int) v.size(); i++) (*this) (i+1) = ((vector<double>&) v) [i]; }
};

class RowVector_to_vector_adaptor : public vector<double>
{
 public:
  RowVector_to_vector_adaptor (const RowVector& v) : vector<double> (v.Ncols())
    { for (int i = 0; i < v.Ncols(); i++) (*this) [i] = ((RowVector&) v) (i+1); }
};

class ColumnVector_to_vector_adaptor : public vector<double>
{
 public:
  ColumnVector_to_vector_adaptor (const ColumnVector& v) : vector<double> (v.Nrows())
    { for (int i = 0; i < v.Nrows(); i++) (*this) [i] = ((ColumnVector&) v) (i+1); }
};

#endif
