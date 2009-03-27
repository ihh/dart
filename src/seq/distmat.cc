#include <numeric>
#include "seq/distmat.h"

double Dist_func::operator() (int i, int j) { return 0; }
Dist_func::~Dist_func() { }

Nearest::Nearest (Dist_func& f, int N, int M, int k) : vector<int> (M)
{
  vector<double> dist (M);
  for (int i = 0; i < M; ++i) dist[i] = f(k,i);
  Schwartzian<double> index_by_distance (dist);
  for (int i = 0; i < M; ++i)
    (*this)[i] = i;
  partial_sort (begin(), begin() + N, end(), index_by_distance);
  erase (begin() + N, end());
}

Dist_matrix::Dist_matrix (int size, Dist_func& f) : array2d<double> (size, size, 0.0)
{
  for (int i = 0; i < size; ++i)
    {
      (*this)(i,i) = 0;
      for (int j = i + 1; j < size; ++j)
	{
	  CTAG(3,DISTMATRIX) << "Estimating distances of sequences " << i << " and " << j << "\n";
	  (*this)(i,j) = (*this)(j,i) = f(i,j);
	}
    }
  CTAG(2,DISTMATRIX) << "Distance matrix:\n" << *this;
}

vector<int> Dist_matrix::similar (int idx, int n_sim) const
{
  vector<double> dist (xsize());
  for (int i = 0; i < xsize(); ++i) dist[i] = (*this)(idx,i);
  Schwartzian<double> index_by_distance (dist);
  vector<int> order (xsize());
  for (int i = 0; i < (int) order.size(); ++i)
    order[i] = i;
  partial_sort (order.begin(), order.begin() + n_sim, order.end(), index_by_distance);
  order.erase (order.begin() + n_sim, order.end());
  return order;
}

int Dist_matrix::most_similar (int idx) const
{
  const vector<int> most = similar (idx, 1);
  return most[0];
}
