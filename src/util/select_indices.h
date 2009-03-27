#ifndef SELECT_INDICES_INCLUDED
#define SELECT_INDICES_INCLUDED

#include <vector>
#include "util/macros.h"

template <class T>
class select_indices_t : public vector<T>
{
 public:
  select_indices_t (const vector<T>& source_vector, const vector<int>& indices) : vector<T> ()
    {
      reserve (indices.size());
      for_const_contents (vector<int>, indices, i) push_back ((vector<T>&) source_vector) [*i];
    }
};

template <class T>
select_indices_t<T>
select_indices (const vector<T>& source_vector, const vector<int>& indices)
{
  return select_indices_t (source_vector, indices);
}

#endif
