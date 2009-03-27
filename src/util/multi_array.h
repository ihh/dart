// multidimensional array template

#ifndef MULTI_ARRAY_INCLUDED
#define MULTI_ARRAY_INCLUDED

#include <stdarg.h>
#include <vector>
#include <functional>
#include "util/counter.h"
#include "util/dexception.h"
#include "util/logfile.h"

template <class T>
class multi_array
{
 public:
  typedef T value_type;
  typedef vector<int> key_type;
  typedef Counter counter_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;

  typedef ptrdiff_t difference_type;
  typedef size_t size_type;

  typedef typename vector<value_type>::iterator iterator;

 private:
  vector<int>         _dimensions;
  vector<value_type>  _data;
  
  counter_type        _begin;
  counter_type        _end;

  static difference_type _calculate_volume (const key_type& x)
    { key_type z (x.size()); return counter_type (x,z,x) - counter_type (z,z,x); }

 public:
  
  const key_type& dim() const { return (const key_type&) _dimensions; }
  int rank() const { return _dimensions.size(); }

  reference operator[](const key_type& k) { return _data[counter_type (k, _begin, _end) - _begin]; }
  const_reference operator[](const key_type& k) const { return _data[counter_type (k, _begin, _end) - _begin]; }

  reference operator()(int x0, ...)
    {
      int index = 0;
  /*
   * the following lines are commented out because gcc doesn't seem to like the va_start macros
   *
      va_list argp;
      va_start (argp, x0);
      for (int i = 0; i < rank(); i++) index += _begin.multiplier()[i] * (va_arg (argp, int) - _begin[i]);
      va_end (argp);
  */
      THROW Standard_exception ("multi_array::operator() not implemented due to bug in gcc");
      return _data[index];
    }

  iterator begin() { return _data.begin(); }
  iterator end() { return _data.end(); }

  counter_type counter_begin() const { return _begin; }
  counter_type counter_end() const { return _end; }

  // swap method
  multi_array<T>& swap (multi_array<T>& m)
  {
    _dimensions.swap(m._dimensions);
    _data.swap(m._data);
    _begin.swap(m._begin);
    _end.swap(m._end);
    return *this;
  }

  // constructors
  multi_array () : _dimensions(), _data(), _begin(), _end() {}

  multi_array (const multi_array& m) :
    _dimensions(m._dimensions), _data(m._data), _begin(m._begin), _end(m._end) {}
  
  multi_array (const key_type& dimensions, const value_type& t) :
    _dimensions(dimensions), _data(_calculate_volume(dimensions), t),
    _begin(counter_type::create_begin(dimensions)), _end(counter_type::create_end(dimensions))
    {}
  
  multi_array (const key_type& dimensions) :
    _dimensions(dimensions), _data(_calculate_volume(dimensions)),
    _begin(counter_type::create_begin(dimensions)), _end(counter_type::create_end(dimensions))
    {}
  
  multi_array (int rank, int length, const value_type& t) :
    _dimensions(vector<int>(rank, length)), _data(_calculate_volume(key_type(rank, length)), t),
    _begin(counter_type::create_begin(vector<int>(rank, length))), _end(counter_type::create_end(vector<int>(rank, length)))
    {}
  
  multi_array (int rank, int length) :
    _dimensions(vector<int>(rank, length)), _data(_calculate_volume(key_type(rank, length))),
    _begin(counter_type::create_begin(vector<int>(rank, length))), _end(counter_type::create_end(vector<int>(rank, length)))
    {}
  
};

#endif
