// multi-digit counter class and iterator adaptors
// counter is little-endian, ie LSD-first.

#ifndef COUNTER_INCLUDED
#define COUNTER_INCLUDED

#include <iostream>
#include <vector>
#include <iterator>
#include "util/dexception.h"
#include "util/logfile.h"

// class for iterating through a tuple of bounded integers, e.g. co-ordinates in a multidimensional array
class Counter : public vector<int>
{
public:
  typedef vector<int>::iterator vector_iterator;
  typedef vector<int>::const_iterator const_vector_iterator;

private:
  vector<int> _begin;
  vector<int> _end;
  vector<int> _multiplier;
  
  void _inc(int i)
    {
      for (; i < (int) size(); i++)
	if (++((*this)[i]) >= _end[i]) (*this)[i] = _begin[i]; else break;
      if (i == (int) size()) *this = _end;
    }

  void _dec(int i)
    {
      for (; i < (int) size(); i++)
	if (--((*this)[i]) < _begin[i]) (*this)[i] = _end[i] - 1; else break;
      if (i == (int) size()) *this = _begin;
    }

 public:
  void reset() { *this = _begin; }

  Counter& swap (Counter& c) { _begin.swap(c._begin); _end.swap(c._end); _multiplier.swap(c._multiplier); vector<int>::swap (c); return *this; }

  Counter& operator++() { _inc(0); return *this; }
  Counter operator++(int) { Counter tmp = *this; ++(*this); return tmp; }

  Counter& operator--() { _dec(0); return *this; }
  Counter operator--(int) { Counter tmp = *this; --(*this); return tmp; }

  Counter& operator+=(ptrdiff_t d)
    {
      for (int i = 0; i < (int) size(); i++)
	{
	  int diff = _end[i] - _begin[i];
	  if (((*this)[i] += d % diff) >= _end[i])
	    {
	      (*this)[i] -= diff;
	      d += diff;
	    }
	  d /= diff;
	}
      return *this;
    }

  Counter& operator-=(ptrdiff_t d)
    {
      for (int i = 0; i < (int) size(); i++)
	{
	  int diff = _end[i] - _begin[i];
	  if (((*this)[i] -= d % diff) < _begin[i])
	    {
	      (*this)[i] += diff;
	      d += diff;
	    }
	  d /= diff;
	}
      return *this;
    }

  Counter operator+(ptrdiff_t d) const { Counter result (*this); result += d; return result; }
  Counter operator-(ptrdiff_t d) const { Counter result (*this); result -= d; return result; }
  
  ptrdiff_t operator-(const vector<int>& c) const
    {
      ptrdiff_t d = 0;
      for (int i = 0; i < (int) size(); i++)
	d += ((*this)[i] - c[i]) * _multiplier[i];
      return d;
    }

  Counter& operator=(const vector<int>& v)
    {
      if (v.size() != size()) THROW Standard_exception ("Counter dimensionality mismatch");
      ((vector<int>&)*this) = v;
      return *this;
    }

  bool operator<(const vector<int>& c2) const { return (*this - c2) < 0; }

  const vector<int>& begin() const { return _begin; }
  const vector<int>& end() const { return _end; }
  const vector<int>& multiplier() const { return _multiplier; }

  vector_iterator vector_begin() { return ((vector<int>&) *this) . begin(); }
  vector_iterator vector_end() { return ((vector<int>&) *this) . end(); }

  const_vector_iterator vector_begin() const { return ((vector<int>&) *this) . begin(); }
  const_vector_iterator vector_end() const { return ((vector<int>&) *this) . end(); }

  Counter() {}

  Counter(const vector<int>& init, const vector<int>& begin, const vector<int>& end) :
    vector<int> (init), _begin(begin), _end(end)
    {
      int mul = 1;
      for (int i = 0; i < (int) size(); i++) { _multiplier.push_back(mul); mul *= _end[i] - _begin[i]; }
    }

  static Counter create_begin (const vector<int>& end) { const vector<int> v (end.size(),0); return Counter (v, v, end); }
  static Counter create_end (const vector<int>& end) { return Counter (end, end, end); }

  static Counter create_begin (int rank, int value) { const vector<int> v (rank, value); return create_begin (v); }
  static Counter create_end (int rank, int value) { const vector<int> v (rank, value); return create_end (v); }

  friend ostream& operator<<(ostream&o, const Counter& c)
    { o << '['; for (int i=0; i<(int) c.size(); i++) { o << c[i]; if (i < (int) c.size()-1) o << ' '; } o << ']'; return o; }

};

#endif
