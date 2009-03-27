#ifndef VECTOR_OUTPUT_INCLUDED
#define VECTOR_OUTPUT_INCLUDED

#include <iostream>
#include <vector>

using namespace std;

template <class T>
ostream& operator<< (ostream&o, const vector<T>&v)
{
  for (int i = 0; i < (int) v.size(); i++) { o << v[i]; if (i < (int) v.size()-1) o << " "; }
  return o;
}

template <class S, class T>
ostream& operator<< (ostream&o, const pair<S,T>& p)
{
  o << "(" << p.first << "," << p.second << ")";
  return o;
}

#endif
