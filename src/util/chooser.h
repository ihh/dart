// ordered multi-digit counter class
// ordered counter is big-endian, ie MSD-first,
// with consecutive elements (x_{n}, x_{n+1}) satisfying x_{n} < x_{n+1}.
//
// useful for iterating through all possible K-selections from a set of N integers.

#ifndef CHOOSER_INCLUDED
#define CHOOSER_INCLUDED

class Chooser : public vector<int>
{
 private:
  int _N;

 public:
  Chooser& operator++()
    {
      int i;
      for (i = size() - 1; i >= 0; i--)
	if (++((*this)[i]) < _N) break;
      while (++i < size()) (*this)[i] = (*this)[i-1] + 1;
      return *this;
    }

  Chooser operator++(int) { Chooser tmp = *this; (*this)++; return tmp; }

  bool done() const { return back() >= _N; }

  Chooser (int N, int K) : vector<int>(K), _N(N)
    {
      for (int i = 0; i < K; i++) (*this)[i] = i;
    }

};

#endif

