// debugging functions

#ifndef CHECK_DERIV_INCLUDED
#define CHECK_DERIV_INCLUDED

#include <math.h>
#include "util/macros.h"

template<class Function1, class Function2>
void check_deriv (Function1 f, Function2 df, double xmin, double xmax, double xstep, ostream& o)
{
  int old_prec = o.precision(3);
  ios::fmtflags old_flags = o.flags();
  ios::fmtflags right_flags = (old_flags & ~ios::left & ~ios::fixed) | ios::right | ios::scientific;
  o.flags(right_flags);
  o.fill(' ');
      
  vector<double> f1v;
  vector<double> dfnumv;
  vector<double> dfanalv;
  double ssq = 0;
  double n = 0;
  double maxerr = 0;
  double maxerrfrac = 0;
  for (double x = xmin; x <= xmax; x += xstep)
    {
      double dx = xstep;
      double f1 = f(x);
      double f2 = f(x+dx);
      double dfnum = (f2 - f1) / dx;  // numeric
      double dfanal = df(x);  // analytic
      f1v.push_back(f1);
      dfnumv.push_back(dfnum);
      dfanalv.push_back(dfanal);
      double err = abs(dfnum-dfanal);
      ssq += err * err;
      n += 1.0;
      if (err > maxerr) maxerr = err;
      if (dfnum != 0) if (err/abs(dfnum) > maxerrfrac) maxerrfrac = err/abs(dfnum);
    }
  if (n > 0.0) ssq = ssq / n;
  // double rms = sqrt(ssq);
  //  o << "rms error = " << rms << "\n";
  o << "max error = " << maxerr << ", max %error = " << ((int) 100*maxerrfrac) << "\n";

  o << "  x: ";
  for (double x = xmin; x <= xmax + xstep/2; x += xstep)
    {
      o.width(10);
      o.fill(' ');
      o << x << " ";
    }
  o << "\n";

  o << "  f: ";
  for_contents (vector<double>, f1v, i)
    {
      o.width(10);
      o.fill(' ');
      o << *i << " ";
    }
  o << "\n";

  o << "num: ";
  for_contents (vector<double>, dfnumv, i)
    {
      o.width(10);
      o.fill(' ');
      o << *i << " ";
    }
  o << "\n";

  o << "anl: ";
  for_contents (vector<double>, dfanalv, i)
    {
      o.width(10);
      o.fill(' ');
      o << *i << " ";
    }
  o << "\n";

  o.flags(old_flags);
  o.precision(old_prec);
}

#endif
