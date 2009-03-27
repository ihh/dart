#ifndef MAXIMISE_INCLUDED
#define MAXIMISE_INCLUDED

// Generic one-dimensional optimisation routines
// These are exclusively used for estimating branch lengths in trees.
// Don't use them for anything else! Lots of magic numbers and tweakery in here.
// You don't want to break the tree-fitting code, now do ya... IH 6/1/06

#include <math.h>

#include <map>
#include <list>
#include <algorithm>

#include "util/macros.h"
#include "util/logfile.h"
#include "util/math_fn.h"

// Templated class to cache the results of function calls on a one-dimensional real domain
//
template<class AdaptableUnaryFunction>
class Function_cache
{
private:
  typedef typename AdaptableUnaryFunction::result_type result_type;
  typedef          double                              argument_type;
  typedef          map <double, result_type*>          Cache;
  typedef typename Cache::value_type                   Cache_point;
  typedef typename Cache::iterator                     Cache_iterator;

  AdaptableUnaryFunction&  f;
  const double             resolution;
  Cache                    cache;

public:
  Function_cache (AdaptableUnaryFunction& f, double resolution = TINY) : f(f), resolution(resolution) { }
  ~Function_cache() { template_for_contents (Cache, cache, c) delete (*c).second; }
  result_type operator() (argument_type x)
    {
      Cache_iterator i = lower_bound (cache.begin(), cache.end(), Cache_point (x*(1+resolution), (result_type*) 0));
      Cache_iterator j = i;
      if (i != cache.begin())
	if ((*--j).first >= x - resolution) return *(*j).second;
      result_type* r = new result_type (f ((typename AdaptableUnaryFunction::argument_type) x));      // cast is just in case f() doesn't take doubles. We always do, though.
      cache.insert (i, Cache_point(x,r));
      return *r;
    }
};

// Templated routine to find a preliminary bracketing interval for the maximum of a function on a one-dimensional real domain
//
// On entry, bracket is [xmin,xmax]
// On exit, bracket is [x1,x2,x3] with  xmin <= x1 < x2 < x3 <= xmax  and  max(f(x1),f(x3)) < f(x2)
// Return code is FALSE if above condition is not satisfied.
//
// Use a Function_cache (above) to avoid duplicating time-consuming function evaluations
// 
//
// IH, 12/31/2005:
// This code is mainly used to maximise time-likelihood functions.
// It's called from lots of places, generally with (xmin=0, xmax=10) for no particular reason.
// Unfortunately, the maxima of these functions tends to be at x<1, which the bracketing routine
// doesn't handle very well.
// To hack around this, I've hard-coded an initial bracketing interval that's very close to the
// start of the range, with x1 just positive enough so that we don't fall afoul of the singularity
// at x=0.
// This is highly unsatisfactory, and is probably due for a major refactoring.
template<class Function>
bool bracket_maximum (Function& f, double& x1, double& x2, double& x3, double xmin, double xmax)
{
  // some magic numbers
  // IH, 4/26/2006: The following numbers are really crude hacks.
  // Bracketing appears to be very sensitive to init_x1 in particular
  // (set it too low and all time estimates come out to be this value).
  // We ***really*** need a better maximise routine.
  const double default_mag = 1.618034;
  const double max_mag = 50.0;
  const double init_x1 = .00001;  // expressed as a fraction of (xmax-xmin); equals 0.0001 if (xmin,xmax)=(0,10)
  const double init_x2 = .050;  // expressed as a fraction of (xmax-xmin); equals 0.50 if (xmin,xmax)=(0,10)

  CTAG(2,MAXIMISE MAXIMISE_BRACKET) << "Scanning for bracketing interval in range [" << xmin << " , " << xmax << "]\n";
  
  if (CTAGGING(-5,MAXIMISE_DEBUG))
    {
      const double n_points = 100.;
      CL << "Debug: " << n_points << " points from " << x1 << " to " << x3 << "\n";
      for (double x = xmin; x < xmax; x += (xmax-xmin)/n_points)
	{
	  const double fx = f(x);
	  CTAG(-5,MAXIMISE_DEBUG) << "f(" << x << ") = " << fx << "\n";
	}
    }

  // start with x1 and x2 near the beginning of the range (since this code is mostly used to maximise time-likelihood functions)
  //
  double f1 = f(x1 = xmin + init_x1 * (xmax - xmin));
  double f2 = f(x2 = xmin + init_x2 * (xmax - xmin));
  if (f2 < f1) { swap(x1,x2); swap(f1,f2); }       // ensure that x1->x2 is uphill

  x3 = x1 < x2 ? (xmax + x2)/2 : (xmin + x1)/2;   // first guess for x3
  double f3 = f(x3);
  int steps = 0;
  while (f2 < f3 && x3-TINY >= xmin && x3+TINY <= xmax)
    {
      ++steps;
      CTAG(0,MAXIMISE_BRACKET) << "Bracketing, step #" << steps << ": f(" << x1 << ")=" << f1 << ", f(" << x2 << ")=" << f2 << ", f(" << x3 << ")=" << f3 << "\n";

      // try a parabolic estimate for x3 first; if that fails, try stepping uphill
      //
      double a = (x3-x1) * (f2-f1);
      double b = (x2-x1) * (f3-f1);
      double d = 2 * sgn(a-b) * max(abs(a-b),TINY);
      double xnew = ((x3-x1) * a - (x2-x1) * b) / d + x1;
      double xlim = max (xmin, min (xmax, x2 + max_mag*(x3-x2)));
      double fnew;
      
      if ((xnew-x2)*(x3-xnew) > 0.0)          // parabolic estimate is between x2 and x3
	{
	  if ((fnew = f(xnew)) > f3) { x1 = x2; f1 = f2; x2 = xnew; f2 = fnew; break; }   // maximum is between xnew and x3
	  else if (fnew < f2) { x3 = xnew; f3 = fnew; break; }       // maximum is between x1 and x3
	  fnew = f(xnew = max (xmin, min (xmax, x3 + default_mag * (x3-x2))));      // parabolic fit didn't work, so increase the step size instead
	}
      else if ((xnew-x3)*(xlim-xnew) > 0.0)   // parabolic estimate is between x3 and xlim
	{
	  if ((fnew = f(xnew)) > f3) { x2 = x3; f2 = f3; x3 = xnew; f3 = fnew; fnew = f(xnew = max (xmin, min (xmax, x3 + default_mag * (x3-x2)))); }    // if parabolic estimate is uphill, keep climbing
	}
      else if ((xnew-xlim)*(xlim-x3) > 0.0)   // parabolic estimate is beyond xlim
	fnew = f(xnew = xlim);
      else
	fnew = f(xnew = max (xmin, min (xmax, x3 + default_mag * (x3-x2))));    // failsafe: step uphill
      
      // eliminate oldest point
      //
      x1 = x2; f1 = f2;
      x2 = x3; f2 = f3;
      x3 = xnew; f3 = fnew;
    }

  CTAG(2,MAXIMISE MAXIMISE_BRACKET) << "Bracketing interval: f(" << x1 << ")=" << f1 << ", f(" << x2 << ")=" << f2 << ", f(" << x3 << ")=" << f3 << "   (found in " << steps << " steps)\n";
  
  return max(f1,f3) < f2;
}


// Templated routine to find the bracketed maximum of a differentiable function on a one-dimensional real domain
//  - based on Brent's method in one dimension (slips into a conservative bisection method if secant extrapolation method starts jumping around)
//
// On entry, bracket is   [x1,x3]   with x2 lying between x1 and x3 and   max(f(x1),f(x3)) < f(x2)
// On exit,   xmax = argmax[f] +/- res*xmax   and   fmax = f(xmax)
// Returns FALSE if maximisation failed
//
template<class Function1, class Function2>
bool brent_deriv (Function1& f, Function2& df, double x1, double x2, double x3, double res, double& xmax, double& fmax)
{
  const int max_iter = 100;

  double f1 = f(x1), f2 = f(x2), f3 = f(x3);

  double a = min(x1,x3);    // [a,b] = bracketing interval
  double b = max(x1,x3);

  double& x = xmax;         // (x,f,df) = best point so far
  double& fx = fmax;
  x = x2;
  fx = f2;
  double dx = df(x);

  double w=x, v=x;          // w = second best point, v = third best
  double fw=fx, fv=fx;
  double dw=dx, dv=dx;

  double last_step = 0.0;

  CTAG(2,MAXIMISE) << "Looking for a maximum between f(" << a << ")=" << (x1<x3?f1:f3) << " and f(" << b << ")=" << (x1>x3?f1:f3) << ",  fractional resolution " << res << "\n";
  bool log_flag = CTAGGING(-1,MAXIMISE_PROGRESS);
  sstring report, move_type;

  int iter;
  for (iter = 1; iter < max_iter; iter++)
    {
      if (log_flag) { report = ""; report << "Iteration #" << iter << ": best f(" << x << ")=" << fx << ", df/dx=" << dx << ", 2nd best x=" << w << ", 3rd best x=" << v; }
      double step = 0;
      double xmedian = 0.5 * (a+b);
      double xres = res*abs(x) + TINY;
      double xres2 = 2*xres;
      if (abs(x-xmedian) + 0.5*(b-a) <= xres2) { if (log_flag) CTAG(-1,MAXIMISE_PROGRESS) << report << "\n"; break; }    // finish here if distance between x & furthest bracket is small enough

      if (abs(last_step) > xres)         // try secant method unless it's giving steps that are too small
	{
	  double wstep = 2*(b-a), vstep = wstep;       // initialise these to be out-of-range
	  if (dw != dx) wstep = (w-x) * dx/(dx-dw);         // extrapolate derivative to zero, using secant method from x to w and v
	  if (dw != dx) vstep = (v-x) * dx/(dx-dv);
	  bool accept_wstep =  x+wstep > a  &&  x+wstep < b  &&  dx*wstep > 0.0;   // make sure we move to an in-bracket point in the same direction as the derivative
	  bool accept_vstep =  x+vstep > a  &&  x+vstep < b  &&  dx*vstep > 0.0;
	  double last_tmp = last_step;
	  last_step = step;
	  if (accept_wstep || accept_vstep)
	    {
	      if (accept_wstep && accept_vstep)                  // take smallest of two eligible steps
		step = abs(wstep) < abs(vstep) ? wstep : vstep;
	      else
		step = accept_wstep ? wstep : vstep;
	      if (abs(step) < 0.5*abs(last_tmp))                 // Brent condition: accept only if steps are decreasing in size by factor sqrt(2) per step, averaged over two steps
		{
		  double xnew = x + step;
		  if (xnew-a < xres2 || b-xnew < xres2)          // if step takes us too close to a bracket, then inch away from the nearest bracket instead
		    {
		      step = sgn(xnew-x) * xres;
		      if (log_flag) move_type = "adjusting";
		    }
		  else if (log_flag) move_type = "interpolating";
		}
	      else                 // conservative method: bisect interval pointed to by derivative
		{
		  step = 0.5 * (last_step = (dx < 0.0 ? a-x : b-x));
		  if (log_flag) move_type = "bisecting";
		}
	    }
	  else
	    {
	      step = 0.5 * (last_step = (dx < 0.0 ? a-x : b-x));
	      if (log_flag) move_type = "bisecting";
	    }
	}
      else
	{
	  step = 0.5 * (last_step = (dx < 0.0 ? a-x : b-x));
	  if (log_flag) move_type = "bisecting";
	}
      
      double xnew, fnew;
      if (abs(step) >= xres)
	{
	  fnew = f(xnew = x+step);
	  if (log_flag) CTAG(-1,MAXIMISE_PROGRESS) << report << " (" << move_type << ")  New f(" << xnew << ")=" << fnew << "\n";
	}
      else
	{
	  fnew = f(xnew = x + sgn(step)*xres);
	  if (log_flag) CTAG(-1,MAXIMISE_PROGRESS) << report << " (inching)  New f(" << xnew << ")=" << fnew << "\n";
	  if (fnew < fx) break;              // finish here if the smallest allowable step in the downhill direction takes us uphill
	}
      double dnew = df(xnew);
      if (fnew > fx)  // changed from ">=", IH 6/1/06
	{
	  if (xnew > x) a=x; else b=x;
	  v = w; fv = fw; dv = dw;
	  w = x; fw = fx; dw = dx;
	  x = xnew; fx = fnew; dx = dnew;
	}
      else
	{
	  if (xnew > x) b=xnew; else a=xnew;
	  if (fnew >= fw || w==x)
	    {
	      v = w; fv = fw; dv = dw;
	      w = xnew; fw = fnew; dw = dnew;
	    }
	  else if (fnew > fv || v==x || v==w)
	    {
	      v = xnew; fv = fnew; dv = dnew;
	    }
	}
    }
  
  if (iter == max_iter)
    CTAG(2,MAXIMISE) << "Giving up after " << iter << " iterations; best point found was f(" << xmax << ")=" << fmax << "\n";
  else
    CTAG(2,MAXIMISE) << "Maximum value f(" << xmax << ")=" << fmax << " found in " << iter << " iterations\n";

  // hacky patch (IH 6/1/06): if f(x1)>=f(x), then set x=x1, to keep branch lengths minimal
  if (f1 >= fmax)
    {
      CTAG(2,MAXIMISE) << "Reverting to equally good but smaller-valued point at f(" << x1 << ")=" << f1 << "\n";
      xmax = x1;
    }
  
  return iter < max_iter;
}


// struct encapsulating Function_cache, bracket_maximum & brent_deriv
template <class F, class D>
struct Cached_function
{
  // typedefs
  typedef F func_type;
  typedef D deriv_type;

  // data
  Function_cache<func_type> func_cache;
  Function_cache<deriv_type> deriv_cache;
  double min_arg, max_arg, arg_res;

  // constructor
  Cached_function (func_type& func, deriv_type& deriv, double min_arg = 0., double max_arg = 10., double arg_res = .01)
    : func_cache (func), deriv_cache (deriv),
      min_arg (min_arg), max_arg (max_arg), arg_res (arg_res)
  { }

  // method to find mode
  double find_max()
  {
    double x1, x2, x3, xbest, fbest;
    bracket_maximum (func_cache, x1, x2, x3, min_arg, max_arg);
    brent_deriv (func_cache, deriv_cache, x1, x2, x3, arg_res, xbest, fbest);
    return xbest;
  }
};

#endif
