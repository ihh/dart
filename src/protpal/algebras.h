/*
 *    This file is part of HMMoC 1.3, a hidden Markov model compiler.
 *    Copyright (C) 2007 by Gerton Lunter, Oxford University.
 *
 *    HMMoC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    HMMOC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HMMoC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
\*/
//
// algebras.h - extended real types
//
// Gerton Lunter, 27/8/04
//

#ifndef _algebras_h_
#define _algebras_h_

#include <math.h>
#include <iostream>
#include <cstdlib>
#include <limits>

using namespace std;

// typedefs
typedef float BFMantissa;
const BFMantissa cBFloatRange = 20282409603651670423947251286016.0;  // 2.03e+31; 2^104
const BFMantissa cBFloatRangeInv = 1.0/cBFloatRange;
// Aaron E. Darling 6/7/7: need to typecast to avoid compiler warnings about imprecise FP representations
const BFMantissa cBFloatRangeSqrt    = (BFMantissa)1.0e+18;          // Value between square root of the exponent, and the exponent
const BFMantissa cBFloatRangeInvSqrt = (BFMantissa)1.0e-18;          // Square of this should still be representable, with full mantissa!
const BFMantissa logcBFloatRange     = log(cBFloatRange);
const int cBFloatDigits              = 7;                 // Number of significant digits for printing (7 for floats, 16 for doubles?)
const int cBFloatInfinity            = 1000000000;        // Tiniest number representable is cBFloatRangeInv ^ BFloatInfinity
const int cBFloatConvTableSize       = 100;               // This includes many zero entries, it makes additions a bit faster
const int cBFloatDoubleConvTableSize = 50;                // Table size for bfloat -> double conversion; cBFloatRange^(-size/2) is double 0
//#define BFLOAT_CHECK_UOFLOW                             // Don't bother with under- and overflow checking.


//
// BFloats: more buoyant floats.
//
// struct{ float + int } is 8 bytes; nice size makes noticable speed difference
//
class BFloat {                   
 public:
  static BFMantissa* aConversionLookup;              // used by addition
  static double* aDoubleConversionLookup;            // used by Value()
  BFMantissa f;
  int e;
 public:
  BFloat(BFMantissa iF, int iE) : f(iF), e(iE) {};
  BFloat() {};
  ~BFloat() {};
  inline double Value() const { 
    if (abs(e) < cBFloatDoubleConvTableSize/2) {
      return (double)f * aDoubleConversionLookup[ e + cBFloatDoubleConvTableSize/2 ];
    } else if (e < cBFloatDoubleConvTableSize/2) {
      return 0.0;
    } else {
      return (double)f * exp((double)e * logcBFloatRange);
    }
  }
  void clear() { f=0; e=-cBFloatInfinity; }
};


//
// dummy class to initialise BFloat lookup table
//
class _BFloatInitialize {
public:
  _BFloatInitialize();
};


//
// implementations of BFloat calculations
//


// Normalization of BFloat result of a single operation
#ifdef BFLOAT_CHECK_UOFLOW
static inline void BFloatNormalise(BFloat& a)
     //#define BFloatNormalise(a)
{\
  if (a.f > cBFloatRangeSqrt) {\
    a.f *= cBFloatRangeInv;\
    a.e++;\
  } else if (a.f < cBFloatRangeInvSqrt) {\
    if (a.f == 0.0) {\
      a.e = -cBFloatInfinity;\
    } else {\
      a.f *= cBFloatRange;\
      a.e--;\
    }\
  }\
  if (a.e > cBFloatInfinity) {\
    cerr << "BFloat: Overflow" << endl;\
    a.e = cBFloatInfinity;\
  } else if (a.e < -cBFloatInfinity) {\
    cerr << "BFloat: Underflow" << endl;\
    a.e = -cBFloatInfinity;\
    a.f = 0.0;\
  }\
};
#else
static inline void BFloatNormDown(BFloat& a) { 
  a.f *= cBFloatRangeInv;
  a.e++;
}
static inline void BFloatNormUp(BFloat& a) { 
  if (a.f == 0.0) {
    a.e = -cBFloatInfinity;
  } else {
    a.f *= cBFloatRange;
    a.e--;
  }
}
static inline void BFloatNormalise(BFloat& a)
     //#define BFloatNormalise(a) 
{
  if (a.f > cBFloatRangeSqrt) {
    BFloatNormDown(a);
  } else if (a.f < cBFloatRangeInvSqrt) {
    BFloatNormUp(a);
  }
};
#endif

static inline void DoubleNormalise(double& f, int& e)
{
  // comparing to 0.0 here fails, because the comparison is done
  // using higher-precision doubles, but the subsequent while-loop
  // uses true doubles, resulting in an infinite loop. (G.L. 3/9/07)
  if (f < std::numeric_limits<double>::min()) {
    if (f < 0.0) cerr << "BFloat: Negative number: " << f << endl;
    f = 0.0; 
    e=-cBFloatInfinity;
  } else {
    while (f > (double)cBFloatRangeSqrt) {
      f *= (double)cBFloatRangeInv;
      e++;
    }
    while (f < (double)cBFloatRangeInvSqrt) {
      f *= (double)cBFloatRange;
      e--;
    }
  }
};

// Logarithm of a BFloat
static inline double bfloat_doublelog( const BFloat& a ) { return a.e*logcBFloatRange+log(a.f); }

// BFloat exp of a double
static inline BFloat bfloat_doubleexp( double iA ) 
{
  int iE = (int)floor( iA / log(cBFloatRange) );
  iA -= iE * log(cBFloatRange);
  BFloat iX( exp(iA), iE );
  BFloatNormalise( iX );
  return iX;
}

// Returns a double value - or underflow/overflow if it does not fit.
static inline double bfloat2double( const BFloat bfloat) { return bfloat.Value(); }

// Simplistic double-to-BFloat conversion - can be slow if 'standard' numbers get very large/small
static inline BFloat double2bfloat( double prob) { 
  if (prob < std::numeric_limits<double>::min()) {
    if (prob < 0.0)
      cerr << "BFloat: Negative number: " << prob << endl;
    return BFloat (0.0, -cBFloatInfinity );
  } else {
    register BFloat a( 0.0, 0 );
    while (prob > cBFloatRangeSqrt) {
      prob *= cBFloatRangeInv;
      a.e++;
    }
    while ((prob < cBFloatRangeInvSqrt)) {
      prob *= cBFloatRange;
      a.e--;
    }
    a.f = prob; 
    return a;
  }
}

static inline BFloat bfloat_pr_product (const BFloat& a, const BFloat& b) 
{ 
  register BFloat sf(a.f*b.f,a.e+b.e); 
  BFloatNormalise(sf); 
  return sf; 
}

static inline BFloat bfloat_pr_double_product (const BFloat& a, double b) 
{ 
  register double mantisse = a.f*b;
  int exponent = a.e;
  DoubleNormalise(mantisse, exponent);
  return BFloat(mantisse, exponent);
}

static inline void bfloat_pr_product_accum( BFloat& a, const BFloat& b) { 
  a.f *= b.f; a.e += b.e; 
  BFloatNormalise( a ); 
}

static inline void bfloat_pr_double_product_accum (BFloat& a, double b) 
{ 
  register double mantisse = a.f*b;
  DoubleNormalise(mantisse, a.e);
  a.f = mantisse;
}

static inline BFloat bfloat_pr_quotient( const BFloat& a, const BFloat& b) 
{ 
  register BFloat sf(a.f/b.f, a.e-b.e); 
  BFloatNormalise(sf); 
  return sf;
}
  
static inline void bfloat_pr_quotient_accum( BFloat& a, const BFloat& b) 
{ 
  a.f /= b.f; 
  a.e -= b.e; 
  BFloatNormalise( a ); 
}

static inline BFloat bfloat_pr_sum(const BFloat& a, const BFloat& b) 
{
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return a;
    else
      return BFloat( a.f + b.f * BFloat::aConversionLookup[ a.e - b.e ], a.e );
  } else {
    if (a.e <= b.e - cBFloatConvTableSize)
      return b;
    else
      return BFloat( b.f + a.f * BFloat::aConversionLookup[ b.e - a.e ], b.e );
  }
}
 
static inline void bfloat_pr_sum_accum( BFloat& a, const BFloat& b) 
{
  if (a.e >= b.e) {
    if (a.e < b.e + cBFloatConvTableSize)
      a.f += b.f * BFloat::aConversionLookup[ a.e - b.e ];
  } else {
    if (a.e > b.e - cBFloatConvTableSize) {
      a.f = b.f + a.f * BFloat::aConversionLookup[ b.e - a.e ];
      a.e = b.e;
    } else {
      a = b;
    }
  }
}

  
static inline bool bfloat_equal( const BFloat& a, const BFloat& b) 
{
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return false;
    else
      return a.f == b.f * BFloat::aConversionLookup[ a.e - b.e ];
  }
  if (a.e <= b.e - cBFloatConvTableSize)
    return false;
  else
    return a.f * BFloat::aConversionLookup[ b.e - a.e ] == b.f;
};
static inline bool bfloat_less( const BFloat& a, const BFloat& b) 
{
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return false;
    else
      return a.f < b.f * BFloat::aConversionLookup[ a.e - b.e ];
  }
  if (a.e <= b.e - cBFloatConvTableSize)
    return true;
  else
    return a.f * BFloat::aConversionLookup[ b.e - a.e ] < b.f;
};

static inline bool bfloat_lessequal( const BFloat& a, const BFloat& b) 
{
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return false;
    else
      return a.f <= b.f * BFloat::aConversionLookup[ a.e - b.e ];
  }
  if (a.e <= b.e - cBFloatConvTableSize)
    return true;
  else
    return a.f * BFloat::aConversionLookup[ b.e - a.e ] <= b.f;
};

static inline ostream& bfloat_print( ostream& out, const BFloat& x ) 
{
  static const double log10 = log(10.0);
  static const double maxmantisse = 10.0 * (1.0 - 0.55 * exp(-cBFloatDigits * log10));
  //out.setf(ios::fixed,ios::floatfield);
  out.precision( cBFloatDigits );
  if (x.e == cBFloatInfinity) {
    out << 1.0 << "e+Inf";
  }
  if (x.e == -cBFloatInfinity) {
    out << 1.0 << "e-Inf";
  } else {
    double iM = (log(x.f) + log(cBFloatRange)*(double)x.e) / log10;
    long iExp = long(floor(iM));
    iM = exp((iM - iExp) * log10);
    if (iM > maxmantisse) {
      iExp += 1;
      iM = 1.0;
    }
    out << iM << ( iExp<0 ? "e" : "e+" ) << iExp;
  }
  //out.setf(ios::fixed,ios::floatfield);  // default
  out.precision( 6 );           // default
  return out;
}


//
// Wrapper to allow BFloats to be used by Algebra template
//
struct BFloatMethods
{
  typedef BFloat Value;
  static inline double to_prob (BFloat iX) { return bfloat2double(iX); }
  static inline BFloat from_prob (double iP) { return double2bfloat(iP); }
  static inline BFloat pmul( BFloat iX, BFloat iY) { return bfloat_pr_product(iX,iY); }
  static inline BFloat pmuldouble( BFloat iX, double iY) { return bfloat_pr_double_product(iX,iY); }
  static inline BFloat pdiv( BFloat iX, BFloat iY) { return bfloat_pr_quotient(iX,iY); }
  static inline BFloat psum( BFloat iX, BFloat iY) { return bfloat_pr_sum(iX,iY); }
  static inline BFloat pdiff( BFloat iX, BFloat iY) { cerr << "Bfloat pdiff: Not implemented." << endl; return BFloat(0,0); }
  static inline BFloat doubleexp( double iX) { return bfloat_doubleexp(iX); }
  static inline double doublelog( BFloat iX) { return bfloat_doublelog(iX); }
  static inline void pmulacc( BFloat& iX, BFloat iY) { bfloat_pr_product_accum(iX,iY); }
  static inline void pmulaccdouble( BFloat& iX, double iY) { bfloat_pr_double_product_accum(iX,iY); }
  static inline void pdivacc( BFloat& iX, BFloat iY) { bfloat_pr_quotient_accum(iX,iY); }
  static inline void psumacc( BFloat& iX, BFloat iY) { bfloat_pr_sum_accum(iX,iY); }
  static inline void pdiffacc( BFloat& iX, BFloat iY) { cerr << "Bfloat pdiffacc: Not implemented." << endl; }
  static inline bool less( BFloat iX, BFloat iY) { return bfloat_less(iX,iY); }
  static inline bool equal( BFloat iX, BFloat iY) { return bfloat_equal(iX,iY); }
  static inline bool lessequal( BFloat iX, BFloat iY) { return bfloat_lessequal(iX,iY); }
  static inline ostream& print( ostream& iOut, BFloat iX ) { return bfloat_print( iOut, iX ); }
};


//
// Simple log-space numbers - don't use, except possibly for Viterbi
//
class Logspace {
  double x;
 public:
  Logspace( double x ) : x(x) {}
  Logspace() {}
  operator double&(){ return x; }
  void clear() {x=-1.0e+300;}
};

inline Logspace logspace_addsmall( Logspace iX, Logspace iY ) {
  if (iX - iY > 36.7) return iX;
  return iX + log(1.0+exp(iY-iX));
}

inline Logspace logspace_add( Logspace iX, Logspace iY ) {
  if (iX>iY) return logspace_addsmall(iX,iY); else return logspace_addsmall(iY,iX);
}

struct LogspaceMethods
{
  typedef Logspace Value;
  static inline double to_prob (Value iX) { return exp(iX); }
  static inline Value from_prob (double iP) { return Value(log(iP)); }
  static inline Value pmul( Value iX, Value iY) { return iX+iY; }
  static inline Value pmuldouble( Value iX, double iY) { return iX+log(iY); }
  static inline Value pdiv( Value iX, Value iY) { return iX-iY; }
  static inline Value psum( Value iX, Value iY) { return logspace_add(iX,iY); }
  static inline Value pdiff( Value iX, Value iY) { cerr << "Logspace pdiff: Not implemented." << endl; return 0.0; }
  static inline Value doubleexp( double iX) { return iX; }
  static inline double doublelog( Value iX) { return iX; }
  static inline void pmulacc( Value& iX, Value iY) { iX+=iY; }
  static inline void pmulaccdouble( Value& iX, double iY) { iX+=log(iY); }
  static inline void pdivacc( Value& iX, Value iY) { iX -= iY; }
  static inline void psumacc( Value& iX, Value iY) { iX = logspace_add(iX,iY); }
  static inline void pdiffacc( Value& iX, Value iY) { cerr << "Logspace pdiffacc: Not implemented." << endl; }
  static inline bool less( Value iX, Value iY) { return iX<iY; }
  static inline bool equal( Value iX, Value iY) { return iX==iY; }
  static inline bool lessequal( Value iX, Value iY) { return iX<=iY; }
  static inline ostream& print( ostream& iOut, Value iX ) { return bfloat_print( iOut, bfloat_doubleexp(iX) ); }
};


//
// Algebra class - Wrapper for overloading all arithmetic operators to use a different algebra.
//
// Gerton Lunter, 19/3/03
// Based on logprob.h by by Ian Holmes.
//

template <class AlgebraMethods>
class Algebra {
public:
  // typedef
  typedef typename AlgebraMethods::Value Value;

  // value
  Value val;

public:
  // constructors
  Algebra() { }  // no initialisation, for speed
  Algebra (double px) : val(from_prob(px)) { }
  Algebra (const Algebra& lx) : val(lx.val) { }
  Algebra (const BFloat v) : val(v) { }

  // fast initialization
  void clear() { val.clear(); }

  // assignment operators
  inline Algebra& operator= (const Algebra& lx) { val = lx.val; return *this; }
  inline Algebra& operator= (double px) { val = from_prob(px); return *this; }

  // arithmetic operators; all combinations of Algebra and double are covered
  inline friend Algebra operator+ (const Algebra& lx, const Algebra& ly) { return from_log (psum (lx.val, ly.val)); }
  inline friend Algebra operator+ (const Algebra& lx, double py) { return from_log (psum (lx.val, from_prob(py))); }
  inline friend Algebra operator+ (double px, const Algebra& ly) { return from_log (psum (from_prob(px), ly.val)); }
  inline Algebra& operator+= (const Algebra& lx) { psumacc (val, lx.val); return *this; }
  inline Algebra& operator+= (double px) { psumacc (val, from_prob(px)); return *this; }

  inline friend Algebra operator- (const Algebra& lx, const Algebra& ly) { return from_log (pdiff (lx.val, ly.val)); }
  inline friend Algebra operator- (const Algebra& lx, double py) { return from_log (pdiff (lx.val, from_prob(py))); }
  inline friend Algebra operator- (double px, const Algebra& ly) { return from_log (pdiff (from_prob(px), ly.val)); }
  inline Algebra& operator-= (const Algebra& lx) { pdiffacc (val, lx.val); return *this; }
  inline Algebra& operator-= (double px) { pdiffacc (val, from_prob(px)); return *this; }

  inline friend Algebra operator* (const Algebra& lx, const Algebra& ly) { return from_log (pmul (lx.val, ly.val)); }
  inline friend Algebra operator* (const Algebra& lx, double py) { return from_log (pmuldouble (lx.val, py)); }
  inline friend Algebra operator* (double px, const Algebra& ly) { return from_log (pmuldouble (ly.val, px)); }
  inline Algebra& operator*= (const Algebra& lx) { pmulacc (val, lx.val); return *this; }
  inline Algebra& operator*= (double px) { pmulaccdouble (val, px); return *this; }

  inline friend Algebra operator/ (const Algebra& lx, const Algebra& ly) { return from_log (pdiv (lx.val, ly.val)); }
  inline friend Algebra operator/ (const Algebra& lx, double py) { return from_log (pdiv (lx.val, from_prob(py))); }
  inline friend Algebra operator/ (double px, const Algebra& ly) { return from_log (pdiv (from_prob(px), ly.val)); }
  inline Algebra& operator/= (const Algebra& lx) { pdivacc (val, lx.val); return *this; }
  inline Algebra& operator/= (double px) { pdivacc (val, from_prob(px)); return *this; }

  // miscellaneous operators
  inline friend double log( const Algebra& lx ) { return doublelog( lx.val ); }
  inline friend Algebra exp( const Algebra& px ) { return doubleexp( to_prob(px) ); }
  
  // increment & decremement
  Algebra& operator++() { *this += 1.; return *this; }
  Algebra operator++(int) { Algebra tmp (*this); ++(*this); return tmp; }

  Algebra& operator--() { *this -= 1.; return *this; }
  Algebra operator--(int) { Algebra tmp (*this); --(*this); return tmp; }

  // relational operators
  inline friend int operator== (const Algebra& lx, const Algebra& ly) { return equal(lx.val, ly.val); }
  inline friend int operator== (const Algebra& lx, const double py) { return equal(lx.val, from_prob(py)); }
  inline friend int operator== (const double px, const Algebra& ly) { return equal(from_prob(px), ly.val); }

  inline friend int operator!= (const Algebra& lx, const Algebra& ly) { return !equal(lx.val, ly.val); }
  inline friend int operator!= (const Algebra& lx, const double py) { return !equal(lx.val, from_prob(py)); }
  inline friend int operator!= (const double px, const Algebra& ly) { return !equal(from_prob(px), ly.val); }

  inline friend int operator< (const Algebra& lx, const Algebra& ly) { return less(lx.val, ly.val); }
  inline friend int operator< (const Algebra& lx, const double py) { return less(lx.val, from_prob(py)); }
  inline friend int operator< (const double px, const Algebra& ly) { return less(from_prob(px), ly.val); }

  inline friend int operator> (const Algebra& lx, const Algebra& ly) { return less(ly.val, lx.val); }
  inline friend int operator> (const Algebra& lx, const double py) { return less(from_prob(py), lx.val); }
  inline friend int operator> (const double px, const Algebra& ly) { return less(ly.val, from_prob(px)); }

  inline friend int operator<= (const Algebra& lx, const Algebra& ly) { return lessequal(lx.val, ly.val); }
  inline friend int operator<= (const Algebra& lx, const double py) { return lessequal( lx.val, from_prob(py) ); }
  inline friend int operator<= (const double px, const Algebra& ly) { return lessequal( from_prob(px), ly.val); }

  inline friend int operator>= (const Algebra& lx, const Algebra& ly) { return lessequal( ly.val, lx.val); }
  inline friend int operator>= (const Algebra& lx, const double py) { return lessequal( from_prob(py), lx.val ); }
  inline friend int operator>= (const double px, const Algebra& ly) { return lessequal( ly.val, from_prob(px) ); }

  // stream operators
  inline friend ostream& operator<< (ostream& out, const Algebra& lx) { return AlgebraMethods::print(out, lx.val); }
  inline friend istream& operator>> (istream& in, const Algebra& lx) { double px; in >> px; lx.val = px; return in; }

  // cast operators
  inline double prob() const { return to_prob (val); }
  inline operator double() const { return to_prob (val); }

private:
  // private AlgebraMethods method wrappers
  static inline double to_prob (Value X) { return AlgebraMethods::to_prob (X); }
  static inline Value from_prob (double P) { return AlgebraMethods::from_prob (P); }
  static inline Value pmul (Value X, Value Y) { return AlgebraMethods::pmul (X, Y); }
  static inline Value pmuldouble (Value X, double Y) { return AlgebraMethods::pmuldouble (X, Y); }
  static inline Value pdiv (Value X, Value Y) { return AlgebraMethods::pdiv( X, Y); }
  static inline Value psum (Value X, Value Y) { return AlgebraMethods::psum (X, Y); }
  static inline Value pdiff (Value X, Value Y) { return AlgebraMethods::pdiff (X, Y); }
  static inline Value doubleexp (double X) { return AlgebraMethods::doubleexp( X ); }
  static inline double doublelog (Value X) { return AlgebraMethods::doublelog( X ); }
  static inline void pmulacc (Value& X, Value Y) { AlgebraMethods::pmulacc (X, Y); }
  static inline void pmulaccdouble (Value& X, double Y) { AlgebraMethods::pmulaccdouble (X, Y); }
  static inline void pdivacc( Value& X, Value Y) { AlgebraMethods::pdivacc( X, Y); }
  static inline void psumacc (Value& X, Value Y) { AlgebraMethods::psumacc (X, Y); }
  static inline void pdiffacc (Value& X, Value Y) { AlgebraMethods::pdiffacc (X, Y); }
  static inline bool less (Value X, Value Y ) { return AlgebraMethods::less( X, Y ); }
  static inline bool equal (Value X, Value Y ) { return AlgebraMethods::equal( X, Y ); }
  static inline bool lessequal( Value X, Value Y ) { return AlgebraMethods::lessequal( X, Y ); }

public:
  // static constructor from logspace value
  static inline Algebra from_log (Value X) { Algebra lx; lx.val = X; return lx; }
};


//
// and these are the the things that we'll use:
//

#define bfloat Algebra<BFloatMethods>

#define logspace Algebra<LogspaceMethods>

#endif
