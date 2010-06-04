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
// algebras.cc - extended real types
//
// Gerton Lunter, 27/8/04
//
//


#include "algebras.h"


BFMantissa *BFloat::aConversionLookup;           // Actual location of the static members of BFloat class
double *BFloat::aDoubleConversionLookup;


_BFloatInitialize _dummyInitializer;             // This initializes aConversionLookup and aDoubleConversionLookup


_BFloatInitialize::_BFloatInitialize() {

  BFloat::aConversionLookup = new BFMantissa[cBFloatConvTableSize];
  BFloat::aDoubleConversionLookup = new double[cBFloatDoubleConvTableSize];

  BFMantissa iBFM = 1.0;
  for (int i = 0; i < cBFloatConvTableSize; i++) {
    BFloat::aConversionLookup[ i ] = iBFM;
    iBFM *= cBFloatRangeInv;
  }

  for (int i = 0; i < cBFloatDoubleConvTableSize; i++) {
    BFloat::aDoubleConversionLookup[ i ] = exp( (i-cBFloatDoubleConvTableSize/2) * logcBFloatRange );
  }

}
