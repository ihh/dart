/*
 *    This file is part of HMMoC 1.2.1, a hidden Markov model compiler.
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
#include <iostream>

#include "dptables.h"

using std::cout;
using std::min;
using std::max;

// Implements a traversal of the dynamic programming table along a diagonal band.
// The band should remain (diagonally) connected even if iWidth==1, and the forward
// and backward iterations must visit the same sites - these two requirements made the
// implementation slightly tricky.
class TwoDBanding : Banding<2> {

 private:

  int iLen0, iLen1, iDiv1;
  int iWidth;
  Position pos;
  int iLastRow;

  int diagonal(int p1) { return (p1 * iLen0 + iLen1/2) / iDiv1; }

 public:

  TwoDBanding( int iLen0, int iLen1, int iWidth ) : iLen0(iLen0), iLen1(iLen1), iDiv1(iLen1 ? iLen1 : 1), iWidth(iWidth) {}

  Position& forwardIterator() {

    pos[0] = pos[1] = 0;
    iLastRow = iLen0;
    return pos;

  }

  Position& backwardIterator() {

    pos[0] = iLen0;
    pos[1] = iLen1;
    iLastRow = 0;
    return pos;

  }

  bool hasNextForward() {

    if (pos[0] < iLen0 && pos[0] < (iWidth-1)/2 + max( diagonal(pos[1]+1)-1, diagonal(pos[1]))) {
      ++pos[0];
      return true;
    }
    if (pos[1]<iLen1) {
      ++pos[1];
      pos[0] = max(0, diagonal(pos[1]) - iWidth/2);
      return true;
    }
    return false;
  }

  bool hasNextBackward() {

    if (pos[0] > 0 && pos[0] > diagonal(pos[1]) - iWidth/2) {
      --pos[0];
      return true;
    }
    if (pos[1]>0) {
      --pos[1];
      pos[0] = min(iLen0, (iWidth-1)/2 + max( diagonal(pos[1]+1)-1, diagonal(pos[1])));
      return true;
    }
    return false;
  }

  bool lastColumnEntry() {

    return pos[0] == iLastRow;

  }

  void warning() {

    cout << "Warning - out of bounds at position (" << pos[0] << "," << pos[1] << ")" << endl;

  }

};


class ThreeDBanding : Banding<3> {

 private:

  int iLen0, iLen1, iLen2, iDiv1, iDiv2;
  int iDiameter;
  Position pos;
  int iLastRow;

  int diagonal0_1(int p1) { return (p1 * iLen0 + iLen1/2) / iDiv1; } //return diagonal position in dim 0 given position in dim 1
  int diagonal1_2(int p2) { return (p2 * iLen1 + iLen2/2) / iDiv2; } //return diagonal position in dim 1 given position in dim 2

 public:

  ThreeDBanding( int iLen0, int iLen1, int iLen2, int iDiameter ) : iLen0(iLen0), iLen1(iLen1), iLen2(iLen2), iDiv1(iLen1 ? iLen1 : 1), iDiv2(iLen2 ? iLen2 : 1), iDiameter(iDiameter) {}

  Position& forwardIterator() {

    pos[0] = pos[1]  = pos[2] = 0;
    iLastRow = iLen0;
    return pos;

  }

  Position& backwardIterator() {

    pos[0] = iLen0;
    pos[1] = iLen1;
    pos[2] = iLen2;
    iLastRow = 0;
   return pos;

  }

  bool hasNextForward() {

    //iterate along 0th dimension... max( diagonalI_J(pos[J]+1)-1, diagonalI_J(pos[J]) assure diagonal connectivity
    if (pos[0] < iLen0 && pos[0] < (iDiameter-1)/2 + max( diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[1]))) { 
      ++pos[0];
      return true;
    }
    //if at the end of 0th dimension, move iterator one in 1st dimension 
    if (pos[1]<iLen1 && pos[1] < (iDiameter-1)/2  + max (diagonal1_2(pos[2]+1)-1, diagonal1_2(pos[2]))) {
      ++pos[1];
      pos[0] = max(0, diagonal0_1(pos[1]) - iDiameter/2); //reset position in dim 0 to beginning of band
      return true;
    }
    //if at the end of the 1st dimension, move iterator in 2nd dimension
    if (pos[2]<iLen2) {
	++pos[2];
	pos[1] = max(0, diagonal1_2(pos[2]) - iDiameter/2); //reset position in dim 1
	pos[0] = max(0, diagonal0_1(pos[1]) - iDiameter/2); //reset position in dim 0
	return true;
    }
    return false;
  }

  bool hasNextBackward() {

    //iterate backwards along 0th dim
    if (pos[0] > 0 && pos[0] > diagonal0_1(pos[1]) - iDiameter/2) {
      --pos[0];

      return true;
    }
    //if at the end of 0th dim, move one back in dim 1
    if (pos[1]>0 && pos[1] > diagonal1_2(pos[2] - iDiameter/2)) {
      --pos[1];
      pos[0] = min(iLen0, (iDiameter-1)/2 + max( diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[1]))); //reset position in dim 0
	 return true;
    }
    //if at the end of dim 1, move one back in dim 2
    if (pos[2]>0){
	--pos[2];
	pos[1] = min(iLen1, (iDiameter-1)/2 + max( diagonal1_2(pos[2]+1)-1, diagonal1_2(pos[2]))); //reset position in dim 1
	pos[0] = min(iLen0, (iDiameter-1)/2 + max( diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[2]))); //reset position in dim 0
	return true;
    }
    return false;
  }

  bool lastColumnEntry() {

    return pos[0] == iLastRow;

  }

  void warning() {

    cout << "Warning - out of bounds at position (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << endl;

  }

};

class FourDBanding : Banding<3> {

 private:

  int iLen0, iLen1, iLen2, iLen3, iDiv1, iDiv2, iDiv3;
  int iDiameter;
  Position pos;
  int iLastRow;

  int diagonal0_1(int p1) { return (p1 * iLen0 + iLen1/2) / iDiv1; } //return diagonal position in dim 0 given position in dim 1
  int diagonal1_2(int p2) { return (p2 * iLen1 + iLen2/2) / iDiv2; } //return diagonal position in dim 1 given position in dim 2
  int diagonal2_3(int p3) { return (p3 * iLen2 + iLen3/2) / iDiv3; } //return diagonal position in dim 2 given position in dim 3

 public:

  FourDBanding( int iLen0, int iLen1, int iLen2, int iLen3, int iDiameter ) : iLen0(iLen0), iLen1(iLen1), iLen2(iLen2), iLen3(iLen3), iDiv1(iLen1 ? iLen1 : 1), iDiv2(iLen2 ? iLen2 : 1), iDiv3(iLen3 ? iLen3 : 1), iDiameter(iDiameter) {}

  Position& forwardIterator() {

    pos[0] = pos[1]  = pos[2] = pos[3] = 0;
    iLastRow = iLen0;
    return pos;

  }

  Position& backwardIterator() {

    pos[0] = iLen0;
    pos[1] = iLen1;
    pos[2] = iLen2;
    pos[3] = iLen3;
    iLastRow = 0;
   return pos;

  }

  bool hasNextForward() {

    //iterate along 0th dimension
    if (pos[0] < iLen0 && pos[0] < (iDiameter-1)/2 + max(diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[1]))) { 
      ++pos[0];
      return true;
    }
    //if at the end of 0th dimension, move iterator one in 1st dimension 
    if (pos[1]<iLen1 && pos[1] < (iDiameter-1)/2  + max (diagonal1_2(pos[2]+1)-1, diagonal1_2(pos[2]))) {
      ++pos[1];
      pos[0] = max(0, diagonal0_1(pos[1]) - iDiameter/2); //reset position in dim 0 to beginning of band
      return true;
    }
    //if at the end of the 1st dimension, move iterator in 2nd dimension
    if (pos[2]<iLen2 && pos[2] < (iDiameter-1)/2 + max (diagonal2_3(pos[3]+1)-1, diagonal2_3(pos[3]))) {
	++pos[2];
	pos[1] = max(0, diagonal1_2(pos[2]) - iDiameter/2); //reset position in dim 1
	pos[0] = max(0, diagonal0_1(pos[1]) - iDiameter/2); //reset position in dim 0
	return true;
    }
    //if at the end of 2nd dim, move in 3rd dim
    if (pos[3]<iLen3) {
	++pos[3];
	pos[2] = max(0, diagonal2_3(pos[3]) - iDiameter/2);
	pos[1] = max(0, diagonal1_2(pos[2]) - iDiameter/2);
	pos[0] = max(0, diagonal0_1(pos[1]) - iDiameter/2); 
	return true;
    }
    return false;
  }

  bool hasNextBackward() {

    //iterate backwards along 0th dim
    if (pos[0] > 0 && pos[0] > diagonal0_1(pos[1]) - iDiameter/2) {
      --pos[0];
      return true;
    }
    //if at the end of 0th dim, move one back in dim 1
    if (pos[1]>0 && pos[1] > diagonal1_2(pos[2] - iDiameter/2)) {
      --pos[1];
      pos[0] = min(iLen0, (iDiameter-1)/2 + max( diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[1]))); //reset position in dim 0
      return true;
    }
    //if at the end of dim 1, move one back in dim 2
    if (pos[2]>0 && pos[2] > diagonal2_3(pos[3] - iDiameter/2)){
	--pos[2];
	pos[1] = min(iLen1, (iDiameter-1)/2 + max( diagonal1_2(pos[2]+1)-1, diagonal1_2(pos[2]))); //reset position in dim 1
	pos[0] = min(iLen0, (iDiameter-1)/2 + max( diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[2]))); //reset position in dim 0
	return true;
    }
    if(pos[3]>0){
	--pos[3];
	pos[2] = min(iLen2, (iDiameter-1)/2 + max( diagonal2_3(pos[3]+1)-1, diagonal2_3(pos[3])));
	pos[1] = min(iLen1, (iDiameter-1)/2 + max( diagonal1_2(pos[2]+1)-1, diagonal1_2(pos[2]))); //reset position in dim 1
	pos[0] = min(iLen0, (iDiameter-1)/2 + max( diagonal0_1(pos[1]+1)-1, diagonal0_1(pos[2]))); //reset position in dim 0
	return true;
    }
    return false;
  }

  bool lastColumnEntry() {

    return pos[0] == iLastRow;

  }

  void warning() {

    cout << "Warning - out of bounds at position (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << endl;

  }

};
