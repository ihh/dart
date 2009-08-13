/////////////////////////////////////////////////////////////////
// SparseMatrix.h
//
// Sparse matrix computations
/////////////////////////////////////////////////////////////////

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>

using namespace std;

const float POSTERIOR_CUTOFF = 0.01;         // minimum posterior probability
                                             // value that is maintained in the
                                             // sparse matrix representation

typedef pair<int,float> PIF;                 // Sparse matrix entry type
                                             //   first --> column
                                             //   second --> value

/////////////////////////////////////////////////////////////////
// SparseMatrix
//
// Class for sparse matrix computations
/////////////////////////////////////////////////////////////////

class SparseMatrix {

  int seq1Length, seq2Length;                     // dimensions of matrix
  VI rowSize;                                     // rowSize[i] = # of cells in row i
  SafeVector<PIF> data;                           // data values
  SafeVector<SafeVector<PIF>::iterator> rowPtrs;  // pointers to the beginning of each row

  VF gapPosteriors;                               // gap posteriors (sum of pairs) for first and second sequences
  int gapPostBase;                                // offset to begining of gapPosteriors

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::SparseMatrix()
  //
  // Private constructor.
  /////////////////////////////////////////////////////////////////

  SparseMatrix (){}

 public:

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::SparseMatrix()
  //
  // Constructor.  Builds a sparse matrix from a posterior matrix.
  // Note that the expected format for the posterior matrix is as
  // a (seq1Length+1) x (seq2Length+1) matrix where the 0th row
  // and 0th column are ignored (they should contain all zeroes).
  /////////////////////////////////////////////////////////////////

  SparseMatrix (int seq1Length, int seq2Length, const VF &posterior) :
  seq1Length (seq1Length), seq2Length (seq2Length), gapPostBase ((seq1Length + 1) * (seq2Length + 1)) {

    int numCells = 0;

    assert (seq1Length > 0);
    assert (seq2Length > 0);

    // calculate memory required; count the number of cells in the
    // posterior matrix above the threshold
    VF::const_iterator postPtr = posterior.begin();
    for (int i = 0; i <= seq1Length; i++){
      for (int j = 0; j <= seq2Length; j++){
        if (*(postPtr++) >= POSTERIOR_CUTOFF){
          assert (i != 0 && j != 0);
          numCells++;
        }
      }
    }

    // allocate memory
    data.resize(numCells);
    rowSize.resize (seq1Length + 1); rowSize[0] = -1;
    rowPtrs.resize (seq1Length + 1); rowPtrs[0] = data.end();
    gapPosteriors.resize(seq1Length + seq2Length + 2);


    // build sparse matrix
    for (int i = 0; i < seq2Length + 1; i++) 
      gapPosteriors[seq1Length + i + 1] = 1;
    postPtr = posterior.begin() + seq2Length + 1;           // note that we're skipping the first row here
    SafeVector<PIF>::iterator dataPtr = data.begin();
    for (int i = 1; i <= seq1Length; i++){
      gapPosteriors[i] = 1;
      postPtr++;                                            // and skipping the first column of each row
      rowPtrs[i] = dataPtr;
      for (int j = 1; j <= seq2Length; j++){
        if (*postPtr >= POSTERIOR_CUTOFF){
          dataPtr->first = j;
          dataPtr->second = *postPtr;
	  gapPosteriors[i] -= *postPtr;
	  gapPosteriors[seq1Length + j + 1] -= *postPtr;
          dataPtr++;
	  if (gapPosteriors[i] < 1e-4)
	    gapPosteriors[i] = 1e-4;
	  if (gapPosteriors[seq1Length + j + 1] < 1e-4)
	    gapPosteriors[seq1Length + j + 1] = 1e-4;
        }
        postPtr++;
      }
      rowSize[i] = dataPtr - rowPtrs[i];
    }
    
    //    for (int i = 0; i < seq1Length + seq2Length + 2; i++)
    //      gapPosteriors[i] = *(postPtr++);
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetRowPtr()
  //
  // Returns the pointer to a particular row in the sparse matrix.
  /////////////////////////////////////////////////////////////////

  SafeVector<PIF>::iterator GetRowPtr (int row) const {
    assert (row >= 1 && row <= seq1Length);
    return rowPtrs[row];
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetValue()
  //
  // Returns value at a particular row, column.
  /////////////////////////////////////////////////////////////////

  float GetValue (int row, int col){
    assert (row >= 1 && row <= seq1Length);
    assert (col >= 1 && col <= seq2Length);
    for (int i = 0; i < rowSize[row]; i++){
      if (rowPtrs[row][i].first == col) 
	return rowPtrs[row][i].second;
    }
    return 0;
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetRowSize()
  //
  // Returns the number of entries in a particular row.
  /////////////////////////////////////////////////////////////////

  int GetRowSize (int row) const {
    assert (row >= 1 && row <= seq1Length);
    return rowSize[row];
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetSeq1Length()
  //
  // Returns the first dimension of the matrix.
  /////////////////////////////////////////////////////////////////

  int GetSeq1Length () const {
    return seq1Length;
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetSeq2Length()
  //
  // Returns the second dimension of the matrix.
  /////////////////////////////////////////////////////////////////

  int GetSeq2Length () const {
    return seq2Length;
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetRowPtr
  //
  // Returns the pointer to a particular row in the sparse matrix.
  /////////////////////////////////////////////////////////////////

  int GetNumCells () const {
    return data.size();
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::Print()
  //
  // Prints out a sparse matrix.
  /////////////////////////////////////////////////////////////////

  void Print (ostream &outfile) const {
    outfile << "Match posteriors:" << endl;
    for (int i = 1; i <= seq1Length; i++){
      outfile << "  " << i << ":";
      for (int j = 0; j < rowSize[i]; j++){
        outfile << " (" << rowPtrs[i][j].first << "," << rowPtrs[i][j].second << ")";
      }
      outfile << endl;
    }
    outfile << "Gap posteriors 0: ";
    for (int i = 1; i <= seq1Length; i++){
      outfile << " (" << i << "," << gapPosteriors[i] << ")";
    }
    outfile << endl << "Gap posteriors 1: ";
    for (int i = 1; i <= seq2Length; i++){
      outfile << " (" << i << "," << gapPosteriors[i + seq1Length + 1] << ")";
    }
    outfile << endl;
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::ComputeTranspose()
  //
  // Returns a new sparse matrix containing the transpose of the
  // current matrix.
  /////////////////////////////////////////////////////////////////

  SparseMatrix *ComputeTranspose () const {

    // create a new sparse matrix
    SparseMatrix *ret = new SparseMatrix();
    int numCells = data.size();

    // gapPostBase must be set properly!
    // (otherwise causes an obscure error with ComputeTranpose()->GetPosterior())
    // RKB -- 2/27/08
    ret->gapPostBase = (*this).gapPostBase;

    ret->seq1Length = seq2Length;
    ret->seq2Length = seq1Length;

    // allocate memory
    ret->data.resize (numCells);
    ret->rowSize.resize (seq2Length + 1); ret->rowSize[0] = -1;
    ret->rowPtrs.resize (seq2Length + 1); ret->rowPtrs[0] = ret->data.end();
    ret->gapPosteriors.resize(seq1Length + seq2Length + 2);

    // compute row sizes
    for (int i = 1; i <= seq2Length; i++) ret->rowSize[i] = 0;
    for (int i = 0; i < numCells; i++)
      ret->rowSize[data[i].first]++;

    // compute row ptrs
    for (int i = 1; i <= seq2Length; i++){
      ret->rowPtrs[i] = (i == 1) ? ret->data.begin() : ret->rowPtrs[i-1] + ret->rowSize[i-1];
    }

    // now fill in data
    SafeVector<SafeVector<PIF>::iterator> currPtrs = ret->rowPtrs;

    for (int i = 1; i <= seq1Length; i++){
      SafeVector<PIF>::iterator row = rowPtrs[i];
      for (int j = 0; j < rowSize[i]; j++){
        currPtrs[row[j].first]->first = i;
        currPtrs[row[j].first]->second = row[j].second;
        currPtrs[row[j].first]++;
      }
    }

    for (int i = 0; i <= seq1Length; i++) {
      ret->gapPosteriors[i + seq2Length + 1] = gapPosteriors[i];
    }
    for (int i = 0; i <= seq2Length; i++) {
      ret->gapPosteriors[i] = gapPosteriors[i + seq1Length + 1];
    }

    return ret;
  }

  /////////////////////////////////////////////////////////////////
  // SparseMatrix::GetPosterior()
  //
  // Return the posterior representation of the sparse matrix.
  /////////////////////////////////////////////////////////////////

  VF *GetPosterior () const {

    // create a new posterior matrix
    VF *posteriorPtr = new VF((seq1Length+1) * (seq2Length+1) + seq1Length + seq2Length + 2,0); assert (posteriorPtr);
    VF &posterior = *posteriorPtr;

    // build the posterior matrix
    for (int i = 0; i < (seq1Length+1) * (seq2Length+1); i++) posterior[i] = 0;
    for (int i = 1; i <= seq1Length; i++){
      VF::iterator postPtr = posterior.begin() + i * (seq2Length+1);
      for (int j = 0; j < rowSize[i]; j++){
        postPtr[rowPtrs[i][j].first] = rowPtrs[i][j].second;
      }
    }
    for (int i = 0; i < seq1Length + seq2Length + 2; i++)
      posterior[gapPostBase + i] = gapPosteriors[i];

    return posteriorPtr;
  }

  float GetGapPosterior(int seqNum, int position) const {
    assert((seqNum == 0 && position <= seq1Length) || (seqNum == 1 && position <= seq2Length));
    return gapPosteriors[position + seqNum * (seq1Length + 1)];
  }

};

#endif
