/*****************************************************************/
// MultiSequenceDag.h
//
// Classes for representing a multiple sequence alignment as a DAG,
// and aligning using sequence-annealing.
/*****************************************************************/


#ifndef MULTISEQUENCEDAG_H
#define MULTISEQUENCEDAG_H

// no annotation character
#define NOANNOT "-"

#include <list>
#include <map>
#include <queue>
#include <iostream>
#include <limits>

#include "MultiSequence.h"
#include "SparseMatrix.h"

using namespace std;
typedef map<int,int> MII;
typedef pair<int,int> PII;
const float INVALID_EDGE = -1e10;

/*****************************************************************/
// Column
//
// A class for storing alignment column information.
/*****************************************************************/

class Column {

  int index;                                // the position of the column in the alignment
  bool visited;                             // indicates if the column has been visited by the dfs 
  MII seqPositions;                         // map of (seq_id,pos) pairs in the column
  Column *mergedInto;                       // the Column this column has merged into
  char consensus;                           // the consensus character for this column

  float matchAccuracy; // expected accuracy of matches for this column
                       // match accuracy defined as sum over (n choose 2) pairs of aligned characters of
                       //                       P(char1 ~ char2) - P(char1 ~ -) - P(char2 ~ -)
                       // It can therefore take a maximum value of (n choose 2) * 2 (factor of 2 b/c counting # chars, not just # pairs).
                       //                   (cf calcTgfWeight())
                       // Note that this value can become negative if the gap factor is < 1!
                       //                   -- RKB 5/1/08

 public:

  /*****************************************************************/
  // Column::Column()
  //
  // Constructor.  Creates a new empty column.
  /*****************************************************************/
    
  //  Column (int pos) : index (pos), visited (false), mergedInto (this), consensus (0) { 
  Column (int pos) : index (pos), visited (false), mergedInto (this), consensus (0), matchAccuracy (0) { 
  }

  /*****************************************************************/
  // Column::operator<()
  //
  // Compares two columns based on their current position.
  /*****************************************************************/

  bool operator< (Column const &c2) {
    return index < c2.index;
  }

  /*****************************************************************/
  // Column::operator==()
  //
  // Compares two columns based on their current position.
  /*****************************************************************/

  bool operator== (Column const &c2) {
    return index == c2.index;
  }

  /*****************************************************************/
  // Column::Marked()
  //
  // Checks if the column has been marked as visited.
  /*****************************************************************/

  bool Marked () {
    return visited;
  }

  /*****************************************************************/
  // Column::Mark()
  //
  // Marks the column as visited.
  /*****************************************************************/

  void Mark () {
    visited = true;
  }

  /*****************************************************************/
  // Column::Unmark()
  //
  // Unmarks the column as visited.
  /*****************************************************************/

  void Unmark () {
    visited = false;
  }

  /*****************************************************************/
  // Column::GetIndex()
  //
  // Returns the current position of the column in the alignment.
  /*****************************************************************/

  int GetIndex () {
    return index;
  }

  /*****************************************************************/
  // Column::SetIndex()
  //
  // Sets the current position of the column in the alignment.
  /*****************************************************************/

  void SetIndex (int idx) {
    index = idx;
  }

  /*****************************************************************/
  // Column::GetMergedInto()
  //
  // Gets a pointer to the column this column has been merged into.
  /*****************************************************************/

  Column* GetMergedInto () {
    return mergedInto;
  }

  /*****************************************************************/
  // Column::SetMergedInto()
  //
  // Sets a pointer to the column this column has been merged into.
  /*****************************************************************/

  void SetMergedInto (Column *newcol) {
    mergedInto = newcol;
  }

  /*****************************************************************/
  // Column::GetSeqPositions()
  //
  // Returns the sequence positions in this column.
  /*****************************************************************/

  const MII & GetSeqPositions () const {
    return seqPositions;
  }

  /*****************************************************************/
  // Column::AddSeqPositions()
  //
  // Adds a sequence positions to this column.
  /*****************************************************************/

  void AddSeqPosition (const PII &seqPos) {
    seqPositions.insert(seqPos);
  }

  void UpdateConsensus(MultiSequence* sequences, int numSequences, bool isNuc = false);

  /*****************************************************************/
  // Column::GetConsensus()
  //
  // Returns the current character consensus
  /*****************************************************************/

  char GetConsensus () {
    return consensus;
  }

  /*****************************************************************/
  // Column::GetMatchAccuracy()
  //
  // Returns the expected match accuracy of the column.
  // -- RKB 5/1/08
  /*****************************************************************/

  float GetMatchAccuracy () {
    return matchAccuracy;
  }


  /*****************************************************************/
  // Column::SetMatchAccuracy()
  //
  // Adds to column expected match accuracy.
  // -- RKB 5/1/08
  /*****************************************************************/

  void SetMatchAccuracy (float delta) {
    matchAccuracy = delta;
  }


  /*****************************************************************/
  // Column::operator<<()
  //
  // Output operator of the column.
  /*****************************************************************/

  friend ostream& operator<<(ostream& os,const Column& col) {
    os << "Column index: " << col.index << endl;
    os << "Visited: " << col.visited << endl;
    os << "Sequence positions: ";
    MII seqs = col.seqPositions;
    os << '[';
    for (MII::iterator iter = seqs.begin(); iter != seqs.end(); iter++)
      os << '(' << iter->first << ',' << iter->second << ')';
    os << ']' << endl;
    return os;
  } 

};

/*****************************************************************/
// Edge
//
// A class for storing information of 
// a candidate edge (match of column-pairs),
// and for calculating the edge weight.
/*****************************************************************/

class Edge {
 public:
  Column* sourceColumn;  
  Column* targetColumn;
  float weight;

  /*****************************************************************/
  // Edge::Edge()
  //
  // Constructor. Creates a new edge with pointers 
  // to two columns and the edge weight.
  /*****************************************************************/

  Edge (Column* c1, Column* c2, float initWeight) : sourceColumn(c1), targetColumn(c2), weight(initWeight) {};

  /*****************************************************************/
  // Edge::calcTgfWeight()
  //
  // Updates weights based on tgf, 
  // and returns the expected accuracy improvement of the edge.
  // (accuracy defined as for Column::matchAccuracy)
  /*****************************************************************/

  float calcTgfWeight (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, bool enableVerbose) {
    float sumPmatch = 0;
    float sumPgap = 0;
    while (sourceColumn->GetMergedInto() != sourceColumn) // Find current source column
      sourceColumn = sourceColumn->GetMergedInto();      
    while (targetColumn->GetMergedInto() != targetColumn) // Find current target column
      targetColumn = targetColumn->GetMergedInto();
    const MII &c1pos = sourceColumn->GetSeqPositions();
    const MII &c2pos = targetColumn->GetSeqPositions();
    for (MII::const_iterator c1posIter = c1pos.begin(); c1posIter != c1pos.end(); c1posIter++) { // loop over seqs i,j and associated seq positions ii,jj in columns joined by edge
      int i = c1posIter->first;   // sequence i
      int ii = c1posIter->second; // position ii
      for (MII::const_iterator c2posIter = c2pos.begin(); c2posIter != c2pos.end(); c2posIter++) {
	int j = c2posIter->first;   // sequence j
	if (i == j)    // Cannot match two columns with the same sequence
	  return INVALID_EDGE;
	int jj = c2posIter->second; // position jj
	SparseMatrix *ijMatrix = sparseMatrices[i][j];
	sumPmatch += ijMatrix->GetValue(ii,jj);
	sumPgap += ijMatrix->GetGapPosterior(0,ii); // so sumPgap really gets "overcounted": for each (seq,position) pair (i,ii) and (j,jj), add P(position ii in seq i is gapped) and P(position jj in seq j is gapped)
	sumPgap += ijMatrix->GetGapPosterior(1,jj); // 
      }
    }
    if (enableVerbose) {
      cerr << "previous weight= " << weight;
      cerr << " sumPmatch= " << sumPmatch << " sumPgap= " << sumPgap << endl;
    }
    weight = 2 * sumPmatch / sumPgap; // like eq (4) in the AMAP paper, but gap factor scaled s.t. value of 1 (instead of 1/2) gives AMA
    if (enableVerbose) {
      cerr << "new weight= " << weight << endl;
    }

    return 2 * sumPmatch - sumPgap;
  }

  /*****************************************************************/
  // Edge::calcMaxStepWeight()
  //
  // Updates weights based on maxstep, 
  // and returns the expected accuracy improvement of the edge. 
  /*****************************************************************/

  float calcMaxStepWeight (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, bool enableVerbose, float gapFactor) {
    float sumPmatch = 0;
    float sumPgap = 0;
    while (sourceColumn->GetMergedInto() != sourceColumn) // Find current source column
      sourceColumn = sourceColumn->GetMergedInto();
    while (targetColumn->GetMergedInto() != targetColumn) // Find current target column
      targetColumn = targetColumn->GetMergedInto();
    const MII &c1pos = sourceColumn->GetSeqPositions();
    const MII &c2pos = targetColumn->GetSeqPositions();
    for (MII::const_iterator c1posIter = c1pos.begin(); c1posIter != c1pos.end(); c1posIter++) {
      int i = c1posIter->first;
      int ii = c1posIter->second;
      for (MII::const_iterator c2posIter = c2pos.begin(); c2posIter != c2pos.end(); c2posIter++) {
	int j = c2posIter->first;
	if (i == j)    // Cannot match two columns with the same sequence
	  return INVALID_EDGE;
	int jj = c2posIter->second;
	SparseMatrix *ijMatrix = sparseMatrices[i][j];
	sumPmatch += ijMatrix->GetValue(ii,jj);
	sumPgap += ijMatrix->GetGapPosterior(0,ii);
	sumPgap += ijMatrix->GetGapPosterior(1,jj);
      }
    }
    if (enableVerbose) {
      cerr << "previous weight= " << weight;
      cerr << " sumPmatch= " << sumPmatch << " sumPgap= " << sumPgap << endl;
    }
    weight = (2 * sumPmatch - gapFactor * sumPgap) / (c1pos.size() * c2pos.size());
    if (enableVerbose) {
      cerr << "new weight= " << weight << endl;
    }
    return 2 * sumPmatch - sumPgap;
  }

  /*****************************************************************/
  // Edge::operator<()
  //
  // Compares two edges based on their weights.
  /*****************************************************************/

  bool operator< (Edge const e2) {
    return weight < e2.weight;
  }

  /*****************************************************************/
  // Edge::ostream&()
  //
  // Output operator for edges.
  /*****************************************************************/

  friend ostream& operator<<(ostream& os,const Edge& edge) {
    os << "sourceColumn: " << endl << *edge.sourceColumn;
    os << "targetColumn: " << endl << *edge.targetColumn;
    os << "weight: " << edge.weight << endl;
    return os;
  }
};


/*****************************************************************/
// greater_index
//
// An empty class. Used only for defining a binary comparison
// operator for column pointers.
/*****************************************************************/

class greater_index : binary_function<Column*, Column*, bool> {
 public:
  bool operator()(Column* x, Column* y) { return *y  < *x; }
};

/*****************************************************************/
// smaller_index
//
// An empty class. Used only for defining a binary comparison
// operator for column pointers.
/*****************************************************************/

class smaller_index : binary_function<Column*, Column*, bool> {
 public:
  bool operator()(Column* x, Column* y) { return *x  < *y; }
};

/*****************************************************************/
// smaller_weight
//
// An empty class. Used only for defining a binary comparison
// operator for edge pointers.
/*****************************************************************/

class smaller_weight : binary_function<Edge*, Edge*, bool> {
 public:
  bool operator()(Edge* x, Edge* y) { return x->weight  < y->weight; }
};

typedef map < pair <int,int>, Column* > MPIIC;

/*****************************************************************/
// MultiSequenceDag
//
// A class for storing a multiple alignment as a DAG,
// and for producing an alignment using sequence-annealing.
/*****************************************************************/

class MultiSequenceDag {
  list<Column*> columns;                // the current columns of the alignment
  list<Column*> oldColumns;             // old columns that have been merged with current columns
  MultiSequence *sequences;             // the sequences to be aligned 
  int numSequences;                     // num of sequences
  int alignLength;                      // num of columns in the alignment
  MPIIC seqPos2colIndex;                // mapping from sequence positions to columns
  float expectedAccuracy;               // the expected accuracy of the alignment
  bool outputForGUI;                    // true if the output is for the AMAP GUI
  float guiStartWeight;                 // the edge weight at which output for GUI starts
  int guiStepSize;                      // the step size for the GUI output
  bool isNuc;                           // true if nucleotides (false for proteins); for GUI output
                                        // -- RKB 4/19/08
 public:
  static const char pepGroup[26];             // peptide groups for GUI coloring (from TEXshade)
  static const bool pepSim[26][26];         // peptide similarities for GUI coloring (from TEXshade)
  static const char nucGroup[26];             // nucleotides
  static const bool nucSim[26][26];            // nucleotide similarities
                                             // -- RKB 4/19/08

 private:
  /*****************************************************************/
  // MultiSequenceDag::init()
  //
  // Initializes the DAG data structures.
  /*****************************************************************/

  void init (bool aligned) {
    if (aligned) {         // initialize DAG to input alignment
      for (int seqNum = 0; seqNum < numSequences; seqNum++) {
	if (alignLength < sequences->GetSequence(seqNum)->GetLength())
	  alignLength = sequences->GetSequence(seqNum)->GetLength();
      }
      for (int col = 1; col <= alignLength; col++) {
	columns.push_back(new Column(col));
      }
      for (int seqNum = 0; seqNum < numSequences; seqNum++) {
	Sequence *seq = sequences->GetSequence(seqNum);
	list<Column*>::iterator colIter = columns.begin();
	for (int col = 1, seqPos = 1; col <= seq->GetLength(); col++, colIter++) {
	  if (seq->GetPosition(col) == '-')
	    continue;
	  PII newSeqPos(seqNum,seqPos++);
	  (*colIter)->AddSeqPosition(newSeqPos);
	  seqPos2colIndex[newSeqPos] = *colIter;
	}
      }
    } else {             // initialize DAG to the null alignment
      alignLength = 0;
      expectedAccuracy = 0;
      for (int seqNum = 0; seqNum < numSequences; seqNum++) {
	Sequence *seq = sequences->GetSequence(seqNum);
	for (int col = 1, seqPos = 1; col <= seq->GetLength(); col++) {
	  if (seq->GetPosition(col) == '-')
	    continue;
	  columns.push_back(new Column(++alignLength));
	  PII newSeqPos(seqNum,seqPos++);
	  (*(--columns.end()))->AddSeqPosition(newSeqPos);
	  seqPos2colIndex[newSeqPos] = *(--columns.end());
	}	
      }
    }
  }


  /*****************************************************************/
  // MultiSequenceDag::DfsF()
  //
  // Implementation of the dfs-f() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  bool DfsF (Column* node, Column* upperBound,vector<Column*> &rForward) {
    node->Mark();
    rForward.push_back(node);
    push_heap(rForward.begin(),rForward.end(),greater_index());
    for (MII::const_iterator posIter = node->GetSeqPositions().begin(); posIter != node->GetSeqPositions().end(); posIter++) {
      PII u = PII(posIter->first,posIter->second + 1);
      if (seqPos2colIndex.find(u) == seqPos2colIndex.end()){  // reached end of the current sequence
	continue;
      }
      Column* w = seqPos2colIndex[u];
      if (*w == *upperBound){
	return true;          // found a cycle
      }
      if (!w->Marked() && *w < *upperBound && DfsF(w, upperBound, rForward))
	return true;          // found a cycle
    }
    return false;           // no cycles found
  }

  /*****************************************************************/
  // MultiSequenceDag::DfsB()
  //
  // Implementation of the dfs-b() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  void DfsB (Column* node, Column* lowerBound,vector<Column*> &rBackward) {
    node->Mark();
    rBackward.push_back(node);
    push_heap(rBackward.begin(),rBackward.end(),greater_index());
    for (MII::const_iterator posIter = node->GetSeqPositions().begin(); posIter != node->GetSeqPositions().end(); posIter++) {
      PII u = PII(posIter->first,posIter->second - 1);
      if (seqPos2colIndex.find(u) == seqPos2colIndex.end())  // reached end of the current sequence
	continue;
      Column*  w = seqPos2colIndex[u];
      if (!w->Marked() && *lowerBound < *w)
	DfsB(w, lowerBound, rBackward);
    }
  }

  /*****************************************************************/
  // MultiSequenceDag::Reorder()
  //
  // Implementation of the reorder() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  void Reorder(vector<Column*> &rForward, vector<Column*> &rBackward) {
    list<int> indexes;

    for (vector<Column*>::iterator rbIter = rBackward.begin(); rbIter != rBackward.end(); rbIter++) {
      indexes.push_back((*rbIter)->GetIndex());
    }
    for (vector<Column*>::iterator rfIter = rForward.begin(); rfIter != rForward.end(); rfIter++) {
      indexes.push_back((*rfIter)->GetIndex());
    }
    indexes.sort();
    list<int>::iterator idxIter = indexes.begin();
    for (unsigned i = 0; i < rBackward.size(); i++) {
      //      if (outputForGUI)
      //	cout << rBackward[0]->GetIndex() << ' ';
      rBackward[0]->SetIndex(*(idxIter++));
      pop_heap(rBackward.begin(),rBackward.end() - i,greater_index());
    }
    for (unsigned i = 0; i < rForward.size(); i++) {
      //      if (outputForGUI)
      //	cout << rForward[0]->GetIndex() << ' ';
      rForward[0]->SetIndex(*(idxIter++));
      pop_heap(rForward.begin(),rForward.end() - i,greater_index());
    }
    columns.sort(smaller_index());           // brute-foce. Can be done more efficiently using splice() in the previous steps.
  }

  /*****************************************************************/
  // MultiSequenceDag::Unmark()
  //
  // Unmarks all columns in the list.
  /*****************************************************************/

  void Unmark(vector<Column*> &rColist) {
    for (vector<Column*>::iterator rcIter = rColist.begin(); rcIter != rColist.end(); rcIter++)
      (*rcIter)->Unmark();
  }

  /*****************************************************************/
  // MultiSequenceDag::Merge()
  //
  // Merges the source and target columns of the edge 
  // into the source column.
  /*****************************************************************/

  //  void Merge (Edge *newEdge) {
  void Merge (Edge *newEdge, float delta) {
    Column* col1 = newEdge->sourceColumn;
    Column* col2 = newEdge->targetColumn;
    MII map1 = col1->GetSeqPositions();
    MII map2 = col2->GetSeqPositions();
    for (MII::iterator iter = map2.begin(); iter != map2.end(); iter++) {
      col1->AddSeqPosition(PII(iter->first,iter->second));
      seqPos2colIndex[PII(iter->first,iter->second)] = col1;
    }
    col2->SetMergedInto(col1);
    if (outputForGUI)
      col1->UpdateConsensus(sequences, numSequences, isNuc);

    // update the expected accuracy for the source column
    col1->SetMatchAccuracy (col1->GetMatchAccuracy() + delta);
  }

 public:
  
  /*****************************************************************/
  // MultiSequenceDag::MultiSequenceDag()
  //
  // Constructor.
  // Initialized alignment DAG. 
  // Uses input alignment if aligned is true.
  /*****************************************************************/

  MultiSequenceDag (MultiSequence *msa, bool aligned, bool forGUI, float start, int step, bool isNuc = false) : 
  sequences (msa), numSequences(msa->GetNumSequences()), outputForGUI(forGUI), guiStartWeight(start), guiStepSize(step), isNuc (isNuc) {
    init(aligned);
  }
  
  /*****************************************************************/
  // MultiSequenceDag::~MultiSequenceDag()
  //
  // Destructor.
  /*****************************************************************/

   ~MultiSequenceDag() {
    for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
      (*colIter)->SetMergedInto(NULL);
    }
    for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
      Column* colPtr = *colIter;
      *colIter = NULL;
      delete colPtr;
    }
    columns.clear();

    for (list<Column*>::iterator colIter = oldColumns.begin(); colIter != oldColumns.end(); colIter++) {
      (*colIter)->SetMergedInto(NULL);
    }
    for (list<Column*>::iterator colIter = oldColumns.begin(); colIter != oldColumns.end(); colIter++) {
      Column* colPtr = *colIter;
      *colIter = NULL;
      delete colPtr;
    }
    oldColumns.clear();
  }

  /*****************************************************************/
  // MultiSequenceDag::AddEdge()
  //
  // Adds a new edge to the alignment DAG.
  // Implementation of the add_edge() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

   //  int AddEdge (Edge *newEdge) {
   int AddEdge (Edge *newEdge, float delta) {
    Column* col1 = newEdge->sourceColumn;
    Column* col2 = newEdge->targetColumn;
    while (col1->GetMergedInto() != col1)  // get current source column
      col1 = col1->GetMergedInto();
    while (col2->GetMergedInto() != col2)  // get current target column
      col2 = col2->GetMergedInto();
    MII colSeqPos1 = col1->GetSeqPositions();
    MII colSeqPos2 = col2->GetSeqPositions();
    for (MII::iterator pos1Iter = colSeqPos1.begin(); pos1Iter != colSeqPos1.end(); pos1Iter++) {
      for (MII::iterator pos2Iter = colSeqPos2.begin(); pos2Iter != colSeqPos2.end(); pos2Iter++){
	if (pos1Iter->first == pos2Iter->first)
	  return 1;                                          // both columns contain positions from the same sequence
      }
    }
    
    vector<Column*> rForward, rBackward;
    Column *lBound, *uBound;
    if (*col1 < *col2) {
      lBound = col1;
      uBound = col2;
    } else {
      lBound = col2;
      uBound = col1;
    }

    if (DfsF(lBound,uBound,rForward)){        // new edge introduces a cycle in the DAG
      for (vector<Column*>::iterator rfIter = rForward.begin(); rfIter != rForward.end(); rfIter++)
	(*rfIter)->Unmark();
      return 2;
    }
    DfsB(uBound,lBound,rBackward);
    
    Unmark(rForward);
    Unmark(rBackward);
    if (rForward.size() == 1) {
      col1 = uBound;
      col2 = lBound;
    } else if (rBackward.size() == 1) {
      col1 = lBound;
      col2 = uBound;
    } else {
      //      if (outputForGUI)
      //	cout << "Reorder ";
      Reorder(rForward, rBackward);
      //      if (outputForGUI)
      //	cout << endl;
    }
    newEdge->sourceColumn = col1;
    newEdge->targetColumn = col2;
    //    if (outputForGUI)
    //      cout << "Merge " << col2->GetIndex() << ' ' << col1->GetIndex() << endl;
    //    Merge(newEdge);
    Merge (newEdge, delta);
    oldColumns.push_back(col2);              // keep pointers to old columns for later memory deallocation
    columns.remove(col2);
    return 0;
  }

  
  /*****************************************************************/
  // MultiSequenceDag::GetColumns()
  //
  // Returns the list of current columns in the the alignment DAG.
  /*****************************************************************/

  list<Column*> *GetColumns() {
    return &columns;
  }
  
  /*****************************************************************/
  // MultiSequenceDag::operator<<()
  //
  // Output operator.
  /*****************************************************************/

  friend ostream& operator<<(ostream& os,const MultiSequenceDag& msaDag) {
    for (list<Column*>::const_iterator iter = msaDag.columns.begin(); iter != msaDag.columns.end(); iter++)
      os << **iter;
    os << "seqPos2colIndex: \n";
    for (MPIIC::const_iterator iter = msaDag.seqPos2colIndex.begin(); 
	 iter != msaDag.seqPos2colIndex.end(); iter++)
      os << '(' << iter->first.first << ',' << iter->first.second << ',' << iter->second->GetIndex() << ')' << endl;
    return os;
  } 

  /*****************************************************************/
  // MultiSequenceDag::GetMultiSequence()
  //
  // Converts the alignment DAG to a standard MSA.
  /*****************************************************************/

  MultiSequence *GetMultiSequence () {
    SafeVector<SafeVector<char>::iterator> oldPtrs(numSequences);
    SafeVector<SafeVector<char> *> newPtrs(numSequences);

    // grab old data
    for (int i = 0; i < numSequences; i++){
      oldPtrs[i] = sequences->GetSequence(i)->GetDataPtr();
    }

    int newLength = columns.size();

    // build new alignments
    for (int i = 0; i < numSequences; i++){
      newPtrs[i] = new SafeVector<char>(); assert (newPtrs[i]);
      newPtrs[i]->push_back ('@');
    }

    // add all needed columns
    for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
      MII colPos = (*colIter)->GetSeqPositions();
      for (int j = 0; j < numSequences; j++) {
	if (colPos.find(j) != colPos.end())
	  newPtrs[j]->push_back(oldPtrs[j][colPos[j]]);
	else
	  newPtrs[j]->push_back('-');
      }
    }

    // wrap sequences in MultiSequence object
    MultiSequence *ret = new MultiSequence();
    for (int i = 0; i < numSequences; i++){
      ret->AddSequence (new Sequence(newPtrs[i], sequences->GetSequence(i)->GetName(), newLength,
                                      sequences->GetSequence(i)->GetSortLabel(), 
				     sequences->GetSequence(i)->GetLabel()));
    }
    return ret;
  }

  /*****************************************************************/
  // MultiSequenceDag::GetMatchAccuracyAnnot()
  //
  // Get the per-column expected match accuracy (formatted as a string for #=GC annotation).
  // -- RKB 5/1/08
  /*****************************************************************/

  sstring GetMatchAccuracyAnnot() {

    sstring annot;
    for (list<Column*>::iterator iter = columns.begin(); iter != columns.end(); iter++) {
      int n = (int) ((*iter)->GetSeqPositions()).size(); // # of non-gap chars in column
      if (n > 1) { // if at least 2 chars in column 
	float accuracy = (*iter)->GetMatchAccuracy() / (n * (n - 1)); // (convert to per-pair accuracy: divide by (n choose 2 = 1/2*n*(n-1)) and multiply by 2 to get num chars)
	char post = (accuracy > 0) ? "0123456789"[ min ((int) floor (accuracy * 10.0), 9) ] : '0'; // normalize to 0-9 (0 if expected accuracy < 0)
	annot << post;
      } else {
	annot << NOANNOT;
      }

    }

    return annot;
  }


  /*****************************************************************/
  // MultiSequenceDag::GetSequences()
  //
  // Returns the current sequences including gaps and color codes
  /*****************************************************************/

  string GetSequences() {
    string s = "";
    string colors = "";
    SafeVector<SafeVector<char>::iterator> oldPtrs(numSequences);

    // grab old data
    for (int i = 0; i < numSequences; i++){
      oldPtrs[i] = sequences->GetSequence(i)->GetDataPtr();
    }

    // add all needed columns
    for (int j = 0; j < numSequences; j++) {
      s += '>' + sequences->GetSequence(j)->GetName() + ' ';      
      colors += '@' + sequences->GetSequence(j)->GetName() + ' ';      
      for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
	MII colPos = (*colIter)->GetSeqPositions();
	if (colPos.find(j) != colPos.end()) {
	  s += oldPtrs[j][colPos[j]];
	  char consensus = (*colIter)->GetConsensus();
	  if (consensus) {
	    if (oldPtrs[j][colPos[j]] == consensus) 
	      colors += '1';
	    else if (oldPtrs[j][colPos[j]] == -consensus) 
	      colors += '2';
	    else if (isNuc && (nucGroup[(int) oldPtrs[j][colPos[j]] - 'A'] == consensus || (consensus >= 'A' && nucSim[(int) consensus - 'A'][(int) oldPtrs[j][colPos[j]] - 'A'])))  // -- RKB 4/19/08
	      colors += '3';
	    else if (!isNuc && (pepGroup[(int) oldPtrs[j][colPos[j]] - 'A'] == consensus || (consensus >= 'A' && pepSim[(int) consensus - 'A'][(int) oldPtrs[j][colPos[j]] - 'A']))) 
	      colors += '3';
	    else
	      colors += '0';
	  } else
	    colors += '0';
	}
	else {
	  s += '-';
	  colors += '0';
	}
      }
      s += '\n';
      colors += '\n';
    }
    return s + colors;
  }


  /*****************************************************************/
  // MultiSequenceDag::GetStockade()
  //
  // Get the Stockade corresponding to this alignment.
  // Mark up with column match accuracy scores.
  /*****************************************************************/

  Stockade GetStockade() {

    MultiSequence* msa = this->GetMultiSequence();
    Stockade stockade = msa->GetStockade();
    delete msa;
    sstring accuracy = GetMatchAccuracyAnnot();
    stockade.align.set_gc_annot ("Match_Accuracy", accuracy);

    return stockade;
  }

  /*****************************************************************/
  // MultiSequenceDag::AlignDag()
  //
  // Produces an alignment using the sequence-annealing method.
  // Input parameters include the posterior probabilities matrices,
  // the gap-factor used in the objective function, 
  // parameter for dynamic/static weights, 
  // tgf/maxweight weight function, and threshold for edge weight.
  /*****************************************************************/

  Stockade AlignDag (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, float gapFactor, 
		 bool enableVerbose, bool enableEdgeReordering, bool useTgf, float edgeWeightThreshold, bool show_intermediates) {

    priority_queue<Edge*, vector<Edge*>, smaller_weight> edges;
    Edge *edge;
    int guiFrame = 0;
    cerr << "Creating candidate edge list" << endl;
    if (outputForGUI && guiStartWeight == numeric_limits<float>::max()) {
      cout << "Weight " << numeric_limits<float>::max() << endl << 
	this->GetSequences() << endl;
    }

    for (int i = 0; i < numSequences; i++) {
      int seq1Length = sequences->GetSequence(i)->GetLength();
      for (int j = i + 1; j < numSequences; j++) {
	SparseMatrix* ijMatrix = sparseMatrices[i][j];
	for (int ii = 1; ii <= seq1Length; ii++) {
	  float pGapii = ijMatrix->GetGapPosterior(0,ii);
	  for (SafeVector<PIF>::iterator rowPtr = ijMatrix->GetRowPtr(ii), 
		 rowEnd = rowPtr + ijMatrix->GetRowSize(ii); rowPtr != rowEnd; rowPtr++) {
	    int jj = rowPtr->first;
	    float pMatch = 2 * rowPtr->second;
	    if (!pMatch)
	      continue;
	    float pGapjj = ijMatrix->GetGapPosterior(1,jj);
	    float weight = useTgf ? 2 * pMatch / (pGapii + pGapjj) : 2* pMatch - gapFactor * (pGapii + pGapjj);
	    if (weight < edgeWeightThreshold || (useTgf && weight < gapFactor))
	      continue;
	    edge = new Edge(seqPos2colIndex[PII(i,ii)], seqPos2colIndex[PII(j,jj)], weight);
	    edges.push(edge);
	  }
	}
      }
    }
    if (enableVerbose)
      cerr << "Adding edges to the DAG" << endl;
    while (!edges.empty()) {
      Edge *edge = edges.top();
      edges.pop();
      float delta = 0;
      if (enableEdgeReordering) { // recalculate edge weight if using dynamic edge weights
	delta = useTgf ? (edge->calcTgfWeight)(sparseMatrices,enableVerbose) : 
	  (edge->calcMaxStepWeight)(sparseMatrices,enableVerbose,gapFactor);
	if (delta == INVALID_EDGE) { // edge is no longer valid
	  delete edge;
	  continue;
	}
	if (edge->weight < edges.top()->weight) {
	  edges.push(edge);
	  if (enableVerbose)
	    cerr << "wrong order" << endl << *edge << endl;
	  continue;
	}
	if (edge->weight < edgeWeightThreshold || (useTgf && edge->weight < gapFactor)) {  // done. Delete remaining edges
	  while (!edges.empty()) {
	    edge = edges.top();
	    edges.pop();
	    delete edge;
	  }
	  break;
	}
      }
      if (enableVerbose)
	cerr << "Adding edge" << endl << *edge;
      //      int result = AddEdge(edge);
      int result = AddEdge (edge, delta);
      if (enableVerbose) {
	if (result)
	  cerr << "Failed to add edge with error code: " << result << endl;
	else
	  cerr << "Edge added successfully" << endl;
      }
      if (!result) { 
	expectedAccuracy += delta;
	if (outputForGUI && guiStartWeight >= edge->weight && !(guiFrame = ++guiFrame % guiStepSize)) {
	  cout << "Weight " << edge->weight << endl << 
	    this->GetSequences() << endl;
	}
	// write intermediates in Stockholm format
	// (but don't mess up GUI output!)
	// RKB -- 3/8/08
	if (show_intermediates && !outputForGUI) {
	  Stockade stockade = this->GetStockade();
	  sstring comment;
	  comment << "Alignment at edge-weight = " << edge->weight;
	  stockade.align.add_gf_annot ("AMAP", comment);
	  comment.clear();
	  comment << "Incremental expected match accuracy improvement = " << delta;
	  stockade.align.add_gf_annot ("AMAP", comment);
	  comment.clear();
	  comment << "Total improvement = " << expectedAccuracy;
	  stockade.align.add_gf_annot ("AMAP", comment);
	  stockade.align.write_Stockholm (cout);

	}
      }

      delete edge;
    }

    return this->GetStockade();
  }


};


                                               // A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
const char MultiSequenceDag::pepGroup[26]     =  {5,0,0,4,4,1,5,3,2,0,3,2,2,7,0,0,7,3,6,6,0,2,1,0,1,0};
const bool MultiSequenceDag::pepSim[26][26]   = {{1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0}, // A
						 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // B
						 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // C
						 {0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}, // D
						 {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}, // E
						 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0}, // F
						 {1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // G
						 {0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // H
						 {0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0}, // I
						 {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // J
						 {0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // K
						 {0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0}, // L
						 {0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0}, // M
						 {0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0}, // N
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}, // O
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}, // P
						 {0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0}, // Q
						 {0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // R
						 {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0}, // S
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0}, // T
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}, // U
						 {0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0}, // V
						 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0}, // W
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}, // X
						 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0}, // Y
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}};// Z


// nucGroup partitions the nucleotides into the purines and pyrimidines,
// {A,G,R} and {C,T,U,Y}
// nucSim says which nucleotides are similar
// (These are defined over the entire alphabet for coding convenience later.)
// -- RKB 4/19/08
                                               // A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
const char MultiSequenceDag::nucGroup[26]     =  {1,0,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,2,2,0,0,0,2,0};
const bool MultiSequenceDag::nucSim[26][26]   = {{1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // A
						 {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // B
						 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0}, // C
						 {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // D
						 {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // E
						 {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // F
						 {1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // G
						 {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // H
						 {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // I
						 {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // J
						 {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // K
						 {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // L
						 {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}, // M
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}, // N
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}, // O
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}, // P
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}, // Q
						 {1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // R
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0}, // S
						 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0}, // T
						 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0}, // U
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}, // V
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, // W
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}, // X
						 {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0}, // Y
						 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}};// Z


/*****************************************************************/
// Column::UpdateConsensus()
//
// Updates the consensus character
/*****************************************************************/
void Column::UpdateConsensus(MultiSequence* sequences, int numSequences, bool isNuc /* = false */) {
  unsigned int charCount[26] = {0}; // formerly pepCount
  unsigned int charGroupCount[10] = {0}; // formerly pepGroupCount

  unsigned int threshold = seqPositions.size() / 2;
  char groupConsensus = 0;
  
  if (!threshold) {
    consensus = 0;
    return;
  }

  for (MII::iterator posIter = seqPositions.begin(); posIter != seqPositions.end(); posIter++) {
    int currChar = sequences->GetSequence(posIter->first)->GetDataPtr()[posIter->second] - 'A';

    // if nucleotides
    // -- RKB 4/19/08
    if (isNuc) {
      if (++charCount[currChar] > threshold) {
	consensus = (char) currChar + 'A';
	if ((int) charCount[currChar] == numSequences)
	  consensus = -consensus;
      }
      if (++charGroupCount[(int) MultiSequenceDag::nucGroup[currChar]] > threshold) {
	groupConsensus = MultiSequenceDag::nucGroup[currChar];
      }
    }

    // if peptides
    else {
      if (++charCount[currChar] > threshold) {
	consensus = (char) currChar + 'A';
	if ((int) charCount[currChar] == numSequences)
	  consensus = -consensus;
      }
      if (++charGroupCount[(int) MultiSequenceDag::pepGroup[currChar]] > threshold) {
	groupConsensus = MultiSequenceDag::pepGroup[currChar];
      }
    }
  }

  if (!consensus) 
    consensus = groupConsensus;
}


#endif
