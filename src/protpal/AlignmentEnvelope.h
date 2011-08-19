#ifndef ALIGNMENT_ENVELOPE_H
#define ALIGNMENT_ENVELOPE_H
#include<iostream>
#include<string>
#include<vector>
#include<map>

#include "util/sstring.h"
#include "seq/biosequence.h"

using namespace std;


class AlignmentEnvelope
{
 public:
  // constructor, initializer
  AlignmentEnvelope(void); 
  void build_index(sstring guide_alignment_filename, sstring gap_char = "-", int sausage_size=0, sstring type="gap-sliding");
  // main data storage unit - nested maps
  map< string, map<int, map<string, pair<int, int> > > > coordinates; 
  pair<int, int> lookup(string seqName, int index, string otherSeqName);
  Stockholm stk; 
  Sequence_database seq_db; 

  void show(ostream& out); 
 private:
  bool is_gap(int, int); 
  bool is_not_gap(int, int); 
  int name2int(sstring);

  int sausage_size; 
  sstring type; 
  map<string, string> alignment;  
  sstring gap_char; 
  // convert from alignment coords to sequence coords
  int sequence_coordinate(int alignment_coordinate, string seqName); 
  int next_match_column(int alignment_coordinate, string direction, string seqName1, string seqName2); 
};

#endif
