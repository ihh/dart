#include<iostream>
#include<string>
#include<vector>
#include<map>

#include "util/sstring.h"
#include "seq/biosequence.h"
#include "protpal/utils.h"
#include "protpal/AlignmentEnvelope.h"

using namespace std;

AlignmentEnvelope::AlignmentEnvelope(void)
{
  // Placeholder constructor
}
void AlignmentEnvelope::build_index(sstring guide_alignment_filename, sstring gap_char_in, int sausage_size)
{
  map<string, string>::iterator seqIter, seqIter2; 
  unsigned int seqSize; 
  pair<int, int> envPair; 
  Alphabet dummy_alpha("Uninitialized", 0); 
  alignment = parse_stockholm(guide_alignment_filename.c_str(), dummy_alpha); 
  gap_char = gap_char_in; 
  
  for (seqIter = alignment.begin(); seqIter != alignment.end(); seqIter++)
    {
      seqSize = (seqIter->second).size(); 
      for (unsigned int idx = 0; idx < seqSize; idx++)
	{
	  if ( stringAt((seqIter->second), idx) != gap_char)
	    {
	      for (seqIter2 = alignment.begin(); seqIter2 != alignment.end(); seqIter2++)
		{
		  if  ( seqIter->first == seqIter2->first ) 
		    continue;
		  else
		    {
		      envPair.first = 
			sequence_coordinate(
			   next_match_column(idx, "backward", seqIter2->first, seqIter->first), seqIter2->first);
		      envPair.second = 
			sequence_coordinate(
			   next_match_column(idx, "forward", seqIter2->first, seqIter->first), seqIter2->first);
		      coordinates[seqIter->first][ sequence_coordinate(idx, seqIter->first) ][ seqIter2->first ] = 
			envPair; 
		    }
		}
	    }
	}
    }
  alignment.clear(); 
}

int AlignmentEnvelope::next_match_column(int alignment_coordinate, string direction, string seqName1, string seqName2)
{
  int matches_found = 0; 
  while ( matches_found < sausage_size )
    {
      if ( stringAt(alignment[seqName1], alignment_coordinate) != gap_char &&
	   stringAt(alignment[seqName2], alignment_coordinate) != gap_char )
	++matches_found; 
      if (direction == "forward")
	++alignment_coordinate; 
      else
	--alignment_coordinate; 
    }
  return alignment_coordinate; 
}
	
int AlignmentEnvelope::sequence_coordinate(int alignment_coordinate, string seqName)
{
  int seqCoord = -1; 
  for(int alignCoord = alignment_coordinate; alignCoord >=0; alignCoord--)
    if ( stringAt(alignment[seqName],alignCoord) != gap_char )
      ++seqCoord; 
  return seqCoord; 
}

pair<int, int> AlignmentEnvelope::lookup(string seqName, int index, string otherSeqName)
{
  bool debug = true; 
  // might need to add some checking here...
  if (debug)
    {
      if ( ! coordinates.count(seqName) )
	std::cerr<<"this sequence name not found: " << seqName << endl; 
      else if ( ! coordinates[seqName].count(index) )
	std::cerr<< "Index " << index  << " in sequence " << seqName << " not found!\n"; 
      else if ( ! coordinates[seqName][index].count(otherSeqName) )
	std::cerr<< "other sequence " << otherSeqName << " in index " << index  << " in sequence " << seqName << " not found!\n"; 
      else 
	return coordinates[seqName][index][otherSeqName]; 
    }
  return coordinates[seqName][index][otherSeqName]; 
}

void AlignmentEnvelope::show(ostream& out)
{
//   for (map<string, string>::iterator seqIter = alignment.begin(); seqIter != alignment.end(); seqIter++)
//     {
//       out << "Displaying homologies for sequence: " << seqIter->first << endl; 
// 	  for (map<string, string>::iterator seqIter2 = alignment.begin(); seqIter2 != alignment.end(); seqIter2++)
// 	    {
// 	      if (seqIter2->first == seqIter->first)
// 		continue; 
// 	      std::cerr << "\tHomologies with sequence " << seqIter2->first << endl; 
// 	      for (unsigned int seqIdx = 0; seqIdx < (seqIter->second).size(); seqIdx++)
// 		{
// 		  std::cerr<< sequence_coordinate(seqIdx, seqIter->first) << " -> " << 
// 		    lookup(seqIter->first, sequence_coordinate(seqIdx, seqIter->first), seqIter2->first) << endl; 
// 		}
// 	    }
//     }
}
	      
	
