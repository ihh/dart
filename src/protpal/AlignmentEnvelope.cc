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
}
void AlignmentEnvelope::build_index(sstring guide_alignment_filename, sstring gap_char_in)
{
  map<string, string>::iterator seqIter, seqIter2; 
  unsigned int seqSize; 
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
		    coordinates[seqIter->first][ sequence_coordinate(idx, seqIter->first) ][ seqIter2->first ] = 
		      sequence_coordinate(idx, seqIter2->first);
		}
	    }
	}
    }
}

int AlignmentEnvelope::sequence_coordinate(int alignment_coordinate, string seqName)
{
  int seqCoord = -1; 
  for(int alignCoord = alignment_coordinate; alignCoord >=0; alignCoord--)
    if ( stringAt(alignment[seqName],alignCoord) != gap_char )
      seqCoord++; 
  return seqCoord; 
}

int AlignmentEnvelope::lookup(string seqName, int index, string otherSeqName)
{
  // might need to add some checking here...
  return coordinates[seqName][index][otherSeqName]; 
}

void AlignmentEnvelope::show(ostream& out)
{
  for (map<string, string>::iterator seqIter = alignment.begin(); seqIter != alignment.end(); seqIter++)
    {
      out << "Displaying homologies for sequence: " << seqIter->first << endl; 
	  for (map<string, string>::iterator seqIter2 = alignment.begin(); seqIter2 != alignment.end(); seqIter2++)
	    {
	      if (seqIter2->first == seqIter->first)
		continue; 
	      std::cerr << "\tHomologies with sequence " << seqIter2->first << endl; 
	      for (unsigned int seqIdx = 0; seqIdx < (seqIter->second).size(); seqIdx++)
		{
		  std::cerr<< sequence_coordinate(seqIdx, seqIter->first) << " -> " << 
		    lookup(seqIter->first, sequence_coordinate(seqIdx, seqIter->first), seqIter2->first) << endl; 
		}
	    }
    }
}
	      
	
