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
  // Default constructor
}
bool AlignmentEnvelope::is_not_gap(int row, int col)
{
  return stk.path(row, col);
}
bool AlignmentEnvelope::is_gap(int row, int col)
{
  return !is_not_gap(row,col); 
}

int AlignmentEnvelope::name2int(sstring name)
{
  // this is temporary...
  for (int i=0; i < stk.rows(); ++i)
    if (name==stk.row_name[i])
      return i; 
  THROWEXPR("No row found with the name: "+name);
}


void AlignmentEnvelope::build_index(sstring guide_alignment_filename, sstring gap_char_in, int sausage_size_in, sstring type_in)
{
  int  seq1, seq2; 
  sstring name1, name2;
  //int seqSize; 
  pair<int, int> envPair; 
  Alphabet dummy_alpha("Uninitialized", 0); 
  Sequence_database seq_db;  // create the sequence database object
  ifstream sequenceFileStream(guide_alignment_filename.c_str());
  // Stockholm object stores the alignment. 
  stk.read_Stockholm(sequenceFileStream, seq_db); // read from file       
  stk.assert_flush(); 

  int rows = stk.rows(); 
  int cols = stk.columns(); 
  gap_char = gap_char_in; 
  type = type_in; 
  sausage_size = sausage_size_in; 
  for (seq1 = 0; seq1 < rows; ++seq1)
    {
      name1 = stk.row_name[seq1]; 
      //      std::cerr<<"\nIndexing from " << seqIter->first << " to: \n"; 
      //seqSize = (seqIter->second).size(); 
      for (seq2 = 0; seq2 < rows; ++seq2)
	{
	  if  ( seq1 == seq2)
	    continue;
	  else
	    {
	      name2 = stk.row_name[seq2]; 
	      //	      std::cerr << " " << seqIter2->first << " "; 
	      for (int idx = 0; idx < cols; idx++)
		{
		  if (is_not_gap(seq1, idx))
		    {
		      if (type == "basic")
			{
			  envPair.first = 
			    sequence_coordinate( idx, name2 ) - sausage_size; 
			  envPair.second = 
			    sequence_coordinate( idx, name2 ) + sausage_size; 
			}
		      else
			{
			  envPair.first = 
			    sequence_coordinate(next_match_column(idx, "backward", name2, name1), name2);
			  envPair.second = 
			    sequence_coordinate(next_match_column(idx, "forward", name2, name1), name2);
			}
		      coordinates[name1][ sequence_coordinate(idx, name1) ][ name2 ] = envPair; 
		      
		    }
		}
	    }
	}
    }
  alignment.clear(); 
  if (coordinates.size() == 0)
    {
      THROWEXPR("\nERROR - no coordinates stored in alignment envelope.\n This is most likely due to an empty or incorrectly formatted guide alignment.  Guide alignment must be in Stockholm format (see www.biowiki.org/StockholmFormat for specifications)\n")
    }
}

int AlignmentEnvelope::next_match_column(int alignment_coordinate, string direction, string seqName1, string seqName2)
{
  bool logging = false; 
  if (logging)
    std::cerr<< "Beginning from alignment coord: " << alignment_coordinate << " "; 
  int matches_found = 0; 
  while ( matches_found < sausage_size )
    {
      if ( is_not_gap(name2int(seqName1), alignment_coordinate) &&
	   is_not_gap(name2int(seqName2), alignment_coordinate))
	//      if ( stringAt(alignment[seqName1], alignment_coordinate) != gap_char &&
	//	   stringAt(alignment[seqName2], alignment_coordinate) != gap_char )
	++matches_found; 
      if (direction == "forward")
	{
	  if ( alignment_coordinate == int(alignment[seqName1].size())-1)
	    return alignment_coordinate; 
	  else
	    ++alignment_coordinate; 
	}
      else
	{
	  if (alignment_coordinate == 0 )
	    return alignment_coordinate; 
	  else
	    --alignment_coordinate; 
	}
      if (logging)
	std::cerr<<" " << alignment_coordinate << " ";
    }
  if (logging)
    std::cerr<< "\nAlignment coordinate " << sausage_size << " match columns " << direction << " between " << 
      seqName1 << " and " << seqName2 << " was " << alignment_coordinate<<endl; 
  return alignment_coordinate; 
}
	
int AlignmentEnvelope::sequence_coordinate(int alignment_coordinate, string seqName)
{
  int seqCoord = -1; 
  for(int alignCoord = alignment_coordinate; alignCoord >=0; alignCoord--)
    if ( is_not_gap(name2int(seqName), alignCoord))
      //    if ( stringAt(alignment[seqName],alignCoord) != gap_char )
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
	      
	
