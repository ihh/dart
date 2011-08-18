#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<ctime>
#include <sys/stat.h>

#include "utils.h"
#include "ecfg/ecfgsexpr.h"
#include "seq/biosequence.h"
#include "util/macros.h"
#include "protpal/MyMap.h"

using namespace std;

string rep(int num, string in)
{
  string out="";
  for (int i=0; i< num; i++) out += in;
  return out;
}

double absoluted(bfloat in)
{
  if ( in<0.0) return -in; 
  else return in;
}


bfloat randomUnit(void)
{
  bfloat retVal; 
  bfloat rm = RAND_MAX;
  bfloat draw = rand();
  retVal = draw/(rm+1.0);
  if (retVal<=1.0)
    return retVal;
  else
    {
      cerr<<"Error - random unit returned a number larger than 1: " << retVal << endl; 
      exit(1);
    }
  return retVal; 
}


int maxIndex(vector<bfloat> &weights)
{
  int maxIdx=0; 
  bfloat max=0; 
  for (unsigned int i=0; i<weights.size(); i++)
	{
	  if (weights[i]>max) { maxIdx =i; max=weights[i];}
	}
  for (unsigned int i=0; i<weights.size(); i++)
	{
	  if (int(i)==maxIdx) continue; 
	  //if (weights[i] == max)std::cerr<<"Warning: non-unique max in maxIndex function\n";
	}
  return maxIdx; 
}


int sample(vector<float> &weights)
{
  vector<float> probs;
  double s = sum(weights);
  unsigned int i;
  float previous = 0 ;
  for (i=0; i<weights.size(); i++)	
	{
	  previous +=  weights[i]/s;
	  probs.push_back(previous);
	}
  
  float draw = randomUnit();
  //  std::cerr<<"number drawn: "<<draw<<endl;
  for (i=0; i<probs.size(); i++)
	{
	  if (draw <= probs[i]) return i;
	}
  return i; 
}

int sample(vector<double> &weights)
{
  vector<double> probs;
  double s = sum(weights);
  unsigned int i;
  float previous = 0 ;
  for (i=0; i<weights.size(); i++)	
	{
	  previous +=  weights[i]/s;
	  probs.push_back(previous);
	}
  
  float draw = randomUnit();
  //  std::cerr<<"number drawn: "<<draw<<endl;
  for (i=0; i<probs.size(); i++)
	{
	  if (draw <= probs[i]) return i;
	}
  return i; 
}


int sample(vector<bfloat> &weights)
{
  vector<bfloat> probs; //not really necessary, but good to be careful
  bfloat s = sum(weights);
  unsigned int i;
  bfloat previous = 0 ;
  for (i=0; i<weights.size(); i++)	
	{
	  previous +=  weights[i]/s;
	  probs.push_back(previous);
	}
  
  bfloat draw = randomUnit();
  //  std::cerr<<"number drawn: "<<draw<<endl;
  for (i=0; i<probs.size(); i++)
    if (draw <= probs[i] || i == probs.size()-1) return i;

  cerr<<"ERROR: reached end of sampling vector without choosing an element! The sampled number was " << draw << endl; 
  cerr<<"The weights vector  was "; displayVector(weights); 
  cerr<<"The probabilities vector  was "; displayVector(probs); 
  exit(1); 
  
  return i; 
}



MyMap<node, bool> merge(MyMap<node, bool>  *map1, MyMap<node, bool>  *map2)
{
  MyMap<node, bool> joinedMap; 
  MyMap<node, bool>::iterator it; 
  for (it = (*map1).begin(); it!=(*map1).end(); it++) joinedMap[(*it).first] = (*it).second;
  for (it = (*map2).begin(); it!=(*map2).end(); it++) joinedMap[(*it).first] = (*it).second;	    
  
  return joinedMap;
}

pair<int, int> merge(pair<int, int> coords1, pair<int, int> coords2)
{
  pair<int, int> joinedCoords; 
  joinedCoords.first = min(coords1.first, coords2.first);
  joinedCoords.second = max(coords1.second, coords2.second);
  return joinedCoords; 
}

double sum(vector<double> in)
{
  double out = 0;
  for (unsigned int i=0; i<in.size(); i++) out += in[i];
  return out;
}

double sum(vector<float> in)
{
  double out = 0;
  for (unsigned int i=0; i<in.size(); i++) out += in[i];
  return out;
}

bfloat sum(vector<bfloat> in)
{
  bfloat out = 0;
  for (unsigned int i=0; i<in.size(); i++) out += in[i];
  return out;
}


double sum(vector<int> in)
{
  double out = 0;
  for (unsigned int i=0; i<in.size(); i++) out += in[i];
  return out;
}

string stringAt(string in, int index)
{
  if (index > int(in.size())-1)
    {
      std::cerr<<"Error: Asking for a substring  beyond the size of the string\n"; 
      std::cerr<<"The  call was position " << index  << " in string " << in << endl; 
      exit(1);
    }
  string output;
  output = in[index];
  return output;
}

int index(string query, string in )
{
  string subStr;
  for (unsigned int i=0; i<in.size(); i++)
	{
	  subStr = in[i];
	  if (subStr == query) return(i);		   
	}
  return(-1);
}


int index(string query, vector<string> in )
{
  string subStr;
  for (unsigned int i=0; i<in.size(); i++)
	{
	  subStr = in[i];
	  if (subStr == query) return(i);		   
	}
  return(-1);
}

int index(sstring query, vector<sstring> in )
{
  sstring subStr;
  for (unsigned int i=0; i<in.size(); i++)
	{
	  subStr = in[i];
	  if (subStr == query) return(i);		   
	}
  return(-1);
}



int index(int query, vector<int> in )
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	if (in[i] == query) return(i);
	}
  return(-1);
}

int index(float query, vector<float> in )
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  if (in[i] == query) return(i);
	}
  return(-1);
}

int index(bfloat query, vector<bfloat> in )
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  if (in[i] == query) return(i);
	}
  return(-1);
}



bool in(string query, string in )
{
  if (index(query,in) == -1) return false;
  else return true;
}

bool in(int query, vector<int> in )
{
  if (index(query,in) == -1) return false;
  else return true;
}

bool in(float query, vector<float> in )
{
  if (index(query,in) == -1) return false;
  else return true;
}

void displayVector(vector<int> in)
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  //  std::cerr<<"\n";
}


void displayVector(vector<string> in)
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  std::cerr<<"\n";
}

void displayVector(vector<double> in)
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  std::cerr<<"\n";
}

void displayVector(vector<bfloat> in)
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  std::cerr<<"\n";
}
  
void displayVector(vector <vector <double> > in)
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  for (unsigned int j=0; j<in[i].size(); j++)
		{
		  std::cerr<<in[i][j]<<"  ";
		}
	  std::cerr<<"\n";
	}
}
void displayVector(vector <vector <int> > in)
{
  for (unsigned int i=0; i<in.size(); i++)
	{
	  for (unsigned int j=0; j<in[i].size(); j++)
		{
		  std::cerr<<in[i][j]<<"  ";
		}
	  std::cerr<<"\n";
	}
}

MyMap<string, string> parse_stockholm(const char* sequenceFileName, Alphabet& alphabet)
{
  MyMap<string, string> sequences;
  ifstream sequenceFileStream(sequenceFileName); 
  if (! sequenceFileStream.is_open())
    {
      std::cerr<<"\nERROR: could not open Stockholm file " << sequenceFileName << endl; 
      exit(1); 
    }
  Sequence_database seq_db;  // create the sequence database object
  Stockholm stk;
  stk.read_Stockholm(sequenceFileStream, seq_db); // read from file
  seq_db.seqs2dsqs (alphabet);   // parse the sequences into tokens using the Alphabet class that you read in with the substitution model
  for_const_contents (list<Named_profile>, seq_db, prof) 
    {
      sequences[prof->name] = (prof->seq); 
      if (sequences[prof->name].size() == 0)
	std::cerr<<"Warning: zero length sequence for name " << prof->name << endl; 
    }
  return sequences; 

}

vector<string> splitWhite(string in)
{
  string buf; // Have a buffer string
  stringstream ss(in); // Insert the string into a stream
  vector<string> tokens; // Create vector to hold our words
  
  while (ss >> buf)
    tokens.push_back(buf);
  return tokens;
}

vector<string> split(string in, string splitChar)
{
  vector<string> out; 
  string tmp; 
  
  for (unsigned int c=0; c<in.size(); c++)
    {
      if (stringAt(in, c)==splitChar)
	{
	  if (tmp == "") continue;
	  else
	    {
	    out.push_back(tmp); 
	    tmp = "";
	    }
	}
      else tmp += stringAt(in, c);
    }
  if (tmp != "") out.push_back(tmp); 
  return out;
}



MyMap<string, string> parse_fasta(const char* sequenceFileName, Alphabet& alphabet)
{
  MyMap<string, string> sequences;
  ifstream sequenceFileStream(sequenceFileName); 
  if (! sequenceFileStream.is_open())
    {
      std::cerr<<"\nERROR: could not open fasta file " << sequenceFileName << endl; 
      exit(1); 
    }
  Sequence_database seq_db;  // create the object
  seq_db.read_FASTA (sequenceFileStream);  // read from file
  seq_db.seqs2dsqs (alphabet);   // parse the sequences into tokens using the Alphabet class that you read in with the substitution model
  for_const_contents (list<Named_profile>, seq_db, prof) 
    {
      sequences[prof->name] = (prof->seq); 
      if (sequences[prof->name].size() == 0)
	std::cerr<<"Warning: zero length sequence for species " << prof->name << endl; 
    }

  return sequences; 
}

void seqDictSize(MyMap<string, string> seqDict)
{
  for (MyMap<string,string>::iterator seqIter =  seqDict.begin(); seqIter!=seqDict.end(); seqIter++)
    std::cerr<< seqIter->first << " " << (seqIter->second).size() << endl;
}




bool FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}

bool bfloat_is_nonzero(bfloat in)
{
  return  (in != in*2.0);
}


// Some codon model related functions. 
bool is_synonymous(string codon1, string codon2, map<sstring,sstring>& codon_table)
{
  return bool(codon_table[codon1] == codon_table[codon2]);
}

bool differ_more_than_one(string codon1, string codon2)
{
  if (codon1.size() == codon2.size() && codon1.size() ==  3)
    {
      int diffs = 0;
      for (int i=0; i<3; i++)
	{
	  if (stringAt(codon1,i) != stringAt(codon2,i))
	    diffs++;
	  if (diffs > 1)
	    return true; 
	}
      return false; 
    }
  else
    {
      cerr <<  "Improperly sized codons: " << codon1 << " "  << codon1.size() << " " << codon1 << " " << codon2.size(); 
  THROWEXPR("Wierd codons"); 
    }
  return true ; 
}

pair<string, string> find_first_difference(string codon1, string codon2)
{
  pair<string, string> out; 
  for (int i=0; i<3; i++)
    {
      if (stringAt(codon1,i) != stringAt(codon2,i))
	{
	  out.first = stringAt(codon1,i);
	  out.second = stringAt(codon2,i);
	  return out; 
	}
    }
  cerr <<"Warning: returning empty pair of differences between codons: "<<codon1 << " " <<codon2<<endl; 
  return out; 
}


bool is_transition(pair<string,string> nucPair)
{
  string nuc1 = nucPair.first; 
  string nuc2 = nucPair.second; 
  if (nuc1 == "a")
    return bool(nuc2=="g");

  if (nuc1 == "c")
    return bool(nuc2=="t");

  if (nuc1 == "g")
    return bool(nuc2=="a");

  if (nuc1 == "t")
    return bool(nuc2=="c");
  
  THROWEXPR("Unknown nucleotides queried for transition: " + nuc1 +" "+nuc2);
  return false; 
}

vector<sstring> all_codons(void)
{
  vector<sstring> out; 
  vector<sstring> dna; 
  string codon; 
  int i,j,k;
  dna.push_back("a");   dna.push_back("c");   dna.push_back("g");   dna.push_back("t"); 
  for (i=0; i<4; i++)
    {
      for (j=0; j<4; j++)
	{
	  for (k=0; k<4; k++)
	    {
	      codon.clear();
	      codon = dna[i] + dna[j] + dna[k];
	      out.push_back(codon); 
	    }
	}
    }
  return out; 
}

map<sstring, sstring> codon_table(void)
{
  map<sstring, sstring> out; 
  // hard-coded translation table.  Not a great solution, but will have to do for now
  out["ctt"]="l";
  out["atg"]="m";
  out["aca"]="t";
  out["acg"]="t";
  out["atc"]="i";
  out["aac"]="n";
  out["ata"]="i";
  out["agg"]="r";
  out["cct"]="p";
  out["ctc"]="l";
  out["agc"]="s";
  out["aag"]="k";
  out["aga"]="r";
  out["cat"]="h";
  out["aat"]="n";
  out["att"]="i";
  out["ctg"]="l";
  out["cta"]="l";
  out["act"]="t";
  out["cac"]="h";
  out["aaa"]="k";
  out["ccg"]="p";
  out["agt"]="s";
  out["cca"]="p";
  out["caa"]="q";
  out["ccc"]="p";
  out["tat"]="y";
  out["ggt"]="g";
  out["tgt"]="c";
  out["cga"]="r";
  out["cag"]="q";
  out["tct"]="s";
  out["gat"]="d";
  out["cgg"]="r";
  out["ttt"]="f";
  out["tgc"]="c";
  out["ggg"]="g";
  out["tag"]="stop";
  out["gga"]="g";
  out["tgg"]="w";
  out["ggc"]="g";
  out["tac"]="y";
  out["ttc"]="f";
  out["tcg"]="s";
  out["tta"]="l";
  out["ttg"]="l";
  out["tcc"]="s";
  out["acc"]="t";
  out["tca"]="s";
  out["gca"]="a";
  out["gta"]="v";
  out["gcc"]="a";
  out["gtc"]="v";
  out["gcg"]="a";
  out["gtg"]="v";
  out["gag"]="e";
  out["gtt"]="v";
  out["gct"]="a";
  out["tga"]="stop";
  out["gac"]="d";
  out["cgt"]="r";
  out["gaa"]="e";
  out["taa"]="stop";
  out["cgc"]="r";
  return out; 
}
