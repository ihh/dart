#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<ctime>

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
  bfloat rm = RAND_MAX;
  bfloat draw = rand();
  return draw/(rm+1.0);
}


int maxIndex(vector<bfloat> &weights)
{
  int maxIdx; 
  bfloat max=0; 
  for (int i=0; i<weights.size(); i++)
	{
	  if (weights[i]>max) { maxIdx =i; max=weights[i];}
	}
  for (int i=0; i<weights.size(); i++)
	{
	  if (i==maxIdx) continue; 
	  //if (weights[i] == max)std::cerr<<"Warning: non-unique max in maxIndex function\n";
	}
  return maxIdx; 
}


int sample(vector<float> &weights)
{
  vector<float> probs;
  double s = sum(weights);
  int i;
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
	  if (draw < probs[i]) return i;
	}
}

int sample(vector<double> &weights)
{
  vector<double> probs;
  double s = sum(weights);
  int i;
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
	  if (draw < probs[i]) return i;
	}
}


int sample(vector<bfloat> &weights)
{
  vector<bfloat> probs; //not really necessary, but good to be careful
  bfloat s = sum(weights);
  int i;
  bfloat previous = 0 ;
  for (i=0; i<weights.size(); i++)	
	{
	  previous +=  weights[i]/s;
	  probs.push_back(previous);
	}
  
  bfloat draw = randomUnit();
  //  std::cerr<<"number drawn: "<<draw<<endl;
  for (i=0; i<probs.size(); i++)
	{
	  if (draw < probs[i]) return i;
	}
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
  for (int i=0; i<in.size(); i++) out += in[i];
  return out;
}

double sum(vector<float> in)
{
  double out = 0;
  for (int i=0; i<in.size(); i++) out += in[i];
  return out;
}

bfloat sum(vector<bfloat> in)
{
  bfloat out = 0;
  for (int i=0; i<in.size(); i++) out += in[i];
  return out;
}


double sum(vector<int> in)
{
  double out = 0;
  for (int i=0; i<in.size(); i++) out += in[i];
  return out;
}

string stringAt(string in, int index)
{
  if (index > in.size()-1)
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
  for (int i=0; i<in.size(); i++)
	{
	  subStr = in[i];
	  if (subStr == query) return(i);		   
	}
  return(-1);
}


int index(string query, vector<string> in )
{
  string subStr;
  for (int i=0; i<in.size(); i++)
	{
	  subStr = in[i];
	  if (subStr == query) return(i);		   
	}
  return(-1);
}

int index(sstring query, vector<sstring> in )
{
  sstring subStr;
  for (int i=0; i<in.size(); i++)
	{
	  subStr = in[i];
	  if (subStr == query) return(i);		   
	}
  return(-1);
}



int index(int query, vector<int> in )
{
  for (int i=0; i<in.size(); i++)
	{
	if (in[i] == query) return(i);
	}
  return(-1);
}

int index(float query, vector<float> in )
{
  for (int i=0; i<in.size(); i++)
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
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  //  std::cerr<<"\n";
}


void displayVector(vector<string> in)
{
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  std::cerr<<"\n";
}

void displayVector(vector<double> in)
{
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  std::cerr<<"\n";
}

void displayVector(vector<bfloat> in)
{
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in[i]<<"  ";
	}
  std::cerr<<"\n";
}
  
void displayVector(vector <vector <double> > in)
{
  for (int i=0; i<in.size(); i++)
	{
	  for (int j=0; j<in[i].size(); j++)
		{
		  std::cerr<<in[i][j]<<"  ";
		}
	  std::cerr<<"\n";
	}
}
void displayVector(vector <vector <int> > in)
{
  for (int i=0; i<in.size(); i++)
	{
	  for (int j=0; j<in[i].size(); j++)
		{
		  std::cerr<<in[i][j]<<"  ";
		}
	  std::cerr<<"\n";
	}
}

MyMap<string, string> parse_stockholm(const char* fileName, Alphabet alphabet)
{
  MyMap<string, string> sequences; 
  string line;
  ifstream seqFile(fileName);

  // Parse stockholm file
  if (seqFile.is_open())
    {
      while (! seqFile.eof() )
	{
	  getline(seqFile,line);
	  //std::cerr<<"Line: "<<line<<endl;
	  if (stringAt(line,0) == "#" || splitWhite(line).size()<2) continue;
	  else
	    {
	      //std::cerr<<"splitting line...\n";
	      sequences[splitWhite(line)[0]] += splitWhite(line)[1]; 
	      //sequences.insert(pair<string, string>(splitWhite(line)[0], sequences[splitWhite(line)[0]] + splitWhite(line)[1]));
	    }
	}
      seqFile.close();
    }
  else 
    {
      std::cerr << "\nERROR: Unable to open Stockholm file: "<<fileName<< "\n"; 
      exit(1); 
    }
  for (MyMap<string, string>::iterator seqIter = sequences.begin(); seqIter!=sequences.end(); seqIter++)
    if ( (seqIter->second).size() == 0)
      std::cerr<<"Warning: 0 length sequence for species: " << seqIter->first << endl; 
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
  
  for (int c=0; c<in.size(); c++)
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



MyMap<string, string> parse_fasta(const char* sequenceFileName, Alphabet alphabet)
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
