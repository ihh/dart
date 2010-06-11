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



map<node, bool> merge(map<node, bool>  *map1, map<node, bool>  *map2)
{
  map<node, bool> joinedMap; 
  map<node, bool>::iterator it; 
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
	  std::cerr<<in.at(i)<<"  ";
	}
  //  std::cerr<<"\n";
}


void displayVector(vector<string> in)
{
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in.at(i)<<"  ";
	}
  std::cerr<<"\n";
}

void displayVector(vector<double> in)
{
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in.at(i)<<"  ";
	}
  std::cerr<<"\n";
}

void displayVector(vector<bfloat> in)
{
  for (int i=0; i<in.size(); i++)
	{
	  std::cerr<<in.at(i)<<"  ";
	}
  std::cerr<<"\n";
}
  
void displayVector(vector <vector <double> > in)
{
  for (int i=0; i<in.size(); i++)
	{
	  for (int j=0; j<in.at(i).size(); j++)
		{
		  std::cerr<<in.at(i).at(j)<<"  ";
		}
	  std::cerr<<"\n";
	}
}
void displayVector(vector <vector <int> > in)
{
  for (int i=0; i<in.size(); i++)
	{
	  for (int j=0; j<in.at(i).size(); j++)
		{
		  std::cerr<<in.at(i).at(j)<<"  ";
		}
	  std::cerr<<"\n";
	}
}

map<string, string> parse_stockholm(const char* fileName)
{
  map<string, string> sequences; 
  string line;
  ifstream seqFile(fileName);

  // Parse stockholm file
  if (seqFile.is_open())
    {
      while (! seqFile.eof() )
	{
	  getline(seqFile,line);
	  //std::cerr<<"Line: "<<line<<endl;
	  if (stringAt(line,0) == "#" || split(line," ").size()<2) continue;
	  else
	    {
	      //std::cerr<<"splitting line...\n";
	      sequences[split(line," ")[0]] += split(line," ")[1];
	    }
	}
      seqFile.close();
    }
  else 
    {
      std::cerr << "Error: Unable to open file: "<<fileName;
    }
  return sequences; 
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



map<string, string> parse_fasta(const char* fileName)
{
  map<string, string> sequences;
  return sequences; 
}
/*
  string current; 
  // Parse fasta file
  if (seqFile.is_open())
	{
	  while (! seqFile.eof() )
		{
		  getline(seqFile,line);
		  if (stringAt(0,line) == "#" || split(line," ").size()<2) continue;
		  else
			{
			  sequences[split(line," ")[0]] = split(line," ")[1];
			}
		}
	  seqFile.close();
	}
  else 
	{
	  std::cerr << "Error: Unable to open file: "<<fileName;
	}
}
*/
