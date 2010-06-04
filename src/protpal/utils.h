#ifndef UTILS_H
#define UTILS_H
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include "algebras.h"
#include "ecfg/ecfgsexpr.h"

using namespace std;

class dart_rate_matrix: public Irrev_EM_matrix
{
 public:
  dart_rate_matrix(void);
};



typedef int state;
typedef int node;
// oops, Ian and I differ in capitalization.  I'll fix this soon.  
typedef int Node; 

char* stockholm_tree(const char*);
map<string, string> parse_stockholm(const char* );
map<string, string> parse_fasta(const char* );
vector<string> split(string,string); 

map<node, bool> merge(map<node, bool>  *map1, map<node, bool>  *map2); 


bfloat randomUnit(void);
double absoluted(bfloat);

int sample(vector<float>&);
int sample(vector<double>&);
int sample(vector<bfloat>&);


int maxIndex(vector<float>&);
int maxIndex(vector<double>&);
int maxIndex(vector<bfloat>&);

string rep(int num, string in); 
string stringAt(string in, int index);

int index(string query, string in );
int index(int query, vector<int> in );
int index(float query, vector<float> in );

bool in(string query, string in );
bool in(int query, vector<int> in );
bool in(float query, vector<float> in );



double sum(vector<double> in);
double sum(vector<float> in);
double sum(vector<int> in);
bfloat sum(vector<bfloat> in);

void displayVector(vector<int> in);

void displayVector(vector<string> in);

void displayVector(vector<double> in);
void displayVector(vector<bfloat> in);
  

void displayVector(vector <vector <double> > in);

void displayVector(vector <vector <int> > in);

#endif
