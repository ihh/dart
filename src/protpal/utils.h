#ifndef UTILS_H
#define UTILS_H
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include <sys/stat.h>

#include "algebras.h"
#include "ecfg/ecfgsexpr.h"
#include "seq/biosequence.h"
#include "util/sstring.h"
#include "protpal/MyMap.h"
using namespace std;

class dart_rate_matrix: public Irrev_EM_matrix
{
 public:
  dart_rate_matrix(void);
};

typedef int state;
typedef int node;
// oops, DART and I differ in capitalization.  I'll fix this soon (maybe). 
typedef int Node; 
typedef vector<bool>  Row;
typedef pair<int,int> Row_pair;
typedef map<Row_pair,Alignment_path> Decomposition;

char* stockholm_tree(const char*);
MyMap<string, string> parse_stockholm(const char*, Alphabet );
MyMap<string, string> parse_fasta(const char*, Alphabet );
vector<string> split(string,string); 
vector<string> splitWhite(string); 

MyMap<node, bool> merge(MyMap<node, bool>  *map1, MyMap<node, bool>  *map2); 
void seqDictSize(MyMap<string, string>);                                                                                                                                     
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
int index(string query, vector<string> in );
int index(sstring query, vector<sstring> in );
int index(int query, vector<int> in );
int index(float query, vector<float> in );
int index(bfloat query, vector<bfloat> in );


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

bool FileExists(string); 

#endif
