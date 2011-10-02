#include<iostream>
#include<string>
#include<fstream>
#include<sstream>

#include "util/unixenv.h"
#include "util/sstring.h"
#include "util/opts_list.h"
#include "util/sexpr.h"
#include "protpal/ReadProfileScore.h"
#include "seq/biosequence.h"
#include "util/dexception.h"
#include "util/sstring.h"
#include "protpal/utils.h"
#include "protpal/algebras.h"

class PlacementLimiter : public map<string, vector<Node> >
{
 public:
  PlacementLimiter(void); 
  void parse_JSON(const char* JSON_file); 
  map<int, string> branches2nodes; 
  map<string, int> nodes2branches; 
  map<string, vector<int> > read_node_map;
  bool is_allowed(string, string); 
  string treeString; 
};


class Placement
{
 public:
  // Top level fns
  Placement (int argc, char* argv[]); 
  void Run(); 

  // Key internal data structures
  ScoreMap scores; 

  // Key input data
  map<int, string> profile_filenames; 
  map<int, string> check_profile_filenames(); 
  map<string, string> reads; 

  // General args
  double minPostProb; 
  stringstream args; 
  sstring readsFileName; 
  sstring savedPosteriorProfilesDir; 
  sstring treeFileName; 
  sstring default_chain_filename; 
  sstring rateMatrixFileName; 
  sstring jsonPlacementsFileName; 
  void write_placement_JSON(ostream&, ScoreMap&); 
  void write_numbered_newick(ostream&, bool quotes=true); 
  int loggingLevel;   

  // Model/tree related 
  Alphabet alphabet; 
  Irrev_EM_matrix rate_matrix; 
  PHYLIP_tree tree; 
  void parse_tree(const char*); 
  vector<string> alphabetVector;

  // Limiting placements
  bool using_limiter; 
  sstring jsonInputFileName; 
  PlacementLimiter placementLimiter; 
};

