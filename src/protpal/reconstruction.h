#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include<iostream>
#include<string>
#include<map>

#include "protpal/utils.h"
#include "protpal/profile.h"
#include "protpal/transducer.h"
#include "tree/phylogeny.h"
#include "ecfg/ecfgsexpr.h"

class Reconstruction
{
 public:
  // Constructor  - not much happens here.
  Reconstruction(void);

  // store the possible cmd line options, and brief descriptions
  map<string, string> options; 

  // Set command line options, get data (sequences, tree) pointed to, etc. 
  void get_cmd_args(int argc, char* argv[]);
  
  // Each leaf node is assigned a sequence. Possibly delete these after making ExactMatch transducers for each
  // leaf.  Hmm, maybe I should have used pointers instead. 
  map<string, string> sequences; 

  // not (yet) using Ian's fancy Alphabet class
  string alphabet; 

  // the rate matrix is from DART - imported on top-level
  sstring rate_matrix_filename;   

  // the phylogenetic tree, stored as a PHYLYP_tree object.  This is a fairly hairy object, with most of its code
  // living in phylogeny.*
  PHYLIP_tree tree; 

  // Each tree node will at some point contain a sequence profile.  The profile at a node  is constructed via its
  // two children, which are viciously deleted after the parent is constructed.  
  map<node, AbsorbingTransducer> profiles; 


  // default naming of internal nodes.  (not yet implemented)
  vector<string> get_node_names(void);
  
  // simulation
  bool simulate; 
  void make_sexpr_file(Alphabet, Irrev_EM_matrix);
  string show_branch(node); 
  void simulate_alignment(Alphabet, Irrev_EM_matrix);
  int rootLength; 
  
  //simulation utility functions
  string sample_pairwise(string, BranchTrans); 
  string sample_root(SingletTrans); 
  
  
  // Reconstruction algorithm parameters
  int num_sampled_paths;
  bool show_alignments; 
  int loggingLevel;
  int envelope_distance; 
  bool leaves_only; 

  
  //model parameters
  double ins_rate;
  double del_rate;   
  double gap_extend; 

 private:
  // Help message
  void display_opts(void);
  
  // loading data 
  bool have_tree;
  bool have_sequences;   
  
  void load_rate_matrix(const string);
  void loadTreeString(const char*);
  void get_stockholm_tree(const char*);  
};
#endif
