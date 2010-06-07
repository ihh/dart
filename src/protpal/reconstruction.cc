#include<iostream>
#include<string>
#include<fstream>
#include<sstream>

#include "protpal/utils.h"
#include "protpal/profile.h"
#include "protpal/reconstruction.h"
#include "tree/phylogeny.h"
#include "util/piper.h"
#include "ecfg/ecfgsexpr.h"
#include "util/unixenv.h"

#define DEFAULT_CHAIN_FILE "data/handalign/prot1.hsm"

using namespace std; 
Reconstruction::Reconstruction(void)
{
  options["-stk"] = "<filename> (Unaligned) sequences in Stockholm file format.  If there is a #=GF NH line, this will be used as the phylogenetic tree, though this can be overridden by the -t option. \n";

  //  options["-fa"] = "<filename> (Unaligned) sequences in fasta file format\n";

  options["-t"] = "<string> Tree in newick file format\n";

  options["-i"] = "<float> Insert rate (default .02) \n";
  ins_rate = .019; 
  options["-d"] = "<float> Delete rate (default .02)\n";    
  del_rate = .02; 
  options["-g"] = "<float> Gap-extend probability (default .4) \n";
  gap_extend = .912; 

  options["-n"] = "<int> Number of paths to sample in traceback (default 10) \n";
  num_sampled_paths = 10;

  options["-e"] = "<int> Maximum allowed distance between aligned leaf characters (default 100)\n";
  envelope_distance = 100; 

  options["-g"] = "<grammar file> DART format chain file to be used for match states absorb/emit likelihoods.  Default is data/handalign/prot1.hsm ";
  rate_matrix_filename << Dart_Unix::get_DARTDIR() << '/' << DEFAULT_CHAIN_FILE; 
  
  num_sampled_paths = 10; 

  options["-sa"] = "<bool> Show sampled alignments for internal nodes (default false) \n";
  show_alignments = false;

  options["-l"] = "<bool> Show only leaves when displaying alignments (default false) \n";
  leaves_only = false; 
 
  options["-log"] = 
	"<int> Show log messages of varying levels of verbosity. (default = level 1 )"
"\n\n\tLevel  0 : Show no log messages.  The alignment is sent to standard out."
  "\n\tLevel  1 : Top-level logging related to profile creating/sampling, etc.  This is default"
  "\n\tLevel  2 : Detailed logging including DP, traceback, etc (dense)"    
  "\n\tLevel  3 : Very detailed logging including transducer composition algorithms. (dense)\n";
  
  loggingLevel = 1;  

  have_tree=false; 
  have_sequences=false; 

}


void Reconstruction::get_cmd_args(int argc, char* argv[])
{
  string error = ""; 
  bool all_reqd_args = true; 
  for (int i = 0; i < argc; i++) 
	{ 
	  if (string(argv[i]) == "-log")
		{
		  const char* logChar = argv[i+1]; 
		  loggingLevel = atoi(logChar);
		  if( loggingLevel>=1)
			std::cout<<"Logging level set as "<< loggingLevel<<endl; 
		}
	}
  for (int i = 0; i < argc; i++) 
	{ 
	  if(loggingLevel>=2) std::cout<<"Examining argument: "<<argv[i]<<endl; 
	  if (string(argv[i]) == "-stk") 
		{
		  const char* sequenceFileName = argv[i + 1];
		  sequences = parse_stockholm(sequenceFileName); 
		  have_sequences = true; 
		  get_stockholm_tree(sequenceFileName); 
		}
	  else if (string(argv[i]) == "-fa") 
		{
		  const char* sequenceFileName = argv[i + 1];
		  sequences = parse_fasta(sequenceFileName); 
		  have_sequences = true; 
		}
	  else if (string(argv[i]) == "-t") 
		{
		  if(loggingLevel >= 2) std::cout<<"loading tree string:"<<argv[i+1]<<endl;; 
		  loadTreeString(argv[i+1]); 
		  have_tree = true; 
		  if(loggingLevel >= 2) std::cout<<"Tree loaded successfully\n";
		}
	  else if (string(argv[i]) == "-n") 
		{		
		  const char* numPaths = argv[i+1]; 
		  num_sampled_paths = atoi(numPaths);
		  if (num_sampled_paths == 0)
			{
			  std::cerr<<"ERROR: You must sample at least one path at each iteration\n";
			  display_opts(); 
			  exit(0); 
			}
		}

	  else if (string(argv[i]) == "-i") 
		{		
		  const char* ins = argv[i+1]; 
		  ins_rate = atof(ins);
		  if (ins_rate < 0.0)
			{
			  std::cerr<<"ERROR: Insert rate must be greater than 0\n";
			  display_opts(); 
			  exit(0); 
			}
		}

	  else if (string(argv[i]) == "-d") 
		{		
		  const char* del = argv[i+1]; 
		  del_rate = atof(del);
		  if (del_rate < 0.0 )
			{
			  std::cerr<<"ERROR: delete rate must be greater than 0\n";
			  display_opts(); 
			  exit(0); 
			}
		}

	  else if (string(argv[i]) == "-g") 
		{		
		  const char* gap = argv[i+1]; 
		  gap_extend = atof(gap);
		  if (gap_extend < 0.0 || gap_extend >1)
			{
			  std::cerr<<"ERROR: Gap extend probability must be within 0 and 1\n";
			  display_opts(); 
			  exit(0); 
			}
		}


	  else if (string(argv[i]) == "-e") 
		{		
		  const char* envDist = argv[i+1]; 
		  envelope_distance = atoi(envDist);
		}

	  else if (string(argv[i]) == "-h") 
		{
		  display_opts(); 
		  exit(0);
		}
	  else if (string(argv[i]) == "-sa") 
		{
		  if (i==argc-1) show_alignments = true; 
		  else if (string(argv[i+1]) == "true" || string(argv[i+1]) == "TRUE" || string(argv[i+1]) == "1")
			show_alignments = true;
		  else
			show_alignments = false; 
		}
	  else if (string(argv[i]) == "-l") 
		{
		  if (i==argc-1) leaves_only = true; 
		  else if (string(argv[i+1]) == "true" || string(argv[i+1]) == "TRUE" || string(argv[i+1]) == "1")
			leaves_only = true;
		  else
			leaves_only = false; 
		}
			  
	  else if(string(argv[i]) == "-g")
		rate_matrix_filename = argv[i+1]; 
	}
  if(!have_sequences)
	{
	  error+= "ERROR: No sequence file could be imported.  Use -stk to specify a sequence file\n";
	  all_reqd_args =false; 
	}
  if (!have_tree)
	{
	  error+= "ERROR: No tree string could be imported.  Use -t  to specify a phylogenetic tree, or include it the stockholm file as a  '#=GF NH' line \n";
	  all_reqd_args =false; 
	}
  
  
  if(! all_reqd_args){std::cerr<<error<<endl;  display_opts(); exit(0); }
	  
}

void Reconstruction::display_opts(void)
{
  map<string, string>::iterator opt_iter;
  std::cout<<"\nProtPal - Progressive alignment and reconstruction of proteins\n";
  std::cout<<"\nUsage: \nprotpal [OPTIONS] > reconstruction.stk\n";
  std::cout<<"\n\n";
  std::cout<<"Command line options:\n\n";

  for (opt_iter = options.begin(); opt_iter != options.end(); opt_iter++)
	{
	  std::cout<<opt_iter->first<<"\t"<<opt_iter->second<<endl;
	}
  std::cout<<"-h\tDisplay help message\n\n";
}

void Reconstruction::loadTreeString(const char* in)
{
  if (loggingLevel>=2) std::cerr<<"Loading tree string: "<<in<<endl; 
  const sstring tree_string = in; 
  istringstream tree_input (tree_string);
  PHYLIP_tree in_tree;
  in_tree.read (tree_input);
  tree = in_tree; 
}

// vector<string> Reconstruction::get_node_names(void)
// {
//   //not yet implemented...
// }  


void Reconstruction::get_stockholm_tree(const char* fileName)
{
  string line;
  ifstream seqFile(fileName);
  string tree_tmp = ""; 
  const char* tree_string; 
  vector<string> splitLine; 
  // Parse stockholm file                                                                                                    
  if (seqFile.is_open())
    {
      while (! seqFile.eof() )
        {
		  getline(seqFile,line);
		  splitLine = split(line, " "); 
		  //std::cerr<<"Line: "<<line<<endl;
          if (splitLine.size()<2) continue;
		  else if (splitLine[1] == "NH")
            {
			  for (int i= 2; i<splitLine.size(); i++) // newick string possibly has spaces!
				tree_tmp += splitLine[i]; 

			  tree_string = tree_tmp.c_str(); 
			  loadTreeString(tree_string); 
			  have_tree = true; 
            }
		}
      seqFile.close();
    }
  else
    {
	  std::cerr << "Error: Unable to open file: "<<fileName;
    }
}

  
