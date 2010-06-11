#include<iostream>
#include<string>
#include<fstream>
#include<sstream>

#include "protpal/utils.h"
#include "protpal/profile.h"
#include "protpal/reconstruction.h"
#include "protpal/transducer.h"
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
  
  //  options["-s"] = "<bool> Rather than aligning, create an sexpr file to pipe to phylcomposer.  Use protpal -s -t <tree> | phylocomposer --simulate | phylcomposter -al <alignment_file> to create a simulated alignment for the tree of interest\n";
  options["-s"] = "<bool> Rather than aligning, simulate a set of (unaligned) sequences according to the models (e.g. transducers, rate matrix, tree) specified. \n";
  simulate = false; 

  options["-i"] = "<float> Insert rate (default .02) \n";
  ins_rate = .019; 
  options["-d"] = "<float> Delete rate (default .02)\n";    
  del_rate = .02; 
  options["-g"] = "<float> Gap-extend probability (default .4) \n";
  gap_extend = .912; 

  options["-r"] = "<int> Root sequence length in simulation.  Default is to sample direclty from singlet transducer.  \n";
  rootLength = -1; 

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
	  else if (string(argv[i]) == "-r") 
		{		
		  const char* rL = argv[i+1]; 
		  rootLength = atoi(rL);
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
	  else if (string(argv[i]) == "-s") 
		{
		  if (i==argc-1) simulate = true; 
		  else if (string(argv[i+1]) == "true" || string(argv[i+1]) == "TRUE" || string(argv[i+1]) == "1")
			simulate = true;
		  else
			simulate = false; 
		}
			  
	  else if(string(argv[i]) == "-g")
		rate_matrix_filename = argv[i+1]; 
	}
  if(!have_sequences && !simulate)
	{
	  error+= "ERROR: No sequence file could be imported.  Use -stk to specify a sequence file\n";
	  all_reqd_args =false; 
	}
  if (!have_tree)
	{
	  error+= "ERROR: No tree string could be imported.  Use -t  to specify a phylogenetic tree, or include it the stockholm file as a  '#=GF NH' line \n";
	  all_reqd_args =false; 
	}
  
  
  if(! all_reqd_args){std::cerr<<"Error: not all required arguments were supplied:\n"<<error<<endl;  display_opts(); exit(0); }
	  
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

  
void Reconstruction::make_sexpr_file(Alphabet alphabet, Irrev_EM_matrix rate_matrix)
{
  map<string,string> prot2pc;
  prot2pc["S"] = "start";
  prot2pc["M"] = "match"; 
  prot2pc["I"] = "insert"; 
  prot2pc["D"] = "delete";
  prot2pc["W"] = "wait";
  prot2pc["E"] = "end";   
  
  std::cout << ";; phylocomposer file for simulating branch-length dependent transducers\n";
  std::cout << "(token ("  ; 
  string tokens = string(alphabet.nondegenerate_chars());
  unsigned int alph_size = tokens.size();
  string::iterator tok, tok2; 
  unsigned int tokIdx, tokIdx2; 
  vector<state>::iterator state1, state2; 
  vector<state> outgoing; 
  string parent, child; 

  for (tok = tokens.begin(); tok != tokens.end(); tok++)
	std::cout << " "<< *tok << " "; 
  std::cout<< ") )\n"; 

  // define pi, the shared equilibrium distribution over tokens
  vector<double> equilibrium = rate_matrix.create_prior();
  for (tokIdx = 0 ; tokIdx != alph_size; tokIdx++)
	std::cout << "(value ((pi " << tokens[tokIdx] << ") " << equilibrium[tokIdx] << "))\n";
  
  
  // Define parameters for, then declare singlet transducer
  SingletTrans R(alphabet, rate_matrix);
  
  // Define parameters in 'value' blocks
  std::cout<<"\n";
  for (state1 = R.states.begin(); state1 != R.states.end(); state1++)
	{
	  if (R.get_state_type(*state1) == "E") continue; 
	  outgoing = R.get_outgoing(*state1); 
	  for (state2 = outgoing.begin(); state2 != outgoing.end(); state2++)
		{
		  std::cout<< "(value (root_"<<R.get_state_name(*state1)<<"_"<<R.get_state_name(*state2)<<" ";
		  std::cout<< R.get_transition_weight(*state1, *state2)<<"))\n";
		}
	} 

  std::cout<<"\n(transducer \n\n\t(name ROOT)\n\n"; 
  
  for (state1 = R.states.begin(); state1 != R.states.end(); state1++)
	{
	  std::cout<<"\t(state (name "<< R.get_state_name(*state1)<< ") (type "<< prot2pc[R.get_state_type(*state1)] << ")";
	  if (R.get_state_type(*state1) == "I")
		std::cout<<" (label pi)"; 
	  std::cout<<")\n";
	}

  // Define transitions, referring back to earlier specified values
  std::cout<<"\n";
  for (state1 = R.states.begin(); state1 != R.states.end(); state1++)
	{
	  if (R.get_state_type(*state1) == "E") continue; 
	  outgoing = R.get_outgoing(*state1); 
	  for (state2 = outgoing.begin(); state2 != outgoing.end(); state2++)
		{
		  std::cout<<"\t(transition (from "<< R.get_state_name(*state1) << ") ";
		  std::cout<<"(to "<< R.get_state_name(*state2) << ") ";
		  std::cout<< "(label root_"<<R.get_state_name(*state1)<<"_"<<R.get_state_name(*state2)<<" ";
		  std::cout<<"))\n";
		}
	} 
  std::cout<<");; end transducer ROOT\n\n";



  for_nodes_post (tree, tree.root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      node treeNode = b.second;
	  if (treeNode != tree.root)
		{
		  parent = string( tree.node_name[b.first].c_str() );
		  child = string( tree.node_name[b.second].c_str() ); 	  

		  // construct a branch transducer with the appropriate branch length.
		  BranchTrans branch(tree.branch_length(b.first, b.second), alphabet, rate_matrix, ins_rate, del_rate, gap_extend); 
		  
		  // define the Q_child match matrix
		  for (state1 = branch.states.begin(); state1 != branch.states.end(); state1++)
			{
			  if ( branch.get_state_type(*state1) != "M" ) 
				continue;
			  for (tokIdx = 0 ; tokIdx != alph_size; tokIdx++)
				{
				  for (tokIdx2 = 0 ; tokIdx2 != alph_size; tokIdx2++)
					{
					  std::cout<< "(value ((Q_"<< child<< " " << tokens[tokIdx] << " "<< tokens[tokIdx2] << ") "; 
					  std::cout<<  branch.get_match_weight(*state1, tokIdx, tokIdx2);
					  std::cout<< "))\n";
					}
				}
			}
		  

		  // Define transition parameters in 'value' blocks
		  std::cout<<"\n";
		  for (state1 = branch.states.begin(); state1 != branch.states.end(); state1++)
			{
			  if (branch.get_state_type(*state1) == "E") continue; 
			  outgoing = branch.get_outgoing(*state1); 
			  for (state2 = outgoing.begin(); state2 != outgoing.end(); state2++)
				{
				  std::cout<< "(value ("<< child <<"_"<<branch.get_state_name(*state1)<<"_"<<branch.get_state_name(*state2)<<" ";
				  std::cout<< branch.get_transition_weight(*state1, *state2)<<"))\n";
				}
			} 


		  // begin transducer declaration
		  std::cout<<"\n(transducer \n\n\t(name "<<child<< ")\n\n"; 

  		  //  declare states, possibly with labels for insert and match states
		  for (state1 = branch.states.begin(); state1 != branch.states.end(); state1++)
			{
			  std::cout<<"\t(state (name "<< branch.get_state_name(*state1)<< ") (type "<< prot2pc[branch.get_state_type(*state1)] << ")";
			  if (branch.get_state_type(*state1) == "I")
				std::cout<<" (label pi)"; 
			  else if  (branch.get_state_type(*state1) == "M")
				std::cout<<" (label Q_"<<child<<")"; 
					 
			  std::cout<<")\n";
			}

		  // Define transitions, referring back to earlier specified values
		  std::cout<<"\n";
		  for (state1 = branch.states.begin(); state1 != branch.states.end(); state1++)
			{
			  if (branch.get_state_type(*state1) == "E") continue; 
			  outgoing = branch.get_outgoing(*state1); 
			  for (state2 = outgoing.begin(); state2 != outgoing.end(); state2++)
				{
				  std::cout<<"\t(transition (from "<< branch.get_state_name(*state1) << ") ";
				  std::cout<<"(to "<< branch.get_state_name(*state2) << ") ";
				  std::cout<< "(label " << child << "_"<<branch.get_state_name(*state1)<<"_"<<branch.get_state_name(*state2)<<" ";
				  std::cout<<"))\n";
				}
			} 
		  std::cout<<");; end transducer " << child << "\n\n";
		  
		}
	}
  vector<node> rootsKids; 
  for_rooted_children(tree, tree.root , child)
	rootsKids.push_back(*child);

  std::cout<<"(branch (name ROOT)\n\t(from SUBROOT) (to root)\n\t(transducer ROOT) \n\t";

  std::cout<< show_branch(rootsKids[0]) ;
  std::cout<< show_branch(rootsKids[1]);
  std::cout<< ");; end phylogenetic tree"; 
}

string Reconstruction::show_branch(node startNode)
{
  string out; 
  vector<node> kids; 
  out += "(branch (name " + string(tree.node_name[startNode].c_str()) + ")\n";

  out += "\t(from " + string(tree.node_name[tree.parent[startNode]].c_str()) + ") ";
  out += "(to " + string(tree.node_name[startNode].c_str()) + ") \n";  

  out += "\t(transducer " + string(tree.node_name[startNode].c_str()) + ")\n";

  if (tree.is_internal(startNode))
	{
	    for_rooted_children(tree, startNode , child)
		  kids.push_back(*child);
		out += show_branch(kids[0]);
		out += show_branch(kids[1]); 
	}		
  out += "\t);; end branch " + string(tree.node_name[startNode].c_str()) +  " \n\n";
  return out; 
}

void Reconstruction::simulate_alignment(Alphabet alphabet, Irrev_EM_matrix rate_matrix)
{
  map<node, string> sequences; 
  string parentSeq, childName;

  SingletTrans R(alphabet, rate_matrix);
  //write the tree in stockolm style
  std::cout<<"#=GF NH\t";
  tree.write(std::cout, 0); 

  for_nodes_pre (tree, tree.root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      node treeNode = b.second;
	  //std::cout<<"visiting node:" << treeNode<<endl; 
	  if (treeNode == tree.root)
		{
		  sequences[treeNode] = sample_root(R); 
		  std::cout<< childName << tree.node_name[treeNode] << "   " << sequences[treeNode] <<endl; 
			}
	  else
		{
		  // construct a branch transducer with the appropriate branch length.
		  BranchTrans branch(tree.branch_length(b.first, b.second), alphabet, rate_matrix, ins_rate, del_rate, gap_extend); 

		  parentSeq = sequences[b.first]; 
		  childName = string( tree.node_name[b.second].c_str() ); 	  
		  
		  sequences[treeNode] = sample_pairwise(parentSeq, branch);
		  std::cout<< childName << "    " << sequences[treeNode] << endl; 
		}
	}
}

string Reconstruction::sample_root(SingletTrans R)
{
  state s = 0;  // 0 is the start state
  string sType;
  vector<state>::iterator sIter; 
  string childSeq; 
  vector<state> outgoing; 
  vector<double> outgoingWeights; 
  int sampled; 

  if (rootLength != -1)
	{
	  for (s=0; s<R.states.size(); s++)
		if (R.get_state_type(s) == "I")
		  outgoingWeights = R.get_emission_distribution(s);

	  for (int sampled=0; sampled< rootLength; sampled++)
		childSeq += R.alphabet[ sample(outgoingWeights) ];
	  return childSeq; 
	}
		  

  while (s != R.states.size()-1) // not the end state
	{
	  outgoing = R.get_outgoing(s); 
	  outgoingWeights.clear();
	  for (sIter = outgoing.begin(); sIter != outgoing.end(); sIter++)
		outgoingWeights.push_back( R.get_transition_weight(s, *sIter) );
	  //The new state is sampled
	  sampled = sample(outgoingWeights)	  ;
	  s = outgoing[sampled]; 
	  sType = R.get_state_type(s); 
	  
	  if  (sType == "I")
		{
		  outgoingWeights = R.get_emission_distribution(s);
		  childSeq += R.alphabet[ sample(outgoingWeights) ];
		}
	}
  return childSeq;
}
  
string Reconstruction::sample_pairwise(string parentSeq, BranchTrans branch)
{
  state s = 0, sNew, endState = branch.states.size()-1;  // 0 is the start state
  string sType;
  vector<state>::iterator sIter; 
  string childSeq; 
  vector<state> outgoing, possible; 
  vector<double> outgoingWeights; 
  int inIdx = 0, incomingCharacter; 

  while (s != endState) // not the end state
	{
	  outgoing.clear();
	  possible = branch.get_outgoing(s); 
	  for (sIter=possible.begin(); sIter!=possible.end(); sIter++)
		{
		  if (inIdx != parentSeq.size())
			{
			  //we can only transition to end state if we've used up all the input characters
			  if (branch.get_state_type(*sIter) != "E") 
				outgoing.push_back(*sIter);
			}
		  else
			//Similarly, after the input characters have been exhausted, we can't go to match or delete
			if (branch.get_state_type(*sIter) != "M" && (branch.get_state_type(*sIter) != "D") ) 
			  outgoing.push_back(*sIter);
		}
	  
	  outgoingWeights.clear(); 
	  for (sIter = outgoing.begin(); sIter != outgoing.end(); sIter++)
		outgoingWeights.push_back( branch.get_transition_weight(s, *sIter) );

	  if (inIdx == parentSeq.size() && branch.get_state_type(s) == "W")
		if ( !in( endState, outgoing) )
		  {
			outgoing.push_back(endState); 
			outgoingWeights.push_back(1); 
		  }
	  
	  
	  //The new state is sampled
	  sNew = outgoing[sample(outgoingWeights)]; 
	  if (branch.get_state_type(s) == "W" && branch.get_state_type(sNew) == "W")
		{
		  std::cerr<<"Wait-wait transition, problem!\n Outgoing size, states, possible: "<<outgoing.size() << " , ";
		  std::cerr<<"outgoing states"; displayVector(outgoing);
		  std::cerr<<"\n";
		  std::cerr<<"outgoingWeights: "; displayVector(outgoingWeights);
		  std::cerr<<"\n";
		  std::cerr<<"possible states"; displayVector(possible);
		  exit(1); 
		}


	  s = sNew; 
	  sType = branch.get_state_type(s); 
	  if (loggingLevel >=2)
		std::cerr<<"Pairwise: new state: "<<s<<"  type: "<<branch.get_state_type(s)<<endl; 
	  
	  if  (sType == "I")
		{
		  outgoingWeights = branch.get_emission_distribution(s);
		  childSeq += branch.alphabet[ sample(outgoingWeights)  ];
		}

	  else if (sType == "M") 
		{
		  incomingCharacter = index(stringAt(parentSeq, inIdx), branch.alphabet); 
		  outgoingWeights.clear(); 
		  for (int charIdx=0; charIdx < branch.alphabet_size; charIdx++)
			outgoingWeights.push_back(branch.get_match_weight(s, incomingCharacter, charIdx) );
		  childSeq += branch.alphabet[sample(outgoingWeights)];
		  inIdx++; 
		}
	  else if ( sType == "D")
		inIdx++;
	}
  return childSeq;
}

