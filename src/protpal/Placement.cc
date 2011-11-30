#include "protpal/Placement.h"
#include "protpal/profile.h"
#include "protpal/json.h"
#include "util/Regexp.h"
#include "util/sexpr.h"
#include "util/unixenv.h"

#define DEFAULT_CHAIN_FILE "data/handalign/hidden.hsm"
#define nullValue ""

Placement::Placement(int argc, char* argv[])
{
  default_chain_filename << Dart_Unix::get_DARTDIR() << '/' << DEFAULT_CHAIN_FILE;
  INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
                  "Place reads\n");
  opts.program_name = "MePal"; 
  opts.short_description = "Phylogenetic placement of sequence reads using ancestral profiles"; 
  opts.syntax = "-r <reads.fa> -spp <directory containing posterior profiles> -tf <TREE_FILE> [OPTIONS] > [PLACEMENT FILE]"; 
  opts.add("l -logging-level", loggingLevel=1, "Logging level"); 

  opts.newline();
  opts.print_title("Input options");
  opts.add("r -reads", readsFileName=nullValue, "Reads to place", false); 
  opts.add("spp -saved-posterior-profiles", savedPosteriorProfilesDir=nullValue, "Directory containing saved posterior profiles (made w/ ProtPal)", false); 
  opts.add("tf -tree-file", treeFileName=nullValue, "Tree file, in Newick format\n", false); 

  opts.print_title("Limiting heuristics");
  opts.add("ji -json-input", jsonInputFileName=nullValue, "Use JSON placement file to constrain read/node pairs considered", false); 
  opts.add("dc -distance-cutoff", distance_cutoff=0.0, "Consider only placements which are within this distance to a placement in the input JSON file. "); 

  
  opts.print_title("Read-profile scoring"); 
  opts.add("b -subst-model", rateMatrixFileName = default_chain_filename, "DART format chain file to be used for profile - read pairHMM.  NB: this should be the same chain used when calling ProtPal to make profiles.");
  opts.add("is -invert-score", invert = false, "Invert read-profile score", false); 
  // Possibly: parameters of pairHMM (indel rate, etc)

  opts.print_title("Output");
  opts.add("json -write-json", jsonPlacementsFileName=nullValue, "Write JSON format summary of placements (as per pplacer JSON spec)", false);
  opts.add("mp -min-postprob", minPostProb=1e-5, "Omit placements whose posterior probability falls below this threshold");

  opts.expect_args = -1; 
  opts.parse_or_die();  
  string errors; 

  SExpr_file ecfg_sexpr_file (rateMatrixFileName.c_str());
  SExpr& ecfg_sexpr = ecfg_sexpr_file.sexpr;
  ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, ecfg_sexpr);

  // Holdover from profile I/O quirk
  vector<sstring> toks = alphabet.tokens();
  for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
    alphabetVector.push_back(string(a->c_str()));

  // Some things which we absolutely need: profiles, tree, and sequences
  // Profiles are assumed to be living in a directory specified by -spp, named according to tree nodes
  if (savedPosteriorProfilesDir == nullValue)
    errors += "\tNo directory for saved profiles was specified.  Profiles are necessary for placement.  Use the -spp option to specify this\n"; 

  // Tree must be provided in a Newick file
  if (treeFileName != nullValue)
    {
      if (!FileExists( string(treeFileName)))
	THROWEXPR("\nERROR: tree file " + treeFileName  + " does not exist. Exiting...\n\n");
      parse_tree(treeFileName.c_str()); 
    }
  else
    errors += "\tNo tree was specified.  Use the -tf option to specify a tree file\n"; 
  
  // Reads must be provided in a FASTA file
  if (readsFileName != nullValue)
    {
      if (!FileExists( string(readsFileName)))
        THROWEXPR("\nERROR: sequence file " + readsFileName  + " does not exist. Exiting...\n\n");
      reads = parse_fasta(readsFileName.c_str(), alphabet); 
    }
  else
    errors += "\tNo reads were specified.  Use the -r option to specify a FASTA file of sequence reads\n"; 
  
  // If any of the three above weren't satisfied, complain and exit.  
  if (errors.size())
    {
      cerr << opts.short_help(); 
      THROWEXPR("\nNot all required arguments were supplied:\n" + errors); 
    }
  
  // If we got a JSON file on input, we can use this to limit which node/read pairs are considered
  // (e.g. only try nodes "close" to those which appear in this file
  if (jsonInputFileName != nullValue)
    {
      placementLimiter.parse_JSON(jsonInputFileName.c_str());
      using_limiter = true; 
      placementLimiter.tree = &tree; 
      if (distance_cutoff < 0.0)
	THROWEXPR("Distance cutoff ( -dc --distance-cutoff ) value must be positive!")
      placementLimiter.distance_cutoff = distance_cutoff; 
      if ( loggingLevel >= 1)
	cerr << "Limiter initialized using file: " << jsonInputFileName << endl;
    }
  
  if (distance_cutoff > 0.0 and jsonInputFileName == nullValue)
    cerr << "Warning: you specified a distance cutoff ( -dc --distance-cutoff ) but no JSON input file, making the distance cutoff meaningless.  Proceeding with exhaustive read/profile scoring...\n"; 

  // Store the arguments to MePal for the JSON string.
  for (int i =0; i< argc; i++)
    {
      args << argv[i];
      if (i != argc-1)
	args << " "; 
    }
}

void Placement::Run()
{
  // Make sure we have a profile for each node - warn the user if some are missing
  profile_filenames = check_profile_filenames(); 
  Read read;
  int readsAtNode; 

  if (loggingLevel >=1 )
    {
      cerr << "Beginning placement operations with the following input data:\n"; 
      cerr << "\t" << reads.size() << " reads read from input.\n"; 
      cerr << "\t" << profile_filenames.size() << " posterior profiles found.\n"; 
      if (using_limiter)
	cerr <<"\tUsing limiter to constrain placements.\n\t\tMax allowed distance to input placement: " << distance_cutoff << "\n"; 
    }
  

  // Move through the nodes, scoring a subset (possibly all) of the reads to the profile
  // associated with that node
  //  for (int profileNode = 0; profileNode < tree.nodes(); profileNode ++)
  //    {
  for_nodes_pre (tree, tree.root, -1, bi)
    {
      const Phylogeny::Branch& b = *bi;
      int profileNode = b.second; 
      //      for_rooted_children(tree, profileNode, child
      //	cerr << tree.branch_length(profileNode, *child) << endl;

      if (profile_filenames[profileNode] == nullValue)
	{
	  cerr << "Warning: no profile found for node: " << tree.node_name[profileNode] << endl; 
	  continue;
	}
      readsAtNode = 0; 
      if (using_limiter)
	if ( !placementLimiter.nodes_with_reads.count( tree.node_name[profileNode] ))
	  continue;

      // This is slightly clunky - and not really needed unless testing on small datasets.  
      // Probably rm it soon.
      bool some_read_possible=false;
      if (using_limiter)
	{
	  for ( map<string, string>::iterator readIter = reads.begin();
		readIter != reads.end(); ++readIter )
	    if (placementLimiter.is_allowed(readIter->first, tree.node_name[profileNode]))
	      { some_read_possible=true; break;}
	  if (! some_read_possible)
	    continue; 
	}
      // Read in the profile from a file
      AbsorbingTransducer ancestralProfile(profile_filenames[profileNode].c_str(),
					   alphabetVector, tree);
      if (ancestralProfile.num_delete_states == 0)
	{
	  cerr << "WARNING: profile " << tree.node_name[profileNode] << " has no delete states.  This makes matching reads to it meaningless, so it is being skipped\n"; 
	  continue;
	}
      
      if (loggingLevel >=1)
	cerr << "Scoring reads to profile at node: " << tree.node_name[profileNode]<<" (" << ancestralProfile.num_delete_states << " states) ... ";
      // The ancestor name is not currently written/parsed in sexpr format. 
      // Fix this later!
      ancestralProfile.name = tree.node_name[profileNode];
      ReadProfileScore profile_scorer(&ancestralProfile,
				      alphabet, rate_matrix);

      // Loop over all reads and possibly score them
      for ( map<string, string>::iterator readIter = reads.begin();
	    readIter != reads.end(); ++readIter )
	{
	  // We may want to avoid this read/node pair based on input data (e.g. prior placements). 
	  if (using_limiter)
	    {
	      if ( ! placementLimiter.is_allowed(readIter->first, tree.node_name[profileNode]) )
		{
		  if (loggingLevel >=3)
		    cerr << "\n\tSkipping read/node pair: " << readIter->first << " " << tree.node_name[profileNode] << endl; 
		  continue;
		}
	    }
	  // Update 'read'
	  read.identifier = readIter->first;
	  read.set( readIter->second );
	  ++readsAtNode; 
	  if (loggingLevel >= 2)
	    cerr<<"\n\tScoring read/profile pair: " << read.identifier << " - " << profile_scorer.name <<";\t";
	  else
	    cerr << "."; 
	  // For testing:
	  //	  profile_scorer.get_score(read, false, //no viterbi                    
	  //				   true ); // do log                         
	  //	  cerr << "Score has been stored (DP function exited normally)\n";    	  
	  //	  exit(0); 

	  if (using_limiter)
	    if ( 0 )// placementLimiter.is_single_possibility(readIter->first))
	      scores[readIter->first][profile_scorer.name] = 1.0; 
	    else
	      // Here goes the atual dynamic programming - the score is deposited in the scores map
	      profile_scorer.score_and_store(read, scores, false, invert);	  
	  else
	    // Awkward repetition
	    profile_scorer.score_and_store(read, scores, false, invert);	  
	  

	  if (loggingLevel >= 2)
	    cerr << "Score: " << scores.at(read.identifier).at(profile_scorer.name) << endl;

	  #ifdef DART_DEBUG
	  bfloat score;
	  
	  if (scores.count(read.identifier))
	    if (scores[read.identifier].count(profile_scorer.name))
	      {
		score = scores[read.identifier][profile_scorer.name];
		if (loggingLevel >=2 )
		  cerr << "Score has been retreived: " << score << endl; 
		if ( ! bfloat_is_nonzero(score))
		  {
		    cerr<<"ERROR: score is non-positive!\n";
		    profile_scorer.clear_DP_matrix();
		    profile_scorer.get_score(read, false, //no viterbi                    
					     true ); // do log                            
		    THROWEXPR("Nonpositive read/profile score");
		  }

	      }
	    else
	      THROWEXPR("Read had no score for profile " + profile_scorer.name + " (this shouldn't happen)"); 
	  #endif
	}
      if (loggingLevel >=1)
	cerr << readsAtNode << " reads scored.\n"; 
    }
  if (loggingLevel >=1 )  
    cerr <<"\n\n";
  // Display scores that were stored, either to file or stdout
  if (jsonPlacementsFileName != nullValue)
    {
      if (loggingLevel >=1 )
	cerr<<"Writing placement JSON to " << jsonPlacementsFileName << endl; 
      ofstream json_file;
      json_file.open(jsonPlacementsFileName.c_str());
      write_placement_JSON(json_file, scores);
      json_file.close(); 
    }
  else 
    write_placement_JSON(cout, scores);    

  
  if (loggingLevel >=1 )
    cerr<<"Read placement finished without errors. \n";
  exit(0);
}
// END PHYLO-PLACEMENT                        


map<int, string> Placement::check_profile_filenames(void)
{
  map<int, string> filenames;
  for (int treeNode = 0; treeNode < tree.nodes(); ++treeNode)
    {
      stringstream candidateName;
      candidateName << savedPosteriorProfilesDir << "/" << tree.node_name[treeNode] << ".sexpr";
      if ( FileExists(candidateName.str()) )
	filenames[treeNode] = candidateName.str();
      else
        {
          cerr<<"Warning: no saved profile for tree node " << tree.node_name[treeNode] << endl;
          filenames[treeNode] = nullValue;
        }
    }
  return filenames;
}
void Placement::parse_tree(const char* fileName)
{
  if (!FileExists( string(fileName)))
    {
      THROWEXPR("Could not parse tree file " + string(fileName) + " since it does not exist.\n"); 
    }
  string line;
  ifstream treeFile(fileName);
  try
    {
      if (treeFile.is_open())
        tree.read (treeFile);
    }
  catch (const Dart_exception& e)
    {
      cerr << "ERROR: input tree was not readable:\n";
      cerr << e.what();
      exit(1);
    }
}


void Placement::write_placement_JSON(ostream& out, ScoreMap& scores)
// Now that we're using a 'real' JSON parser/writer, I might do away with this business below:
// (which is an "ungainly mix of code and formatting text" if I ever saw one!!)
{
  out <<"{\n";
  out <<"\"tree\":   ";
  // Write a tree string with branches labeled (rather than nodes)                                
  write_numbered_newick(out, true);
  out <<",\n";
  out << "\"placements\": [\n";
  bfloat totalScore;
  for (ScoreMap::iterator readIter = scores.begin(); readIter != scores.end();
       ++readIter)
    {
      totalScore=0.0;
      for (map<string, bfloat>::iterator nodeIter = (readIter->second).begin();
           nodeIter != (readIter->second).end(); ++nodeIter)
        totalScore += nodeIter->second;
      if (! bfloat_is_nonzero(totalScore) )
	THROWEXPR("Sum-over-nodes score for read is zero:" + readIter->first); 
      if (bfloat_is_nan(totalScore) )
	THROWEXPR("Sum-over-nodes score for read is NaN:" + readIter->first); 
      //      out << "\t{\"p\":\n\t[";                                                            
      out << "\t{\"p\": [\n\t";
      for (map<string, bfloat>::iterator nodeIter = (readIter->second).begin();
	   nodeIter != (readIter->second).end(); ++nodeIter)
        {
          if (nodeIter->second/totalScore < minPostProb && nodeIter != --(readIter->second).end())
            continue;
          node nodeIdx = tree.find_node((nodeIter->first).c_str());
          //      if (nodeIter != readIter->second.begin())                                       
          out << "   [" << nodeIdx << ", " << nodeIter->second << ", " <<  nodeIter->second/totalScore << ", " << "0.000008" << ", " << "0.01";

          if (nodeIter != --(readIter->second).end())
            out << "],\n\t";
          else
	    out << "]],\n\t";
        }
      out << " \"n\": [\"" << readIter->first << "\"] }";
      if (readIter != --scores.end())
        out << ",\n\n";
      else
        out << "\n";
    }
  out << "\t],\n\n";

  out << "\"metadata\": {\"invocation\": \"" << args.str() << "\"},\n";
  out << "\"version\": 1,\n";
  out << "\"fields\":   [\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\", \
\"pendant_length\" ]\n";
  out << "}\n";
}


void Placement::write_numbered_newick(ostream& out, bool quotes)
{
  if (quotes)
    out <<"\"";
  tree.write(out, -1, -1 , false, // no newline                                                   
             true // add branch numberings as per pplacer JSON spec                               
             );
  if (quotes)
    out <<"\"";
}

void PlacementLimiter::parse_JSON(const char* JSON_file)
{
  Json::Value root; 
  Json::Reader reader; 
  ifstream jsonFile; 
  string readName; 
  unsigned int idx, edge_field=1e5; 
  
  // Open the file and attempt to parse it.
  jsonFile.open(JSON_file); 
  bool parsingSuccessful = reader.parse(jsonFile, root); 
  if (!parsingSuccessful)
    {
      cerr << "Error in parsing JSON file:\n" << reader.getFormattedErrorMessages(); 
      THROWEXPR("JSON parsing error"); 
    }
  
  // Get the tree and extract the node/branch names
  // This is fairly hacky - originally tried w/ Regexps, but troublesome
  const Json::Value tree = root["tree"]; 
  vector<sstring> node_vector =  sstring((tree.asString())).split("(),;"); 
  int branch; 
  sstring node; 
  vector<sstring> splitName;
  Regexp re("\\{[0-9]+\\}"); 
  for (vector<sstring>::iterator nIter=node_vector.begin(); nIter!=node_vector.end(); ++nIter)
    {
      if (! re.Match(nIter->c_str()))
	THROWEXPR("This node does not appear to have a pplacer-style branch index label like this:\n nodeName:branchLength[branchIndex]\nrendering the JSON data unusable. The string was:\n" + *nIter);
      if (index( ":", *nIter) != -1)
	{
	  // if there is no node name this cannot continue
	    splitName  = nIter->split(":"); 
	    if (splitName.size() < 2)
	      THROWEXPR("Internal nodes must have names, but this one seems not to " + *nIter + "   (This makes mapping between nodes very tricky, and this is not yet supported");
	     
	    node = nIter->split(":")[0]; 
	    branch = atoi(remove_from_string( nIter->split(":")[1].split("{")[1], "}").c_str()); 
	    branches2nodes[branch] = node; 
	    nodes2branches[node] = branch; 
	  }
	else
	  {
	    node = nIter->split("{")[0]; 
	    branch = atoi(remove_from_string( nIter->split("{")[1], "}").c_str());
	    branches2nodes[branch] = node; 
	    nodes2branches[node] = branch; 
	  }
      }
    
//     for( map<int,string>::iterator it=branches2nodes.begin(); it!=branches2nodes.end(); it++)
//       cerr << it->first << "\t" << it->second << endl; 
//     exit(0); 
    const Json::Value fields = root["fields"]; 
    for (idx = 0; idx< fields.size(); ++idx)
      if (fields[idx] == "edge_num" )
	edge_field = idx; 
    
    if (edge_field == 1e5)
      THROWEXPR("No edge field found in JSON file");
    
    const Json::Value placements = root["placements"]; 
    for (idx=0; idx<placements.size(); ++idx)
      {
	readName = remove_from_string(placements[idx]["n"][0].asString(), "\""); 
	const Json::Value probs = placements[idx]["p"]; 
	for (Json::ValueConstIterator pIter = probs.begin(); pIter!=probs.end(); ++pIter)
	  {
	    branch = int( (*pIter)[edge_field].asInt());
	    read_node_map[ readName ].push_back( branch );
	    nodes_with_reads[branches2nodes[branch]] = true; 
	  }
      }    
//     for (map<string, vector<int> >::iterator placeIter=read_node_map.begin();
// 	 placeIter!=read_node_map.end(); placeIter++)
//       {
// 	cerr << "Read " << placeIter->first << " can map to : "; 
// 	for (vector<int>::iterator nIter=placeIter->second.begin(); nIter!=placeIter->second.end(); ++nIter)
// 	  cerr << *nIter << " "; 
// 	cerr <<endl;
//       }
    cerr << nodes_with_reads.size() << " nodes have reads in limiter\n"; 
}

bool PlacementLimiter::is_single_possibility(string readName)
{
  return bool(read_node_map.at(readName).size() == 1);
}  

bool PlacementLimiter::is_allowed(string readName,  string nodeName)
{
  if (distance_cutoff == 0.0)
    {
      if ( !nodes2branches.count(nodeName))
	{
	  cerr << "Warning: tree node has no entry in branch list!\n"; 
	  return true; 
	}
      else
	return bool(index(nodes2branches[nodeName], read_node_map[readName]) != -1);
    }
  else 
    return is_within_distance(nodeName, distance_cutoff, read_node_map[readName]);
}


bool PlacementLimiter::is_within_distance(string nodeName, double distance_cutoff, 
					  vector<int> input_branches)
{
  // Is the node named nodeName within distance of one of the nodes in input_branches?
  Node node1 = tree->find_node_or_die(nodeName), node2; 
  for (vector<int>::iterator branch = input_branches.begin(); 
       branch!= input_branches.end(); ++branch)
    {
      node2 = tree->find_node_or_die(branches2nodes[*branch]);
      if (tree_distance(node1, node2) < distance_cutoff)
	return true; 
    }      
  return false; 
}

double PlacementLimiter::tree_distance(Node node1, Node node2)
{
  double totalDist=0.0; 
  Phylogeny::Node_vector path, parents = tree->find_parents(0,-1); 
  path = tree->node_path(node1, node2, parents); 
  for (unsigned int i=0; i<path.size()-1; ++i)
    {
      totalDist += tree->branch_length( path[i], path[i+1] );
    }	
  //  cerr << "Total distance between nodes: " <<  tree->node_name[node1] << " , ";
  //  cerr << tree->node_name[node2] << " " << totalDist << endl; 
  return totalDist; 
}
  
  

PlacementLimiter::PlacementLimiter(void)
{
  // Null constructor
}
