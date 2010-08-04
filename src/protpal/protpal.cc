#include<iostream>
#include<math.h>
#include<ctime>

#include "protpal/profile.h"
#include "protpal/reconstruction.h"
#include "protpal/CompositePath.h"
#include "protpal/exactMatch.h"
#include "protpal/transducer.h"
#include "protpal/Q.h"
#include "protpal/utils.h"
#include "ecfg/ecfgsexpr.h"
#include "tree/phylogeny.h"
#include "ecfg/ecfgmain.h"
#include "seq/alignment.h"
#include "seq/biosequence.h"

// Main reconstruction program. 

int main(int argc, char* argv[])
{
  // A few utility variables
  node treeNode; 
  vector<Node> children;
  vector<double> branchLengths;   
  vector<string> node_names;
  double verySmall = 0.001; //proxy for zero-length branches
  double branch_length; 

  // create main reconstruction object
  Reconstruction reconstruction(argc, argv);
  
  srand((unsigned)time(0));

  // clunky, sorry
  for (unsigned int i=0; i< reconstruction.tree.nodes(); i++)
	node_names.push_back(string(reconstruction.tree.node_name[i]));


  // Initialize the DART-style rate matrix.
  SExpr_file ecfg_sexpr_file (reconstruction.rate_matrix_filename.c_str());
  SExpr& ecfg_sexpr = ecfg_sexpr_file.sexpr;
  Irrev_EM_matrix rate_matrix(1,1);
  Alphabet alphabet ("uninitialized", 1);
  ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, ecfg_sexpr);
  
  // We need the alphabet to parse sequences via DART's machinery...that's why this is out here away from other
  // data parsing stuff. 
  reconstruction.parse_sequences(alphabet);

  
  // These  transducers remain the same throughout the traversal, so we can initialize them
  // once and for all and leave them.  
  SingletTrans R(alphabet, rate_matrix);
  SplittingTrans Upsilon;

  // If running a simulation was requested, do this instead of reconstruction
  if (reconstruction.simulate)
	{
	  reconstruction.simulate_alignment(alphabet, rate_matrix); 
	  exit(0); 
	}

  
  // Otherwise, let the ancestral reconstruction begin!  
  //  Initialize the exact-match transducers at leaf nodes - can be thought of as a trivial reconstruction
  if(reconstruction.loggingLevel>=1)
	std::cerr<<"\n";
  vector<Node> leaves = reconstruction.tree.leaf_vector(); 
  for (unsigned int i=0; i<leaves.size(); i++)
	{
	  treeNode = leaves[i]; 
	  if(reconstruction.loggingLevel>=1)
		std::cerr<<"Making exact-match transducer for: "<<reconstruction.tree.node_name[treeNode]<<endl;
	  ExactMatch leaf(
			  reconstruction.sequences[reconstruction.tree.node_name[treeNode]], // sequence
			  leaves[i], //tree index
			  alphabet // sequence alphabet
			  );
	  // Then, make an absorbing transducer from this, and place it in the profiles map.  
	  AbsorbingTransducer leafAbsorb(&leaf); 
	  reconstruction.profiles[treeNode] = leafAbsorb; 	  
 	}

  // Postorder traversal over internal nodes.  
  // For each internal node:
  //    1. Instantiate its profile using its two child profiles, delete children
  //    2. Fill DP matrix
  //    3. Sample from DP matrix
  //    4. Instantiate an absorbing profile via the sampled profile

  for_nodes_post (reconstruction.tree, reconstruction.tree.root, -1, bi)
	{
	  const Phylogeny::Branch& b = *bi;
	  treeNode = b.second; 
	  if (reconstruction.tree.is_leaf(treeNode)) continue; 
	  if (reconstruction.loggingLevel>=1)
		std::cerr<<"\nBuilding sequence profile for: "<<reconstruction.tree.node_name[treeNode]<<endl; 
	  children.clear(); 
	  branchLengths.clear(); 
	  for_rooted_children(reconstruction.tree, treeNode, child)
		{
		  children.push_back(*child); 
		  branch_length = max(verySmall, reconstruction.tree.branch_length(treeNode ,*child));
		  if (branch_length != reconstruction.tree.branch_length(treeNode ,*child))
		    {
		      std::cerr<<"NB: branch length "<< reconstruction.tree.branch_length(treeNode ,*child);
		      std::cerr<<" rounded up to "<< branch_length<<endl; 
		    }
		  branchLengths.push_back(branch_length); 
		}

	  // Instantiate the Q transducer object and its prerequisites.  
	  // The two branch transducers must be re-created at each iteration so that the branch lengths
	  // and the depending parameters are correct.  Thus, Q's transitions must be re-built.  
	  // Finally, marginalize Q's null states (e.g. IMDD)
	  
	  BranchTrans B_l(branchLengths[0], alphabet, rate_matrix, 
			    reconstruction.ins_rate, reconstruction.del_rate, reconstruction.gap_extend);
	  B_l.name ="Left branch";
	  BranchTrans B_r(branchLengths[1], alphabet, rate_matrix,
			    reconstruction.ins_rate, reconstruction.del_rate, reconstruction.gap_extend);
	  B_r.name ="Right branch";

	  QTransducer Q(R, //singlet
					B_l, //left-branch
					B_r, // right-branch
					Upsilon, //splitting 
					alphabet //sequence alphabet
					);
	  Q.marginalizeNullStates();
	  // Create the new profile at node "treeNode", using the profiles at the children of treeNode:
	  Profile profile(treeNode, // node on the tree where this profile sits
					  reconstruction.profiles[children[0]], // left absorbing transducer
					  reconstruction.profiles[children[1]], // right absorbing transducer
					  Q // Q transducer object - holds R, left branch, right branch, etc
					  );

	  // Somewhat clunky: allow the profile to see the node-name mapping that the reconstruction.tree object
	  // holds, for multiple alignment display.  
	  profile.leaves = leaves; 
	  profile.node_names = node_names; 
	  profile.envelope_distance = reconstruction.envelope_distance; 
	  profile.max_sampled_externals = reconstruction.max_sampled_externals; 

	  // Fill the Z matrix via the forward-like algorithm- the only argument is logging level
	  if(reconstruction.loggingLevel>=1)
		std::cerr<<"\tFilling forward dynamic programming matrix..."; 
	  profile.fill_DP(reconstruction.loggingLevel, reconstruction.estimate_params);

	  if(reconstruction.loggingLevel>=1)
		std::cerr<<"done.\n\t\tSubalignment likelihood: "<<-log(profile.forward_prob)/log(2)<<" bits\n"; 

	  // For non-root nodes:
	  if (treeNode != reconstruction.tree.root) 
		{
		  // Sample a traceback through the Z matrix.  Relatively quick compared to DP
		  // This step also stores the set of states, and the associated transition information (weights, connectivity)
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"\tSampling "<<reconstruction.num_sampled_paths<<" alignments from DP matrix..."; 
		  profile.sample_DP(
							reconstruction.num_sampled_paths, // number of paths
							reconstruction.loggingLevel, // debugging log messages ?
							reconstruction.show_alignments, // show sampled alignments ?
							reconstruction.leaves_only // only show leaves
							); 

		  // clear the DP matrix - this really only clears the associativity, not the actual objects therein
		  profile.clear_DP();
		  if(reconstruction.loggingLevel>=1)
			{
			  std::cerr<<"done.  ";
			  std::cerr<<"\n\tResulting DAG has "<< profile.num_sampled_externals<< " absorbing states." << endl; 
			}
		  
		  // Transform the (null-in, null-out) transducer into an absorbing transducer:
		  // Remove R-states, modify transition probabilities, sum over null states, index remaining delete states
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"\tTransforming sampled profile DAG into an absorbing transducer...";
		  AbsorbingTransducer absorbTrans(&profile);
		  absorbTrans.test_transitions(); 
		  reconstruction.profiles[treeNode] = absorbTrans; 
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"done.\n";
		}
	  
	  else	  // When reaching the root: 
	    {
	      
	      string alignString = profile.sample_DP(
						     1, // sample only the viterbi path
						     0, // debugging log messages
						     true, // show the final alignment
						     reconstruction.leaves_only // show only the leaf alignment
						     ); 
	      // There is no way that this is the easiest way to do this, but oh well:
	      stringstream treeStream; 
	      reconstruction.tree.write_Stockholm(treeStream);
	      alignString = treeStream.str() + alignString; 
	      Sequence_database db; 
	      istringstream stockStream(alignString);
	      
	      Stockholm stk(1,1);
	      stk.read_Stockholm(stockStream,db); 
	      // if requested, show what was inserted/deleted on each branch  (written to file)
	      if (reconstruction.indel_filename != "None")
		reconstruction.show_indels(stk);

	      if (reconstruction.leaves_only)
		stk.write_Stockholm(std::cout);
	      else
		{
		  ECFG_main ecfg; 
		  ecfg.ancrec_CYK_MAP = true; 
		  if (reconstruction.ancrec_postprob)
		    {
		      ecfg.ancrec_postprob=true;
		      ecfg.min_ancrec_postprob = reconstruction.min_ancrec_postprob; 
		    }
		  
		  SExpr_file grammar_sexpr_file (reconstruction.grammar_filename.c_str()); 
		  SExpr& grammar_ecfg_sexpr = grammar_sexpr_file.sexpr;
		  
		  if (reconstruction.loggingLevel >=1) 
		    std::cerr<< "\tReconstruction ancestral characters conditional on ML indel history..."; 
		  Stockholm annotated = ecfg.run_alignment_annotation(stk, grammar_ecfg_sexpr); 
		  if (reconstruction.loggingLevel >=1) 
		    {
		      std::cerr<<"Done.\n";
		      std::cerr<<"\tDisplaying Viterbi alignment\n\n"; 
		    }
		  if (reconstruction.xrate_output || reconstruction.ancrec_postprob)
		    annotated.write_Stockholm(std::cout);
		  else
		    {
		      Phonebook::iterator seq; 
		      int nameSize, maxNameLength = 0; 
		      sstring sequence; 
		      std::cout<<treeStream.str(); 
		      for (seq = annotated.row_index.begin(); seq!=annotated.row_index.end(); seq++)
			    {
			      nameSize = seq->first.size();
			      maxNameLength = max( maxNameLength, nameSize );
			    }
		      for (seq = annotated.row_index.begin(); seq!=annotated.row_index.end(); seq++)
			{
			  if ( reconstruction.tree.is_leaf(index(seq->first, reconstruction.tree.node_name) ) )
			    sequence =  annotated.get_row_as_string(seq->second); 
			  else 
			    sequence = annotated.gr_annot[seq->first]["ancrec_CYK_MAP"];
			  if (reconstruction.fasta_output)
			    std::cout<< ">" << seq->first << "\n" << sequence << endl; 
			  else
			    std::cout<< seq->first << rep(maxNameLength-seq->first.size()+4," ") << sequence << endl; 
			}
		    }
		}
	      if (reconstruction.estimate_params)
		{
		  if (reconstruction.loggingLevel >= 1)
		    {
		      std::cerr<<"\nBeginning parameter estimation";
		      std::cerr<<"\n\tFilling backward dynamic programming matrix..."; 
		    }
		  profile.fill_backward_DP(reconstruction.loggingLevel); 
		  if (reconstruction.loggingLevel >= 1)
		    std::cerr<<"Done. \n\tAssembling expected transition event counts...";

		  int numStates = profile.backward_states.size(); 
		  for (int testIdx = 0 ; testIdx<15; testIdx++)
		    {
		      M_id m = profile.backward_states[testIdx]; 
		      for (int testIdx2 =0 ; testIdx2 < profile.outgoing[profile.backward_states[testIdx].toVector()].size(); testIdx2++)
			{
			  M_id mPrime = profile.outgoing[profile.backward_states[testIdx].toVector()][testIdx2];
			  std::cerr<<"Left state: " << m.q_state<< "\n";
			  m.display(Q); 
			  std::cerr<<"Right state: "<< mPrime.q_state<< "\n";
			  mPrime.display(Q);
			  CompositePath path(m, mPrime, reconstruction.tree, .5, reconstruction.profiles, Q); 
			  std::cerr<<"Composite path created, length:" << path.path.size() << endl; 
			  path.explode(true); 
			  std::cerr<<"After expansion, path has length " << path.path.size() << endl;
			}
		    }
		  if (reconstruction.loggingLevel >= 1)
		    std::cerr<<"Done\n";
		}
	    }
	}
  return(0);
}

