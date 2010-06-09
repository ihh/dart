#include<iostream>
#include<math.h>

#include "protpal/profile.h"
#include "protpal/reconstruction.h"
#include "protpal/exactMatch.h"
#include "protpal/transducer.h"
#include "protpal/Q.h"
#include "protpal/utils.h"
#include "ecfg/ecfgsexpr.h"
#include "tree/phylogeny.h"

// Main reconstruction program. 

int main(int argc, char* argv[])
{
  // A few utility variables
  node treeNode; 
  vector<Node> children;
  vector<double> branchLengths;   
  vector<string> node_names;
  string alphabet_string; 
  double verySmall = 0.00001; //proxy for zero-length branches
  double branch_length; 

  // create main reconstruction object
  Reconstruction reconstruction;
  
  // Initialize options / arguments from cmd line.  Tree, input sequences, etc. 
  reconstruction.get_cmd_args(argc, argv);
  //get_node_names assigns reasonable names to internal nodes (e.g. all descendent leaves, joined by "_"
  //  node_names = reconstruction.get_node_names(); 
  for (int i=0; i< reconstruction.tree.node_name.size(); i++)
	node_names.push_back(string(reconstruction.tree.node_name[i]));

  // Initialize the DART-style rate matrix. 

  SExpr_file ecfg_sexpr_file (reconstruction.rate_matrix_filename.c_str());
  SExpr& ecfg_sexpr = ecfg_sexpr_file.sexpr;
  Irrev_EM_matrix rate_matrix(1,1);
  Alphabet alphabet ("uninitialized", 1);
  ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, ecfg_sexpr);
  alphabet_string = string(alphabet.nondegenerate_chars());
  
  // These  transducers remain the same throughout the traversal, so we can initialize them
  // once and for all and leave them.  
  SingletTrans R(alphabet, rate_matrix);
  SplittingTrans Upsilon;


  //  Initialize the exact-match transducers at leaf nodes
  vector<Node> leaves = reconstruction.tree.leaf_vector(); 
  for (int i=0; i<leaves.size(); i++)
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
	  if(reconstruction.tree.is_leaf(treeNode)) continue; 
	  if(reconstruction.loggingLevel>=1)
		std::cerr<<"\nBuilding sequence profile for: "<<reconstruction.tree.node_name[treeNode]<<endl; 
	  children.clear(); 
	  branchLengths.clear(); 
	  for_rooted_children(reconstruction.tree, treeNode, child)
		{
		  children.push_back(*child); 
		  branch_length = max(verySmall, reconstruction.tree.branch_length(treeNode ,*child));
		  if (branch_length != reconstruction.tree.branch_length(treeNode ,*child))
		    {
		      std::cerr<<"Warning: branch length "<< reconstruction.tree.branch_length(treeNode ,*child);
		      std::cerr<<" rounded up to "<< branch_length<<endl; 
		    }
		  branchLengths.push_back(branch_length); 
		}

	  // Instantiate the Q transducer object and its prerequisites.  
	  // Thetwo branch transducers must be re-created at each iteration so that the branch lengths
	  // and the resulting parameters are correct.  Thus, Q's transitions must be re-built.  
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

	  reconstruction.profiles.erase( children[0] );
	  reconstruction.profiles.erase( children[1] ); 

	  // Fill the Z matrix via the forward-like algorithm- the only argument is logging
	  if(reconstruction.loggingLevel>=1)
		std::cerr<<"\tFilling dynamic programming matrix..."; 
	  profile.fill_DP(reconstruction.loggingLevel);
	  if(reconstruction.loggingLevel>=1)
		std::cerr<<"done. Sum-over-alignments likelihood: "<<-log(profile.forward_prob)/log(2)<<" bits\n"; 

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
							reconstruction.leaves_only
							); 
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"done.\n";
		  
		  // Transform the (null-in, null-out) transducer into an absorbing transducer:
		  // Remove R-states, modify transition probabilities, sum over null states, index remaining delete states
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"\tTransforming sampled profile into an absorbing transducer...";
		  AbsorbingTransducer absorbTrans(&profile);
		  absorbTrans.test_transitions(); 
		  reconstruction.profiles[treeNode] = absorbTrans; 
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"done.\n";
		}
	  
	  // When reaching the root: 
	  else
		{
		  if(reconstruction.loggingLevel>=1)
			std::cerr<<"\tDisplaying Viterbi alignment\n\n"; 
		  std::cout<<"#=GF NH\t";
		  reconstruction.tree.write(std::cout, 0); 
		  profile.sample_DP(
							1, // sample only the viterbi path
							0, // debugging log messages
							true, // show the final alignment
							reconstruction.leaves_only // show only the leaf alignment
							); 
		}
	}
  return(0);
}

