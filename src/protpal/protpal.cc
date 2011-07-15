#include<iostream>
#include<fstream>
#include<stack>
#include<math.h>
#include<ctime>
#include<time.h>

#include "protpal/profile.h"
#include "protpal/reconstruction.h"
#include "protpal/exactMatch.h"
#include "protpal/transducer.h"
#include "protpal/Q.h"
#include "protpal/utils.h"
#include "protpal/ReadProfileScore.h"

#include "protpal/AlignmentSampler.h"
#include "protpal/MyMap.h"
// #include "protpal/CompositePath.h"


#include "ecfg/ecfgsexpr.h"
#include "tree/phylogeny.h"
#include "ecfg/ecfgmain.h"
#include "seq/alignment.h"
#include "seq/biosequence.h"
#include "util/dexception.h"
#include "util/rnd.h"
#include "util/sexpr.h"

// Main reconstruction program. 

int main(int argc, char* argv[])
{
  time_t start,end;
  time (&start);
  // A few utility variables
  node treeNode; 
  vector<Node> children;
  vector<double> branchLengths;   
  vector<string> node_names;
  double verySmall = 0.001; //proxy for zero-length branches
  double branch_length; 
  string alignString; 
  // create main reconstruction object
  Reconstruction reconstruction(argc, argv);
  
  // srand((unsigned)time(0));

  // clunky, sorry
  for (int i=0; i< reconstruction.tree.nodes(); i++)
	node_names.push_back(string(reconstruction.tree.node_name[i]));

  // Initialize the DART-style rate matrix.
  SExpr_file ecfg_sexpr_file (reconstruction.rate_matrix_filename.c_str());
  SExpr& ecfg_sexpr = ecfg_sexpr_file.sexpr;
  Irrev_EM_matrix rate_matrix(1,1);
  Alphabet alphabet ("uninitialized", 1);
  ECFG_builder::init_chain_and_alphabet (alphabet, rate_matrix, ecfg_sexpr);
  // yeccch
  vector<string> alphabetVector;
  vector<sstring> toks = alphabet.tokens(); 
  for (vector<sstring>::iterator a=toks.begin(); a!=toks.end(); a++)
    alphabetVector.push_back(string(a->c_str()));

  // We need the alphabet to parse sequences via DART's machinery...that's why this is out here away from other
  // data parsing stuff. 
  reconstruction.parse_sequences(alphabet);
  // Same with this.
  if (reconstruction.estimate_root_insert)
    {
      reconstruction.root_insert_prob = reconstruction.get_root_ins_estimate();
      if (reconstruction.loggingLevel >=1 )
	std::cerr<< "Root insert probability estimated as: " << reconstruction.root_insert_prob << endl; 
    }
  
  // These  transducers remain the same throughout the traversal, so we can initialize them
  // once and for all and leave them.  
  SingletTrans R(alphabet, rate_matrix, reconstruction.root_insert_prob);
  SplittingTrans Upsilon;

  // If running a simulation was requested, do this instead of reconstruction
  if (reconstruction.simulate)
	{
	  reconstruction.simulate_alignment(alphabet, rate_matrix); 
	  exit(0); 
	}
  // If generating a phylocomposer input file was requested, do this instead of reconstruction
  if ( reconstruction.phylocomposer_filename != "None" )
    {
      ofstream phylocomp_out;
      phylocomp_out.open( reconstruction.phylocomposer_filename.c_str() );
      if (reconstruction.loggingLevel >= 1)
	std::cerr<<"Generating phylocomposer file..."; 
      reconstruction.make_sexpr_file(alphabet, rate_matrix, phylocomp_out);
      if (reconstruction.loggingLevel >= 1)
	std::cerr<<"done.\n"; 
      exit(0); 
    }
  
  // If using an input alignment was requested, do this and bypass reconstruction
  // This is a bit ambiguous now - an 'input alignment' is used for indel counting and other statistics, whereas a 
  // "guide alignment" is used to constrain DP in aligning sequences.  
  else if (reconstruction.input_alignment)
    {
      int seqLength = -1;
      alignString = "";
      for (treeNode = 0; treeNode < reconstruction.tree.nodes(); treeNode++)
	if (reconstruction.sequences.count(reconstruction.tree.node_name[treeNode]))
	  {
	    alignString += reconstruction.tree.node_name[treeNode] + "   " + reconstruction.sequences[reconstruction.tree.node_name[treeNode]] + "\n";
	    if (seqLength != -1)
	      {
		if ((int) reconstruction.sequences[reconstruction.tree.node_name[treeNode]].size() != seqLength)
		  {
		    std::cerr<<"Error: input alignment's sequences are not all of the same length.  Do not use -ia with unaligned sequences\n"; 
		    exit(1); 
		  }
	      }
	    else 
	      {

		seqLength = reconstruction.sequences[reconstruction.tree.node_name[treeNode]].size();
	      }
	  }
    }
  else
    {
      // Otherwise, let the "progressive transducer-profile-based" ancestral reconstruction begin!  
      //  Initialize the exact-match transducers at leaf nodes - can be thought of as a trivial reconstruction
      if(reconstruction.loggingLevel>=1)
	std::cerr<<"\n";
      vector<Node> leaves = reconstruction.tree.leaf_vector(); 
      if(reconstruction.loggingLevel>=1)
	std::cerr<<"Making exact-match transducers for leaf nodes...";

      for (unsigned int i=0; i<leaves.size(); i++)
	{
	  treeNode = leaves[i]; 
	  string sequence = reconstruction.sequences[reconstruction.tree.node_name[treeNode]]; // sequence
	  if (sequence.size() == 0 )
	    {
	      std::cerr<<"\n\tERROR: the leaf node " << reconstruction.tree.node_name[treeNode] << " has no sequence associated with it!\n";
	      std::cerr<<"This is usually caused by a tree having slightly different names than the sequences in the input sequences. \n";

	      std::cerr<<"Sequence names in input file: \n";
	      for (MyMap<string, string>::iterator seqIter = reconstruction.sequences.begin(); seqIter != reconstruction.sequences.end(); seqIter++)
		if ((seqIter->second).size() > 0)
		  std::cerr<<"\t" << seqIter->first << "\t" << split(seqIter->first, string("|"))[0] << endl; 
		    
	      std::cerr<<"Sequence names in tree leaves: \n";
	      for (unsigned int j=0; j<leaves.size(); j++)
		std::cerr<<"\t" << reconstruction.tree.node_name[leaves[j]] << endl; 
	      exit(1); 
	    }
	  ExactMatch leaf(
			  reconstruction.sequences[reconstruction.tree.node_name[treeNode]], // sequence
			  leaves[i], //tree index
			  alphabet // sequence alphabet
			  );
	  // Then, make an absorbing transducer from this, and place it in the profiles map.  
	  AbsorbingTransducer leafAbsorb(&leaf); 
	  reconstruction.profiles[treeNode] = leafAbsorb; 	  
 	}
      if(reconstruction.loggingLevel>=1)
	std::cerr<<"Done\n"; 
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
	      branch_length = min(reconstruction.max_branch_length, max(verySmall, reconstruction.tree.branch_length(treeNode ,*child)));
	      if (branch_length != reconstruction.tree.branch_length(treeNode ,*child))
		{
		  std::cerr<<"\tNB: branch length "<< reconstruction.tree.branch_length(treeNode ,*child);
		  std::cerr<<" rounded to "<< branch_length<< endl; 
		}
	      branchLengths.push_back(branch_length); 
	    }

	  // Instantiate the Q transducer object and its prerequisites.  
	  // The two branch transducers must be re-created at each iteration so that the branch lengths
	  // and the depending parameters are correct.  Thus, Q's transitions must be re-built.  
	  // Finally, marginalize Q's null states (e.g. IMDD)
	  

	  BranchTrans B_l(branchLengths[0], alphabet, rate_matrix, 
			  reconstruction.ins_rate, reconstruction.del_rate, reconstruction.gap_extend, reconstruction.sub_rate);
	  BranchTrans B_r(branchLengths[1], alphabet, rate_matrix,
			  reconstruction.ins_rate, reconstruction.del_rate, reconstruction.gap_extend, reconstruction.sub_rate);
	  // Boy this is awkward.  Can't seem to declare things within an if-else, so I'll declare it outside, then possibly overwrite it
	  // hmm...
	  if ( reconstruction.mixture_2 )
	    {
	      BranchTrans L_mix_2(branchLengths[0], alphabet, rate_matrix, 
			      reconstruction.ins_rate, reconstruction.del_rate, 
				  reconstruction.gap_extend, 
				  1-reconstruction.mix_prior_2, 
				  reconstruction.gap_extend_2,
				  reconstruction.mix_prior_2);

	      BranchTrans R_mix_2(branchLengths[1], alphabet, rate_matrix,
			      reconstruction.ins_rate, reconstruction.del_rate, 
				  reconstruction.gap_extend, 
				  1-reconstruction.mix_prior_2, 
				  reconstruction.gap_extend_2,
				  reconstruction.mix_prior_2);
	      B_l = L_mix_2; 
	      B_r = R_mix_2; 
	    }
	  else if ( reconstruction.mixture_3 )
	    {
	      BranchTrans L_mix_3(branchLengths[0], alphabet, rate_matrix, 
				  reconstruction.ins_rate, reconstruction.del_rate, 
				  reconstruction.gap_extend, 
				  1-reconstruction.mix_prior_2 - reconstruction.mix_prior_3, 
				  reconstruction.gap_extend_2,
				  reconstruction.mix_prior_2, 
				  reconstruction.gap_extend_3,
				  reconstruction.mix_prior_3);

	      BranchTrans R_mix_3(branchLengths[1], alphabet, rate_matrix,
			      reconstruction.ins_rate, reconstruction.del_rate, 
				  reconstruction.gap_extend, 
				  1-reconstruction.mix_prior_2 - reconstruction.mix_prior_3, 
				  reconstruction.gap_extend_2,
				  reconstruction.mix_prior_2, 
				  reconstruction.gap_extend_3,
				  reconstruction.mix_prior_3);
	      B_l = L_mix_3; 
	      B_r = R_mix_3; 
	    }
	  
	  B_l.name ="Left branch";
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
	  reconstruction.profiles.erase(children[0]); 
	  reconstruction.profiles.erase(children[1]); 
	  if (reconstruction.loggingLevel>=1)
	    {
	      std::cerr<<"\tLeft profile has: "<<profile.left_profile.num_delete_states<< "  absorbing states.\n"; 
	      std::cerr<<"\tRight profile has: "<<profile.right_profile.num_delete_states<< "  absorbing states.\n"; 
	    }

	  // Somewhat clunky: allow the profile to see the node-name mapping that the reconstruction.tree object
	  // holds, for multiple alignment display.  
	  profile.leaves = leaves; 
	  profile.node_names = node_names; 
	  profile.envelope_distance = reconstruction.envelope_distance; 
	  profile.max_sampled_externals = reconstruction.max_sampled_externals; 
	  // If we got a guide alignment on input, let the profile know about it so it can constrain DP breadth
	  if (reconstruction.guide_alignment_filename != "None")
	    {
	      profile.use_guide_alignment = true; 
	      profile.envelope = &reconstruction.envelope; 
	    }
	  else
	    profile.use_guide_alignment = false; 
	  // Fill the Z matrix via the forward-like algorithm- the only argument is logging level
	  if(reconstruction.loggingLevel>=1)
	    std::cerr<<"\tFilling forward dynamic programming matrix..."; 
	  if (treeNode != reconstruction.tree.root || reconstruction.root_profile_filename == "None")
	    profile.fill_DP(reconstruction.loggingLevel, false); // do not store incoming/outgoing information.  
	  else if (reconstruction.root_profile_filename != "None")
	    {
	      if(reconstruction.loggingLevel>=1)
		std::cerr<<" (logging state connectivity) "; 
	      profile.fill_DP(reconstruction.loggingLevel, true); // do not store incoming/outgoing information.  
	      if(reconstruction.loggingLevel>=1)
		std::cerr<<"\n\tFilling backward DP matrix... "; 	      
	      profile.fill_backward_DP(reconstruction.loggingLevel); 
	    }
	  if(reconstruction.loggingLevel>=1)
	    {
	      std::cerr<<"done.\n\t\tSubalignment likelihood: "<<-log(profile.forward_prob)/log(2)<<" bits\n"; 
	      std::cerr<<"\tProfile's DP matrix had " << profile.DP_size() << " cells\n"; 
	      std::cerr<<"\tEncountered " << profile.num_zero_states << " zero-likelihood cells\n"; 
	      if (profile.num_discarded_states > 0)
		std::cerr<<"\tEnvelope allowed for discarding " << profile.num_discarded_states << " candidate states during DP recursion\n"; 
	    }
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
				reconstruction.leaves_only, // only show leaves
				reconstruction.viterbi // sample viterbi path (on the last time through - only really applies to null states)
				); 
	      //	      profile.show_DOT(cout, reconstruction.tree.node_name[treeNode]); 


	      // clear the DP matrix - this really only clears the associativity, not the actual objects therein
	      // UPDATE - with some scope trickery, this should *actually* clear the DP entries
	      // profile.clear_DP();

	      //	      if (reconstruction.estimate_params)
	      //		reconstruction.pre_summed_profiles[treeNode] = profile; 
	      if(reconstruction.loggingLevel>=1)
		{
		  std::cerr<<"done.  \n";
		  std::cerr<<"\tResulting DAG has "<< profile.num_sampled_externals<< " absorbing states." << endl; 

		}
		  
	      // Transform the (null-in, null-out) transducer into an absorbing transducer:
	      // Remove R-states, modify transition probabilities, sum over null states, index remaining delete states
	      if(reconstruction.loggingLevel>=1)
		std::cerr<<"\tTransforming sampled profile DAG into an absorbing transducer...";
	      AbsorbingTransducer absorbTrans(&profile);
	      absorbTrans.test_transitions(); 
	      reconstruction.profiles[treeNode] = absorbTrans; 
	      // Dot code is a good way to investigate a strangely-behaving transducer
	      // (e.g. if the above transitions test fails)
	      //absorbTrans.show_DOT(cout, reconstruction.tree.node_name[treeNode]+"_absorbing"); 


	      // TESTING CAPABILITY  TO READ/WRITE ABSORBING TRANSDUCERS
// 	      ofstream saved_profile;
// 	      saved_profile.open( "saved_profile.sexpr" );
// 	      absorbTrans.write_profile(saved_profile); 
// 	      cerr<<"Wrote profile to file: saved_profile.sexpr\n";
// 	      saved_profile.close();
// 	      AbsorbingTransducer absorbTrans2("saved_profile.sexpr",
// 					       alphabetVector, reconstruction.tree);
	      
// 	      CHECK EQUALITY IN READ/WRITTEN TRANSDUCER
// 		absorbTrans2.test_equality(absorbTrans,true, true); 
// 	      DONE READ/WRITE TESTING


	      // TESTING / BENCHMARKING READ-PROFILE SCORING (E.G PPLACER-LIKE FUNCTIONALITY)

	      time_t readStart,readEnd;
	      time (&readStart);

	      
	      // Writing to hmmoc file - not yet fully functional, but provides a 
	      // significant speedup (~10X)
	      //testScore.HMMoC_adapter("testHMMoC.xml"); 
	      //exit(0); 


	      unsigned int numReads = 1000;
	      for (unsigned int i=0; i<numReads; i++)
		{
		  ReadProfileScore testScore(&absorbTrans, alphabet, rate_matrix); 
		  Read read; 
		  read.set("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); 
		  read.identifier="oscarsRead"; 
		  cout<<"Analyzing read number: " << i << endl; 
		  testScore.readSize=140;
		  testScore.fill_DP_matrix(read, cout, false, false, true);
		  //testScore.score_and_print(read, cout, false); 
		}
	      time(&readEnd); 
	      std::cout<< "Minutes used to map "<<  numReads << " reads to a profile: " << difftime (readEnd, readStart)/60.0 <<endl; 
	      exit(0); 
	      // END TESTING / BENCHMARKING READ-PROFILE SCORING

	      if(reconstruction.loggingLevel>=1)
		std::cerr<<"done.\n";
	    }
	  else
	    {
	      if (reconstruction.root_profile_filename != "None")
		{
		  if(reconstruction.loggingLevel>=1)
		    std::cerr<<"\tSampling "<<reconstruction.num_root_alignments<<" alignments from root alignment..."; 
		  profile.sample_DP(
				    reconstruction.num_root_alignments, // number of paths
				    reconstruction.loggingLevel, // debugging log messages ?
				    reconstruction.show_alignments, // show sampled alignments ?
				    reconstruction.leaves_only, // only show leaves
				    reconstruction.viterbi // sample viterbi path (on the last time through - only really applies to null states)
				    ); 
		  if(reconstruction.loggingLevel>=1)
		    std::cerr<<"\tTransforming sampled root profile DAG into an absorbing transducer...";
		  state_path viterbi_path;
		  if (reconstruction.root_viterbi_path)
		    viterbi_path = profile.sample_DP(
						       1, // sample only one path
						       0, // debugging log messages
						       false, // don't show alignments
						       false, // leaves only
						       true // sample the viterbi path
						       );
		  AbsorbingTransducer absorbTrans(&profile);
		  // ********* Write to file ********* 
		  ofstream saved_profile;
		  saved_profile.open(reconstruction.root_profile_filename.c_str() );
		  absorbTrans.write_profile(saved_profile, viterbi_path); 
		  cerr<<"Wrote profile to file: root_profile.sexpr\n";
		  saved_profile.close();
		}
	      
		  
	      std::cout<<"#=GF alignment_likelihood "<<-log(profile.forward_prob)/log(2) << endl; 
	      ofstream db_file;
	      state_path path = profile.sample_DP(
						  1, // sample only one path
						  0, // debugging log messages
						  false, // don't show alignments
						  false, // leaves only
						  reconstruction.viterbi // sample the viterbi path
						  );
	      alignString = profile.show_alignment( path, reconstruction.leaves_only); 
	      if (reconstruction.num_root_alignments > 1 && reconstruction.indel_filename != "None")
		{
		  if(reconstruction.loggingLevel>=1)
		    std::cerr<<"\nSampling " << reconstruction.num_root_alignments << " alignments at root level..."; 
		  double tot_ins=0, tot_ins_ext=0, tot_del=0, tot_del_ext=0; 

		  if (reconstruction.db_filename != "None")
		    {		  

		      db_file.open (reconstruction.db_filename.c_str());
		      reconstruction.tree.write_Stockholm(db_file);
		      //alignString = treeStream.str()
		    }
		  for (int samples = 1; samples<= reconstruction.num_root_alignments; samples++)
		    {
		      state_path path = profile.sample_DP(
							  1, // sample only one path
							  0, // debugging log messages
							  false, // don't show alignments
							  false, // leaves only
							  false // don't sample the viterbi path!
							  );
		      if (reconstruction.db_filename != "None")
			db_file << profile.show_alignment(path, reconstruction.leaves_only);
			  
		      MyMap<string, string> alignment = profile.alignment_map(path, false); 
		      IndelCounter indels(alignment, &reconstruction.tree); 
		      indels.gather_indel_info(false); 

		      tot_ins += indels.avg_insert_rate; 
		      if (indels.avg_insert_rate > 0.0)
			tot_ins_ext += indels.avg_insert_ext; 
		      tot_del += indels.avg_delete_rate; 
		      if (indels.avg_delete_rate > 0.0)
			tot_del_ext += indels.avg_delete_ext; 
		      
		      // do a simple average over alignments' inferred param estimates
		      if (samples == reconstruction.num_root_alignments)
			{
			  indels.avg_insert_rate = tot_ins / double(reconstruction.num_root_alignments); 
			  indels.avg_insert_ext = tot_ins_ext / double(reconstruction.num_root_alignments); 
			  indels.avg_delete_rate = tot_del / double(reconstruction.num_root_alignments); 
			  indels.avg_delete_ext = tot_del_ext / double(reconstruction.num_root_alignments); 

			  // write to file
			  ofstream indel_file;
			  indel_file.open (reconstruction.indel_filename.c_str());
			  indel_file << "### Indel information averaged over " << reconstruction.num_root_alignments << " alignments ###"; 
			  indels.display_indel_info(indel_file, false);
			  indel_file.close();

			}
		    }
		  if (reconstruction.db_filename != "None")
		    db_file.close();
		}
	    }
	}
    }
  if(reconstruction.loggingLevel>=1)
    std::cerr<<"done.\n";

  
  if (reconstruction.loggingLevel >= 1)
    std::cerr<<"\nFinished with alignment construction, now post-processing for ancestral characters and indels\n"; 
  // There's now an "alignment" created, either from input or via protpal
  // Convert it to a stockholm object - so we can do character recon w/ Xrate
  // There is no way that this is the easiest way to do this, but oh well:
  stringstream treeStream;
  reconstruction.tree.write_Stockholm(treeStream);
  alignString = treeStream.str() + alignString;

  istringstream stockStream(alignString);
  Sequence_database db; 
  Stockholm stk(1,1);
  if (reconstruction.input_alignment)
    stk.read_Stockholm(stockStream,db, 0,"_"); 
  else
    stk.read_Stockholm(stockStream,db);
  
  if (reconstruction.leaves_only)
    stk.write_Stockholm(std::cout);
  else // display all characters
    {
      ECFG_main ecfg; 
      ecfg.ancrec_CYK_MAP = true; 
      if (reconstruction.ancrec_postprob)
	{
	  ecfg.ancrec_postprob=true;
	  ecfg.min_ancrec_postprob = reconstruction.min_ancrec_postprob; 
	}

      SExpr_file gap_grammar_sexpr_file (reconstruction.gap_grammar_filename.c_str()); 
      SExpr_file grammar_sexpr_file (reconstruction.grammar_filename.c_str()); 

      if (reconstruction.input_alignment)
	{
	  grammar_sexpr_file = gap_grammar_sexpr_file; 
	  ecfg.gap_chars = sstring("_"); 
	}
      
      SExpr& grammar_ecfg_sexpr = grammar_sexpr_file.sexpr;

      // Get the alignment and grammar into the ecfg and get it ready for operations 
      ecfg.read_alignments(stk); 
      ecfg.read_grammars(&grammar_ecfg_sexpr); 
      ecfg.convert_sequences(); 

      // Train the Xrate grammar/chain, if requested
      if (reconstruction.train_grammar)
	{
	  if (reconstruction.loggingLevel >= 1)
	    {
	      std::cerr<<"\nTraining grammar used for character reconstruction, using starting point: \n\t ";
	      if (reconstruction.input_alignment)
		std::cerr << reconstruction.gap_grammar_filename << endl;
	      else
		std::cerr << reconstruction.grammar_filename << endl;
	    }
	  ecfg.train = "/dev/null";
	  ecfg.train_grammars();
	  ecfg.delete_trainers(); 
	}

      // 'Annotate' the alignment, using the current grammar
      if (reconstruction.loggingLevel >=1) 
	std::cerr<< "\tReconstructing ancestral characters conditional on ML indel history..."; 
      ecfg.annotate_alignments(); 

      Stockholm annotated = *(ecfg.stock_db.align.begin()); 

      if (reconstruction.loggingLevel >=1) 
	{
	  std::cerr<<"Done.\n";
	  std::cerr<<"\tDisplaying full ancestral alignment\n\n"; 
	}
      

      time(&end); 
      std::cout<< "#=GF TIME_MINUTES " << difftime (end, start)/60.0 <<endl; 
      if (reconstruction.xrate_output || reconstruction.ancrec_postprob)
	annotated.write_Stockholm(std::cout);
      else // display in bare-bones/ non-Xrate format
	{
	  Phonebook::iterator seq; 
	  int nameSize, maxNameLength = 0; 
	  sstring sequence; 
	  reconstruction.tree.write_Stockholm(std::cout);
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
		  
      // if requested, show what was inserted/deleted on the tree  (written to file)
      if (reconstruction.indel_filename != "None" && reconstruction.num_root_alignments == 1)
	{
	  if (reconstruction.loggingLevel >=1)
	    std::cerr<<"\nWriting indel information to file: " << reconstruction.indel_filename << endl; 
	  ofstream indel_file;
	  indel_file.open (reconstruction.indel_filename.c_str());
	  IndelCounter indels(annotated, &reconstruction.tree); 
	  indels.gather_indel_info(false); 
	  indels.display_indel_info(indel_file, reconstruction.per_branch);
	  indel_file.close();
	}
  
    }
  
  if (reconstruction.loggingLevel >= 1)
    {
      std::cerr<<"\nProtPal reconstruction completed without errors. \n";
    }
  
  
  return(0);
}

