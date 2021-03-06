#include "handel/transducer.h"
#include "handel/rivas.h"
#include "handel/tkftrans.h"
#include "handel/alitrans.h"
#include "handel/recorder.h"
#include "handel/initalign.h"

#include "util/rnd.h"
#include "util/logfile.h"
#include "util/unixenv.h"
#include "ecfg/ecfgsexpr.h"

// default prefix for auto-assigned node names
#define DEFAULT_NODE_PREFIX "Node"

// paths to substitution models
#define JUKES_CANTOR_CHAIN_PATH "/data/handalign/jc.hsm"
#define PROT_CHAIN_PATH "/data/handalign/prot.hsm"

#define DEFAULT_CHAIN_PATH PROT_CHAIN_PATH

// branch length for initial guesstimate alignment
#define INIT_BRANCH_LEN 0.5


// main
int main(int argc, char* argv[])
{
  // initialise handler
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <Stockholm alignment or FASTA sequences file>", "MCMC sampler over multiple alignments, phylogenies, and evolutionary parameters");

  // dummy options for Handel_base::anneal()
  const double kT_start = 1.;
  const double kT_end = 1.;
  const bool sample_internal_sequences = false;
  bool refine;
  int refine_period_per_node;

  // option variables
  Tree_shuffler shuffler;
  sstring first_node;
  int annealing_steps_per_node;
  sstring annealing_steps_total;
  bool use_Redelings_Suchard;
  bool use_time_dependent_banding;
  double banding_coefficient;
  bool use_centroid_band;
  int centroid_band_width;
  bool force_binary;
  sstring exec_filename;
  bool factor_indels;
  sstring mcmc_sample_filename;
  bool use_best;

  // mutation model parameters
  double mean_seq_len;
  double delete_rate;
  double gap_len;
  sstring long_gap_cpt_weight_str;
  double long_gap_len;
  bool tkf91;

  // proposal fn variables
  double indel_proposal_variance; 
  
  // set up path to subst model
  sstring subst_model_filename, subst_model_help_string;
  subst_model_filename << Dart_Unix::get_DARTDIR() << DEFAULT_CHAIN_PATH;
  subst_model_help_string << "xrate-format substitution rate matrix (default is $" << DARTDIR_ENV_VAR << DEFAULT_CHAIN_PATH << ')';

  // set up alternate subst models
  sstring prot_model_alias, prot_model_help_string;
  prot_model_alias << "-m " << Dart_Unix::get_DARTDIR() << PROT_CHAIN_PATH;
  prot_model_help_string << "use alternate rate matrix ($" << DARTDIR_ENV_VAR << PROT_CHAIN_PATH << ')';

  sstring jukes_cantor_model_alias, jukes_cantor_model_help_string;
  jukes_cantor_model_alias << "-m " << Dart_Unix::get_DARTDIR() << JUKES_CANTOR_CHAIN_PATH;
  jukes_cantor_model_help_string << "use alternate rate matrix ($" << DARTDIR_ENV_VAR << JUKES_CANTOR_CHAIN_PATH << ')';

  // misc parameters
  sstring dotfile_dir;
  sstring comp_dir;
  sstring transducer_dump_file;

  // add options to handler
  opts.newline();
  opts.print_title ("Mutation parameters");

  opts.add ("l -seq-len", mean_seq_len = 100, "expected sequence length at equilibrium");
  opts.add ("d -delete-rate", delete_rate = .01, "deletion rate");
  opts.add ("g -gap-len", gap_len = 4, "expected deletion length");

  opts.newline();
  opts.add ("lw -long-gap-weight", long_gap_cpt_weight_str = "", 0);
  opts.print ("-lw,--long-gap-weight <real>\tmixture component weight for long gaps (default is zero)\n");
  opts.add ("lg -long-gap-len", long_gap_len = 9, "expected long deletion length");

  opts.newline();
  opts.add ("m -subst-model", subst_model_filename, subst_model_help_string.c_str(), false);
  opts.add ("p -prot-model", prot_model_alias.c_str(), prot_model_help_string.c_str(), false);
  opts.add ("jc -jukes-cantor-model", jukes_cantor_model_alias.c_str(), jukes_cantor_model_help_string.c_str(), false);
  opts.add ("tkf91", tkf91 = false, "use TKF91 model instead of long indel model");

  opts.newline();
  opts.print_title ("Control of MCMC sampler");

  opts.add ("s -samples", annealing_steps_per_node = 1, "number of MCMC sampling steps per node of the tree");
  opts.add ("ts -total-samples", annealing_steps_total = "", "specify *total* (not per-node) number of MCMC sampling steps", false);
  opts.add ("xl -exec-like", exec_filename = "", "use external alignment-likelihood executable for Metropolis-Hastings sampling", false);
  Likelihood_executable::add_help (&opts);
  opts.add ("fi -factor-indels", factor_indels = false, "multiply external alignment-likelihood by transducer indel-likelihood");
  opts.add ("fn -first-node", first_node = "", "start sampling at specified node (useful for debugging)", false);
  opts.add ("fb -force-binary", force_binary = true,    "force binary tree");
  opts.add ("rs -redsuch", use_Redelings_Suchard = true,    "use Redelings-Suchard proposal scheme when sampling");
  opts.add ("ipv -indel-proposal-variance", indel_proposal_variance = 0.01,    "Variance of proposal function for indel rate moves");

  opts.newline();
  opts.print_title ("MCMC move rates");

  opts.add ("bf -branch-freq", shuffler.branch_realign_rate = 1, "relative rate of branch-sampling: realigning parent-child pairs by DP. This resamples the alignment, keeping the tree fixed");
  opts.add ("nf -node-freq", shuffler.node_realign_rate = 1, "relative rate of node-sampling: realigning grandparent-sibling triplets by DP. This resamples the alignment, keeping the tree fixed");
  opts.add ("ff -flip-freq", shuffler.branch_swap_rate = 1, "relative rate of branch-flipping: exchanging a node & its niece, then re-aligning by DP. This is the only move that resamples tree topology");
  opts.add ("sf -slide-freq", shuffler.node_slide_rate = 1, "relative rate of node-sliding moves that preserve topology & total branch length. This resamples tree branch lengths, but not tree topology");
  opts.add ("lf -length-freq", shuffler.branch_scale_rate = 1, "relative rate of sampling branch lengths. This does not change alignment or tree topology, only tree branch lengths");
  opts.add ("ipf -indel-param-freq", shuffler.indel_param_sampling_rate = 1, "relative rate of sampling indel parameters. This does not change tree or alignment");
  opts.add ("spf -subst-param-freq", shuffler.subst_param_sampling_rate = 1, "relative rate of sampling substitution parameters. This does not change tree or alignment");

  opts.newline();
  opts.print_title ("Shorthands for individual moves");

  opts.add ("branch", "-ts 1 -bf 1 -nf 0 -sf 0 -lf 0 -ff 0 -ipf 0 -spf 0", "resample a random branch");
  opts.add ("node", "-ts 1 -bf 0 -nf 1 -sf 0 -lf 0 -ff 0 -ipf 0 -spf 0", "resample a random node's outgoing branches");
  opts.add ("slide", "-ts 1 -bf 0 -nf 0 -sf 1 -lf 0 -ff 0 -ipf 0 -spf 0", "slide a node along a random branch");
  opts.add ("length", "-ts 1 -bf 0 -nf 0 -sf 0 -lf 1 -ff 0 -ipf 0 -spf 0", "resample a random branch's length");
  opts.add ("flip", "-ts 1 -bf 0 -nf 0 -sf 0 -lf 0 -ff 1 -ipf 0 -spf 0", "flip a random niece-node pair");

  opts.newline();
  opts.print_title ("MCMC sampling transcripts");

  opts.add ("dt -dump-transducer", transducer_dump_file = "", "save phylocomposer file before doing any sampling", false);
  opts.add ("c -composition", comp_dir = "", "save phylocomposer files to this directory during sampling", false);
  opts.add ("dot -dotfile", dotfile_dir = "", "save composite transducer dotfiles to this directory", false);
  opts.add ("af -align-file", mcmc_sample_filename = "", "write Stockholm alignments to file during sampling", false);

  opts.newline();
  opts.print_title ("Stochastic search mode");

  opts.newline();
  opts.add ("ub -use-best",    use_best = false, "report best (most probable) alignment, rather than final sampled one");
  opts.add ("r -refine",     refine = false,  "periodically, do iterative refinement to optimize alignment (does not use Redelings-Suchard)");
  opts.add ("rp -refine-period",    refine_period_per_node = 10, "number of sampling steps between each refinement");

  // HMMoC adapter

  HMMoC_adapter_options hmmoc_opts;
  hmmoc_opts.init_opts_list (opts, false);  // don't offer the option to dump the HMMoC file

  opts.add ("tb -time-band", use_time_dependent_banding = false, "use time-dependent banding, if & when it gives a tighter band", false);
  opts.add ("bc -band-coeff", banding_coefficient = DEFAULT_BANDING_COEFFICIENT, "time-dependent banding coefficient: determines width of band");
  opts.add ("cb -centroid-band", use_centroid_band = false, "use centroid banding (Redelings-Suchard MCMC only)", false);
  opts.add ("cbw -centroid-band-width", centroid_band_width = DEFAULT_CENTROID_BAND_WIDTH, "width of centroid band");

  opts.newline();
  opts.newline();
  opts.print ("More documentation available at http://biowiki.org/HandAlign\n");
  
  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // issue warning if hmmoc adapter not in use
#if defined(HMMOC_INCLUDED) && HMMOC_INCLUDED
      if (!hmmoc_opts.try_to_use_hmmoc_adapter)
	CLOGERR << "\nWarning: with the HMMoC adapter turned off, MCMC moves that use dynamic programming may be slow.\n\n";
#else
      if (!hmmoc_opts.try_to_use_hmmoc_adapter)
	CLOGERR << "\nWarning: without HMMoC installed, MCMC moves that use dynamic programming may be slow.\n\n";
#endif

      // display initial logging messages
      Score_fns::describe_scoring_scheme (CLOG(8));

      // check params
      if (mean_seq_len <= 0)
	THROWEXPR ("Expected equilibrium sequence length must be positive");

      if (gap_len <= 1 || long_gap_len <= 1)
	THROWEXPR ("Expected gap lengths must be greater than 1");

      if (refine && !use_best && mcmc_sample_filename.empty())
	CLOGERR << "When using the --refine option, it's usually a good idea either to turn on\n"
		<< "--use-best or to specify an alignment output file with --align-file,\n"
		<< "otherwise the best refined alignments will be discarded.\n";

      // create Transducer_alignment
      Transducer_alignment_with_subst_model* trans_ptr = 0;
      const double gamma = 1 - 1 / (mean_seq_len + 1);
      if (tkf91)
	trans_ptr = new TKF91_transducer_factory (gamma * delete_rate, delete_rate);
      else if (long_gap_cpt_weight_str.size())
	{
	  vector<double> cpt_weight, cpt_gap_extend;
	  const double long_gap_cpt_weight = long_gap_cpt_weight_str.to_double();

	  cpt_weight.push_back (1. - long_gap_cpt_weight);
	  cpt_gap_extend.push_back (1. - 1. / gap_len);

	  cpt_weight.push_back (long_gap_cpt_weight);
	  cpt_gap_extend.push_back (1. - 1. / long_gap_len);

	  trans_ptr = new Convex_transducer_factory (gamma, delete_rate, cpt_weight, cpt_gap_extend);
	}
      else
	trans_ptr = new Affine_transducer_factory (gamma, delete_rate, 1. - 1. / gap_len);

      Transducer_alignment_with_subst_model& trans (*trans_ptr);

      // populate Transducer_alignment
      trans.hmmoc_opts = hmmoc_opts;
      trans.use_Redelings_Suchard = use_Redelings_Suchard;
      trans.use_banding_coefficient = use_time_dependent_banding;
      trans.banding_coefficient = banding_coefficient;
      trans.use_centroid_band = use_centroid_band;
      trans.centroid_band_width = centroid_band_width;
      trans.dotfile_recorder.directory = dotfile_dir;
      trans.composition_recorder.directory = comp_dir;
      trans.indel_proposal_variance = indel_proposal_variance; 
      
      // read alphabet & substitution matrix
      SExpr_file param_chain_alph_sexpr (subst_model_filename.c_str());
      Alphabet alph;
      ECFG_builder::init_param_chain_alphabet (alph, trans.ems, trans.subst_pscores, trans.subst_pcounts, trans.subst_mutable_pgroups, param_chain_alph_sexpr.sexpr);
      trans.update_seq_scores();

      // read in starting alignment
      Sequence_database seqs;
      Stockholm stock;
      Stockade stockade;
      Stockholm_database stock_db;

      ifstream align_file (opts.args[0].c_str());
      if (!align_file) THROW String_exception ("Alignment file not found: ", opts.args[0].c_str());
      CLOG(7) << "Reading initial alignment\n";

      const bool is_Stockholm_input = stock_db.read_Stockholm_or_FASTA (align_file, seqs);
      if (is_Stockholm_input)
	{
	  if (stock_db.size() != 1)
	    THROWEXPR ("Error - expected only one alignment file on the input");
	  stock = stock_db.align.front();
	}
      else
	{
	  // FASTA input, so need to estimate initial alignment
	  FASTA_sequence_database fasta_seq_db (seqs, &alph);
	  stockade = Stockade_initializer::align (fasta_seq_db, alph, trans, INIT_BRANCH_LEN);
	  stock = stockade.align;
	  stock.clear_annot();
	}

      // convert sequences, discard wild sequences
      seqs.seqs2scores (alph);
      for_contents (vector<Named_profile>, stockade.np, np)
	np->seq2score (alph);

      stock.discard_wild_sequences (alph);

      // get tree, convert to binary
      Stockholm_tree tree (stock, false);
      const bool tree_parsed = tree.nodes() > 0;
      if (!tree_parsed)
	{
	  Stockade_initializer::build_tree (stock, trans.submat_factory());
	  tree = Stockholm_tree (stock, true);
	}
      if (force_binary)
	tree.force_binary();

      // in a postorder traversal of tree nodes: (i) assign names to unlabeled nodes and (ii) fill in any missing rows in the alignment
      set<int> new_rows;
      for_rooted_nodes_post (tree, b)
	{
	  const Phylogeny::Node node = (*b).second;
	  // assign node names
	  if (tree.node_name[node].size() == 0)
	    {
	      sstring new_node_name, prefix;
	      prefix << DEFAULT_NODE_PREFIX << node;
	      for (int n = 0; true; ++n)
		{
		  new_node_name = prefix;
		  if (n > 0)
		    new_node_name << '.' << n;
		  if (find (tree.node_name.begin(), tree.node_name.end(), new_node_name) == tree.node_name.end())
		    break;
		}
	      tree.node_name[node] = new_node_name;
	      CTAG(5,HANDALIGN) << "Node #" << node << " of tree lacks a name; auto-assigning unique name '" << new_node_name << "'\n";
	    }

	  // if alignment lacks a row with this name, then create a row with as many gaps as possible
	  if (stock.row_index.find (tree.node_name[node]) == stock.row_index.end())
	    {
	      if (tree.is_leaf (node))
		CLOGERR << "\nWarning: the alignment is missing row '" << tree.node_name[node] << "',\n"
			<< "although it's present in the tree. This program doesn't do reconstruction\n"
			<< "of leaf nodes. An all-gap row will be added to the alignment, which may\n"
			<< "not be the effect you wanted.\n\n"
			<< "If the effect you wanted was a hold-one-out cross-validation test,\n"
			<< "i.e. deleting a row from an existing alignment then 'reconstructing' it,\n"
			<< "a better bet would be to re-root the tree at the missing node (or its parent).\n\n";
		  
	      CTAG(5,HANDALIGN) << "Creating a minimally-gapped alignment row for tree node '" << tree.node_name[node] << "'\n";
	      const int new_row = stock.add_row (stock.path.create_empty_row(), tree.node_name[node], (Named_profile*) 0);
	      new_rows.insert (new_row);
	      for (int col = 0; col < stock.columns(); ++col)
		for_rooted_children (tree, node, child)
		  if (stock.path (stock.row_index[tree.node_name[*child]], col))
		    {
		      stock.path[new_row][col] = 1;
		      break;
		    }
	    }
	}

      // if new rows added, prune unnecessary residues
      if (new_rows.size())
	{
	  // in a preorder traversal, trim unnecessary residues from newly-added rows
	  for_rooted_nodes_pre (tree, b)
	    {
	      const Phylogeny::Node node = (*b).second;
	      const int row = stock.row_index[tree.node_name[node]];
	      if (new_rows.find(row) != new_rows.end())
		{
		  const int parent = tree.parent[node];
		  const int parent_row = parent < 0 ? -1 : stock.row_index[tree.node_name[parent]];

		  vector<int> child_rows;
		  for_rooted_children (tree, node, child)
		    child_rows.push_back (stock.row_index[tree.node_name[*child]]);

		  for (int col = 0; col < stock.columns(); ++col)
		    if (stock.path[row][col])
		      if (parent_row < 0 ? true : !stock.path[parent_row][col])
			{
			  int n_children = 0;
			  for_const_contents (vector<int>, child_rows, child_row)
			    if (stock.path[*child_row][col])
			      ++n_children;
			  if (n_children < 2)
			    stock.path[row][col] = 0;
			}
		}
	    }

	  // print inferred gap structure
	  CTAG(6,HANDALIGN) << "Initial guess at indel structure of ancestors:\n";
	  stock.write_Stockholm (CL);
	}

      // initialise tree & alignment
      trans.set_tree (tree);
      trans.set_alignment (stock);
      trans.build_maps_from_names();

      if (transducer_dump_file)
	{
	  ofstream transducer_dump_stream (transducer_dump_file.c_str());
	  trans.dump (transducer_dump_stream);
	}

      // get actual number of sampling steps
      const int annealing_steps = annealing_steps_total.size() ? annealing_steps_total.to_int() : (annealing_steps_per_node * tree.nodes());
      const int refine_period = refine_period_per_node * tree.nodes();

      // set up Metropolis-Hastings sampling
      Likelihood_executable like_exec (trans.alphabet(), exec_filename.c_str());
      if (exec_filename.size())
	  trans.target_loglike = &like_exec;
      trans.factor_indels_into_target_loglike = factor_indels;
      ofstream* mcmc_sample_stream = 0;
      if (mcmc_sample_filename.size())
	{
	  mcmc_sample_stream = new ofstream (mcmc_sample_filename.c_str());
	  trans.sample_stream = mcmc_sample_stream;
	}

      // do sampling/refinement phase
      shuffler.set_tree (trans.tree);
      if (first_node.size())
	shuffler.cue (first_node);
      vector<int> alignment_scores;
      const bool refine_node_triplets = !use_Redelings_Suchard;   // do not refine node triplets if Redelings-Suchard is being used (memory-saving measure)
      trans.anneal (kT_start, kT_end, annealing_steps, shuffler, alignment_scores, sample_internal_sequences, use_best, refine, refine_period, refine_node_triplets);

      // log final results
      CLOG(6) << "Sampling phase complete\n";
      if (CLOGGING(3))
	{
	  CL << "Final alignment:\n";
	  trans.write_Stockholm_with_score (CL, Stockholm_alignment_type_final, annealing_steps);
	}

      // log score breakdown
      if (CTAGGING(2,BREAKDOWN))
	trans.show_node_score_breakdown (CL, trans.tree.root);

      // output final alignment, with score
      trans.write_Stockholm_with_score (cout, Stockholm_alignment_type_final, annealing_steps);
      if (mcmc_sample_stream)
	trans.write_Stockholm_with_score (*mcmc_sample_stream, Stockholm_alignment_type_final, annealing_steps);

      // close Metropolis-Hastings sampling stream
      if (mcmc_sample_stream)
	{
	  mcmc_sample_stream->close();
	  delete mcmc_sample_stream;
	}

      // cleanup
      delete trans_ptr;
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
