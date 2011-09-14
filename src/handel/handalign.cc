#include "handel/transducer.h"
#include "handel/rivas.h"
#include "handel/tkftrans.h"
#include "handel/alitrans.h"
#include "handel/recorder.h"

#include "util/rnd.h"
#include "util/logfile.h"
#include "util/unixenv.h"
#include "ecfg/ecfgsexpr.h"

// default prefix for auto-assigned node names
#define DEFAULT_NODE_PREFIX "Node"

// paths to substitution models
#define JUKES_CANTOR_CHAIN_PATH "/data/handalign/jc.hsm"
#define PROT_CHAIN_PATH "/data/handalign/prot.hsm"

#define DEFAULT_CHAIN_PATH JUKES_CANTOR_CHAIN_PATH

int main(int argc, char* argv[])
{
  // initialise handler
  INIT_OPTS_LIST (opts, argc, argv, 1, "[options] <alignment>", "sample multiple alignments using a single-event 'long indel' transducer composition");

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
  sstring mcmc_sample_filename;
  bool use_best;

  // mutation model parameters
  double mean_seq_len;
  double delete_rate;
  double gap_len;
  sstring long_gap_cpt_weight_str;
  double long_gap_len;
  bool tkf91;

  // set up path to subst model
  sstring subst_model_filename, subst_model_help_string;
  subst_model_filename << Dart_Unix::get_DARTDIR() << DEFAULT_CHAIN_PATH;
  subst_model_help_string << "xrate-format substitution rate matrix (default is $" << DARTDIR_ENV_VAR << DEFAULT_CHAIN_PATH << ')';

  // set up alternate subst model
  sstring prot_model_alias, prot_model_help_string;
  prot_model_alias << "-m " << Dart_Unix::get_DARTDIR() << PROT_CHAIN_PATH;
  prot_model_help_string << "use alternate rate matrix ($" << DARTDIR_ENV_VAR << PROT_CHAIN_PATH << ')';

  // misc parameters
  sstring dotfile_dir;
  sstring comp_dir;
  sstring transducer_dump_file;

  // add options to handler
  opts.newline();
  opts.print_title ("Mutation parameters");

  opts.add ("l -seq-len", mean_seq_len = 100, "expected sequence length at equilibrium");
  opts.add ("d -delete-rate", delete_rate = .01, "deletion rate");
  opts.add ("g -gap-len", gap_len = 5, "expected deletion length");

  opts.newline();
  opts.add ("lw -long-gap-weight", long_gap_cpt_weight_str = "", 0);
  opts.print ("-lw,--long-gap-weight <real>\tmixture component weight for long gaps (default is zero)\n");
  opts.add ("lg -long-gap-len", long_gap_len = 10, "expected long deletion length");

  opts.newline();
  opts.add ("m -subst-model", subst_model_filename, subst_model_help_string.c_str(), false);
  opts.add ("p -prot-model", prot_model_alias.c_str(), prot_model_help_string.c_str(), false);
  opts.add ("tkf91", tkf91 = false, "use TKF91 model instead of long indel model");

  opts.newline();
  opts.print_title ("Control of MCMC sampler");

  opts.add ("s -samples", annealing_steps_per_node = 1, "number of MCMC sampling steps per node of the tree");
  opts.add ("ts -total-samples", annealing_steps_total = "", "specify *total* (not per-node) number of MCMC sampling steps", false);
  opts.add ("xl -execlike", exec_filename = "", "use external alignment-likelihood executable for Metropolis-Hastings sampling", false);
  Likelihood_executable::add_help (&opts);
  opts.add ("fn -first-node", first_node = "", "start sampling at specified node (useful for debugging)", false);
  opts.add ("fb -force-binary", force_binary = true,    "force binary tree");
  opts.add ("rs -redsuch", use_Redelings_Suchard = false,    "use Redelings-Suchard proposal scheme when sampling");

  opts.newline();
  opts.add ("ub -use-best",    use_best = false, "report best (most probable) alignment, rather than final sampled one");
  opts.add ("r -refine",     refine = false,  "periodically, do iterative refinement to optimize alignment (does not use Redelings-Suchard)");
  opts.add ("rp -refine-period",    refine_period_per_node = 10, "number of sampling steps between each refinement");

  opts.newline();
  opts.print_title ("MCMC move rates");

  opts.add ("bf -branch-freq", shuffler.branch_realign_rate = 1, "relative rate of branch-sampling: realigning adjacent sequences by DP");
  opts.add ("nf -node-freq", shuffler.node_realign_rate = 1, "relative rate of node-sampling: adding/removing unobserved residues by DP");
  opts.add ("ff -flip-freq", shuffler.branch_swap_rate = 0, "relative rate of branch-flipping: exchanging a node & its niece, then re-aligning by DP");
  opts.add ("sf -slide-freq", shuffler.node_slide_rate = 0, "relative rate of node-sliding moves that preserve topology & total branch length");
  opts.add ("lf -length-freq", shuffler.branch_scale_rate = 0, "relative rate of sampling branch lengths");
  opts.add ("pf -param-freq", shuffler.indel_param_sampling_rate = 0, "relative rate of sampling indel parameters");

  opts.newline();
  opts.print_title ("Shorthands for individual moves");

  opts.add ("branch", "-ts 1 -bf 1 -nf 0 -sf 0 -lf 0 -ff 0", "resample a random branch");
  opts.add ("node", "-ts 1 -bf 0 -nf 1 -sf 0 -lf 0 -ff 0", "resample a random node's outgoing branches");
  opts.add ("slide", "-ts 1 -bf 0 -nf 0 -sf 1 -lf 0 -ff 0", "slide a node along a random branch");
  opts.add ("length", "-ts 1 -bf 0 -nf 0 -sf 0 -lf 1 -ff 0", "resample a random branch's length");
  opts.add ("flip", "-ts 1 -bf 0 -nf 0 -sf 0 -lf 0 -ff 1", "flip a random niece-node pair");

  opts.newline();
  opts.print_title ("MCMC sampling transcripts");

  opts.add ("dt -dump-transducer", transducer_dump_file = "", "save phylocomposer file before doing any sampling", false);
  opts.add ("c -composition", comp_dir = "", "save phylocomposer files to this directory during sampling", false);
  opts.add ("dot -dotfile", dotfile_dir = "", "save composite transducer dotfiles to this directory", false);
  opts.add ("af -align-file", mcmc_sample_filename = "", "write Stockholm alignments to file during sampling", false);

  // HMMoC adapter

  HMMoC_adapter_options hmmoc_opts;
  hmmoc_opts.init_opts_list (opts, false);  // don't offer the option to dump the HMMoC file

  opts.add ("tb -time-band", use_time_dependent_banding = false, "use time-dependent banding, if & when it gives a tighter band", false);
  opts.add ("bc -band-coeff", banding_coefficient = DEFAULT_BANDING_COEFFICIENT, "time-dependent banding coefficient: determines width of band");
  opts.add ("cb -centroid-band", use_centroid_band = false, "use centroid banding (Redelings-Suchard MCMC only)", false);
  opts.add ("cbw -centroid-band-width", centroid_band_width = DEFAULT_CENTROID_BAND_WIDTH, "width of centroid band");

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
      // display initial logging messages
      Score_fns::describe_scoring_scheme (CLOG(8));

      // check params
      if (mean_seq_len <= 0 || gap_len <= 0 || long_gap_len <= 0)
	THROWEXPR ("All expected lengths (equilibrium sequence length & gap lengths) must be positive");

      if (refine && !use_best && mcmc_sample_filename.empty())
	CLOGERR << "When using the --refine option, it's usually a good idea either to turn on\n"
		<< "--use-best or to specify an alignment output file with --align-file,\n"
		<< "otherwise the best refined alignments will be discarded.\n";

      // read in starting alignment
      Sequence_database seqs;
      Stockholm stock;

      ifstream align_file (opts.args[0].c_str());
      if (!align_file) THROW String_exception ("Alignment file not found: ", opts.args[0].c_str());
      CLOG(7) << "Reading initial alignment\n";

      stock.read_Stockholm (align_file, seqs);

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
	  cpt_gap_extend.push_back (1. - 1. / (gap_len + 1.));

	  cpt_weight.push_back (long_gap_cpt_weight);
	  cpt_gap_extend.push_back (1. - 1. / (long_gap_len + 1.));

	  trans_ptr = new Convex_transducer_factory (gamma, delete_rate, cpt_weight, cpt_gap_extend);
	}
      else
	trans_ptr = new Affine_transducer_factory (gamma, delete_rate, 1. - 1. / (gap_len + 1.));
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

      // read alphabet & substitution matrix
      SExpr_file param_sexpr (subst_model_filename.c_str());
      Alphabet alph;
      ECFG_builder::init_chain_and_alphabet (alph, trans.subst_model, param_sexpr.sexpr);
      trans.update_seq_scores();

      // convert sequences, discard wild sequences
      seqs.seqs2scores (alph);
      stock.discard_wild_sequences (alph);

      // TODO: if no tree, estimate by neighbor-joining (postpone this for a while: can use xrate at first)

      // get tree, convert to binary
      Stockholm_tree tree (stock, true);
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
