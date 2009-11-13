#include "handel/transducer.h"
#include "handel/rivas.h"
#include "handel/multiwaydp.h"

#include "util/rnd.h"
#include "util/logfile.h"
#include "ecfg/ecfgsexpr.h"

// default prefix for auto-assigned node names
#define DEFAULT_NODE_PREFIX "Node"

int main(int argc, char* argv[])
{
  // initialise handler
  INIT_OPTS_LIST (opts, argc, argv, 2, "[options] <tree> <sequences>", "multiple alignment using a 'long indel' transducer composition & posterior decoding");

  // mutation model parameters
  double insert_rate;
  double delete_rate;
  double gap_extend_prob;
  sstring subst_model_filename;

  // add options to handler
  opts.newline();
  opts.print_title ("Mutation parameters");

  opts.add ("i -insert-rate", insert_rate = .1, "insertion rate");
  opts.add ("d -delete-rate", delete_rate = .11, "deletion rate");
  opts.add ("x -gap-extend", gap_extend_prob = .8, "gap extension probability");
  opts.add ("m -subst-model", subst_model_filename = "t/jc.hsm", "substitution rate matrix (xgram-like format)");

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
      if (delete_rate <= insert_rate)
	THROWEXPR ("Insert rate (" << insert_rate << ") must be less than delete rate (" << delete_rate << ")");

      // read in tree, sequences
      PHYLIP_tree tree;
      Sequence_database seqs;

      ifstream tree_file (opts.args[0].c_str());
      if (!tree_file) THROW String_exception ("Tree file not found: ", opts.args[0].c_str());
      CLOG(7) << "Reading tree\n";
      tree.read (tree_file);
      tree.force_binary();

      ifstream seq_file (opts.args[1].c_str());
      if (!seq_file) THROW String_exception ("Sequence file not found: ", opts.args[1].c_str());
      CLOG(7) << "Reading sequences\n";
      seqs.read_FASTA (seq_file);

      // create Transducer_alignment
      Affine_transducer_factory trans (insert_rate, delete_rate, gap_extend_prob);

      // read alphabet & substitution matrix
      SExpr_file param_sexpr (subst_model_filename.c_str());
      Alphabet alph ("uninitialized", 1);
      ECFG_builder::init_chain_and_alphabet (alph, trans.subst_model, param_sexpr.sexpr);
      trans.update_seq_scores();

      // convert sequences, discard wild sequences
      seqs.seqs2scores (alph);

      // re-order the tree nodes the way ETree likes it
      ETree etree (tree.nodes());
      PHYLIP_tree sorted_tree;

      vector<int> tree2etree (tree.nodes()), etree2tree;
      tree2etree[tree.root] = 0;
      etree2tree.push_back (tree.root);

      sorted_tree.add_node();
      sorted_tree.node_name.push_back (tree.node_specifier (tree.root));
      sorted_tree.root = 0;

      vector<Pair_transducer_scores> pair_trans_sc (etree.nodes(), Pair_transducer_scores (0));  // indexed by ETree node index
      pair_trans_sc[0] = trans.prior_pair_trans_sc();   // root-->subroot branch of ETree

      for_rooted_branches_pre (tree, b)  // this iteration visits nodes in depth-first preorder, as required
	{
	  const int dad = (*b).first, kid = (*b).second;
	  const double t = (*b).length;

	  tree2etree[kid] = etree2tree.size();
	  etree2tree.push_back (kid);

	  etree.parent[tree2etree[kid]] = tree2etree[dad];
	  CTAG(-1,TRANSDUCER) << "Creating ETree: adding branch from " << tree2etree[dad] << " to " << tree2etree[kid] << "\n";

	  pair_trans_sc[tree2etree[kid]] = trans.branch_pair_trans_sc (t);

	  sorted_tree.add_node (tree2etree[dad], t);
	  sorted_tree.node_name.push_back (tree.node_specifier (kid));
	}
      sorted_tree.rebuild_parents();

      if (CTAGGING(0,TRANSDUCER TRANSDUCER_ETREE))
	{
	  CL << "tree2etree: (" << tree2etree << ")\netree2tree: (" << etree2tree << ")\n";
	  CL << "sorted_tree: ";
	  sorted_tree.write (CL);
	  CL << "\n";
	}

      // build a jointly-normalised multi-transducer from the ETree
      EHMM_transducer_scores ehmm_trans_sc (etree, pair_trans_sc);

      // make lists of leaf, internal nodes
      vector<int> leaves, internals;
      for (int n = 0; n < sorted_tree.nodes(); ++n)
	if (sorted_tree.is_leaf(n))
	  leaves.push_back(n);
	else
	  internals.push_back(n);

      // figure out which multi-transducer states only emit to internal nodes
      vector<int> ehmm_null_states, ehmm_emit_states;
      ehmm_trans_sc.get_emit_and_null_states (leaves, ehmm_emit_states, ehmm_null_states);

      // make a smaller transducer with the null states trimmed out
      Eliminated_EHMM_transducer_scores compact_elim_trans_sc (ehmm_trans_sc, ehmm_null_states);

      // sort list of observed sequences by PHYLIP tree node name
      Sequence_database_index seqdb_index (seqs);
      vector<Score_profile*> seq (tree.nodes(), (Score_profile*) 0);
      for_const_contents (vector<int>, leaves, ln)
	{
	  const sstring& leaf_name = sorted_tree.node_specifier (*ln);
	  Score_profile* prof = &seqdb_index.name2profile(leaf_name)->prof_sc;
	  seq[*ln] = prof;
	  trans.align.row_name.push_back (leaf_name);
	  trans.align.prof.push_back (prof);
	}

      // initialize forward, backward & optacc DP matrices
      Transducer_backward_matrix back;
      back.trans_sc = &compact_elim_trans_sc;
      back.seq = seq;

      // allocate & fill forward & backward DP matrices
      back.alloc();
      back.fill();
      back.init_sumpairs_reward();

      // allocate & fill optimal accuracy DP matrix; get traceback alignment
      back.alloc_optacc();
      back.fill_optacc();
      trans.align.path = back.optacc_traceback();

      // print alignment
      cout << "# STOCKHOLM 1.0\n";
      trans.align.write_MUL (cout, alph);
      cout << "//\n";

    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
