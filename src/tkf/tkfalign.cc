#include "tkf/tkfsequence.h"
#include "tkf/tkfopts.h"
#include "tkf/tkfdata.h"
#include "util/rnd.h"
#include "util/logfile.h"

int main(int argc, char* argv[])
{
  TKF_opts opts (argc, argv, 1);
  opts.short_description = "align sequences using a guide tree under the TKF model";
  opts.syntax = "[options] <treefile> <sequences>";

  bool    force_binary;
  sstring tree_save_file;
  int    annealing_steps;
  double kT_start;
  double kT_end;

  Tree_shuffler shuffler;

  int    refine_period;
  bool   use_best;
  bool   viterbi_progressive;
  bool   refine;
  bool   sample_internal_sequences;
  sstring benchmark_filename;
  sstring exec_filename;
  sstring mcmc_sample_filename;

  opts.newline();
  opts.print_title ("Deprecated options");

  opts.add ("-vprog",      viterbi_progressive = 1,            "do initial progressive alignment greedily, i.e. Viterbi");
  opts.add ("-ktstart",    kT_start = 1,                       "initial value of kT for annealing phase");
  opts.add ("-ktend",      kT_end = 1,                         "final value of kT for annealing phase");
  opts.add ("-benchmark", benchmark_filename = "",         "reference alignment for benchmarking", 0);

  opts.newline();
  opts.print_title ("Input format options");

  opts.add ("fb -force-binary", force_binary = 1,    "force tree to be binary, including root");
  opts.add ("st -save-tree",    tree_save_file = "", "save initial (potentially binary-ised) tree to file", 0);
  opts.add_align_detection();

  opts.newline();
  opts.print_title ("MCMC alignment sampling options");

  opts.add ("s -samples", annealing_steps = 0, "number of MCMC sampling steps per node of the tree");
  opts.add ("xl -execlike", exec_filename = "", "use external alignment-likelihood executable for Metropolis-Hastings sampling", 0);
  Likelihood_executable::add_help (&opts);

  opts.add ("af -align-file", mcmc_sample_filename = "", "write Stockholm alignments to file during sampling", 0);  opts.add ("i -internal",   sample_internal_sequences = 0,"sample sequences at internal nodes instead of using Felsenstein wildcards");

  opts.newline();
  opts.add ("ub -use-best",    use_best = 1, "report best alignment, rather than final sampled one");
  opts.add ("r -refine",     refine = 0,                         "periodically, do iterative refinement to optimize alignment");
  opts.add ("rp -refine-period",    refine_period = 100,                "number of sampling steps between each refinement");

  opts.newline();
  opts.add ("bf -branch-freq", shuffler.branch_realign_rate = 1, "relative rate of branch-sampling: realigning adjacent sequences by DP");
  opts.add ("nf -node-freq", shuffler.node_realign_rate = 1, "relative rate of node-sampling: adding/removing unobserved residues by DP");
  opts.add ("sf -slide-freq", shuffler.node_slide_rate = 0, "relative rate of sliding moves that preserve topology & total branch length");
  opts.add ("lf -length-freq", shuffler.branch_scale_rate = 0, "relative rate of sampling branch lengths");
  opts.add ("ff -flip-freq", shuffler.branch_swap_rate = 0, "relative rate of branch-flipping: exchanging a node & its niece over an ungapped branch");

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 2) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
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
      if (refine && !use_best)
	CLOGERR << "When using the --refine option, it's usually a good idea to turn on --use-best, otherwise the best refined alignments may be discarded.\n";

      // read in sequences; autodetect aligned/unaligned and DNA/protein
      Sequence_database seqs;
      Alignment         align;

      ifstream seq_file (opts.args[1].c_str());
      if (!seq_file) THROW String_exception ("Sequence file not found: ", opts.args[1].c_str());
      
      bool prealigned = opts.detect_aligned(seq_file);
      if (prealigned)
	{
	  CLOG(7) << "Reading alignment\n";
	  align.read_MUL (seq_file, seqs);
	}
      else
	{
	  CLOG(7) << "Reading sequences\n";
	  seqs.read_FASTA (seq_file);
	}

      opts.detect_alphabet (seqs);

      // set up TKF_align structure; read in tree
      TKF_params params = opts.params();
      TKF_align  tkf (params);

      seqs.seqs2scores (tkf.alphabet());       // decode sequence data

      ifstream tree_file (opts.args[0].c_str());
      if (!tree_file) THROW String_exception ("Tree file not found: ", opts.args[0].c_str());
      CLOG(7) << "Reading tree\n";
      tkf.read_PHYLIP (tree_file);

      if (force_binary) { tkf.tree.force_binary(); tkf.reset_maps(); tkf.tree_changed(); }
      if (tree_save_file != "") { ofstream tree_save_stream (tree_save_file.c_str()); tkf.tree.write (tree_save_stream); }

      // attach sequences to TKF_align structure; do preliminary alignment phase
      if (prealigned)
	{
	  align.discard_wild_sequences (tkf.alphabet());
	  tkf.set_alignment(align);
	  CLOG(7) << "Mapping alignment to tree\n";
	  tkf.build_maps_from_names();
	  tkf.optimise_missing_nodes (sample_internal_sequences);
	}
      else
	{
	  tkf.make_empty_alignment();
	  CLOG(7) << "Mapping sequences to tree\n";
	  tkf.attach_sequences(seqs);
	  CLOG(6) << "Doing first-pass progressive alignment\n";
	  if (viterbi_progressive) tkf.viterbi_progressive_alignment (sample_internal_sequences);
	  else tkf.sample_progressive_alignment (kT_start, sample_internal_sequences);
	  CLOG(6) << "Preliminary alignment phase complete\n";
	}

      // load benchmark alignment
      if (benchmark_filename != "")
	{
	  tkf.read_benchmark_alignment (benchmark_filename.c_str(), seqs);
	  tkf.log_benchmark_results ("Initial alignment");
	}

      // set up Metropolis-Hastings sampling
      Likelihood_executable like_exec (tkf.alphabet(), exec_filename.c_str());
      if (exec_filename.size())
	  tkf.target_loglike = &like_exec;
      ofstream* mcmc_sample_stream = 0;
      if (mcmc_sample_filename.size()) {
	  mcmc_sample_stream = new ofstream (mcmc_sample_filename.c_str());
	  tkf.sample_stream = mcmc_sample_stream;
      }

      // do sampling/refinement phase
      shuffler.set_tree (tkf.tree);
      vector<int> alignment_scores;
      const Score final_score = tkf.anneal (kT_start, kT_end, annealing_steps * tkf.tree.nodes(), shuffler, alignment_scores, sample_internal_sequences, use_best, refine, refine_period * tkf.tree.nodes());

      // close Metropolis-Hastings sampling stream
      if (mcmc_sample_stream) {
	  mcmc_sample_stream->close();
	  delete mcmc_sample_stream;
      }

      // log final results
      CLOG(6) << "Sampling phase complete\n";
      if (CLOGGING(3)) { CL << "Final alignment:\n"; tkf.write_Stockholm (CL, tkf.alphabet()); }
      
      tkf.log_benchmark_results ("Final alignment");

      if (CTAGGING(2,BREAKDOWN))
	tkf.show_node_score_breakdown (CL, tkf.tree.root);

      // output final results
      const Score null_score = tkf.null_score();
      const Score odds_ratio_score = ScorePMul (final_score, -null_score);

      cout << Stockholm_header;
      tkf.tree.write_Stockholm (cout);
      tkf.align.write_MUL (cout, tkf.alphabet());
      cout << Stockholm_file_annotation << ' ' << Stockholm_comment_tag << ' ';
      cout << "LgP(alignment|tree)= " << Score2Bits (final_score) << " bits; ";
      cout << "LgP(sequences|unrelated)= " << Score2Bits (null_score) << " bits; ";
      cout << "LgOddsRatio= " << Score2Bits (odds_ratio_score) << " bits\n";
      cout << Stockholm_file_annotation << ' ' << Stockholm_bit_score_tag << ' ' << Score2Bits (odds_ratio_score)  << '\n';
      cout << Stockholm_alignment_separator << '\n';
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
