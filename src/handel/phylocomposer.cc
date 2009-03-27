#include "handel/movement.h"
#include "handel/recorder.h"
#include "handel/hmmoc_adapter.h"

int main(int argc, char* argv[])
{
  // initialise options handler
  INIT_OPTS_LIST (opts, argc, argv, -1, "[options] [<transducer file>]", "perform transducer composition on a phylogenetic tree");

  // logging options
  opts.add ("w -warn", "-log TRANSDUCER_WARNINGS", "log transducer-related warning messages");
  opts.newline();

  // initialize Handel_movement
  Handel_movement movement;

  // display options

  opts.print_title ("Transducer display options");

  opts.add ("c -composite", movement.composite = true, "display composite transducer");
  opts.add ("cd -compdot", movement.dotfile = "", "save composite transducer as graphviz dotfile", false);
  opts.add ("cp -composite-path", movement.show_constrained_composite_paths = false, "show composite paths for constrained branches");
  opts.add ("a -acyclic", movement.acyclic = false, "display acyclic composite transducer, from which \"unobserved\" states are eliminated");
  opts.add ("q -quiet", movement.quiet = false, "don't display branch transducers");
  opts.add ("qq -hush", "--quiet --nocomposite", "don't display any transducers");
  opts.add ("al -alignment", movement.stockfile = "", "save input sequences & paths as a Stockholm alignment", false);

  // DP options

  opts.newline();
  opts.print_title ("Algorithm options");

  opts.add ("sim -simulate", movement.simulate = false, "simulate sequences and state paths");
  opts.add ("ml vit -viterbi", movement.viterbi = false, "report Viterbi likelihood and ML state path");
  opts.add ("f -forward", "--nforward 0", "report Forward likelihood", false);
  opts.add ("s -sample", "--nforward 1", "report Forward likelihood; sample one state path");
  opts.add ("nf -nforward", movement.nforward = -1, "report Forward likelihood; sample multiple state paths", false);
  opts.add ("e -expect", movement.want_expected_counts = false, "report expected counts from Forward-Backward & peeling");
  opts.add ("o -optacc", movement.optacc = false, "report optimal accuracy alignment path");
  opts.add ("np -normprof", movement.normalize_peeled_profiles = true, "normalize bit-profiles in peeled composition");
  opts.add ("rsp -redsuch-propose", movement.propose_redelings_suchard_move = false, "propose & compute Hastings terms for Redelings-Suchard MCMC move");
  opts.add ("rsi -redsuch-inverse", movement.evaluate_redelings_suchard_inverse_move = false, "compute Hastings terms for Redelings-Suchard inverse MCMC move");
  opts.add ("rs -redsuch", "--redsuch-propose --redsuch-inverse", "propose & compute Hastings ratio for Redelings-Suchard MCMC move");

  // HMMoC options

  movement.hmmoc_opts.init_opts_list (opts);

  // notes

  opts.newline();
  opts.newline();
  opts.print ("A detailed users' guide is available at the following URL:\n");
  opts.print ("http://biowiki.org/PhyloComposer\n");

  // main program

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
      // read in the Transducer_funcs and the ETree
      movement.read_composition (opts.args);

      // run
      movement.dump_composition (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}
