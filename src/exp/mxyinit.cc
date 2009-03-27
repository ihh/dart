#include "util/logfile.h"
#include "util/rnd.h"
#include "exp/matrix_array.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  Opts_list opts (argc, argv);
  opts.short_description = "randomly seed a hidden substitution matrix array with sequence, fixed clique & variable clique classes\n";
  opts.syntax = "[options]";
  Rnd::add_opts (opts);

  opts.newline();
  Log_stream::add_opts (opts);

  int C;
  int F;
  int S;
  double tres;
  double seq_dev;
  double fixed_dev;
  double prior_dev;
  double intra_min;
  double intra_dev;
  double inter_min;
  double inter_dev;

  opts.newline();
  opts.add ("c classes -classes", C = 1, "\t\tnumber of variable classes");
  opts.add ("f fixedclasses -fixedclasses", F = 1, "\tnumber of fixed classes");
  opts.add ("s seqclasses -seqclasses", S = 1, "\tnumber of sequence classes");
  opts.add ("tr tres -timeres", tres = .01, "\t\tfractional resolution of separation times");

  opts.newline();
  opts.add ("seqdev -seqdev",   seq_dev = .1,   "deviance of sequence class prior probability");
  opts.add ("fixdev -fixdev",   fixed_dev = .1,   "deviance of fixed residue class prior probability");
  opts.add ("priordev -priordev", prior_dev = .1, "deviance of variable residue class prior probability");
  opts.add ("intramin -intramin", intra_min = .05, "minimum intra-class substitution rate");
  opts.add ("intradev -intradev", intra_dev = .01, "deviance of intra-class substitution rate");
  opts.add ("intermin -intermin", inter_min = .005, "minimum inter-class substitution rate");
  opts.add ("interdev -interdev", inter_dev = .001, "deviance of inter-class substitution rate");

  // parse the command line
  //
  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 0) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << opts.short_help();
      CLOGERR << e.what();
      exit(1);
    }

  try
    {
      // initialise alphabet & database
      const Alphabet& alphabet = Protein_alphabet;
      Sequence_database dummy_seq_db;
      const Tree_alignment_database dummy_align_db (dummy_seq_db);

      // build the matrix
      const int max_fork = 1;
      Matrix_array matrix_array (S, F, C, alphabet.size(), max_fork, &dummy_align_db, tres);

      // randomise it
      matrix_array.randomise (seq_dev, fixed_dev, prior_dev, intra_min, intra_dev, inter_min, inter_dev);

      // output it
      matrix_array.write (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
