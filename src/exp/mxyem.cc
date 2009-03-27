#include "util/logfile.h"
#include "util/rnd.h"
#include "exp/matrix_array.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  // initialise the options parser
  //
  Opts_list opts (argc, argv);
  opts.short_description = "train a hidden substitution matrix array on an alignment database using EM\n";
  opts.syntax = "[options] <initial MXY file> <index file>";
  Rnd::add_opts (opts);

  opts.newline();
  Log_stream::add_opts (opts);

  int C;
  int F;
  int S;
  double tres;
  int forgive;
  bool rind;
  bool intra;
  bool inter;
  int max_fork;

  opts.newline();
  opts.add ("c classes -classes", C = 1, "\t\tnumber of classes");
  opts.add ("f fixedclasses -fixedclasses", F = 1, "\tnumber of fixed classes");
  opts.add ("s seqclasses -seqclasses", S = 1, "\tnumber of sequence classes");
  opts.add ("r rind -rind", rind = 0, "\t\t\tuse RIND-constrained M-step");
  opts.add ("intra -intra", intra = 1, "\t\t\ttrain inter-class substitution rates");
  opts.add ("inter -inter", inter = 1, "\t\t\ttrain inter-class substitution rates");
  opts.add ("tr tres -timeres", tres = .01, "\t\tfractional resolution of separation times");
  opts.add ("fg forgive -forgive", forgive = 20, "\t\tnumber of bad EM rounds to forgive");
  opts.add ("fk fork -fork", max_fork = 2, "\t\tnumber of processes to fork");

  // parse the command line
  //
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
      // get args
      const char* initial_mxy_filename = opts.args[0].c_str();
      const char* index_filename = opts.args[1].c_str();

      // initialise alphabet & database
      const Alphabet& base_alphabet = Protein_alphabet;
      Sequence_database seq_db;
      Tree_alignment_database align_db (seq_db, index_filename);

      // initialise array
      Matrix_array matrix_array (1, 1, 1, base_alphabet.size(), max_fork, &align_db, tres);
      
      // read in array
      ifstream mxy_file (initial_mxy_filename);
      if (!mxy_file) THROWEXPR ("MXY file '" << initial_mxy_filename << "' not found");
      matrix_array.read (mxy_file);

      // digitise sequences
      matrix_array.init_alphabet (base_alphabet);
      const Alphabet& hidden_alphabet = matrix_array.alphabet();
      seq_db.seqs2scores (hidden_alphabet);

      // train matrix
      const Loge log_likelihood = matrix_array.iterate_quick_EM (rind, intra, inter, forgive);

      // show log-likelihood
      CLOG(7) << "Final log-likelihood = " << log_likelihood << "\n";

      // output matrix
      matrix_array.write (cout);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
