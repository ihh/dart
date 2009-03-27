#include "empath/p53.h"
#include "util/vector_output.h"

int main (int argc, char** argv)
{
  Opts_list opts (argc, argv);
  opts.short_description = "P53 binding site searcher";
  opts.syntax = "[options] <training set> <database>";

  opts.newline();
  Log_stream::add_opts (opts);

  int unit_size;
  int max_spacer;
  double min_bits;
  double pseud_mul;

  opts.newline();
  opts.add ("unit",    unit_size = 5,      "length of quarter-unit");
  opts.add ("space",   max_spacer = 10,    "max distance between units");
  opts.add ("minbits", min_bits = 0,       "score threshold (in bits) for matches");
  opts.add ("mpseud",  pseud_mul = .1,     "multiplier to go from null model probabilities to pseudocounts");

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
      sstring training_file = opts.args[0];
      sstring search_file = opts.args[1];

      // read in sequences
      FASTA_sequence_database training_db (training_file.c_str());
      FASTA_sequence_database search_db (search_file.c_str());

      const Alphabet& alphabet = DNA_alphabet;

      // create P53_model and Local_trainer
      PScore pscore;
      P53_model model (unit_size, max_spacer, pscore);
      Local_trainer trainer (model, 0);  // no masking

      // set up model params
      model.optimise_null_model (search_db, pseud_mul);

      // train
      trainer.global_train (training_db);

      // search
      GFF_list hits;
      trainer.set_bit_threshold (min_bits);
      trainer.local_search (search_db, hits);

      // output
      cout << hits;
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
