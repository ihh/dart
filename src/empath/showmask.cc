#include "seq/biosequence.h"
#include "seq/gff.h"

int main (int argc, char** argv)
{
  Opts_list opts (argc, argv);
  opts.short_description = "show a probabilistic mask for a sequence database";
  opts.syntax = "[options] <sequence file> <mask file>";

  opts.newline();
  Log_stream::add_opts (opts);

  Prob mask_pr;
  Prob show_pr;

  opts.newline();
  opts.add ("m -mask",     mask_pr = 0,   "probability threshold below which residues are masked");
  opts.add ("s -show",     show_pr = .9,  "probability threshold below 10% probability increments are indicated");

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
      sstring seq_file = opts.args[0];
      sstring mask_file = opts.args[1];

      // read in sequences
      FASTA_sequence_database seq_db (seq_file.c_str(), 0, Profile_flags_enum::NONE);
      const Alphabet& alphabet = seq_db.alphabet();

      // create & load mask
      seq_db.reset_metascores(1);
      ifstream mask_stream (mask_file.c_str());
      GFF_list gff_list;
      mask_stream >> gff_list;
      gff_list.apply_mask (seq_db.index, 0, 1);

      // output seq_db
      const Score mask_sc = Prob2Score (mask_pr);
      const char mask_char = alphabet.unknown_char();
      const int columns = 50;  // output width
      for_const_contents (Sequence_database, seq_db, np)
	{
	  cout << ">" << np->name << " " << np->cruft << "\n";
	  for (int pos = 0; pos < np->size(); pos += columns)
	    {
	      for (int i = pos; i < min (np->size(), pos + columns); ++i)
		cout << (np->meta_sc[0][i] < mask_sc ? mask_char : np->seq[i]);
	      cout << "\n";
	      if (show_pr > 0)
		{
		  for (int i = pos; i < min (np->size(), pos + columns); ++i)
		    {
		      const Prob p = Score2Prob (np->meta_sc[0][i]);
		      cout << (char) (p <= show_pr ? ('0' + min ((int) (p * 10), 9)) : ' ');
		    }
		  cout << "\n";
		}
	    }
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
