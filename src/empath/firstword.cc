#include "util/vector_output.h"
#include "empath/ungapped.h"
#include "empath/trainer.h"
#include "empath/multihit.h"
#include "empath/repeater.h"
#include "empath/shortlist.h"

#define DEFAULT_MASK  "EMPMASK"
#define MOTIF_FEATURE "MOTIF"

int main (int argc, char** argv)
{
  Opts_list opts (argc, argv);
  opts.short_description = "search initial words of an miRNA dataset against a query UTR dataset";
  opts.syntax = "[options] <miRNA> <UTR>";
  opts.expect_args = 2;

  opts.newline();
  Log_stream::add_opts (opts);

  int word_len;
  int max_fork;

  opts.newline();
  opts.add ("wl -wordlen", word_len = 8, "word length");
  opts.add ("fk -fork", max_fork = 2, "\tnumber of processes to fork");

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
      sstring miRNA_file = opts.args[0];
      sstring UTR_file = opts.args[1];

      // read in sequences
      FASTA_sequence_database miRNA_db (miRNA_file.c_str(), 0, Profile_flags_enum::DSQ);
      FASTA_sequence_database UTR_db (UTR_file.c_str(), 0, Profile_flags_enum::DSQ);
      const Alphabet& alphabet = RNA_alphabet;
      const int alph_sz = alphabet.size();
      if (miRNA_db.alphabet().size() != alph_sz || UTR_db.alphabet().size() != alph_sz)
	THROWEXPR ("Expected RNA sequences");

      // initialise mask
      UTR_db.reset_metascores(1);

      // create motif
      PScores pscore;
      Ungapped_model ungapped (word_len, word_len, 0, pscore);

      // initialise scores
      pscore[ungapped.null_emit] = vector<Score> (alph_sz, (Score) 0);
      pscore[ungapped.null_extend] = vector<Score> (2, (Score) 0);
      pscore[ungapped.skip_start[0]] = 0;
      pscore[ungapped.skip_end[0]] = 0;

      // create Local_trainer
      Single_fast_matrix_factory dp_factory;
      Local_trainer local_trainer (ungapped, dp_factory);
      local_trainer.local_hmm.allow_multiple_hits();

      // cycle through miRNA database
      for (int n_miRNA = 0; n_miRNA < miRNA_db.size(); ++n_miRNA)
	{
	  // set up HMM for initial word of miRNA
	  const Named_profile& np_miRNA = *miRNA_db.index.profile[n_miRNA];
	  if (np_miRNA.size() < word_len) continue;
	  for (int pos = 0; pos < word_len; ++pos)
	    {
	      vector<Score> match_sc (alph_sz, (Score) 0);
	      const int sym = alphabet.complement (np_miRNA.dsq[word_len - 1 - pos]);
	      match_sc[sym] = Prob2Score(2);  // set match score to 1 bit for now
	      pscore[ungapped.match_emit[pos]] = match_sc;
	    }
	  
	  // search using HMM
	  GFF_list gff_results;
	  local_trainer.local_search (UTR_db, gff_results, opts.program_name.c_str(), np_miRNA.name.c_str());
	  cout << gff_results;  // dump GFF to cout if nowhere else to report results
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}

