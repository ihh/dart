#include "util/vector_output.h"
#include "empath/bindpaircfg.h"
#include "scfg/paircfgdp.h"

int main (int argc, char** argv)
{
  Opts_list opts (argc, argv);
  opts.short_description = "search two RNA datasets for complementarity";
  opts.syntax = "[options] <query> <target>";
  opts.expect_args = 2;

  opts.newline();
  Log_stream::add_opts (opts);

  int max_fork;
  double min_energy;

  opts.newline();
  opts.add ("f -fork", max_fork = 2, "\tnumber of processes to fork");
  opts.add ("e -energy", min_energy = 11.0, "\tminimum free energy increase");

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
      sstring query_file = opts.args[0];
      sstring target_file = opts.args[1];

      // read in sequences
      FASTA_sequence_database query_db (query_file.c_str(), 0, Profile_flags_enum::DSQ);
      FASTA_sequence_database target_db (target_file.c_str(), 0, Profile_flags_enum::DSQ);
      const Alphabet& alphabet = RNA_alphabet;
      const int alph_sz = alphabet.size();
      if (query_db.alphabet().size() != alph_sz || target_db.alphabet().size() != alph_sz)
	THROWEXPR ("Expected RNA sequences");

      // create Pair_CFG
      Energy21 energy;
      Bind_pair_CFG cfg (energy);

      // cycle through query database
      for (int nq = 0; nq < query_db.size(); ++nq)
	{
	  const Named_profile& query_np = *query_db.index.profile[nq];
	  const Digitized_biosequence& query_dsq = query_np.dsq;
	  Fold_envelope query_env;
	  query_env.initialise_3prime (query_dsq.size());
	  // cycle through target database
	  for (int nt = 0; nt < target_db.size(); ++nt)
	    {
	      Named_profile& target_np = *target_db.index.profile[nt];
	      Digitized_biosequence& target_dsq = target_np.dsq;
	      Fold_envelope target_env;
	      target_env.initialise_5prime (target_dsq.size());
	      // save old DSQ
	      const Digitized_biosequence orig_target_dsq = target_dsq;
	      // search using CFG
	      while (1)
		{
		  Pair_CYK_matrix cyk_matrix (query_dsq, target_dsq, query_env, target_env, cfg);
		  const vector<int> state_path = cyk_matrix.traceback();
		  const Pair_CFG_parse_tree parse_tree = cfg.parse (state_path);
		  const Pair_CFG_scores::XYMask xymask = cfg.get_mask (parse_tree, query_dsq, target_dsq, cfg.mask_states);
		  GFF gff = cfg.ymask2gff (xymask, query_np, target_np);
		  gff.source = opts.program_name;
		  gff.score = Score2Bits (cfg.path_score (parse_tree, query_dsq, target_dsq));
		  gff.strand = GFF::Minus;
		  if (gff.score < min_energy) break;
		  cout << gff;
		  // hacky randomise-masking
		  for (int y = gff.start-1; y <= gff.end-1; ++y)
		    target_dsq[y] = Rnd::rnd_int (alph_sz);
		}
	      // restore old DSQ
	      target_dsq = orig_target_dsq;
	    }
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}

