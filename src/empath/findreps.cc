#include "util/vector_output.h"
#include "empath/trirep.h"
#include "empath/trainer.h"

#define DEFAULT_MASK   "REPMASK"
#define DEFAULT_MODEL  "REPMODEL"
#define REPEAT_FEATURE "Microsatellite"

int main (int argc, char** argv)
{
  Opts_list opts (argc, argv);
  opts.short_description = "exhaustive EM-based short repeat finder";
  opts.syntax = "[options] <sequence file>";

  opts.newline();
  opts.print ("Logfile options\n");
  opts.print ("---------------\n");
  Log_stream::add_opts (opts);

  int min_reps;
  int shape;
  int rep_len;
  int invrep_len;
  bool low_complexity;
  double mismatch_prob;
  double insert_prob;
  double delete_prob;
  // double bit_threshold;
  bool train;
  int forgive;
  bool mask;
  Prob report;
  sstring inmask;
  sstring outmask;
  sstring inmodel;
  sstring outmodel;
  sstring outdist;

  opts.newline();
  opts.print ("Training options\n");
  opts.print ("----------------\n");
  opts.add ("t -train",     train = 0,                 "\toptimise repeat lengths & frequencies by EM");
  opts.add ("fg -forgive",   forgive = 3,               "number of bad EM rounds to forgive");

  opts.newline();
  opts.print ("Masking options\n");
  opts.print ("---------------\n");
  opts.add ("m -mask",      mask = 1,                  "\tmask sequences");
  opts.add ("r -report",    report = .01,              "posterior probability threshold for reported hits");

  opts.newline();
  opts.print ("Modeling options\n");
  opts.print ("----------------\n");
  opts.add ("rl -replen",    rep_len = 1,               "\tmax length of repeat motifs");
  opts.add ("il -invreplen", invrep_len = 0,            "\tmax length of inverted repeat motifs");
  opts.add ("lc -lowcomplexity",   low_complexity = 1,        "\tinclude 'low-complexity' regions in model");
  opts.add ("mr -minreps",   min_reps = 5,              "\tmin number of copies of repeat");
  opts.add ("s -shape",     shape = 5,                 "\tshape of copy distribution, i.e. number of states used to model");
  opts.add ("mp -mismatchprob",       mismatch_prob = .01,       "mismatch probability");
  opts.add ("ip -insertprob",       insert_prob = .1,          "insert probability");
  opts.add ("dp -deleteprob",       delete_prob = .01,         "delete probability");
  //  opts.add ("hitbits",   bit_threshold = 0,         "hit threshold in bits");

  opts.newline();
  opts.print ("Loading/saving\n");
  opts.print ("--------------\n");
  opts.add ("lk -loadmask",  inmask = "",               "load initial mask from file", 0);
  opts.add ("sk -savemask",  outmask = DEFAULT_MASK,    "save final mask to file");
  opts.add ("lm -loadmodel", inmodel = "",              "load initial model scores from file", 0);
  opts.add ("sm -savemodel", outmodel = DEFAULT_MODEL,  "save final model scores to file");
  opts.add ("sd -savedist",  outdist = "",              "save repeat length distribution to file", 0);

  try
    {
      if (!opts.parse()) { CLOGERR << opts.short_help(); exit(1); }
      if (opts.args.size() != 1) { CLOGERR << opts.short_help(); CLOGERR << "Wrong number of arguments\n\n"; exit(1); }
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

      // read in sequences
      FASTA_sequence_database seq_db (seq_file.c_str(), 0, Profile_flags_enum::DSQ);
      const Alphabet& alphabet = seq_db.alphabet();

      // initialise mask
      seq_db.reset_metascores(1);
      if (inmask != "")
	{
	  ifstream mask_stream (inmask.c_str());
	  GFF_list gff_list;
	  mask_stream >> gff_list;
	  gff_list.apply_mask (seq_db.index, 0, 1);
	}

      // create Trirep_model
      Trirep_set trirep_set (rep_len, invrep_len, alphabet);
      if (low_complexity) trirep_set.add_degenerates();
      PScores pscore;
      Trirep_model trirep_model (trirep_set, pscore, min_reps, shape, mismatch_prob, insert_prob, delete_prob);
      Single_fast_matrix_factory dp_factory;

      // initialise null model & priors
      const double pseud_mul = 1.0;
      trirep_model.optimise_null_model (seq_db, pseud_mul);

      // create Local_trainer
      Local_trainer local_trainer (trirep_model, dp_factory);
      
      // load model, or use default seed
      if (inmodel == "")
	local_trainer.prior.initialise();
      else
	{
	  ifstream inmodel_stream (inmodel.c_str());
	  trirep_model.pscore.read (inmodel_stream);
	}

      // training
      if (train)
	{
	  // do EM until convergence
	  PScores best_pscore = trirep_model.pscore;
	  Loge best_loglike = -InfinityLoge;
	  int dec_count = 0;
	  for (int em_iter = 1; 1; ++em_iter)  // EM outer loop
	    {
	      const PScores old_pscore = trirep_model.pscore;
	      const Loge loglike = local_trainer.local_EM (seq_db);  // do one round of EM
	      CLOG(6) << "EM iteration #" << em_iter << ": score " << Nats2Bits(loglike) << " bits\n";
	      if (loglike > best_loglike)
		{
		  best_loglike = loglike;  // record hiscores
		  best_pscore = old_pscore;
		  dec_count = 0;
		  // save model
		  if (outmodel != "")
		    {
		      ofstream outmodel_stream (outmodel.c_str());
		      best_pscore.write (outmodel_stream);
		    }
		}
	      else if (++dec_count > forgive)
	        {
		  CLOG(6) << "Likelihood decreased for " << dec_count << "th time; stopping EM";
		  break;
		}
	    }
	  local_trainer.model.pscore = best_pscore;
	}

      if (outdist != "")
	{
	  ofstream outdist_stream (outdist.c_str());
	  trirep_model.write_repdist (outdist_stream);
	}

      // mask repeats, save mask
      if (mask)
	{
	  GFF_list gff_list;
	  local_trainer.get_local_mask (seq_db, gff_list, report, opts.program_name.c_str(), REPEAT_FEATURE);
	  ofstream mask_stream (outmask.c_str());
	  mask_stream << gff_list;
	}
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
