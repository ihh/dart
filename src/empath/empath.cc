#include "util/vector_output.h"
#include "empath/ungapped.h"
#include "empath/trainer.h"
#include "empath/multihit.h"
#include "empath/repeater.h"
#include "empath/shortlist.h"
#include "seq/suffix.h"
#include "empath/exprhit.h"

#define DEFAULT_MASK  "EMPMASK"
#define MOTIF_FEATURE "MOTIF"

int main (int argc, char** argv)
{
  Opts_list opts (argc, argv);
  opts.short_description = "EM patterning heuristic, simple ungapped motif";
  opts.syntax = "[options] <sequence file>";

  opts.newline();
  opts.print ("Logfile options\n");
  opts.print ("---------------\n");
  Log_stream::add_opts (opts);

  int n_motifs;
  int min_motif_len;
  int max_motif_len;
  bool revcomp;
  int min_copies;
  int max_copies;
  int min_seqs;
  int max_seqs;
  int min_word_reps;
  int repeated_word_len;
  int selfsim_shortlist_len;
  int local_shortlist_len;
  int expr_neighbors;
  int expr_shortlist_len;
  int null_order;
  sstring null_file;
  double pseud_mul;
  sstring expr_file;
  double alpha;
  double lambda;
  double min_inc;
  int forgive;
  int max_fork;

  sstring inmask;
  sstring outmask;
  sstring motif_file_prefix;
  sstring gff_file_prefix;

  opts.newline();
  opts.print ("Sequence modeling options\n");
  opts.print ("-------------------------\n");
  opts.add ("n -nmotifs",  n_motifs = 3,       "\ttotal number of motifs to find");
  opts.add ("minl -minlen",     min_motif_len = 5,  "\t} motif length { minimum");
  opts.add ("maxl -maxlen",     max_motif_len = 12, "\t}              { maximum");
  opts.add ("rc -revcomp", revcomp = 1,        "\t\tallow reverse-complemented hits to motif");
  opts.add ("minc -mincopies",     min_copies = 1,     "} number of copies of   { minimum");
  opts.add ("maxc -maxcopies",    max_copies = 5,      "}   motif per sequence  { maximum");
  opts.add ("mins -minseqs",   min_seqs = 1,       "\t} number of matching     { minimum");
  opts.add ("maxs -maxseqs",   max_seqs = 20,      "\t}  sequences in database { maximum");
  opts.add ("no -nullorder", null_order = 6, "\torder of null model");
  opts.add ("nf -nullfile", null_file = "", "\ttrain null model from separate sequence file, instead of target dataset", FALSE);
  opts.add ("pw -pseudoweight",     pseud_mul = .1,     "multiplier to go from residue frequencies to pseudocounts");

  opts.newline();
  opts.print ("Expression modeling options\n");
  opts.print ("---------------------------\n");
  opts.add ("x expr -exprfile",       expr_file = "",  "read expression data from file", FALSE);
  opts.add ("efh -exprformathelp",      &Expr_table::format_help, "\t\tdescribe expression file format");
  opts.add ("al -alpha",      alpha = .5,   "\tnormal hyperparameter for variance; higher for stronger priors");
  opts.add ("la -lambda",     lambda = 1,   "\tgamma hyperparameter for mean; higher for stronger priors");

  opts.newline();
  opts.print ("Masking options\n");
  opts.print ("---------------\n");
  opts.add ("lm -loadmask",  inmask = "",               "\tload initial mask from file", FALSE);
  opts.add ("sm -savemask",  outmask = DEFAULT_MASK,    "\tsave final mask to file");

  opts.newline();
  opts.print ("Seed drill options\n");
  opts.print ("------------------\n");
  opts.add ("swr -seedwordreps", min_word_reps = 1, "minimum number of repetitions of a word required to seed on it");
  opts.add ("swl -seedwordlen", repeated_word_len = 5, "length of repeated seed words");
  opts.add ("ssl -selfsimlen", selfsim_shortlist_len = 100, "\tnumber of seeds to keep from self-similarity shortlist; 0 to unlimit");
  opts.add ("lml -localmodellen", local_shortlist_len = 10,   "number of seeds to keep from local-modeling shortlist; 0 to unlimit");
  opts.add ("nx -nearexp", expr_neighbors = 100,   "\tnumber of sequences with similar expression profiles to use for database-nearest shortlist");
  opts.add ("nxl -nearexplen", expr_shortlist_len = 0,   "\tnumber of seeds to keep from database-nearest shortlist");

  opts.newline();
  opts.print ("Output options\n");
  opts.print ("--------------\n");
  opts.add ("swm -saveweightmatrix",  motif_file_prefix = "",     "save motifs to files with this prefix", FALSE);
  opts.add ("sgff -savegff",    gff_file_prefix = "",       "\tsave GFF hits to files with this prefix", FALSE);

  opts.newline();
  opts.print ("Miscellaneous options\n");
  opts.print ("---------------------\n");
  opts.add ("mi -mininc",    min_inc = .001,     "\tminimum fractional log-likelihood increase to continue EM");
  opts.add ("fg -forgive",   forgive = 0,        "\tnumber of bad EM rounds to forgive");
  opts.add ("fk -fork",      max_fork = 2,       "\tnumber of processes to fork");

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
      // get filename
      sstring seq_file = opts.args[0];

      // read in sequences
      FASTA_sequence_database seq_db (seq_file.c_str(), 0, Profile_flags_enum::DSQ);
      if (min_seqs > seq_db.size()) THROWEXPR ("--minseqs parameter is greater than size of database");
      if (min_seqs > max_seqs) THROWEXPR ("--minseqs parameter is greater than --maxseqs");
      const Alphabet& alphabet = seq_db.alphabet();

      // metascores: M0=M1-M2 is the mask the HMM uses, M1 is the actual mask, M2 is the null model
      seq_db.reset_metascores (3);

      // initialise mask metascores
      if (inmask != "")
	{
	  ifstream mask_stream (inmask.c_str());
	  GFF_list gff_list;
	  mask_stream >> gff_list;
	  gff_list.apply_mask (seq_db.index, 1, TRUE);
	}
      list<GFF_list> motif_masks;

      // initialise null model metascores
      Suffix_tree suff_tree (FALSE);
      if (null_file.size())
	{
	  CTAG(6,EMPATH) << "Loading null model database\n";
	  FASTA_sequence_database null_db (null_file.c_str(), 0, Profile_flags_enum::DSQ);
	  suff_tree.build (null_db.index, null_order + 1, null_order + 1, 0, 1);
	}
      else
	suff_tree.build (seq_db.index, null_order + 1, null_order + 1, 0, 1);
      suff_tree.set_null_metascores (seq_db.index, 2, null_order, alphabet.size());

      // create motif
      Single_fast_matrix_factory dp_factory;
      PScores pscore;
      Ungapped_model ungapped (min_motif_len, max_motif_len, revcomp, pscore, TRUE);

      // get residue frequencies (i.e. zeroth-order null model)
      ungapped.optimise_null_model (seq_db, pseud_mul);

      // create repeater
      Repeater repeater (ungapped, min_copies, max_copies);
      
      // create Local_trainer's
      Local_trainer local_trainer (ungapped, dp_factory);
      Local_trainer repeat_trainer (repeater, dp_factory);
      local_trainer.local_hmm.allow_multiple_hits();  // allow multiple hits for local_trainer but not for repeat_trainer

      // read in expression data, if specified
      Expr_table expr_tab;
      if (expr_file.size()) expr_tab.read (expr_file.c_str());

      // create Normal-Gamma model
      Expr_table::Stats expr_stats = expr_tab.experiment_stats();
      Normal_Gamma_expr_model ng_expr (expr_stats, alpha, lambda);
      
      // create Hitters
      Sequence_hitter seq_hitter (repeat_trainer);
      Expr_hitter expr_hitter (seq_hitter, ng_expr, expr_tab);
      Hitter* hitter = 0;
      if (expr_file.size())
	hitter = &expr_hitter;
      else
	hitter = &seq_hitter;

      // create multihit stuff
      Range_transitions range_transitions (min_seqs, max_seqs);
      Multihitter multihit_model (range_transitions, *hitter, pscore, max_fork, 1.0);

      // set up the Seed_drill
      Seed_drill drill;

      // create self-similarity comparison shortlist
      const Prob selfsim_match_prob = .9;
      drill.add_fixed_length_shortlist (new Self_similarity_shortlist (seq_db, max_fork, alphabet, min_motif_len, (max_motif_len - min_motif_len) / 2, selfsim_match_prob, revcomp, ungapped.mask_metascore_idx), selfsim_shortlist_len);

      // create local-modeling shortlist
      drill.add_fixed_length_shortlist (new Local_shortlist (seq_db, max_fork, repeat_trainer), local_shortlist_len);

      // create expression shortlist, if expression data are available
      if (expr_file.size())
	{
	  Expr_dist_func expr_dist_func (expr_tab, expr_stats);
	  drill.add_fixed_length_shortlist (new Multihit_shortlist (seq_db, multihit_model, expr_dist_func, expr_neighbors), expr_shortlist_len);
	}

      // create global-modeling shortlist
      drill.add_fixed_length_shortlist (new Multihit_shortlist (seq_db, multihit_model), 0);

      // seed the drill
/*
      const int CG12758 = seq_db.index.name2profile_index.lookup (sstring ("CG12758"));
      const int CG14477 = seq_db.index.name2profile_index.lookup (sstring ("CG14477"));
      drill.front().seed.push_back (Seq_pos (CG12758, 1151));
      drill.front().seed.push_back (Seq_pos (CG14477, 39));
*/
      if (min_word_reps > 1)
	drill.front().seed_locally_repeated (max_motif_len, repeated_word_len, min_word_reps);
      else
	drill.front().seed_everywhere (max_motif_len);
      
      // find n_motifs motifs
      for (int n_motif = 0; n_motif < n_motifs; ++n_motif)
	{
	  CTAG(7,EMPATH) << "Preparing to estimate motif #" << n_motif+1 << "\n";

	  // update metascores
	  for_contents (Sequence_database, seq_db, np)
	    {
	      const int seqlen = np->size();
	      Metascore& effective_mask (np->meta_sc[0]);
	      const Metascore& real_mask (np->meta_sc[1]);
	      const Metascore& null_mask (np->meta_sc[2]);
	      for (int pos = 0; pos < seqlen; ++pos)
		effective_mask[pos] = ScorePMul (real_mask[pos], -null_mask[pos]);
	    }

	  // run the drill and get the best seed
	  drill.refresh();
	  const Seq_pos& best_seed = drill.best_seed();
	  if (CTAGGING(6,EMPATH_SEED)) CL << "Best seed is at sequence '" << seq_db.index.name[best_seed.seq] << "' position " << best_seed.pos << "\n";
	  
	  // do EM until convergence
	  repeat_trainer.local_seed (*seq_db.index.profile[best_seed.seq], best_seed.pos);   // seed at the best spot
	  multihit_model.iterate_multihit_EM (seq_db.index.profile, min_inc, forgive);
	  
	  // display the results
	  if (CTAGGING(6,EMPATH_VARS)) repeat_trainer.display_vars (CL);
	  sstring seqlogo = local_trainer.seqlogo();  // create seqlogo text first, because seqlogo() call corrupts logging flag
	  if (CTAGGING(4,EMPATH_SEQLOGO)) CL << seqlogo;
	  
	  // create a consensus string
	  Score_profile consensus_prof_sc = ungapped.prof_sc();
	  Biosequence consensus_seq;
	  alphabet.score2seq (consensus_prof_sc, consensus_seq);

	  // dump motif to file
	  if (motif_file_prefix != "")
	    {
	      sstring motif_file;
	      motif_file << motif_file_prefix << n_motif+1;
	      ofstream motif_stream (motif_file.c_str());
	      motif_stream << "# Consensus " << consensus_seq << "\n";
	      repeater.instruct_turtle (motif_stream);
	    }
	  
	  // search
	  sstring feature;
	  feature << MOTIF_FEATURE << n_motif+1 << '(' << consensus_seq << ')';
	  GFF_list gff_results;
	  local_trainer.local_search (seq_db, gff_results, opts.program_name.c_str(), feature.c_str());
	  if (gff_file_prefix != "")
	    {
	      sstring gff_file;
	      gff_file << gff_file_prefix << n_motif+1;
	      ofstream gff_stream (gff_file.c_str());
	      gff_stream << gff_results;
	    }
	  else
	    cout << gff_results;  // dump GFF to cout if nowhere else to report results
	  
	  // mask hits
	  motif_masks.push_back (GFF_list());
	  GFF_list& gff_mask = motif_masks.back();
	  const Prob mask_report_pr = .01;
	  local_trainer.get_local_mask (seq_db, gff_mask, mask_report_pr, opts.program_name.c_str(), feature.c_str());
	  gff_mask.apply_mask (seq_db.index, 1, TRUE);

	  // save masks
	  ofstream mask_stream (outmask.c_str());
	  for_const_contents (list<GFF_list>, motif_masks, gff_list)
	    mask_stream << *gff_list;
	}
      // delete shortlists
      drill.delete_shortlists();
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
}
