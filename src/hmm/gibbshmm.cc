#include <stdlib.h>
#include "util/opts_list.h"
#include "hmm/singlealign.h"
#include "util/logfile.h"

// Random motif class
//
struct Random_motif : Single_HMM_scores
{
  int motif_len;

  Random_motif (int motif_len, int pad_len, double motif_specificity, const Alphabet* alphabet = 0)
    : Single_HMM_scores (motif_len + 2, alphabet), motif_len (motif_len)
    {
      int pad_extend_score = Prob2Score (1 - 1.0 / (double) (pad_len + 1));
      int pad_end_score    = Prob2Score (1.0 / (double) (pad_len + 1));
      for (int s = 0; s < motif_len + 2; ++s)
	{
	  start[s] = end[s] = -InfinityScore;
	  for (int d = 0; d < motif_len + 2; ++d) transition(s,d) = -InfinityScore;
	  if (s == 0) { start[s] = 0; transition(s,s) = pad_extend_score; transition(s,s+1) = pad_end_score; }
	  else if (s == motif_len+1) { transition(s,s) = pad_extend_score; end[s] = pad_end_score; }
	  else transition(s,s+1) = 0;
	  if (s == 0 || s == motif_len+1)
	    emit[s] = vector<int> (4, Prob2Score (0.25));
	  else
	    {
	      vector<double> dirichlet_prior (4, 1.0);
	      dirichlet_prior [Rnd::rnd_int(4)] = 1.0 / motif_specificity;                       // ??? is this really how Dirichlets are meant to work ?!?
	      emit[s] = Prob2Score_vector (Rnd::sample_dirichlet (dirichlet_prior));
	    }
	}
    }

  void show_alignment (ostream& o, const Single_HMM_alignment& align)
    {
      vector<int>     state_sequence (motif_len);
      sstring         label_sequence (motif_len);
      Local_alignment local;
      sstring         column_labels;

      for (int s = 0; s < motif_len; ++s)
	{
	  state_sequence[s] = s + 1;
	  label_sequence[s] = s % 35 < 9 ? s % 35 + '1' : s % 35 + 'A' - 9;
	}
      align.make_local_alignment (state_sequence, label_sequence, local, column_labels);
      local.write_MUL (o, *alphabet, &column_labels);
    }
};

int main(int argc, char* argv[])
{
  Opts_list opts (argc, argv, CLOGERR);
  opts.short_description = "simulated Gibbs sampling test";
  opts.syntax = "[options]";
  Rnd::add_opts (opts);
  opts.newline();
  Log_stream::add_opts (opts);

  int    seqs;
  int    motif_len;
  double motif_specificity;
  int    pad_len;
  int    cycles;
  int    iterations;
  double kT_start;
  double kT_end;
  bool   deja_vu;

  opts.newline();
  opts.add ("seqs",       seqs = 100,                  "\tnumber of sequences");
  opts.add ("mlen",       motif_len = 6,               "\tmotif length");
  opts.add ("mspec",      motif_specificity = 100,     "\tmotif specificity: principal Dirichlet pseudocount for emission vectors");
  opts.add ("padlen",     pad_len = 200,               "\tpadding length");
  opts.add ("cycles",     cycles = 10,                 "\tnumber of re-seeding cycles");
  opts.add ("iter",       iterations = 10,             "\tnumber of iterations (annealing/sampling steps per sequence) per re-seeding cycle");
  opts.add ("ktstart",    kT_start = 1,                "initial value of kT for annealing");
  opts.add ("ktend",      kT_end = 1,                  "\tfinal value of kT for annealing");
  opts.add ("dejavu",     deja_vu = 1,                 "stop if an already-visited local maximum is found");

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
      // create the original motif

      const Alphabet&   alphabet = DNA_alphabet;
      Random_motif      motif (motif_len, pad_len, motif_specificity, &alphabet);
      
      if (CLOGGING(7)) { CL << "Original motif:\n"; motif.show(CL); }

      // create the original alignment

      Single_HMM_alignment true_align (motif);
      Sequence_database    db;
      for (int r = 0; r < seqs; ++r)
	{
	  sstring name;
	  name << "Row#" << r;
	  true_align.emit_row (name, db);
	}

      if (CLOGGING(7)) { CL << "Original alignment:\n"; motif.show_alignment (CL, true_align); }

      // set up other persistent structures

      Single_HMM_mask mask (motif, 1);
      mask.clear_all_transitions();
      mask.clear_all_start();
      mask.clear_all_end();
      mask.clear_emit (0);
      mask.clear_emit (motif_len+1);
      
      Single_HMM_counts pseudocounts (motif);
      for (int s = 0; s < motif.states(); ++s)
	for (int i = 0; i < alphabet.size(); ++i)
	  pseudocounts.emit[s][i] = 1;
      
      // start the re-seeding cycles
      
      Random_motif          best_hmm (motif_len, pad_len, 1.0, &alphabet);
      Single_HMM_alignment  best_align (best_hmm);
      int                   best_score = -InfinityScore;
      int                   best_cycle = 0;

      for (int cycle = 1; cycle <= cycles; ++cycle)
	{
	  // set up the test HMM & alignment
	  
	  Random_motif          hmm (motif_len, pad_len, 1.0, &alphabet);
	  Single_HMM_alignment  align (hmm);
	  
	  for (int r = 0; r < seqs; ++r)
	    align.add_row (true_align.row_name[r], true_align.prof[r]);
	  
	  Single_HMM_alignment  cycle_best_align (hmm);
	  int                   cycle_best_score = -InfinityScore;

	  // go into a Gibbs sampling frenzy
	  
	  for (int iter = 0; iter < iterations; ++iter)
	    {
	      double kT = kT_start + (kT_end - kT_start) * ((double) iter) / (double) iterations;
	      vector<int> schedule (seqs);
	      for (int i = 0; i < seqs; ++i) schedule[i] = i;
	      for (int i = 0; i < seqs; ++i)
		{
		  int k = Rnd::rnd_int(seqs-i) + i;
		  swap (schedule[i], schedule[k]);
		  int row = schedule[i];

		  if (CLOGGING(6))
		    {
		      CL << "Re-seeding cycle #" << cycle << ", iteration #" << iter+1 << ", sample #" << i+1 << " (row #" << row << ")";
		      if (iter > 0 || i > 0)
			{
			  CL << ", ";
			  if (iter == 0 && i < seqs-1) CL << "partial ";
			  CL << "score = " << Score2Bits (align.total_path_score()) << " bits";
			}
		      if (best_score > -InfinityScore) CL << ", best score = " << Score2Bits (best_score) << " bits";
		      if (cycle_best_score > -InfinityScore) CL << ", best this cycle = " << Score2Bits (cycle_best_score) << " bits";
		      CL << "\n";
		    }

		  align.Gibbs_sample_row (row, mask, pseudocounts, 1);
		  if (CLOGGING(4)) { CL << "New alignment:\n"; hmm.show_alignment (CL, align); }

		  bool got_all_rows = iter > 0 || i == seqs-1;
		  if (got_all_rows)
		    {
		      int score = align.total_path_score();
		      CLOG(5) << "Alignment score = " << Score2Bits (score);
		      if (score > cycle_best_score) { cycle_best_align = align; cycle_best_score = score; CL << " ***** best score this cycle *****\n"; }
		      CL << "\n";
		    }
		  
		  if (CLOGGING(5)) { CL << "New motif:\n"; hmm.show(CL); }
		}
	    }
	  
	  // refine alignment

	  CLOG(6) << "Refining best alignment\n";
	  if (CLOGGING(5)) hmm.show_alignment (CL, cycle_best_align);
	  
	  int row = 0;
	  int rows_since_last_update = 0;
	  while (rows_since_last_update < seqs)
	    if (cycle_best_align.Gibbs_optimise_row (row = (row + 1) % seqs, mask, pseudocounts, 1))
	      {
		CLOG(6) << "Improved row #" << row << "; new score = " << Score2Bits (cycle_best_align.total_path_score()) << " bits\n";
		rows_since_last_update = 0;
	      }
	    else
	      {
		CLOG(6) << "No improvement for row #" << row << "\n";
		++rows_since_last_update;
	      }
	  CLOG(6) << "No more greedy improvements to be found; stopping refinement.\n";
	  if (CLOGGING(5)) hmm.show_alignment (CL, cycle_best_align);
	  
	  cycle_best_align.optimise_HMM (mask, pseudocounts);
	  cycle_best_score = cycle_best_align.total_path_score();

	  if (deja_vu && cycle > 0)
	    if (cycle_best_align == best_align)
	      {
		CLOG(6) << "Deja vu! This alignment has been seen before; stopping.\n";
		break;
	      }

	  if (cycle_best_score > best_score) { best_hmm = hmm; best_align = cycle_best_align; best_score = cycle_best_score; CLOG(6) << "New best alignment\n"; }
	}
      
      CLOGOUT << "Best alignment:\n";
      best_hmm.show_alignment (CL, best_align);
      CL << "Best alignment score = " << Score2Bits(best_score) << ", accuracy = " << best_align.accuracy (true_align) << "\n";
      CL << "Best motif:\n";
      best_hmm.show(CL);
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }
  
  return 0;
}


