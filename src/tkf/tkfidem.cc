#include "tkf/tkfcoeff.h"
#include "tkf/tkfsequence.h"
#include "tkf/tkfopts.h"
#include "tkf/tkfdata.h"

#include "util/logfile.h"
#include "util/rnd.h"
#include "util/vector_output.h"

int main(int argc, char* argv[])
{
  TKF_opts opts (argc, argv, 1);
  opts.syntax = "[options] <index file>";
  opts.short_description = "estimate TKF model indel rates using EM";
  opts.expect_args = 1;

  double dt;
  bool rev;
  int forgive;
  bool naive;

  opts.newline();
  opts.print ("TKF EM parameters (Thorne et al, 1991)\n");
  opts.print ("--------------------------------------\n");
  opts.add ("d -dt", dt = .001, "time discretisation interval for numerical integration");
  opts.add ("r -rev", rev = TRUE, "count each branch in both forward & reverse directions");
  opts.add ("f -forgive", forgive = 0, "number of bad rounds of EM to forgive");
  opts.add ("naive", naive = FALSE, "use 'naive' E-step, i.e. simple indel-counting rather than probabilistic integration");

  opts.parse_or_die();
  try
    {
      // get args
      const char* index_filename = opts.args[0].c_str();

      // load database
      Sequence_database seq_db;
      Tree_alignment_database align_db (seq_db, index_filename);

      // digitise sequences
      const Alphabet& alph = seq_db.detect_alphabet();
      opts.use_pam = &alph == &Protein_alphabet;  // ensures that TKF_params.submat_factory has the correct alphabet
      seq_db.seqs2scores (alph);

      // find maximum branch length
      double max_len = 0;
      for_const_contents (list<Tree_alignment>, align_db.tree_align_list, tree_align)
	{
	  for_rooted_branches_post (tree_align->tree, b)
	    if ((*b).length > max_len)
	      max_len = (*b).length;
	}

      // get initial TKF_params
      TKF_params params = opts.params();

      // construct TKF_coeff
      TKF_coeff coeff (dt, max_len);

      // main loop
      int iter = 0;
      int dec = 0;
      Score best_sc = -InfinityScore;
      while (1)
	{
	  // increment EM iteration count
	  ++iter;
	  CLOG(7) << "Beginning EM iteration #" << iter << "\n";

	  // initialise TKF_coeff
	  coeff.initialise (params);

	  // accumulate counts
	  TKF_counts em_counts;
	  Score total_sc = 0;
	  for (int i = 0; i < align_db.size(); ++i)
	    {
	      // print name of alignment
	      CLOG(7) << "Alignment: " << align_db.name[i] << "\n";
	      const Tree_alignment& tree_align = *align_db.tree_align[i];
	      if (CLOGGING(3)) tree_align.align.write_MUL (CL, alph);

	      // create TKF_align object
	      TKF_align tkf (params);
	      tkf.set_tree (tree_align.tree);
	      tkf.set_alignment (tree_align.align);
	      tkf.build_maps_from_names();
	      tkf.optimise_missing_nodes (FALSE);   // estimate sequences at internal nodes

	      // add log-likelihood to total
	      const Score align_sc = tkf.alignment_path_score();
	      CLOG(5) << "Alignment score: " << Score2Bits(align_sc) << " bits\n";
	      ScorePMulAcc (total_sc, align_sc);

	      // loop over branches
	      for_rooted_branches_post (tkf.tree, b)
		if ((*b).length >= dt*3)  // minimum branch length...
		  {
		    CLOG(6) << "Branch: " << tkf.tree.branch_specifier (*b) << "\n";
		
		    // get pairwise alignment for this branch
		    Pairwise_path path = tkf.subpath (*b, FALSE);
		    path.erase_empty_columns();
		    if (CLOGGING(4))
		      {
			vector<sstring> row_name (2);
			row_name[0] = tkf.tree.node_specifier ((*b).first);
			row_name[1] = tkf.tree.node_specifier ((*b).second);
			path.show (CL, row_name);
		      }

		    // get transition counts for this pairwise alignment
		    TKF_transition_counts tcounts;
		    tcounts.add_path_counts (path);
		    if (rev)
		      {
			// when time-reversing a TKF pairwise alignment, we must reverse the sequences too, so D->I transition counts are preserved
			path.swap_parent_child();
			for (int r = 0; r < 2; ++r)
			  for (int c = 0; c < path.columns(); ++c)
			    {
			      // STL's swap() doesn't work here with gcc 4.0, hence neither does reverse()... hey ho...
			      const bool old_val = path[r][c];
			      path[r][c] = path[r][path.columns()-1-c];
			      path[r][path.columns()-1-c] = old_val;
			    }
			tcounts.add_path_counts (path);
		      }
		    tcounts.show (CLOG(5));

		    // accumulate EM counts for this set of transition counts & branch length
		    if (naive)
		      coeff.accumulate_naive_counts (tcounts, (*b).length, em_counts);
		    else
		      coeff.accumulate_counts (tcounts, (*b).length, em_counts);
		  }
	    }

	  // output total score
	  CLOG(7) << "EM iteration #" << iter << ", total score: " << Score2Bits(total_sc) << " bits, previous best: " << Score2Bits(best_sc) << " bits\n";
	  if (total_sc > best_sc) { best_sc = total_sc; dec = 0; }
	  else if (++dec > forgive)
	    {
	      CLOG(7) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
	      break;
	    }
	  
	  // output EM counts
	  em_counts.show (CLOG(8));

	  // set new params
	  em_counts.set_params (params);
	}

      // output final params
      CLOG(8) << "Final lambda= " << params.lambda << " mu= " << params.mu << "\n";
      cout << "lambda: " << params.lambda << "\n";
      cout << "mu:     " << params.mu << "\n";
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      exit(1);
    }

  return 0;
}
