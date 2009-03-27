#include <math.h>
#include "scfg/paircfg.h"
#include "util/logfile.h"
#include "util/Regexp.h"
#include "util/vector_output.h"

// This is not a great test. It tries to be generic, by generating lots of
// pseudorandom data, and I've come to the conclusion that succinct, specific
// hand-crafted tests (such as the test parse tree at the very end of this file)
// are generally better.

// However, operating on the assumption that a flawed test is better than no test,
// here goes:

// (Updated 10/3/07: a concrete example of why random tests are bad...
//  the make_single and make_pair functions apparently contain bugs that mean
//  some randomly generated HMMs are not ergodic. This causes later tests to
//  hang indefinitely. And this only happens on Linux for some godforsaken reason.
//  One wasted hour. Not sure I've fixed it very well either.)

// repeatable pseudo-random number generator

unsigned int X;
void pseudseed (unsigned int x) { X = x; }
int pseudrand (int xmod)
{
  // got this algorithm off the web somewhere. ask me if i care how it works
  const int m = 1 << 31;
  const int a = 16807;
  const int q = 127773;
  const int r = 2836;
  int X_new = a * (X % q) - r * (X / q);
  if (X_new < 0) X_new += m;
  X = (unsigned int) X_new;
  return ((int) (X & 0x7fffffff)) % xmod;
}

Score random_score (Score sc) { return pseudrand(sc) - sc / 2; }

// HMM & sequence factory methods

Single_HMM_scores make_single (int null_states, int emit_states, int transitions, Score sc)
{
  const int total_states = null_states + emit_states;
  if (transitions < total_states) THROWEXPR ("Too few transitions");
  Single_HMM_scores hmm (total_states, CFG_alphabet);
  for (int i = 0; i < null_states; ++i) hmm.state_type[i] = Single_state_typing::Null;
  for (int i = null_states; i < total_states; ++i) hmm.state_type[i] = Single_state_typing::Emit;
  for (int i = 0; i < total_states; ++i)
    {
      hmm.start[i] = random_score(sc);
      hmm.end[i] = random_score(sc);
    }
  hmm.start_to_end() = random_score(sc);
  for (int t = 0; t < transitions; ++t)
    {
      int imax = emit_states == 0 ? null_states - 1 : total_states;
      int i = t < imax ? t : pseudrand (imax);  // ensure ergodicity
      int j = pseudrand (total_states - min (i + 1, null_states)) + min (i + 1, null_states);
      hmm.transition(i,j) = random_score(sc);
    }

  // the "ensure ergodicity" line above doesn't seem to work (causing crashes during tests) and i can't be arsed to fix it right now. so i'm just going to ensure ergodicity by blunt force here. -- IH, 10/3/07
  for (int i = 0; i < total_states - 1; ++i)
    hmm.transition (i, pseudrand(total_states-i) + i) = 0;
  if (total_states > 0)
    hmm.transition(total_states-1,total_states-1) = 0;

  for (int i = null_states; i < null_states + emit_states; ++i)
    for (int c = 0; c < CFG_alphabet_size; ++c)
      hmm.emit[i][c] = random_score(sc);
  return hmm;
}

Pair_HMM_scores make_pair (int x, int y, int xy, int transitions, Score sc)
{
  const int total_states = x + y + xy;
  if (transitions < total_states) THROWEXPR ("Too few transitions");
  Pair_HMM_scores hmm (total_states, &CFG_alphabet);
  for (int i = 0; i < x; ++i) hmm.state_type[i] = Pair_HMM_scores::EmitX;
  for (int i = x; i < x + y; ++i) hmm.state_type[i] = Pair_HMM_scores::EmitY;
  for (int i = x + y; i < x + y + xy; ++i) hmm.state_type[i] = Pair_HMM_scores::EmitXY;
  for (int i = 0; i < total_states; ++i)
    {
      hmm.start[i] = random_score(sc);
      hmm.end[i] = random_score(sc);
    }
  hmm.start_to_end() = random_score(sc);
  for (int t = 0; t < transitions; ++t)
    {
      int i = t < total_states ? t : pseudrand (total_states);  // ensure ergodicity
      int j = pseudrand (total_states);
      hmm.transition(i,j) = random_score(sc);
    }

  // the "ensure ergodicity" line above doesn't seem to work in make_single (causing crashes during tests) and i can't be arsed to fix it right now. so i'm just going to ensure ergodicity by blunt force here. -- IH, 10/3/07
  for (int i = 0; i < total_states - 1; ++i)
    hmm.transition (i, pseudrand(total_states-i) + i) = 0;
  hmm.transition(total_states-1,total_states-1) = 0;

  for (int i = 0; i < x + y; ++i)
    {
      hmm.single_emit[i] = vector<Score> (CFG_alphabet_size);
      for (int c = 0; c < CFG_alphabet_size; ++c)
	hmm.single_emit[i][c] = random_score(sc);
    }
  for (int i = x + y; i < x + y + xy; ++i)
    {
      hmm.pair_emit[i] = array2d<Score> (CFG_alphabet_size, CFG_alphabet_size);
      for (int cx = 0; cx < CFG_alphabet_size; ++cx)
	for (int cy = 0; cy < CFG_alphabet_size; ++cy)
	  hmm.pair_emit[i](cx,cy) = random_score(sc);
    }
  return hmm;
}

Digitized_biosequence make_seq (int seqlen)
{
  Digitized_biosequence dsq;
  for (int i = 0; i < seqlen; ++i) dsq.push_back (pseudrand (CFG_alphabet_size));
  return dsq;
}

// comparison test methods

void test_single_emit (int null_states, int emit_states, int transitions, Score sc, int max_pathlen, int reps, bool reject_inf = 1)
{
  int n_inf = 0;
  double sum_hmmsc = 0;
  double sum_hmmsc_sq = 0;
  for (int rep = 0; rep < reps; ++rep)
    {
      CLOG(5) << "Starting rep " << rep << "\n";

      const Single_HMM_scores hmm = make_single (null_states, emit_states, transitions, sc);
      CLOG(5) << "Built Single_HMM_scores\n";

      const Pair_CFG_scores cfg (hmm);
      CLOG(5) << "Built Pair_CFG_scores\n";

      if (!cfg.test_valid()) THROWEXPR ("CFG invalid");
      CLOG(5) << "Pair_CFG_scores is valid\n";
      
      Pair_CFG_counts test_cfg_counts (cfg);
      CLOG(5) << "Built Pair_CFG_counts\n";

      const int pathlen = pseudrand (max_pathlen + 1);

      vector<int> global_path;
      global_path.push_back (HMM_state_enum::Start);
      Named_profile np;
      Score test_transition_score = 0;
      Score test_emit_score = 0;

      CLOG(5) << "Looping i from 0 to " << pathlen << "\n";
      for (int i = 0; i < pathlen; ++i)
	{
	  if (emit_states == 0 && global_path.back() == null_states - 1) break;  // don't get trapped in all-null models
	  const int prev = global_path.back();
	  CLOG(5) << "In loop: i=" << i << ", prev = " << prev << "\n";
	  int s;
	  do s = pseudrand (hmm.states());
	  while (reject_inf && hmm.transition (prev, s) == -InfinityScore);
	  CLOG(5) << "Found state " << s << "\n";
	  ScorePMulAcc (test_transition_score, hmm.transition (prev, s));
	  test_cfg_counts.transition (prev, s) += 1.0;
	  global_path.push_back (s);
	  if (hmm.state_type[s] == Single_state_typing::Emit)
	    {
	      const int c = pseudrand (CFG_alphabet_size);
	      ScorePMulAcc (test_emit_score, hmm.emit[s][c]);
	      test_cfg_counts.emit[s][c] += 1.0;
	      np.dsq.push_back (c);
	    }
	}
      const int prev = global_path.back();
      ScorePMulAcc (test_transition_score, hmm.transition (prev, HMM_state_enum::End));
      test_cfg_counts.transition (prev, HMM_state_enum::End) += 1.0;
      global_path.push_back (HMM_state_enum::End);

      const Score test_score = ScorePMul (test_transition_score, test_emit_score);

      const vector<int> local_path = Transition_methods::make_local_path (global_path);
      CLOG(5) << "Built local_path\n";

      np.dsq2score (CFG_alphabet);
      CLOG(5) << "Built DSQ\n";

      const Pair_CFG_parse_tree parse_tree = cfg.parse (global_path);
      CLOG(5) << "Built Pair_CFG_parse_tree\n";

      const Digitized_biosequence null_dsq;

      if (!cfg.test_branch_coords_consistent (parse_tree)) THROWEXPR ("Branch coords inconsistent");
      CLOG(5) << "Branch coords consistent\n";

      if (!parse_tree.test_global (np.dsq, null_dsq)) THROWEXPR ("Not global");
      CLOG(5) << "Global\n";

      const Score hmm_score = hmm.path_score (local_path, np);
      CLOG(5) << "Got hmm_score\n";

      const Score cfg_score = cfg.path_score (parse_tree, np.dsq, null_dsq);
      CLOG(5) << "Got cfg_score\n";

      if (cfg_score != hmm_score || cfg_score != test_score)
	{
	  cerr << "HMM:\n";
	  hmm.show (cerr);
	  cerr << "CFG:\n";
	  cfg.show (cerr);
	  np.dsq2seq (CFG_alphabet);
	  cerr << "Sequence:\n" << np.seq << "\n";
	  cerr << "State path: (" << local_path << ")\n";
	  cerr << "HMM score " << hmm_score << " = " << hmm.path_transition_score (local_path) << "(trans) + " << hmm.path_emit_score (local_path, np) << "(emit)\n";
	  cerr << "CFG score " << cfg_score << " = " << cfg.path_transition_score (parse_tree) << "(trans) + " << cfg.path_emit_score (parse_tree, np.dsq, null_dsq) << "(emit)\n";
	  cerr << "Test score " << test_score << " = " << test_transition_score << "(trans) + " << test_emit_score << "(emit)\n";
	  THROWEXPR ("CFG score looks bad");
	}

      if (hmm_score < -sc * max_pathlen)
	++n_inf;
      else
	{
	  sum_hmmsc += hmm_score;
	  sum_hmmsc_sq += pow ((double) hmm_score, 2.);
	}

      Single_HMM_counts hmm_counts (hmm);
      CLOG(5) << "Built hmm_counts\n";

      hmm_counts.add_counts_from_state_path (hmm, np, local_path);
      CLOG(5) << "Called hmm_counts.add_counts_from_state_path\n";

      Pair_CFG_counts cfg_counts (cfg);
      CLOG(5) << "Built cfg_counts\n";

      cfg_counts.add_counts_from_parse_tree (cfg, parse_tree, np.dsq, null_dsq);
      CLOG(5) << "Called cfg_counts.add_counts_from_parse_tree\n";
      
      sstring count_err;

      for (int i = 0; i < hmm.states(); ++i)
	for (int j = 0; j < hmm.states(); ++j)
	  if (abs (cfg_counts.transition (i, j) - test_cfg_counts.transition (i, j)) > TINY
	      || abs (cfg_counts.transition (i, j) - hmm_counts.transition (i, j)) > TINY)
	    count_err << " t(" << i << "," << j << ")";

      for (int s = null_states; s < null_states + emit_states; ++s)
	for (int c = 0; c < CFG_alphabet_size; ++c)
	  if (abs (cfg_counts.emit[s][c] - test_cfg_counts.emit[s][c]) > TINY
	      || abs (cfg_counts.emit[s][c] - hmm_counts.emit[s][c]) > TINY)
	    count_err << " emit[" << s << "][" << c << "]";
      
      if (count_err.size())
	{
	  cerr << "HMM counts:\n";
	  hmm_counts.show (cerr);
	  cerr << "CFG counts:\n";
	  cfg_counts.show (cerr);
	  cerr << "Test CFG counts:\n";
	  test_cfg_counts.show (cerr);
	  THROWEXPR ("Counts don't match:" << count_err);
	}

      CLOG(5) << "Finished rep " << rep << "\n";
    }
  cerr << "(checked " << reps << " random Single HMMs with ";
  cerr << null_states << " null, " << emit_states << " emit, ";
  cerr << transitions << " transitions, |sc|<" << sc;
  cerr << ", max pathlen " << max_pathlen;
  if (n_inf) cerr << ", " << n_inf << " neginf paths";
  if (n_inf < reps)
    {
      const double mean_hmmsc = sum_hmmsc / (double) (reps - n_inf);
      const double sd_hmmsc = sqrt (sum_hmmsc_sq / (double) (reps - n_inf) - mean_hmmsc * mean_hmmsc);
      cerr << ", ";
      if (!reject_inf) cerr << "finite ";
      cerr << "path score " << (int) mean_hmmsc << " +/- " << (int) sd_hmmsc;
    }
  cerr << ")\n";
}

void test_pair_emit (int x, int y, int xy, int transitions, Score sc, int max_pathlen, int reps, bool reject_inf = 1)
{
  int n_inf = 0;
  double sum_hmmsc = 0;
  double sum_hmmsc_sq = 0;
  for (int rep = 0; rep < reps; ++rep)
    {
      CLOG(5) << "Starting rep " << rep << "\n";

      const Pair_HMM_scores hmm = make_pair (x, y, xy, transitions, sc);
      CLOG(5) << "Built Pair_HMM_scores\n";

      const Pair_CFG_scores cfg (hmm);
      CLOG(5) << "Built Pair_CFG_scores\n";

      if (!cfg.test_valid()) THROWEXPR ("CFG invalid");
      CLOG(5) << "Pair_CFG_scores is valid\n";
      
      Pair_CFG_counts test_cfg_counts (cfg);
      CLOG(5) << "Built Pair_CFG_counts\n";

      const int pathlen = pseudrand (max_pathlen + 1);

      vector<int> global_path;
      global_path.push_back (HMM_state_enum::Start);

      Score test_transition_score = 0;
      Score test_emit_score = 0;

      if (pathlen)
	{
	  global_path.push_back (pseudrand (hmm.states()));
	  ScorePMulAcc (test_transition_score, hmm.start[global_path.back()]);
	  test_cfg_counts.start[global_path.back()] += 1.0;
	}

      Named_profile npx;
      Named_profile npy;
      for (int i = 0; i < pathlen; ++i)
	{
	  const int prev = global_path.back();
	  switch (hmm.state_type[prev])
	    {
	    case Pair_HMM_scores::EmitX:
	      {
		const int c = pseudrand (CFG_alphabet_size);
		ScorePMulAcc (test_emit_score, hmm.single_emit[prev][c]);
		test_cfg_counts.emit[prev][c] += 1.0;
		npx.dsq.push_back (c);
		break;
	      }
	    case Pair_HMM_scores::EmitY:
	      {
		const int c = pseudrand (CFG_alphabet_size);
		ScorePMulAcc (test_emit_score, hmm.single_emit[prev][c]);
		test_cfg_counts.emit[prev][c] += 1.0;
		npy.dsq.push_back (c);
		break;
	      }
	    case Pair_HMM_scores::EmitXY:
	      {
		const int cx = pseudrand (CFG_alphabet_size);
		const int cy = pseudrand (CFG_alphabet_size);
		ScorePMulAcc (test_emit_score, hmm.pair_emit[prev](cx,cy));
		test_cfg_counts.emit[prev][cx * cfg.emit_xl_mul(Pair_CFG_scores::EmitXLYL)
					  +cy * cfg.emit_yl_mul(Pair_CFG_scores::EmitXLYL)] += 1.0;
		npx.dsq.push_back (cx);
		npy.dsq.push_back (cy);
		break;
	      }
	    default:
	      THROWEXPR ("State type unknown");
	      break;
	    }
	  if (i < pathlen - 1)
	    {
	      int s;
	      do s = pseudrand (hmm.states());
	      while (reject_inf && hmm.transition (prev, s) == -InfinityScore);
	      ScorePMulAcc (test_transition_score, hmm.transition (prev, s));
	      test_cfg_counts.transition (prev, s) += 1.0;
	      global_path.push_back (s);
	    }
	}
      const int prev = global_path.back();
      if (prev < 0)
	{
	  ScorePMulAcc (test_transition_score, hmm.start_to_end());
	  test_cfg_counts.start_to_end() += 1.0;
	}
      else
	{
	  ScorePMulAcc (test_transition_score, hmm.end[prev]);
	  test_cfg_counts.end[prev] += 1.0;
	}
      global_path.push_back (HMM_state_enum::End);

      const Score test_score = ScorePMul (test_transition_score, test_emit_score);

      const vector<int> local_path = Transition_methods::make_local_path (global_path);
      CLOG(5) << "Built make_local_path\n";

      const Pair_CFG_parse_tree parse_tree = cfg.parse (global_path);
      CLOG(5) << "Built Pair_CFG_parse_tree\n";

      const Digitized_biosequence null_dsq;

      if (!cfg.test_branch_coords_consistent (parse_tree)) THROWEXPR ("Branch coords inconsistent");
      CLOG(5) << "Branch coords consistent\n";

      if (!parse_tree.test_global (npx.dsq, npy.dsq)) THROWEXPR ("Not global");
      CLOG(5) << "Global\n";

      npx.dsq2score (CFG_alphabet);
      npy.dsq2score (CFG_alphabet);
      CLOG(5) << "Built DSQs\n";

      const Score hmm_score = hmm.path_score (local_path, npx.prof_sc, npy.prof_sc);
      CLOG(5) << "Got HMM path score\n";

      const Score cfg_score = cfg.path_score (parse_tree, npx.dsq, npy.dsq);
      CLOG(5) << "Got CFG path score\n";

      if (cfg_score != hmm_score || cfg_score != test_score)
	{
	  cerr << "HMM:\n";
	  hmm.show (cerr);
	  cerr << "CFG:\n";
	  cfg.show (cerr);
	  npx.dsq2seq (CFG_alphabet);
	  npy.dsq2seq (CFG_alphabet);
	  cerr << "Sequence X:\n" << npx.seq << "\n";
	  cerr << "Sequence Y:\n" << npy.seq << "\n";
	  cerr << "State path: (" << local_path << ")\n";
	  cerr << "HMM score " << hmm_score << " = " << hmm.path_transition_score (local_path) << "(trans) + " << hmm.path_emit_score (local_path, npx.prof_sc, npy.prof_sc) << "(emit)\n";
	  cerr << "CFG score " << cfg_score << " = " << cfg.path_transition_score (parse_tree) << "(trans) + " << cfg.path_emit_score (parse_tree, npx.dsq, npy.dsq) << "(emit)\n";
	  cerr << "Test score " << test_score << " = " << test_transition_score << "(trans) + " << test_emit_score << "(emit)\n";
	  THROWEXPR ("CFG score looks bad");
	}

      if (hmm_score < -sc * max_pathlen)
	++n_inf;
      else
	{
	  sum_hmmsc += hmm_score;
	  sum_hmmsc_sq += pow ((double) hmm_score, 2.);
	}

      Pair_HMM_counts hmm_counts (hmm);
      CLOG(5) << "Built Pair_HMM_counts\n";

      hmm_counts.add_counts_from_state_path (hmm, npx.prof_sc, npy.prof_sc, local_path);
      CLOG(5) << "Called Pair_HMM_counts::add_counts_from_state_path\n";

      Pair_CFG_counts cfg_counts (cfg);
      CLOG(5) << "Built Pair_CFG_counts\n";

      cfg_counts.add_counts_from_parse_tree (cfg, parse_tree, npx.dsq, npy.dsq);
      CLOG(5) << "Called Pair_CFG_counts::add_counts_from_parse_tree\n";
      
      sstring count_err;

      for (int i = 0; i < hmm.states(); ++i)
	for (int j = 0; j < hmm.states(); ++j)
	  if (abs (cfg_counts.transition (i, j) - test_cfg_counts.transition (i, j)) > TINY
	      || abs (cfg_counts.transition (i, j) - hmm_counts.transition (i, j)) > TINY)
	    count_err << " t(" << i << "," << j << ")";

      for (int s = 0; s < hmm.states(); ++s)
	if (hmm.state_type[s] == Pair_HMM_scores::EmitXY)
	  for (int cx = 0; cx < CFG_alphabet_size; ++cx)
	    for (int cy = 0; cy < CFG_alphabet_size; ++cy)
	      {
		const int emit_idx =
		  cx * cfg.emit_xl_mul(Pair_CFG_scores::EmitXLYL) +
		  cy * cfg.emit_yl_mul(Pair_CFG_scores::EmitXLYL);
		
		if (abs (cfg_counts.emit[s][emit_idx] - test_cfg_counts.emit[s][emit_idx]) > TINY
		    || abs (cfg_counts.emit[s][emit_idx] - hmm_counts.pair_emit[s](cx,cy)) > TINY)
		  count_err << " emit[" << s << "](" << cx << "," << cy << ")";
	      }
	else
	  for (int c = 0; c < CFG_alphabet_size; ++c)
	    if (abs (cfg_counts.emit[s][c] - test_cfg_counts.emit[s][c]) > TINY
		|| abs (cfg_counts.emit[s][c] - hmm_counts.single_emit[s][c]) > TINY)
	      count_err << " emit[" << s << "][" << c << "]";
      
      if (count_err.size())
	{
	  cerr << "HMM counts:\n";
	  hmm_counts.show (cerr);
	  cerr << "CFG counts:\n";
	  cfg_counts.show (cerr);
	  cerr << "Test CFG counts:\n";
	  test_cfg_counts.show (cerr);
	  THROWEXPR ("Counts don't match:" << count_err);
	}

      CLOG(5) << "Finished rep " << rep << "\n";
    }
  cerr << "(checked " << reps << " random Pair HMMs with ";
  cerr << x << " x, " << y << " y, " << xy << " xy, ";
  cerr << transitions << " transitions, |sc|<" << sc;
  cerr << ", max pathlen " << max_pathlen;
  if (n_inf) cerr << ", " << n_inf << " neginf paths";
  if (n_inf < reps)
    {
      const double mean_hmmsc = sum_hmmsc / (double) (reps - n_inf);
      const double sd_hmmsc = sqrt (sum_hmmsc_sq / (double) (reps - n_inf) - mean_hmmsc * mean_hmmsc);
      cerr << ", ";
      if (!reject_inf) cerr << "finite ";
      cerr << "path score " << (int) mean_hmmsc << " +/- " << (int) sd_hmmsc;
    }
  cerr << ")\n";
}

// main

int main (int argc, char** argv)
{
  pseudseed (27815386);  // seed is irrelevant really

  INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
		  "test pairwise SCFGs");

  opts.parse_or_die();

  try
    {
      // create a test parse tree that uses all 18 state types

      cerr << "(testing a simple pair CFG parse)\n";

      Biosequence xseq = "tacgagggatcgagctagcgtatgctcgta";
      // xfold:           >>>>....><><..>><><<.><><.<<<<
      Biosequence yseq = "cacagggctatcctatcggctctacgcctg";
      // yfold:           >>.>..><><..>.><><<...><><..<<

      // NB "abcd" means YL=a, XL=b, XR=c, YR=d
      int state_array[44] = { /* root */ -1, 15, 15, 7, 11, 16,          // Start, tcga, aatt, cc-g, g-cc, Bifurc
			      /* root->left */ 0, 13, 5, 9, 17,          // Null,  aat-, gg--, g-c-, BifurcRevY
			      /* root->left->left */ 0, 14, 6, 1, 16,    // Null,  -cgg, -t-a, g---, Bifurc
			      /* root->left->left->left */ 0, 15, -2,    // Null,  aatt, End
			      /* root->left->left->right */ 0, 15, -2,   // Null,  ccgg, End
			      /* root->left->right */ 0, 12, 10, 3, 17,  // Null,  -gc-, --tt, c--g, BifurcRevY
			      /* root->left->right->left */ 0, 15, -2,   // Null,  ttaa, End
			      /* root->left->right->right */ 0, 15, -2,  // Null,  ggcc, End
			      /* root->right */ 0, 4, 8, 2, 16,          // Null,  -c--, --c-, ---t, Bifurc
			      /* root->right->left */ 0, 15, -2,         // Null,  atat, End
			      /* root->right->right */ 0, 15, -2 };      // Null,  gcgc, End
      
      const Score test_trans_sc =
	9915 + 1515 + 1507 +  711 + 1116
	+   13 + 1305 +  509 +  917
	+   14 + 1406 +  601 +  116
	+   15 + 1599
	+   15 + 1599
	+   12 + 1210 + 1003 +  317
	+   15 + 1599
	+   15 + 1599
	+    4 +  408 +  802 +  216
	+   15 + 1599
	+   15 + 1599;
      
      const Score test_emit_sc =
	4231 + 31144 + 12293 + 23922  // tcga, aatt, cc-g, g-cc
	+ 11149 + 13399 +  3929       // aat-, gg--, g-c-
	+ 29233 +  9491 +  3999       // -cgg, -t-a, g---
	+ 31144                       // aatt
	+ 32233                       // ccgg
	+  9329 + 29944 +  2993       // -gc-, --tt, c--g
	+ 34411                       // ttaa
	+ 33322                       // ggcc
	+  9299 +  9929 +  9994       // -c--, --c-, ---t
	+  1414                       // atat
	+  3232;                      // gcgc

      vector<int> state_path (state_array, state_array + 44);
      
      Digitized_biosequence xdsq;
      Digitized_biosequence ydsq;

      DNA_alphabet.seq2dsq (xseq, xdsq);
      DNA_alphabet.seq2dsq (yseq, ydsq);

      Pair_CFG_scores cfg (18);
      cfg.start_to_end() = 9999;
      for (int i = 0; i < 18; ++i)
	{
	  cfg.start[i] = 9900 + i;
	  cfg.end[i] = i * 100 + 99;
	  cfg.init_emit (i, (Pair_CFG_scores::State_type) (i < 17 ? i : 32));
	  const Pair_CFG_scores::State_type t = cfg.state_type[i];
	  if (Pair_CFG_scores::is_emit_type (t))
	    for (int cxl = 0; cxl < 4; ++cxl)
	      for (int cxr = 0; cxr < 4; ++cxr)
		for (int cyl = 0; cyl < 4; ++cyl)
		  for (int cyr = 0; cyr < 4; ++cyr)
		    if (cxl == 3 - cxr && cyl == 3 - cyr)  // no mispairings allowed
		      {
			const int emit_idx =
			  cxl * cfg.emit_xl_mul (t) +
			  cxr * cfg.emit_xr_mul (t) +
			  cyl * cfg.emit_yl_mul (t) +
			  cyr * cfg.emit_yr_mul (t);
			int emit_code = 0;
			emit_code += ((t & Pair_CFG_scores::EmitXL) ? (cxl + 1) : 9) * 1000;
			emit_code += ((t & Pair_CFG_scores::EmitYL) ? (cyl + 1) : 9) * 100;
			emit_code += ((t & Pair_CFG_scores::EmitYR) ? (cyr + 1) : 9) * 10;
			emit_code += ((t & Pair_CFG_scores::EmitXR) ? (cxr + 1) : 9) * 1;
			if ((t & Pair_CFG_scores::EmitXL) && (t & Pair_CFG_scores::EmitYL) && cxl == cyl) emit_code += 10000;
			if ((t & Pair_CFG_scores::EmitXR) && (t & Pair_CFG_scores::EmitYR) && cxr == cyr) emit_code += 20000;
			cfg.emit[i][emit_idx] = emit_code;
		      }
	  for (int j = 0; j < 18; ++j)
	    cfg.transition (i, j) = i * 100 + j;
	}
      cfg.bifurc[16] = cfg.bifurc[17] = Pair_CFG_scores::Bifurcation (0, 0);

      const Pair_CFG_parse_tree parse_tree = cfg.parse (state_path);
      const vector<int> test_path = cfg.unparse (parse_tree);
      if (test_path != state_path) THROWEXPR ("state path mangled");

      if (!cfg.test_valid()) THROWEXPR ("CFG invalid");
      if (!cfg.test_branch_coords_consistent (parse_tree)) THROWEXPR ("Branch coords inconsistent");
      if (!parse_tree.test_global (xdsq, ydsq)) THROWEXPR ("Not global");

      const Score trans_sc = cfg.path_transition_score (parse_tree);
      const Score emit_sc = cfg.path_emit_score (parse_tree, xdsq, ydsq);

      if (trans_sc != test_trans_sc) THROWEXPR ("CFG path transition score " << trans_sc << " != expected score " << test_trans_sc);
      if (emit_sc != test_emit_sc) THROWEXPR ("CFG path emit score " << emit_sc << " != expected score " << test_emit_sc);

      // generate paths through Single HMMs & compare path scores

      cerr << "(generating some pseudorandom HMMs & comparing path scores to CFG)\n";

      test_single_emit (0, 1, 1, Bits2Score(1), 100, 50);   // one emit state
      test_single_emit (0, 5, 4*5, Bits2Score(3), 100, 50);  // several emit states
      test_single_emit (5, 10, 15*15, Bits2Score(6), 20, 50);  // several emit & null states
      test_single_emit (5, 10, 15*15, Bits2Score(6), 20, 10, 0);  // allow zero-probability paths

      test_single_emit (10, 0, 10*10, Bits2Score(20), 100, 10);  // no emit states
      test_single_emit (0, 0, 0, Bits2Score(3), 0, 1);  // zero path length, no states at all
      test_single_emit (3, 3, 3*3, Bits2Score(3), 0, 1);  // zero path length, but with states

      // generate paths through Pair HMMs & compare path scores

      test_pair_emit (1, 0, 0, 1, Bits2Score(1), 100, 10);   // one EmitX state
      test_pair_emit (0, 1, 0, 1, Bits2Score(1), 100, 10);   // one EmitY state
      test_pair_emit (0, 0, 1, 1, Bits2Score(1), 100, 10);   // one EmitXY state

      test_pair_emit (1, 1, 1, 3*3*2, Bits2Score(5), 200, 20);   // one state of each type
      test_pair_emit (5, 5, 5, 15*15*2, Bits2Score(10), 300, 20);   // five states of each type

    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << "ok\n";
  return 0;
}
