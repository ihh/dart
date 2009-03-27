#include <math.h>
#include <map>
#include "scfg/paircfgdp.h"
#include "hmm/singlehmm.h"
#include "util/logfile.h"
#include "util/Regexp.h"
#include "util/vector_output.h"

const double tol = 0.1;  // % error tolerance

// score comparison that allows for approx infinities
bool scores_equal (Score a, Score b)
{
  const double a_inf = 100.0 * ((double) (a + InfinityScore)) / (double) InfinityScore;
  const double b_inf = 100.0 * ((double) (b + InfinityScore)) / (double) InfinityScore;
  if (a_inf < tol) return b_inf < tol;
  return a == b;
}

// state path comparison that looks at transitions only
bool same_transition_usage (const vector<int>& path1, const vector<int>& path2)
{
  map<pair<int,int>,int> trans1;
  map<pair<int,int>,int> trans2;
  for (int i = 1; i < (int) path1.size(); ++i) ++trans1 [pair<int,int> (path1[i-1], path1[i])];
  for (int i = 1; i < (int) path2.size(); ++i) ++trans2 [pair<int,int> (path2[i-1], path2[i])];
  return trans1 == trans2;
}

// comparison of Pair_CFG_counts

void test_counts_equal (Prob a, Prob b, const char* counts_type)
{
  if (abs(a-b) >= .01 && abs((a-b)/(b==0?1:b)) >= .01)
    THROWEXPR (counts_type << " counts don't match: a=" << a << " b=" << b);
}

void test_cfg_counts_equal (const Pair_CFG_counts& c1, const Pair_CFG_counts& c2)
{
  if (!c1.same_state_types(c2)) THROWEXPR ("State types don't match");
  test_counts_equal (c1.start_to_end(), c2.start_to_end(), "Start-to-end");
  for (int s = 0; s < c1.states(); ++s)
    {
      test_counts_equal (c1.start[s], c2.start[s], "Start");
      test_counts_equal (c1.end[s], c2.end[s], "End");
      for (int d = 0; d < c1.states(); ++d)
	test_counts_equal (c1.transition(s,d), c2.transition(s,d), "Transition");
      for (int c = 0; c < (int) c1.emit[s].size(); ++c)
	test_counts_equal (c1.emit[s][c], c2.emit[s][c], "Emit");
    }
}

// repeatable pseudo-random number generator
unsigned int X = 27815386;  // seed is irrelevant really
void showseed() { cerr << "pseudseed(" << X << ");\n"; }
void pseudseed (unsigned int x)
{
  //  if (X != x) cerr << "(changed seed!)\n";
  X = x;
}
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

void test_single_emit (int null_states, int emit_states, int transitions, Score sc, int max_seqlen, int reps, bool reject_inf = 1)
{
  //  showseed();
  cerr << "(checking " << reps << " random Single HMMs with ";
  cerr << null_states << " null, " << emit_states << " emit, ";
  cerr << transitions << " transitions, |sc|<" << sc;
  cerr << ", max seqlen " << max_seqlen;
  cerr << ")\n";
  vector<double> pe_log;
  for (int rep = 0; rep < reps; ++rep)
    {
      // create HMM & CFG
      
      const Single_HMM_scores hmm = make_single (null_states, emit_states, transitions, sc);
      const Pair_CFG_scores cfg (hmm);

      Named_profile npx;
      Named_profile npy;

      npx.dsq = make_seq (pseudrand (max_seqlen + 1));
      npx.dsq2score (CFG_alphabet);

      Fold_envelope xenv;
      Fold_envelope yenv;

      xenv.initialise_full (npx.dsq.size());
      yenv.initialise_full (npy.dsq.size());

      // compare forward/inside
      
      Single_fast_forward_matrix forward (hmm, npx);
      const Score hmm_forward_score = forward.final_score;

      Pair_inside_matrix inside (npx, npy, xenv, yenv, cfg);
      const Score cfg_inside_score = inside.final_score;

      double percent_error = 100 * abs ((double) (cfg_inside_score - hmm_forward_score) / (double) hmm_forward_score);
      if (percent_error > tol)
	{
	  cerr << "HMM:\n";
	  hmm.show (cerr);
	  cerr << "CFG:\n";
	  cfg.show (cerr);
	  npx.dsq2seq (CFG_alphabet);
	  cerr << "Sequence:\n" << npx.seq << "\n";
	  cerr << "HMM forward matrix:\n";
	  forward.show (cerr);
	  cerr << "CFG inside matrix:\n";
	  inside.show (cerr);
	  
	  sstring tol_err;
	  tol_err << "CFG score " << cfg_inside_score << " != HMM score " << hmm_forward_score << ", off by ";
	  tol_err.precision (4);
	  tol_err << percent_error << "%";
	  THROWEXPR (tol_err);
	}
      else if (percent_error > 0.0)
	pe_log.push_back (percent_error);

      Single_fast_Viterbi_matrix viterbi (hmm, npx);
      const Score hmm_viterbi_score = viterbi.final_score;
      
      Pair_CYK_matrix cyk (npx, npy, xenv, yenv, cfg);
      const Score cfg_cyk_score = cyk.final_score;

      if (!scores_equal (hmm_viterbi_score, cfg_cyk_score))
	{
	  cerr << "HMM:\n";
	  hmm.show (cerr);
	  cerr << "CFG:\n";
	  cfg.show (cerr);
	  npx.dsq2seq (CFG_alphabet);
	  cerr << "Sequence:\n" << npx.seq << "\n";
	  cerr << "HMM Viterbi matrix:\n";
	  viterbi.show (cerr);
	  cerr << "CFG CYK matrix:\n";
	  cyk.show (cerr);
	  THROWEXPR ("CFG CYK score " << cfg_cyk_score << " != HMM Viterbi score " << hmm_viterbi_score);
	}

      // compare Viterbi/CYK
      
      const vector<int> viterbi_traceback = viterbi.optimal_state_path();
      const vector<int> global_viterbi_traceback = Transition_methods::make_global_path (viterbi_traceback);
      
      const vector<int> cyk_traceback = cyk.traceback();
      
      if (cyk_traceback != global_viterbi_traceback)
	{
	  cerr << "Warning: tracebacks don't match\n";
	  cerr << "HMM Viterbi traceback: (" << global_viterbi_traceback << ")\n";
	  cerr << "CFG CYK traceback:     (" << cyk_traceback << ")\n";
	  if (same_transition_usage (global_viterbi_traceback, cyk_traceback))
	    cerr << "Same transition usage; won't look too closely\n";
	  else
	    {
	      cerr << "HMM:\n";
	      hmm.show (cerr);
	      cerr << "CFG:\n";
	      cfg.show (cerr);
	      npx.dsq2seq (CFG_alphabet);
	      cerr << "Sequence:\n" << npx.seq << "\n";
	      cerr << "HMM Viterbi matrix:\n";
	      viterbi.show (cerr);
	      cerr << "CFG CYK matrix:\n";
	      cyk.show (cerr);
	      THROWEXPR ("Tracebacks don't match");
	    }
	}

      const Score viterbi_traceback_score = hmm.path_score (viterbi_traceback, npx);
      if (!scores_equal (viterbi_traceback_score, hmm_viterbi_score))
	THROWEXPR ("HMM traceback score " << viterbi_traceback_score << " != matrix score " << hmm_viterbi_score);
      
      const Pair_CFG_parse_tree cyk_parse_tree = cfg.parse (cyk_traceback);

      if (!cfg.test_branch_coords_consistent (cyk_parse_tree)) THROWEXPR ("CYK traceback: branch coords inconsistent");
      if (!cyk_parse_tree.test_global (npx.dsq, npy.dsq)) THROWEXPR ("CYK traceback: not global");

      const Score cyk_parse_tree_score = cfg.path_score (cyk_parse_tree, npx.dsq, npy.dsq);
      if (!scores_equal (cyk_parse_tree_score, cfg_cyk_score))
	THROWEXPR ("CFG traceback score " << cyk_parse_tree_score << " != matrix score " << cfg_cyk_score);

      // compare forward-backward / inside-outside
      
      const Single_forward_backward_matrix fwd_back (hmm, npx);
      const Pair_CFG_counts fwd_back_count (fwd_back.counts);
      const Pair_inside_outside_matrix in_out (npx, npy, xenv, yenv, cfg);
      //      cout << "Forward-backward counts:\n";
      //      fwd_back_count.show(cout);
      //      cout << "Inside-outside matrix:\n";
      //      in_out.show(cout);
      test_cfg_counts_equal (in_out.count, fwd_back_count);
    }
  if (pe_log.size())
    {
      vector<int> count;
      vector<double> level;
      do
	{
	  count.push_back (0);
	  level.push_back (tol / pow (10., (double) count.size()));
	  for_const_contents (vector<double>, pe_log, pe) if (*pe >= level.back()) ++count.back();
	}
      while (count.back() < (int) pe_log.size());
      cerr << "(inside/forward score deltas: ";
      for (int i = count.size() - 1; i >= 0; --i)
	{
	  if (count[i] == 0) break;
	  if (i < (int) count.size() - 1) cerr << ", ";
	  cerr << count[i] << " over " << level[i] << "%";
	}
      cerr << ")\n";
    }
}

void test_pair_emit (int x, int y, int xy, int transitions, Score sc, int max_xlen, int max_ylen, int reps, bool reject_inf = 1)
{
  //  showseed();
  cerr << "(checking " << reps << " random Pair HMMs with ";
  cerr << x << " x, " << y << " y, " << xy << " xy, ";
  cerr << transitions << " transitions, |sc|<" << sc;
  cerr << ", max xlen " << max_xlen;
  cerr << ", max ylen " << max_ylen;
  cerr << ")\n";
  vector<double> pe_log;
  for (int rep = 0; rep < reps; ++rep)
    {
      // create HMM & CFG

      const Pair_HMM_scores hmm = make_pair (x, y, xy, transitions, sc);
      const Pair_CFG_scores cfg (hmm);
      
      Named_profile npx;
      Named_profile npy;

      npx.dsq = make_seq (pseudrand (max_xlen + 1));
      npy.dsq = make_seq (x==0 && y==0 ? npx.dsq.size() : pseudrand (max_ylen + 1));

      npx.dsq2score (CFG_alphabet);
      npy.dsq2score (CFG_alphabet);

      Fold_envelope xenv;
      Fold_envelope yenv;

      xenv.initialise_full (npx.dsq.size());
      yenv.initialise_full (npy.dsq.size());

      // compare forward/inside

      Pair_forward_DP_matrix forward (hmm, npx.prof_sc, npy.prof_sc);
      const Score hmm_forward_score = forward.final_score;

      Pair_inside_matrix inside (npx, npy, xenv, yenv, cfg);
      const Score cfg_inside_score = inside.final_score;

      double percent_error = 100 * abs ((double) (cfg_inside_score - hmm_forward_score) / (double) hmm_forward_score);
      if (percent_error > tol)
	{
	  cerr << "HMM:\n";
	  hmm.show (cerr);
	  cerr << "CFG:\n";
	  cfg.show (cerr);
	  npx.dsq2seq (CFG_alphabet);
	  npy.dsq2seq (CFG_alphabet);
	  cerr << "Sequence X:\n" << npx.seq << "\n";
	  cerr << "Sequence Y:\n" << npy.seq << "\n";
	  cerr << "HMM forward matrix:\n";
	  forward.show (cerr);
	  cerr << "CFG inside matrix:\n";
	  inside.show (cerr);
	  
	  sstring tol_err;
	  tol_err << "CFG inside score " << cfg_inside_score << " != HMM forward score " << hmm_forward_score << ", off by ";
	  tol_err.precision (4);
	  tol_err << percent_error << "%";
	  THROWEXPR (tol_err);
	}
      else if (percent_error > 0.0)
	pe_log.push_back (percent_error);

      // compare Viterbi/CYK

      Pair_Viterbi_DP_matrix viterbi (hmm, npx.prof_sc, npy.prof_sc);
      const Score hmm_viterbi_score = viterbi.final_score;
      
      Pair_CYK_matrix cyk (npx, npy, xenv, yenv, cfg);
      const Score cfg_cyk_score = cyk.final_score;

      if (!scores_equal (hmm_viterbi_score, cfg_cyk_score))
	{
	  cerr << "HMM:\n";
	  hmm.show (cerr);
	  cerr << "CFG:\n";
	  cfg.show (cerr);
	  npx.dsq2seq (CFG_alphabet);
	  npy.dsq2seq (CFG_alphabet);
	  cerr << "Sequence X:\n" << npx.seq << "\n";
	  cerr << "Sequence Y:\n" << npy.seq << "\n";
	  cerr << "HMM Viterbi matrix:\n";
	  viterbi.show (cerr);
	  cerr << "CFG CYK matrix:\n";
	  cyk.show (cerr);
	  THROWEXPR ("CFG CYK score " << cfg_cyk_score << " != HMM Viterbi score " << hmm_viterbi_score);
	}

      if (hmm_viterbi_score == -InfinityScore)
	{
	  cerr << "(oops! score was -infinity; skipping traceback)\n";
	  continue;
	}
      
      const vector<int> viterbi_traceback = viterbi.optimal_state_path();
      const vector<int> global_viterbi_traceback = Transition_methods::make_global_path (viterbi_traceback);
      
      const vector<int> cyk_traceback = cyk.traceback();
      
      if (cyk_traceback != global_viterbi_traceback)
	{
	  cerr << "Warning: tracebacks don't match\n";
	  cerr << "HMM Viterbi traceback: (" << global_viterbi_traceback << ")\n";
	  cerr << "CFG CYK traceback      (" << cyk_traceback << ")\n";
	  if (same_transition_usage (global_viterbi_traceback, cyk_traceback))
	    cerr << "Same transitions, though; we'll let it pass\n";
	  else
	    {
	      cerr << "HMM:\n";
	      hmm.show (cerr);
	      cerr << "CFG:\n";
	      cfg.show (cerr);
	      npx.dsq2seq (CFG_alphabet);
	      npy.dsq2seq (CFG_alphabet);
	      cerr << "Sequence X:\n" << npx.seq << "\n";
	      cerr << "Sequence Y:\n" << npy.seq << "\n";
	      cerr << "HMM Viterbi matrix:\n";
	      viterbi.show (cerr);
	      cerr << "CFG CYK matrix:\n";
	      cyk.show (cerr);
	      THROWEXPR ("Tracebacks don't match");
	    }
	}

      const Score viterbi_traceback_score = hmm.path_score (viterbi_traceback, npx.prof_sc, npy.prof_sc);
      if (!scores_equal (viterbi_traceback_score, hmm_viterbi_score))
	THROWEXPR ("HMM traceback score " << viterbi_traceback_score << " != matrix score " << hmm_viterbi_score);
      
      const Pair_CFG_parse_tree cyk_parse_tree = cfg.parse (cyk_traceback);

      if (!cfg.test_branch_coords_consistent (cyk_parse_tree)) THROWEXPR ("CYK traceback: branch coords inconsistent");
      if (!cyk_parse_tree.test_global (npx.dsq, npy.dsq)) THROWEXPR ("CYK traceback: not global");

      const Score cyk_parse_tree_score = cfg.path_score (cyk_parse_tree, npx.dsq, npy.dsq);
      if (!scores_equal (cyk_parse_tree_score, cfg_cyk_score))
	THROWEXPR ("CFG traceback score " << cyk_parse_tree_score << " != matrix score " << cfg_cyk_score);

      // compare forward-backward / inside-outside
      
      const Pair_forward_backward_DP_matrix fwd_back (hmm, npx.prof_sc, npy.prof_sc);
      const Pair_CFG_counts fwd_back_count (fwd_back.counts);
      const Pair_inside_outside_matrix in_out (npx, npy, xenv, yenv, cfg);
      //      cout << "Forward-backward counts:\n";
      //      fwd_back_count.show(cout);
      //      cout << "Inside-outside matrix:\n";
      //      in_out.show(cout);
      test_cfg_counts_equal (in_out.count, fwd_back_count);
    }
  if (pe_log.size())
    {
      vector<int> count;
      vector<double> level;
      do
	{
	  count.push_back (0);
	  level.push_back (tol / pow (10., (double) count.size()));
	  for_const_contents (vector<double>, pe_log, pe) if (*pe >= level.back()) ++count.back();
	}
      while (count.back() < (int) pe_log.size());
      cerr << "(inside/forward score deltas: ";
      for (int i = count.size() - 1; i >= 0; --i)
	{
	  if (count[i] == 0) break;
	  if (i < (int) count.size() - 1) cerr << ", ";
	  cerr << count[i] << " over " << level[i] << "%";
	}
      cerr << ")\n";
    }
}

// main

int main (int argc, char** argv)
{
  pseudseed (27815386);

  INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
		  "test pairwise SCFG dynamic programming");

  opts.parse_or_die();

  try
    {
      // test CYK deterministically

      Named_profile npx;
      Named_profile npy;

      Biosequence& xseq = npx.seq;
      Biosequence& yseq = npy.seq;

      xseq          = "tacgagggatcgagctagcgtatgctcgta";
      sstring xfold = "<<<<....<><>..<<><>>.<><>.>>>>";
      yseq          = "cacagggctatcctatcggctctacggctg";
      sstring yfold = "<<.<.<<><>.><.<><>>.>.<><>..>>";

      // NB "abcd" means YL=a, XL=b, XR=c, YR=d
      int constrained_state_array[44] = { /* root */ -1, 15, 15, 7, 11, 16,          // Start, tcga, aatt, cc-g, g-cc, Bifurc
					  /* root->left */ 0, 13, 5, 9, 17,          // Null,  aat-, gg--, g-c-, BifurcRevY
					  /* root->left->left */ 0, 14, 6, 1, 16,    // Null,  -cgg, -t-a, g---, Bifurc
					  /* root->left->left->left */ 0, 15, -2,    // Null,  aatt, End
					  /* root->left->left->right */ 0, 15, -2,   // Null,  ccgg, End
					  /* root->left->right */ 0, 12, 10, 3, 17,  // Null,  -gc-, --tt, c--g, BifurcRevY
					  /* root->left->right->left */ 0, 15, -2,   // Null,  ttaa, End
					  /* root->left->right->right */ 0, 15, -2,  // Null,  ggcc, End
					  /* root->right */ 0, 4, 8, 2, 16,          // Null,  -c--, --g-, ---t, Bifurc
					  /* root->right->left */ 0, 15, -2,         // Null,  atat, End
					  /* root->right->right */ 0, 15, -2 };      // Null,  gcgc, End
      
      const Score P = -1000;
      const Score constrained_emit_sc =
	P + 0 + 2*P + 2*P  // tcga, aatt, cc-g, g-cc
	+ 2*P + 4*P + 4*P  // aat-, gg--, g-c-
	+ 2*P + 4*P + 6*P  // -cgg, -t-a, g---
	+ 0                // aatt
	+ 0                // ccgg
	+ 4*P + 4*P + 4*P  // -gc-, --tt, c--g
	+ 0                // ttaa
	+ 0                // ggcc
	+ 6*P + 6*P + 6*P  // -c--, --g-, ---t
	+ P                // atat
	+ P;               // gcgc

      int unconstrained_state_array[43] = { /* root */ -1, 15, 15, 7, 11, 16,          // Start, tcga, aatt, cc-g, g-cc, Bifurc
					    /* root->left */ 0, 13, 5, 9, 17,          // Null,  aat-, gg--, g-c-, BifurcRevY
					    /* root->left->left */ 0, 14, 6, 1, 16,    // Null,  -cgg, -t-a, g---, Bifurc
					    /* root->left->left->left */ 0, 15, -2,    // Null,  aatt, End
					    /* root->left->left->right */ 0, 15, -2,   // Null,  ccgg, End
					    /* root->left->right */ 0, 12, 10, 3, 17,  // Null,  -gc-, --tt, c--g, BifurcRevY
					    /* root->left->right->left */ 0, 15, -2,   // Null,  ttaa, End
					    /* root->left->right->right */ 0, 15, -2,  // Null,  ggcc, End
					    /* root->right */ 0, 12, 2, 16,            // Null,  -cg-, ---t, Bifurc
					    /* root->right->left */ 0, 15, -2,         // Null,  atat, End
					    /* root->right->right */ 0, 15, -2 };      // Null,  gcgc, End

      const Score unconstrained_emit_sc =
	P + 0 + 2*P + 2*P  // tcga, aatt, cc-g, g-cc
	+ 2*P + 4*P + 4*P  // aat-, gg--, g-c-
	+ 2*P + 4*P + 6*P  // -cgg, -t-a, g---
	+ 0                // aatt
	+ 0                // ccgg
	+ 4*P + 4*P + 4*P  // -gc-, --tt, c--g
	+ 0                // ttaa
	+ 0                // ggcc
	+ 4*P + 6*P        // -cg-, ---t
	+ P                // atat
	+ P;               // gcgc

      cerr << " (emit Scores: constrained " << constrained_emit_sc << ", unconstrained " << unconstrained_emit_sc << ")\n";

      vector<int> constrained_state_path (constrained_state_array, constrained_state_array + 44);
      vector<int> unconstrained_state_path (unconstrained_state_array, unconstrained_state_array + 43);
      
      Digitized_biosequence& xdsq = npx.dsq;
      Digitized_biosequence& ydsq = npy.dsq;

      DNA_alphabet.seq2dsq (xseq, xdsq);
      DNA_alphabet.seq2dsq (yseq, ydsq);

      Pair_CFG_scores cfg (18);
      cfg.start_to_end() = 9999;
      for (int i = 0; i < 18; ++i)
	{
	  cfg.start[i] = -9999;
	  cfg.end[i] = -9999;
	  cfg.init_emit (i, (Pair_CFG_scores::State_type) (i < 17 ? i : 32));
	  const Pair_CFG_scores::State_type t = cfg.state_type[i];
	  if (Pair_CFG_scores::is_emit_type (t))
	    for (int cxl = 0; cxl < 4; ++cxl)
	      for (int cxr = 0; cxr < 4; ++cxr)
		for (int cyl = 0; cyl < 4; ++cyl)
		  for (int cyr = 0; cyr < 4; ++cyr)
		    {
		      const int emit_idx =
			cxl * cfg.emit_xl_mul (t) +
			cxr * cfg.emit_xr_mul (t) +
			cyl * cfg.emit_yl_mul (t) +
			cyr * cfg.emit_yr_mul (t);
		      Score emit_sc = 0;
		      // eliminate mispairings & mismatches
		      const bool xlr = (t & Pair_CFG_scores::EmitXL) && (t & Pair_CFG_scores::EmitXR);
		      const bool ylr = (t & Pair_CFG_scores::EmitYL) && (t & Pair_CFG_scores::EmitYR);
		      const bool xyl = (t & Pair_CFG_scores::EmitXL) && (t & Pair_CFG_scores::EmitYL);
		      const bool xyr = (t & Pair_CFG_scores::EmitXR) && (t & Pair_CFG_scores::EmitYR);
		      const bool xlyr = (t & Pair_CFG_scores::EmitXL) && (t & Pair_CFG_scores::EmitYR);
		      const bool ylxr = (t & Pair_CFG_scores::EmitYL) && (t & Pair_CFG_scores::EmitXR);
		      bool mismatch = 0;
		      if (xlr && cxl != 3 - cxr) mismatch = 1;
		      if (ylr && cyl != 3 - cyr) mismatch = 1;
		      if (xyl && !(xlr && ylr) && cxl != cyl) mismatch = 1;
		      if (xyr && !(xlr && ylr) && cxr != cyr) mismatch = 1;
		      if (xlyr && !(xlr || ylr) && cxl != 3 - cyr) mismatch = 1;
		      if (ylxr && !(xlr || ylr) && cyl != 3 - cxr) mismatch = 1;
		      if (mismatch)
			emit_sc = -InfinityScore;
		      else
			{
			  int penalty = 0;
			  for (int l = 0; l < 4; ++l) if (!(t & (1 << l))) penalty += 2;
			  if ((xyl && cxl != cyl) || (xyr && cxr != cyr)) ++penalty;
			  emit_sc = -1000 * penalty;
			}
		      cfg.emit[i][emit_idx] = emit_sc;
		    }
	  if (i < 16)
	    for (int j = 0; j < 18; ++j)
	      if (i > 0 || (j > 0 && j < 16))  // no direct transitions between null states
		cfg.transition (i, j) = -9999;
	}
      cfg.bifurc[16] = cfg.bifurc[17] = Pair_CFG_scores::Bifurcation (0, 0);

      // set up transitions to favour constrained path
      cfg.start[15] = 0;
      cfg.transition (15, 15) = 0;
      cfg.transition (15,  7) = 0;
      cfg.transition ( 7, 11) = 0;
      cfg.transition (11, 16) = 0;
      cfg.transition ( 0, 13) = 0;
      cfg.transition (13,  5) = 0;
      cfg.transition ( 5,  9) = 0;
      cfg.transition ( 9, 17) = 0;
      cfg.transition ( 0, 14) = 0;
      cfg.transition (14,  6) = 0;
      cfg.transition ( 6,  1) = 0;
      cfg.transition ( 1, 16) = 0;
      cfg.transition ( 0, 15) = 0;
      cfg.end[15] = 0;
      cfg.transition ( 0, 12) = 0;
      cfg.transition (12, 10) = 0;
      cfg.transition (10,  3) = 0;
      cfg.transition ( 3, 17) = 0;
      cfg.transition ( 0,  4) = 0;
      cfg.transition ( 4,  8) = 0;
      cfg.transition ( 8,  2) = 0;
      cfg.transition ( 2, 16) = 0;

      // add some extra transitions to make things different for unconstrained path (& generally confuse things)
      cfg.transition (12,  2) = 0;  // this lets the unconstrained CYK algorithm pair up a yseq CG that the constrained algorithm can't
      cfg.transition ( 3, 16) = 0;  // to check that BifurcRevY really is appealing
      cfg.transition (15, 15) = 0;  // for the easy test

      {
	cerr << "(co-folding '" << xseq << "' and '" << yseq << "' using full envelopes)\n";

#ifndef DART_DEBUG  /* this test is kinda slow, so leave it out during debugging runs */

	Fold_envelope xenv;
	Fold_envelope yenv;
	
	xenv.initialise_full (xdsq.size());
	yenv.initialise_full (ydsq.size());
	
	cerr << "(testing CYK algorithm... please be patient)\n";

	Pair_CYK_matrix cyk (npx, npy, xenv, yenv, cfg);
	const Score cyk_score = cyk.final_score;
	
	const vector<int> cyk_traceback = cyk.traceback();
	const Pair_CFG_parse_tree cyk_tree = cfg.parse (cyk_traceback);
	
	const Score cyk_tree_score = cfg.path_score (cyk_tree, xdsq, ydsq);
	if (!scores_equal (cyk_tree_score, cyk_score))
	  THROWEXPR ("CFG traceback score " << cyk_tree_score << " != matrix score " << cyk_score);
	
	if (!cfg.test_branch_coords_consistent (cyk_tree)) THROWEXPR ("Branch coords inconsistent");
	if (!cyk_tree.test_global (xdsq, ydsq)) THROWEXPR ("Not global");

	if (cyk_traceback != unconstrained_state_path)
	  {
	    cerr << "cyk_score=" << cyk_score << "\n";
	    cerr << "traceback (" << cyk_traceback << ")\n";
	    THROWEXPR ("Unconstrained CYK path != test path");
	  }

#else  /* defined DART_DEBUG */
	cerr << "(...skipping this test 'cos DART_DEBUG was #define'd...)\n";
#endif /* DART_DEBUG */
      }

      {
	cerr << "(now trying using envelopes '" << xfold << "' and '" << yfold << "')\n";
	
	Fold_envelope xenv;
	Fold_envelope yenv;
	
	xenv.initialise_from_fold_string (xfold);
	yenv.initialise_from_fold_string (yfold);

	cerr << "(testing CYK algorithm...)\n";

	Pair_CYK_matrix cyk (npx, npy, xenv, yenv, cfg);
	const Score cyk_score = cyk.final_score;
	
	const vector<int> cyk_traceback = cyk.traceback();
	const Pair_CFG_parse_tree cyk_tree = cfg.parse (cyk_traceback);
	
	const Score cyk_tree_score = cfg.path_score (cyk_tree, xdsq, ydsq);
	if (!scores_equal (cyk_tree_score, cyk_score))
	  THROWEXPR ("CFG traceback score " << cyk_tree_score << " != matrix score " << cyk_score);
	
	if (!cfg.test_branch_coords_consistent (cyk_tree)) THROWEXPR ("Branch coords inconsistent");
	if (!cyk_tree.test_global (xdsq, ydsq)) THROWEXPR ("Not global");

	if (cyk_traceback != constrained_state_path)
	  {
	    cerr << "cyk_score=" << cyk_score << "\n";
	    cerr << "traceback (" << cyk_traceback << ")\n";
	    THROWEXPR ("Constrained CYK path != test path");
	  }

	cerr << "(testing inside-outside algorithm...)\n";

	Named_profile easyx;
	Named_profile easyy;

	Biosequence& easyxseq = easyx.seq;
	Biosequence& easyyseq = easyy.seq;

	easyxseq          = "atat";
	easyyseq          = "cgcg";
	sstring easyfold  = "<><>";
	cerr << "(trying with '" << easyxseq << "' and '" << easyyseq << "' first)\n";

	Digitized_biosequence& easyxdsq = easyx.dsq;
	Digitized_biosequence& easyydsq = easyy.dsq;
	
	DNA_alphabet.seq2dsq (easyxseq, easyxdsq);
	DNA_alphabet.seq2dsq (easyyseq, easyydsq);

	Fold_envelope easyenv;
	easyenv.initialise_from_fold_string (easyfold);
	
	Pair_inside_outside_matrix in_out (easyx, easyy, easyenv, easyenv, cfg);
	//	in_out.show (cerr);

	Pair_CFG_counts easycount (cfg);
	easycount.start[16] = easycount.start[17] = 0.5;
	easycount.transition (0, 15) = 2;
	easycount.emit[15] [0 + 3*4 + 1*16 + 2*64] = 2;
	easycount.end[15] = 2;

	test_cfg_counts_equal (in_out.count, easycount);

	cerr << "(now trying the longer sequences)\n";

	Pair_CFG_counts cfg_count (cfg);
	
	cfg_count.emit_by_chars (15, "tacg") = 1;
	cfg_count.emit_by_chars (15, "atat") = 2;
	cfg_count.emit_by_chars (7,  "cgc")  = 1;
	cfg_count.emit_by_chars (11, "gcc")  = 1;
	cfg_count.emit_by_chars (13, "aat")  = 1;
	cfg_count.emit_by_chars (5,  "gg")   = 1;
	cfg_count.emit_by_chars (9,  "gc")   = 1;
	cfg_count.emit_by_chars (14, "gcg")  = 1;
	cfg_count.emit_by_chars (6,  "at")   = 1;
	cfg_count.emit_by_chars (1,  "g")    = 1;
	cfg_count.emit_by_chars (15, "cgcg") = 1;
	cfg_count.emit_by_chars (12, "gc")   = 1;
	cfg_count.emit_by_chars (10, "tt")   = 1;
	cfg_count.emit_by_chars (3,  "cg")   = 1;
	cfg_count.emit_by_chars (15, "tata") = .8;  // NB extra transition
	cfg_count.emit_by_chars (15, "gcgc") = .8;  // ---see comment for
	cfg_count.emit_by_chars (15, "tagc") = .2;  // transition(3,16) and
	cfg_count.emit_by_chars (15, "gcta") = .2;  // transition(3,17) below
	cfg_count.emit_by_chars (4,  "c")    = 1;
	cfg_count.emit_by_chars (8,  "g")    = 1;
	cfg_count.emit_by_chars (2,  "t")    = 1;
	cfg_count.emit_by_chars (15, "atta") = 1;
	cfg_count.emit_by_chars (15, "gccg") = 1;
	
	cfg_count.start[15] = 1;
	cfg_count.transition (15, 15) = 1;
	cfg_count.transition (15,  7) = 1;
	cfg_count.transition ( 7, 11) = 1;
	cfg_count.transition (11, 16) = 1;
	cfg_count.transition ( 0, 13) = 1;
	cfg_count.transition (13,  5) = 1;
	cfg_count.transition ( 5,  9) = 1;
	cfg_count.transition ( 9, 17) = 1;
	cfg_count.transition ( 0, 14) = 1;
	cfg_count.transition (14,  6) = 1;
	cfg_count.transition ( 6,  1) = 1;
	cfg_count.transition ( 1, 16) = 1;
	cfg_count.transition ( 0, 15) = 6;
	cfg_count.end[15] = 6;
	cfg_count.transition ( 0, 12) = 1;
	cfg_count.transition (12, 10) = 1;
	cfg_count.transition (10,  3) = 1;
	cfg_count.transition ( 3, 17) = .8;  // NB extra transition
	cfg_count.transition ( 3, 16) = .2;  // probability ratio is (1*1) / (.5*.5) due to mismatches
	cfg_count.transition ( 0,  4) = 1;
	cfg_count.transition ( 4,  8) = 1;
	cfg_count.transition ( 8,  2) = 1;
	cfg_count.transition ( 2, 16) = 1;

	Pair_inside_outside_matrix in_out2 (npx, npy, xenv, yenv, cfg);
	//	in_out2.count.show (cerr);
	test_cfg_counts_equal (in_out2.count, cfg_count);
	
      }

      // generate paths through Single HMMs & compare path scores

      cerr << "(generating some pseudorandom HMMs & comparing forward/Viterbi/fwdback to CFG inside/CYK/inout scores & counts)\n";
      cerr << "(score delta tolerance " << tol << "%)\n";

      pseudseed(27815386);
      test_single_emit (0, 1, 1, Bits2Score(1), 2, 1);   // one emit state, short sequence

      pseudseed(681143914);
      test_single_emit (0, 1, 1, Bits2Score(1), 100, 10);   // one emit state, longer sequences

      pseudseed(1599793035);
      test_single_emit (0, 5, 4*5, Bits2Score(3), 100, 10);  // several emit states

      pseudseed(867130349);
      test_single_emit (5, 10, 15*15, Bits2Score(6), 20, 10);  // several emit & null states

      pseudseed(1645288491);
      test_single_emit (5, 10, 15*15, Bits2Score(6), 20, 10, 0);  // allow zero-probability paths

      pseudseed(40122798);
      test_single_emit (0, 0, 0, Bits2Score(3), 0, 1);  // zero path length, no states at all
      
      pseudseed(221266094);
      test_single_emit (3, 3, 3*3, Bits2Score(3), 0, 1);  // zero path length, but with states

      // generate paths through Pair HMMs & compare path scores

      pseudseed(1427370996);
      test_pair_emit (1, 0, 0, 1, Bits2Score(1), 20, 0, 5);   // one EmitX state

      pseudseed(1917438203);
      test_pair_emit (0, 1, 0, 1, Bits2Score(1), 0, 20, 5);   // one EmitY state
      
      pseudseed(1701969387);
      test_pair_emit (0, 0, 1, 1, Bits2Score(1), 20, 20, 5);   // one EmitXY state

      pseudseed(795225163);
      test_pair_emit (1, 1, 1, 3*3*2, Bits2Score(5), 20, 20, 5);   // one state of each type

      pseudseed(311239244);
      test_pair_emit (5, 5, 5, 15*15*2, Bits2Score(10), 20, 20, 5);   // five states of each type

      cerr << "(all done-- what a superb piece of DP!)\n";
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
