#include "util/score.h"
#include "scfg/foldenv.h"
#include "scfg/paircfg.h"
#include "scfg/paircfgdp.h"
#include "indiegram/tripletscfgdp.h"
#include "indiegram/postprobs.h"

const double tol = 0.1;  // % error tolerance

// score comparison that allows for approx infinities
bool scores_equal (Score a, Score b)
{
  const double a_inf = 100.0 * ((double) (a + InfinityScore)) / (double) InfinityScore;
  const double b_inf = 100.0 * ((double) (b + InfinityScore)) / (double) InfinityScore;
  if (a_inf < tol) return b_inf < tol;
  return a == b;
}

// test DP against a Pair SCFG
void test_against_paircfg()
{
  
  cerr << "(testing against a Pair_CFG)\n";
  cerr << "(initializing Pair_CFG_scores object)\n";
  // this is just the Pair_CFG_scores from scfg/testpaircfgdp.cc with the BifurcRevY state (17) removed
  int num_states = 16 + 1; // 15 Emit, 1 Null and 1 Bifurc
  Pair_CFG_scores paircfg (num_states);
  paircfg.start_to_end() = 9999;
  for (int i = 0; i < num_states; ++i)
    {
      paircfg.start[i] = -9999;
      paircfg.end[i] = -9999;
      paircfg.init_emit (i, (Pair_CFG_scores::State_type) i);
      const Pair_CFG_scores::State_type t = paircfg.state_type[i];
      if (Pair_CFG_scores::is_emit_type (t))
	for (int cxl = 0; cxl < 4; ++cxl)
	  for (int cxr = 0; cxr < 4; ++cxr)
	    for (int cyl = 0; cyl < 4; ++cyl)
	      for (int cyr = 0; cyr < 4; ++cyr)
		{
		  const int emit_idx =
		    cxl * paircfg.emit_xl_mul (t) +
		    cxr * paircfg.emit_xr_mul (t) +
		    cyl * paircfg.emit_yl_mul (t) +
		    cyr * paircfg.emit_yr_mul (t);
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
		  paircfg.emit[i][emit_idx] = emit_sc;
		}
      if (i < 16)
	for (int j = 0; j < num_states; ++j)
	  if (i > 0 || (j > 0 && j < 16))  // no direct transitions between null states
	    paircfg.transition (i, j) = -9999;
    }

  // init bifurc states
  paircfg.bifurc[16] = Pair_CFG_scores::Bifurcation (0, 0);

  // set up transitions to favour constrained path
  paircfg.start[15] = 0;
  paircfg.transition (15, 15) = 0;
  paircfg.transition (15,  7) = 0;
  paircfg.transition ( 7, 11) = 0;
  paircfg.transition (11, 16) = 0;
  paircfg.transition ( 0, 13) = 0;
  paircfg.transition (13,  5) = 0;
  paircfg.transition ( 5,  9) = 0;
  paircfg.transition ( 0, 14) = 0;
  paircfg.transition (14,  6) = 0;
  paircfg.transition ( 6,  1) = 0;
  paircfg.transition ( 1, 16) = 0;
  paircfg.transition ( 0, 15) = 0;
  paircfg.end[15] = 0;
  paircfg.transition ( 0, 12) = 0;
  paircfg.transition (12, 10) = 0;
  paircfg.transition (10,  3) = 0;
  paircfg.transition ( 0,  4) = 0;
  paircfg.transition ( 4,  8) = 0;
  paircfg.transition ( 8,  2) = 0;
  paircfg.transition ( 2, 16) = 0;

  // now initialize a Triplet_SCFG from this Pair_CFG_scores
  cerr << "(initializing Triplet_SCFG from the Pair_CFG_scores object)\n";
  Triplet_SCFG triplet (paircfg);

  // set up the sequences
  Biosequence xseq = "accagatggt";
  Biosequence yseq = "acgtgaaatt";
  Biosequence zseq = "";

  // set up Named_profile's
  Named_profile np_x, np_y, np_z;
  np_x.name = "X"; np_y.name = "Y";
  np_x.seq = xseq; np_y.seq = yseq;
  DNA_alphabet.seq2dsq (np_x.seq, np_x.dsq); DNA_alphabet.seq2dsq (np_y.seq, np_y.dsq);

  // show sequences and foldstrings
  cout << "X = " << np_x.seq << "\n";
  cout << "Y = " << np_y.seq << "\n";

  // initialize fold envelopes
  Fold_envelope xenv, yenv, zenv;
  xenv.initialise_full (np_x.size());
  yenv.initialise_full (np_y.size());
  zenv.initialise_full (np_z.size());

  // test Inside
  cerr << "(testing unconstrained Inside)\n";
  Pair_inside_matrix pair_inside (np_x, np_y, xenv, yenv, paircfg, false, true);
  Triplet_inside_matrix triplet_inside (triplet, np_x, np_y, np_z, xenv, yenv, zenv, true);
  if (!scores_equal (pair_inside.final_score, triplet_inside.final_sc))
    THROWEXPR ("pair Inside score " << pair_inside.final_score << " != triplet Inside score " << triplet_inside.final_sc);
  cerr << "Inside scores match (" << pair_inside.final_score << ").\n";

  // test CYK
  cerr << "(testing unconstrained CYK)\n";
  Pair_CYK_matrix pair_CYK (np_x, np_y, xenv, yenv, paircfg, false, true);
  Triplet_CYK_matrix triplet_CYK (triplet, np_x, np_y, np_z, xenv, yenv, zenv, true);
  if (!scores_equal (pair_CYK.final_score, triplet_CYK.final_sc))
    THROWEXPR ("pair CYK score " << pair_CYK.final_score << " != triplet CYK score " << triplet_CYK.final_sc);
  cerr << "CYK scores match (" << pair_CYK.final_score << ").\n";

  // test CYK traceback
  cerr << "(testing unconstrained CYK traceback)\n";
  vector<int> pair_traceback = pair_CYK.traceback();
  vector<int> triplet_traceback = triplet_CYK.traceback();
  if (pair_traceback != triplet_traceback)
    {
      THROWEXPR ("Tracebacks don't match:\n"
		 << "pair traceback:\n(" << pair_traceback << ")\n"
		 << "triplet traceback:\n(" << triplet_traceback << ")\n");
    }
  else
    cerr << "Tracebacks match (" << pair_traceback << ").\n";

  // test Outside (compare some random cells; this is a very poor test)
  cerr << "(testing posterior probabilities and estimated counts)\n";
  Pair_outside_matrix pair_outside (pair_inside, true);
  Triplet_outside_matrix triplet_outside (triplet_inside, true);
  bool outside_agrees = true;
  for (int i = 0; i < pair_inside.cfg.states(); ++i)
    {
      int subseq_idx_x = 1;
      int subseq_idx_y = 4;
      cerr << "outside.read_cell (" << i << "," << xenv.subseq[subseq_idx_x].terser_desc() << "," << yenv.subseq[subseq_idx_y].terser_desc() << "):  "
	   << pair_outside.read_cell (i, subseq_idx_x, subseq_idx_y) << " <--> "
	   << triplet_outside.read_cell (i, subseq_idx_x, subseq_idx_y, 0) << "\n";
      if (!scores_equal (pair_outside.read_cell (i, subseq_idx_x, subseq_idx_y), triplet_outside.read_cell (i, subseq_idx_x, subseq_idx_y, 0)))
	  outside_agrees = false;
    }
  if (!outside_agrees)
    THROWEXPR ("Outside scores disagree.");

  // to do: compare posterior probs/counts here
  // test post probs
  Triplet_inside_outside_matrix triplet_in_out (triplet, np_x, np_y, np_z, xenv, yenv, zenv, true);


  cerr << "\n\n";

}

// test CYK deterministically (reproduce a best path through a scfg hacked by hand)
// Do it for both constrained and unconstrained fold envelopes.
void test_deterministically()
{

  cerr << "(testing a simple triplet SCFG parse)\n";

  // Create a toy grammar.
  // It has "stem" emit states which emit only WC-paired nucleotides as well as "loop" emit states.
  int num_null_states = 1;
  int num_emit_states = 63;
  int num_bifurc_states = 1;
  int num_states = num_null_states + num_emit_states + num_bifurc_states;

  Triplet_SCFG scfg (num_states);
  scfg.transition_scores.start_to_end() = 9999;

  // initialize states: 1 state for each state type
  for (int i = 0; i < num_states; ++i) {

    // initialize state
    const State_type t = static_cast<State_type> (i);
    scfg.init_emit (i, t); // 1 state for each state type, so let 'state' = 'type'

    // set transition scores unfavorably (we'll make some favorable later)
    scfg.transition_scores.start[i] = -999;
    scfg.transition_scores.end[i] = -999;
    for (int j = 0; j < num_states; ++j)
      scfg.transition_scores.transition (i, j) = -999;

    // if an emit state, initialize scores for emissions
    if (SCFG_state_typing::is_emit_type (t))
      {
	for (int xlc = 0; xlc < SCFG_alphabet_size; ++xlc)
	  for (int xrc = 0; xrc < SCFG_alphabet_size; ++xrc)
	    for (int ylc = 0; ylc < SCFG_alphabet_size; ++ylc)
	      for (int yrc = 0; yrc < SCFG_alphabet_size; ++yrc)
		for (int zlc = 0; zlc < SCFG_alphabet_size; ++zlc)
		  for (int zrc = 0; zrc < SCFG_alphabet_size; ++zrc)
		    {
		      Score emit_sc = 0;
		      // Now we need to assign a score to each possible emission from every emit state
		      // We have designed our test sequences and folds so that paired emissions are always WC-basepaired,
		      // and so it's sufficient to assign scores only to basepaired paired emissions.
		      // Note that this does properly assign scores to single (unpaired) emissions as well.

		      // if paired emission
		      const bool xpaired = (t & SCFG_state_typing::EmitXL) && (t & SCFG_state_typing::EmitXR);
		      const bool ypaired = (t & SCFG_state_typing::EmitYL) && (t & SCFG_state_typing::EmitYR);
		      const bool zpaired = (t & SCFG_state_typing::EmitZL) && (t & SCFG_state_typing::EmitZR);
		      // if unpaired emissions (to left or right) to > 1 sequence
		      const bool xyl = (t & SCFG_state_typing::EmitXL) && (t & SCFG_state_typing::EmitYL);
		      const bool xzl = (t & SCFG_state_typing::EmitXL) && (t & SCFG_state_typing::EmitZL);
		      const bool yzl = (t & SCFG_state_typing::EmitYL) && (t & SCFG_state_typing::EmitZL);
		      const bool xyr = (t & SCFG_state_typing::EmitXR) && (t & SCFG_state_typing::EmitYR);
		      const bool xzr = (t & SCFG_state_typing::EmitXR) && (t & SCFG_state_typing::EmitZR);
		      const bool yzr = (t & SCFG_state_typing::EmitYR) && (t & SCFG_state_typing::EmitZR);

		      // disallowed paired emissions which aren't WC-paired
		      bool notwc = 0;
		      if (xpaired && (xlc+xrc != 3)) notwc = 1;
		      if (ypaired && (ylc+yrc != 3)) notwc = 1;
		      if (zpaired && (zlc+zrc != 3)) notwc = 1;
		      if (notwc)
			emit_sc = -InfinityScore;

		      // disallow mismatches for paired emissions
		      bool mismatch = 0;
		      if ((xpaired && ypaired) && (xlc != ylc || xrc != yrc)) mismatch = 1;
		      if ((xpaired && zpaired) && (xlc != zlc || xrc != zrc)) mismatch = 1;
		      if ((ypaired && zpaired) && (ylc != zlc || yrc != zrc)) mismatch = 1;

		      // disallow mismatches for single emissions
		      if (xyl && !(xpaired || ypaired) && (xlc != ylc)) mismatch = 1;
		      if (xzl && !(xpaired || zpaired) && (xlc != zlc)) mismatch = 1;
		      if (yzl && !(ypaired || zpaired) && (ylc != zlc)) mismatch = 1;
		      if (xyr && !(xpaired || ypaired) && (xrc != yrc)) mismatch = 1;
		      if (xzr && !(xpaired || zpaired) && (xrc != zrc)) mismatch = 1;
		      if (yzr && !(ypaired || zpaired) && (yrc != zrc)) mismatch = 1;

		      if (mismatch)
			emit_sc = -InfinityScore;

		      // randomized emission scores for the rest
		      else
			{
			  // penalty for mismatches
			  int penalty = 0;
			  if (xyl && (xlc != ylc)) ++penalty;
			  if (xzl && (xlc != zlc)) ++penalty;
			  if (yzl && (ylc != zlc)) ++penalty;
			  if (xyr && (xrc != yrc)) ++penalty;
			  if (xzr && (xrc != zrc)) ++penalty;
			  if (yzr && (yrc != zrc)) ++penalty;
			  // prefer paired emissions over unpaired emissions
			  if (xpaired || ypaired || zpaired)
			    emit_sc = -100 * penalty;
			  else
			    emit_sc = -1000 * penalty;
			}

		      // now set the score
		      const int emit_idx = scfg.emit_idx (t, xlc, xrc, ylc, yrc, zlc, zrc);
		      scfg.emit[i][emit_idx] = emit_sc;
		    }
      }

  }
  // deal with the bifurcation state
  scfg.init_bifurc (64, 0, 0);

  // set some transition scores by hand to favor the parse we want
  // (see state_array[], defined below)
  scfg.transition_scores.start[63] = 0; // root
  scfg.transition_scores.transition (63,63) = 0;
  scfg.transition_scores.transition (63,32) = 0;
  scfg.transition_scores.transition (32,32) = 0;
  scfg.transition_scores.transition (32,64) = 0;
  scfg.transition_scores.transition (0,55) = 0; // root->left
  scfg.transition_scores.transition (55,55) = 0;
  scfg.transition_scores.transition (55,63) = 0;
  scfg.transition_scores.end[63] = 0;
  scfg.transition_scores.transition (0,12) = 0; // root->right
  scfg.transition_scores.transition (12,12) = 0;
  scfg.transition_scores.end[12] = 0;


  // set up the sequences and foldstrings for the constrained test
  Biosequence xseq = "aaggaaattttt";
  sstring xfold    = "<<..<<<>>>>>";
  Biosequence yseq = "aaaaattttt";
  sstring yfold    = "<<<<<>>>>>";
  Biosequence zseq = "aaatggaatttt";
  sstring zfold    = "<<<>..<<>>>>";

  // set up Named_profile's
  Named_profile np_x, np_y, np_z;
  np_x.name = "X"; np_y.name = "Y"; np_z.name = "Z"; 
  np_x.seq = xseq; np_y.seq = yseq; np_z.seq = zseq;
  DNA_alphabet.seq2dsq (np_x.seq, np_x.dsq); DNA_alphabet.seq2dsq (np_y.seq, np_y.dsq); DNA_alphabet.seq2dsq (np_z.seq, np_z.dsq);

  // show sequences and foldstrings
  cout << "X = " << np_x.seq << "\n";
  cout << "    " << xfold << "\n";
  cout << "Y = " << np_y.seq << "\n";
  cout << "    " << yfold << "\n";
  cout << "Z = " << np_z.seq << "\n";
  cout << "    " << zfold << "\n";

  // this is the state path which we want cyk to return
  int state_array[15] = { /* root */ -1, 63, 63, 32, 32, 64,   // Start, (tct aga), (aaa ttt), (ccc g-g), (
			  /* root->left */ 0, 55, 55, 63, -2,  //
			  /* root->right */ 0, 12, 12, -2 };
  vector<int> state_path (state_array, state_array + 15);
  // get the corresponding parse tree and alignment
  Triplet_SCFG_parse_tree parse_tree = scfg.parse (state_path);
  Triplet_SCFG_alignment alignment = parse_tree.alignment (scfg.state_type, np_x, np_y, np_z);


  // The below is for by-hand score calculation, nothing more.
  bool dump_scores = false;
  set<int> states_used;
  for (int i = 0; i < (int) state_path.size(); ++i)
    {
      if (states_used.find (state_path[i]) != states_used.end()) continue;
      states_used.insert (state_path[i]);
    }
  // transitions
  if (dump_scores)
    {
      // display the relevant transition scores
      cout << "\nTransition scores:\n";
      for (int i = 0; i < (int) state_path.size()-1; ++i)
	{
	  int s = state_path[i];
	  int sp = state_path[i+1];

	  if (s == -2)
	    continue;

	  Score score;
	  if (SCFG_state_typing::is_bifurc_type (scfg.state_type[s]))
	    score = 0;
	  else if (s == -1)
	    {
	      if (sp == -2)
		score = scfg.transition_scores.start_to_end();		    
	      else
		score = scfg.transition_scores.start[sp];
	    }
	  else
	    {
	      if (sp == -2)
		score = scfg.transition_scores.end[s];
	      else
		score = scfg.transition_scores.transition (s, sp);
	    }
	  cout << s << " (" << SCFG_state_typing::state_type_string (scfg.state_type[s]) << ") -> " << sp << " (" << SCFG_state_typing::state_type_string (scfg.state_type[sp]) << ") = " << score << "\n";
	}
    }
  // emissions
  dump_scores = false;
  if (dump_scores)
    {
      // display the relevant emission scores
      cout << "\nEmission scores:\n";
      for (set<int>::iterator iter = states_used.begin(); iter != states_used.end(); ++iter)
	{
	  int s = *iter;
	  State_type t = scfg.state_type[s];
	  cout << "\nstate " << s << " (type " << SCFG_state_typing::state_type_string (t) << "):\n";
	  // for all possible emissions x of state s
	  for (int x = 0; x < (int) scfg.emit[s].size(); x++)
	    {
	      cout << scfg.emit_hash_to_string (s, x) << " = ";
	      scfg.show_element (scfg.emit[s][x], cout);
	      cout << ", ";
	    }
	}
      cout << "\n\n";
    }

  /*******************************************
   * Test constrained CYK (with foldstrings) *
   *******************************************/
  {
    Fold_envelope xenv, yenv, zenv;
    xenv.initialise_from_fold_string (xfold);
    yenv.initialise_from_fold_string (yfold);
    zenv.initialise_from_fold_string (zfold);

    cerr << "(testing constrained CYK)\n";
    Triplet_CYK_matrix cyk (scfg, np_x, np_y, np_z, xenv, yenv, zenv, true);
    Score cyk_sc = cyk.final_sc;

    // do traceback and get traceback parse tree
    vector<int> cyk_traceback = cyk.traceback();
    Triplet_SCFG_parse_tree cyk_tree = scfg.parse (cyk_traceback);
    // misc checks on cyk parse tree
    if (!cyk_tree.test_connections()) { cyk_tree.show(cout, &(scfg.state_type)); THROWEXPR ("Parse tree not properly connected"); }
    if (!cyk_tree.test_global (np_x.dsq, np_y.dsq, np_z.dsq)) THROWEXPR ("Parse tree not global");

    // sanity check: compare cyk score and cyk parse tree score to test traceback for self-consistency
    Score cyk_tree_sc = scfg.path_score (cyk_tree, np_x.dsq, np_y.dsq, np_z.dsq);
    if (!scores_equal (cyk_tree_sc, cyk_sc))
      THROWEXPR ("SCFG traceback score " << cyk_tree_sc << " != CYK matrix score " << cyk_sc);

    // test cyk scores (check against by-hand calculation)
    Score trans_sc = (0 + 0 + 0 + 0 + 0) + (0 + 0 + 0 + 0) + (0 + 0 + 0);
    Score test_trans_sc = scfg.path_transition_score (cyk_tree);
    cout << "\ntransition score = " << trans_sc << "\n";
    cout << "computed transition score = " << test_trans_sc << "\n";
    if (trans_sc != test_trans_sc) THROWEXPR ("SCFG path transition score " << trans_sc << " != expected score " << test_trans_sc);

    Score emit_sc = (0 + 0 + 0 + 0) + (-200 + -200 + 0) + (0 + 0);
    Score test_emit_sc = scfg.path_emit_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq);
    cout << "\nemit score = " << emit_sc << "\n";
    cout << "computed emit score = " << test_emit_sc << "\n";
    if (emit_sc != test_emit_sc) THROWEXPR ("SCFG path emit score " << emit_sc << " != expected score " << test_emit_sc);

    Score total_sc = ScorePMul (trans_sc, emit_sc);
    Score test_total_sc = scfg.path_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq);
    cout << "\ntotal score = " << total_sc << "\n";
    cout << "computed total score = " << test_total_sc << "\n";
    if (total_sc != test_total_sc) THROWEXPR ("SCFG path total score " << total_sc << " != expected score " << test_total_sc);


    // for comparison, show the parse tree we want to produce as well as the cyk tree
    cout << "\nParse tree:\n";
    parse_tree.show (cout, &(scfg.state_type));
    cout << "\nParse tree found by CYK:\n";
    cyk_tree.show (cout, &(scfg.state_type));
    // and alignment
    // of course could also get this as cyk.alignment();
    Triplet_SCFG_alignment cyk_alignment = cyk_tree.alignment (scfg.state_type, np_x, np_y, np_z);
    cyk_alignment.score = cyk_sc; // store score also
    cout << "\nAlignment:\n";
    alignment.show (cout);
    cout << "\nAlignment due to CYK:\n";
    cyk_alignment.show (cout);
  }


  /*****************************************************
   * Test unconstrained CYK (with full fold envelopes) *
   * (the grammar is hacked so that the unconstrained test should produce the desired parse tree)
   *****************************************************/
  {
    Fold_envelope xenv, yenv, zenv;
    xenv.initialise_full (np_x.size());
    yenv.initialise_full (np_y.size());
    zenv.initialise_full (np_z.size());

    cerr << "(testing unconstrained CYK)\n";
    Triplet_CYK_matrix cyk (scfg, np_x, np_y, np_z, xenv, yenv, zenv, true);
    Score cyk_sc = cyk.final_sc;

    // do traceback and get traceback parse tree
    vector<int> cyk_traceback = cyk.traceback();
    Triplet_SCFG_parse_tree cyk_tree = scfg.parse (cyk_traceback);
    // misc checks on cyk parse tree
    if (!cyk_tree.test_connections()) { cyk_tree.show(cout, &(scfg.state_type)); THROWEXPR ("Parse tree not properly connected"); }
    if (!cyk_tree.test_global (np_x.dsq, np_y.dsq, np_z.dsq)) THROWEXPR ("Parse tree not global");

    // sanity check: compare cyk score and cyk parse tree score to test traceback for self-consistency
    Score cyk_tree_sc = scfg.path_score (cyk_tree, np_x.dsq, np_y.dsq, np_z.dsq);
    if (!scores_equal (cyk_tree_sc, cyk_sc))
      THROWEXPR ("SCFG traceback score " << cyk_tree_sc << " != CYK matrix score " << cyk_sc);
    cerr << "cyk_sc = " << cyk_sc << " (" << Score2Bits (cyk_sc) << "bits) \n";


    // test cyk scores (check against by-hand calculation)
    Score trans_sc = (0 + 0 + 0 + 0 + 0) + (0 + 0 + 0 + 0) + (0 + 0 + 0);
    Score test_trans_sc = scfg.path_transition_score (cyk_tree);
    cout << "\ntransition score = " << trans_sc << "\n";
    cout << "computed transition score = " << test_trans_sc << "\n";
    if (trans_sc != test_trans_sc) THROWEXPR ("SCFG path transition score " << trans_sc << " != expected score " << test_trans_sc);

    Score emit_sc = (0 + 0 + 0 + 0) + (-200 + -200 + 0) + (0 + 0);
    Score test_emit_sc = scfg.path_emit_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq);
    cout << "\nemit score = " << emit_sc << "\n";
    cout << "computed emit score = " << test_emit_sc << "\n";
    if (emit_sc != test_emit_sc) THROWEXPR ("SCFG path emit score " << emit_sc << " != expected score " << test_emit_sc);

    Score total_sc = ScorePMul (trans_sc, emit_sc);
    Score test_total_sc = scfg.path_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq);
    cout << "\ntotal score = " << total_sc << "\n";
    cout << "computed total score = " << test_total_sc << "\n";
    if (total_sc != test_total_sc) THROWEXPR ("SCFG path total score " << total_sc << " != expected score " << test_total_sc);

    // for comparison, show the parse tree we want to produce as well as the cyk tree
    cout << "\nParse tree:\n";
    parse_tree.show (cout, &(scfg.state_type));
    cout << "\nParse tree found by CYK:\n";
    cyk_tree.show (cout, &(scfg.state_type));
    // and alignment
    // of course could also get this as cyk.alignment();
    Triplet_SCFG_alignment cyk_alignment = cyk_tree.alignment (scfg.state_type, np_x, np_y, np_z);
    cyk_alignment.score = cyk_sc;
    cout << "\nAlignment:\n";
    alignment.show (cout);
    cout << "\nAlignment due to CYK:\n";
    cyk_alignment.show (cout);	

  }


  /*****************************************************
   * Test unconstrained Inside and Outside (with full fold envelopes) *
   *****************************************************/
  {
    Fold_envelope xenv, yenv, zenv;
    xenv.initialise_full (np_x.size());
    yenv.initialise_full (np_y.size());
    zenv.initialise_full (np_z.size());

    cerr << "(testing unconstrained Inside)\n";
    Triplet_inside_matrix inside (scfg, np_x, np_y, np_z, xenv, yenv, zenv, true);
    Score inside_sc = inside.final_sc;
    cerr << "Inside score = " << inside_sc << "\n";

    Triplet_outside_matrix outside (inside, true);
    Score outside_sc = 0;

    // loop over X start positions
    for (int startpos_x = np_x.size(); startpos_x >= 0; --startpos_x)
      {
	if (xenv.by_start[startpos_x].size())
	  {
	    for_const_contents (vector<int>, xenv.by_start[startpos_x], xsi)
	      {
		const int subseq_idx_x = *xsi;
		const Subseq& subseq_x = xenv.subseq[subseq_idx_x];
		if (subseq_x.len != 0)
		  continue;

		// loop over Y start positions and corresponding subseqs
		for (int startpos_y = np_y.size(); startpos_y >= 0; --startpos_y)
		  {
		    if (yenv.by_start[startpos_y].size())
		      {
			for_const_contents (vector<int>, yenv.by_start[startpos_y], ysi)
			  {
			    const int subseq_idx_y = *ysi;
			    const Subseq& subseq_y = yenv.subseq[subseq_idx_y];
			    if (subseq_y.len != 0)
			      continue;

			    // loop over Z start positions and corresponding subseqs
			    for (int startpos_z = np_z.size(); startpos_z >= 0; --startpos_z)
			      {
				if (zenv.by_start[startpos_z].size())
				  {
				    for_const_contents (vector<int>, zenv.by_start[startpos_z], zsi)
				      {
					const int subseq_idx_z = *zsi;
					const Subseq& subseq_z = zenv.subseq[subseq_idx_z];
					if (subseq_z.len != 0)
					  continue;

					for (int i = 0; i <= inside.scfg.num_states(); ++i)
					  ScorePSumAcc (outside_sc, ScorePMul (outside.read_cell (i, subseq_idx_x, subseq_idx_y, subseq_idx_z), inside.scfg.transition_scores.end[i]));
					cerr << "sum_i cell(i," << subseq_idx_x << "," << subseq_idx_y << "," << subseq_idx_z << ") * end[i] = " << outside_sc << "\n";

				      }
				  }
			      } // z loop

			  }
		      }
		  } // y loop

	      }
	  }
      }

    cerr << "outside score = " << outside_sc << "\n";

  }

}



int main (int argc, char** argv)
{

  try
    {

      cout << "test out\n";
      cerr << "test err\n";

      test_against_paircfg();

      test_deterministically();





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
