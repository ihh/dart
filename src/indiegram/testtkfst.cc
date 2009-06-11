#include "scfg/foldenv.h"
#include "indiegram/tkfst.h"
#include "indiegram/tripletscfgdp.h"
#include "util/score.h"

const double tol = 0.1;  // % error tolerance

// score comparison that allows for approx infinities
bool scores_equal (Score a, Score b)
{
  const double a_inf = 100.0 * ((double) (a + InfinityScore)) / (double) InfinityScore;
  const double b_inf = 100.0 * ((double) (b + InfinityScore)) / (double) InfinityScore;
  if (a_inf < tol) return b_inf < tol;
  return a == b;
}

int main (int argc, char** argv)
{

  try
    {
      // test CYK deterministically

      cout << "test out\n";
      cerr << "test err\n";
      cerr << "(testing a simple triplet SCFG parse)\n";

      // create a TKFST_Triplet_SCFG
      TKFST_Triplet_SCFG scfg = TKFST_Triplet_SCFG (1, 1, 1);

      // show probs
      scfg.dump_params (cerr);
      scfg.dump_trans_probs (cerr);
      scfg.dump_emit_probs (cerr);

      // set up the sequences and foldstrings for the constrained test
      Biosequence xseq = "aaccaaaaggtt";
      sstring xfold    = "<<<<....>>>>";
      Biosequence yseq = "agtcaataggct";
      sstring yfold    = "<<<<....>>>>";
      Biosequence zseq = "aaccaataggtt";
      sstring zfold    = "<<<<....>>>>";

      // set up Named_profile's
      Named_profile np_x, np_y, np_z;
      np_x.name = "X"; np_y.name = "Y"; np_z.name = "Z"; 
      np_x.seq = xseq; np_y.seq = yseq; np_z.seq = zseq;
      DNA_alphabet.seq2dsq (np_x.seq, np_x.dsq); DNA_alphabet.seq2dsq (np_y.seq, np_y.dsq); DNA_alphabet.seq2dsq (np_z.seq, np_z.dsq);

      // this is the state path which we want cyk to return
      const int num_states = 16;
      int state_array[num_states] = { /* root */ -1, 35, 205, // Start = -1, L_L_L_L = 35, BiSL_BmSL_BmSL_BmSL = 205
				      /* left */ 48, 90, 90, 90, 90, 35, // S_S_S_S = 48, IS_MS_MS_MS = 90, IS_MS_MS_MS = 90, IS_MS_MS_MS = 90, IS_MS_MS_MS = 90, L_L_L_L = 35, 
				      95, 95, 95, 95, -2, // IL_ML_ML_ML= 95, IL_ML_ML_ML= 95, IL_ML_ML_ML= 95, IL_ML_ML_ML= 95, End = -2
				      /* right */ 35, -2 };// L_L_L_L = 35, End = -2

      vector<int> state_path (state_array, state_array + num_states);
      // get the corresponding parse tree and alignment
      Triplet_SCFG_parse_tree parse_tree = scfg.parse (state_path);
      Triplet_SCFG_alignment alignment = parse_tree.alignment (scfg.state_type, scfg.state_type_ancestral_map, np_x, np_y, np_z);

      // The below is for by-hand score calculation, nothing more.
      bool dump_scores = true;
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
      // end bit for by-hand calculation
      
      /*****************************************************
       * Do unconstrained CYK (with full fold envelopes) *
       *****************************************************/
      {
	Fold_envelope xenv, yenv, zenv;
	xenv.initialise_full (np_x.size());
	yenv.initialise_full (np_y.size());
	zenv.initialise_full (np_z.size());

	cerr << "(testing unconstrained CYK)\n";
	Triplet_CYK_matrix cyk (scfg, np_x, np_y, np_z, xenv, yenv, zenv, true);

	// do traceback and get traceback parse tree
	vector<int> cyk_traceback = cyk.traceback();
	Triplet_SCFG_parse_tree cyk_tree = scfg.parse (cyk_traceback);
	// misc checks on cyk parse tree
	if (!cyk_tree.test_connections()) { cyk_tree.show(cout, &(scfg.state_type)); THROWEXPR ("Parse tree not properly connected"); }
	if (!cyk_tree.test_global (np_x.dsq, np_y.dsq, np_z.dsq)) THROWEXPR ("Parse tree not global");

	// show desired scores
	cout << "Desired scores:\n";
	cout << "transition score = " << scfg.path_transition_score (parse_tree) << "\n";
	cout << "emit score = " << scfg.path_emit_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq) << "\n";
	cout << "total score = " << scfg.path_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq) << "\n";

	// show cyk scores
	cout << "\nCYK scores:\n";
	cout << "cyk transition score = " << scfg.path_transition_score (cyk_tree) << "\n";
	cout << "cyk emit score = " << scfg.path_emit_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq) << "\n";
	cerr << "cyk total score = " << cyk.final_sc << " (" << Score2Bits (cyk.final_sc) << " bits) \n";

	// show the desired and cyk parse trees
	cout << "\nDesired parse tree:\n";
	parse_tree.show (cout, &(scfg.state_type));
	cout << "\nCYK parse tree:\n";
	cyk_tree.show (cout, &(scfg.state_type));

	// and alignments
	// of course could also get this as cyk.alignment();
	Triplet_SCFG_alignment cyk_alignment = cyk_tree.alignment (scfg.state_type, scfg.state_type_ancestral_map, np_x, np_y, np_z);
	cyk_alignment.score = cyk.final_sc;
	cout << "\nDesired alignment:\n";
	alignment.show (cout);
	cout << "\nCYK alignment:\n";
	cyk_alignment.show (cout);	

	// check that CYK found a parse at least as good as my by-hand one
	if (scfg.path_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq) > cyk.final_sc)
	  THROWEXPR ("CYK failed: the by-hand 'desired' parse tree scores better than the parse tree found by CYK.\n");
      }


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
