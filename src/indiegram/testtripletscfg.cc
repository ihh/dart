#include "indiegram/tripletscfg.h"

int main (int argc, char** argv)
{

  try
    {
      // create a test parse tree
      cout << "test out\n";
      cerr << "test err\n";
      cerr << "(testing a simple triplet SCFG parse)";

      // Create a toy grammar.
      // It has "stem" emit states which emit only WC-paired nucleotides and "loop" emit states which emit to 
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
	cout << "Created state " << i << " of type " << t << " (" << SCFG_state_typing::state_type_string (t) << ")\n";

	// set transition scores
	scfg.transition_scores.start[i] = 9900 + i;
	scfg.transition_scores.end[i] = i * 100 + 99;
	for (int j = 0; j < num_states; ++j)
	  scfg.transition_scores.transition (i, j) = i * 100 + j;

	//cout << "state " << i << " (type " << t << ")\n";
	// if an emit state, initialize scores for emissions
	if (SCFG_state_typing::is_emit_type (t))
	  {
	    for (int xlc = 0; xlc < SCFG_alphabet_size; ++xlc)
	      for (int xrc = 0; xrc < SCFG_alphabet_size; ++xrc)
		for (int ylc = 0; ylc < SCFG_alphabet_size; ++ylc)
		  for (int yrc = 0; yrc < SCFG_alphabet_size; ++yrc)
		    for (int zlc = 0; zlc < SCFG_alphabet_size; ++zlc)
		      for (int zrc = 0; zrc < SCFG_alphabet_size; ++zrc)
			// Now we need to assign a unique score to each possible emission from every emit state
			// We have designed our test sequences so that paired emissions are always WC-basepaired,
			// and so it's sufficient to assign scores only to basepaired paired emissions.
			// Note that this does properly assign scores to single (unpaired) emissions as well.
			if (xlc + xrc == 3 && ylc + yrc == 3 && zlc + zrc == 3) // if basepaired
			  {
 			    const int emit_idx = scfg.emit_idx (t, xlc, xrc, ylc, yrc, zlc, zrc);
			    int emit_code = 0; // the associated score
			    emit_code += ((t & SCFG_state_typing::EmitXL) ? (xlc + 1) : 9) * 100000;
			    emit_code += ((t & SCFG_state_typing::EmitYL) ? (ylc + 1) : 9) * 10000;
			    emit_code += ((t & SCFG_state_typing::EmitZL) ? (zlc + 1) : 9) * 1000;
			    emit_code += ((t & SCFG_state_typing::EmitZR) ? (zrc + 1) : 9) * 100;
			    emit_code += ((t & SCFG_state_typing::EmitYR) ? (yrc + 1) : 9) * 10;
			    emit_code += ((t & SCFG_state_typing::EmitXR) ? (xrc + 1) : 9) * 1;
			    scfg.emit[i][emit_idx] = emit_code;
			  }
	  }

      }

      // now deal with the bifurcation state
      scfg.init_bifurc (64, 0, 0);
      cout << "Created state " << 64 << " of type " << scfg.state_type[64] << " (" << SCFG_state_typing::state_type_string (scfg.state_type[64]) << ")\n";

      // set up the sequences
      Biosequence xseq = "aaggaaattttt";
      // xfold:           <<..<<<>>>>>
      Biosequence yseq = "aaaaattttt";
      // yfold:           <<<<<>>>>>
      Biosequence zseq = "aaatggaatttt";
      // zfold:           <<<>..<<>>>>

      // set up Named_profile's
      Named_profile np_x, np_y, np_z;
      np_x.name = "X"; np_y.name = "Y"; np_z.name = "Z"; 
      np_x.seq = xseq; np_y.seq = yseq; np_z.seq = zseq;
      DNA_alphabet.seq2dsq (np_x.seq, np_x.dsq); DNA_alphabet.seq2dsq (np_y.seq, np_y.dsq); DNA_alphabet.seq2dsq (np_z.seq, np_z.dsq);

      int state_array[15] = { /* root */ -1, 63, 63, 32, 32, 64,   // Start, (tct aga), (aaa ttt), (ccc g-g), (
			      /* root->left */ 0, 55, 55, 63, -2,  //
			      /* root->right */ 0, 12, 12, -2 };
      vector<int> state_path (state_array, state_array + 15);

      // scfg.show (cout);

      // parse tree
      cout << "State array used to build parse tree via scfg.parse (state_path):\n";
      for (int i = 0; i < (int) state_path.size(); ++i)
	cout << state_path[i] << " ";
      Triplet_SCFG_parse_tree parse_tree = scfg.parse (state_path);
      cout << "\nParse tree:\n";
      parse_tree.show (cout, &(scfg.state_type));

      // alignment
      Triplet_SCFG_alignment alignment = parse_tree.alignment (scfg.state_type, np_x, np_y, np_z);
      cout << "\nAlignment:\n";
      alignment.show (cout);

      // scoring
      cout << "\nScoring:\n";

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
      dump_scores = true;
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


      // test scoring: this implicitly tests the parse tree, etc. as well
      Score trans_score = (9963 + 6363 + 6332 + 3232 + 3264) + (0 + 55 + 5555 + 5563 + 6399) + (0 + 12 + 1212 + 1299); // note 0 score associated with 64 -> 0 (bifurcation -> null)
      Score test_trans_score = scfg.path_transition_score (parse_tree);
      cout << "\ntransition score = " << trans_score << "\n";
      cout << "computed transition score = " << test_trans_score << "\n";
      if (trans_score != test_trans_score) THROWEXPR ("SCFG path transition score " << trans_score << " != expected score " << test_trans_score);

      Score emit_score = (111444 + 111444 + 399999 + 399999) + (119344 + 119344 + 111444) + (991499 + 991499);
      Score test_emit_score = scfg.path_emit_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq);
      cout << "\nemit score = " << emit_score << "\n";
      cout << "computed emit score = " << test_emit_score << "\n";
      if (emit_score != test_emit_score) THROWEXPR ("SCFG path emit score " << emit_score << " != expected score " << test_emit_score);

      Score total_score = ScorePMul (trans_score, emit_score);
      Score test_total_score = scfg.path_score (parse_tree, np_x.dsq, np_y.dsq, np_z.dsq);
      cout << "\ntotal score = " << total_score << "\n";
      cout << "computed total score = " << test_total_score << "\n";
      if (total_score != test_total_score) THROWEXPR ("SCFG path total score " << total_score << " != expected score " << test_total_score);




      // display sequences and folds
//      cout << "np_x: " << alignment.np_x.name << ": " << alignment.np_x.seq << "\n";
//      cout << "foldstring_x: " << alignment.foldstring_x << "\n";
//      cout << "foldstring_ungapped_x(): " << alignment.foldstring_ungapped_x() << "\n";
//      cout << "np_y: " << alignment.np_y.name << ": " << alignment.np_y.seq << "\n";
//      cout << "foldstring_y: " << alignment.foldstring_y << "\n";
//      cout << "foldstring_ungapped_y(): " << alignment.foldstring_ungapped_y() << "\n";
//      cout << "np_z: " << alignment.np_z.name << ": " << alignment.np_z.seq << "\n";
//      cout << "foldstring_z: " << alignment.foldstring_z << "\n";
//      cout << "foldstring_ungapped_z(): " << alignment.foldstring_ungapped_z() << "\n";




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
