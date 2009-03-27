#include <math.h>
#include "scfg/paircfgdp.h"

// test_transition function
#define SCORE_TOLERANCE 3
void test_transition (const Pair_CFG_scores& cfg, int src, int dest, Prob p, bool& ok)
{
  const Score sc_expect = Prob2Score(p);
  const Score sc_actual = cfg.transition(src,dest);
  const Score diff = sc_actual - sc_expect;
  if (abs(diff) <= SCORE_TOLERANCE)
      CTAG(1,TESTELIM) << "Passed transition test: transition from " << src << " to " << dest << " was " << sc_actual << ", should be " << sc_expect << "\n";
  else
    {
      CTAG(1,TESTELIM) << "Failed transition test: transition from " << src << " to " << dest << " was " << sc_actual << ", should be " << sc_expect << "\n";
      ok = false;
    }
}

sstring path2string (const vector<int>& path, const vector<sstring>& state_name, bool include_numeric_states = true)
{
  vector<sstring> str;
  for_const_contents (vector<int>, path, s)
    if (*s == Grammar_state_enum::Start)
      str.push_back (sstring ("Start"));
    else if (*s == Grammar_state_enum::End)
      str.push_back (sstring ("End"));
    else if (*s < 0)
      str.push_back (sstring ("Undef"));
    else
      str.push_back (state_name[*s]);
  sstring result;
  result << str;
  if (include_numeric_states)
    result << " (" << path << ')';
  return result;
}

void test_path (const Pair_CFG_scores& cfg, const vector<int>& path, const char* expect, bool& ok)
{
  const sstring str = path2string (path, cfg.state_name, false);
  if (str == expect)
      CTAG(1,TESTELIM) << "Passed path test: path was " << str << ", as expected\n";
  else
    {
      CTAG(1,TESTELIM) << "Failed transition test: path was " << str << ", expected " << expect << '\n';
      ok = false;
    }
}


// test #1: set up a simple SCFG with null states & bifurcations,
// and verify that the eliminated SCFG looks as it should.
void do_elim_test (bool& ok)
{
  // B  -> L L
  // L  -> L    (.9)
  //     | L2   (.1)
  // L2 -> E    (.5)
  //     | B    (.25)
  //     | End  (.25)
  // E  -> x E'
  // E' -> End

  // Types: bifurc(B=0), null(L=1,L2=2), emit(E)

  // Let qX be the null subtree probability for X. Thus
  // qB  = qL*qL
  // qL  = qL2
  // qL2 = .25*qB + .25

  // 16*qB = (qB+1)^2
  //     0 = qB^2 - 14*qB + 1

  //    qB ~= 0.0718
  //    qL ~= 0.268

  // effective transitions:
  // t(B->L) = 2*qL = .536  = 2^-0.900   [nonempty parse tree only]
  // t(B->end) = qB = .0718 = 2^-3.800      [empty parse tree only]

  // t(L->L) = t(B->L) * t(L->B) = .25*.536 = .134 = 2^-2.900
  //  (marginalizing out this self-looping null transition...)

  // t(L->E) = .5 / (1 - .134) = .5774 = 2^-.792
  // t(L->nonemptyB) = .25 / (1 - .134) = .2887 = 2^-1.792
  // t(L->End) = .25 + .25*.0718 = .26795 = 2^-1.900

  // t(Start->nonemptyB) = 1 / (1 - .134) = 1.15473 = 2^.208
  // t(Start->E) = t(B->L) * t(L->E) = .536 * .5774 = .3095 = 2^-1.692

  // create SCFG
  Pair_CFG_scores a (4);
  a.state_name[0] = "B";
  a.state_name[1] = "L";
  a.state_name[2] = "L2";
  a.state_name[3] = "E";

  a.init_bifurc (0, 1, 1);    // B  -> L L
  a.state_type[1] = Pair_CFG_state_type_enum::Null;  // L
  a.state_type[2] = Pair_CFG_state_type_enum::Null;  // L2
  a.init_emit (3, Pair_CFG_state_type_enum::EmitXL, (Score) 0);  // E

  const int s = Grammar_state_enum::Start;
  const int e = Grammar_state_enum::End;

  a.transition(s,0) = Prob2Score(1);    // Start -> B

  a.transition(1,1) = Prob2Score(.9);    // L  -> L
  a.transition(1,2) = Prob2Score(.1);    // L  -> L2

  a.transition(2,3) = Prob2Score(.5);   // L2 -> x E  (.5)
  a.transition(2,0) = Prob2Score(.25);  // L2 -> B    (.25)
  a.transition(2,e) = Prob2Score(.25);  // L2 -> end  (.25)

  a.transition(3,e) = Prob2Score(1);    // E  -> end

  // do null state elimination
  vector<Prob> p_empty;
  Concrete_transition_probs tp_orig(1,0.), tp_elim(1,0.);
  vector<int> b2a_state, c2b_state, bifurc_prec_states;

  const Pair_CFG_scores b = a.add_null_bifurc_transitions (p_empty, b2a_state, bifurc_prec_states);
  const Pair_CFG_scores c = b.eliminate_null_states (tp_orig, tp_elim, bifurc_prec_states, c2b_state);

  // log
  CTAG(6,TESTELIM) << "a.state_name: " << a.state_name << '\n';
  a.show (CL);

  CL << "p_empty: " << p_empty << '\n';
  CL << "b2a: " << b2a_state << "\nc2b: " << c2b_state << '\n';
  CL << "bifurc_prec_states: " << bifurc_prec_states << '\n';

  CL << "b.state_name: " << b.state_name << '\n';
  b.show (CL);

  CL << "c.state_name: " << c.state_name << '\n';
  c.show (CL);

  // do tests
  test_transition (c, s, 2, 1.15473, ok);  // t(Start->nonemptyB) = 1 / (1 - .134) = 1.15473 = 2^.208
  test_transition (c, s, 1, .3095, ok);  // t(Start->E) = t(B->L) * t(L->E) = .536 * .5774 = .3095 = 2^-1.692

  test_transition (c, 0, 1, .5774, ok);  // t(L->E) = .5 / (1 - .134) = .5774 = 2^-.792
  test_transition (c, 0, 2, .2887, ok);  // t(L->nonemptyB) = .25 / (1 - .134) = .2887 = 2^-1.792
  test_transition (c, 0, e, .26795, ok);  // t(L->End) = .25 + .25*.0718 = .26795 = 2^-1.900
}

// test #2: construct an SCFG that is guaranteed to contain exactly one null subtree.
// verify (a) state path from traceback to eliminated SCFG and (b) state path with restored null subtree.
void do_restore_test (bool& ok)
{
  // B  -> L R
  // L  -> N N
  // N  -> End
  // R  -> F N
  // F  -> E
  // E  -> x E'
  // E' -> End

  // create SCFG
  Pair_CFG_scores a (6);
  a.state_name[0] = "B";
  a.state_name[1] = "L";
  a.state_name[2] = "N";
  a.state_name[3] = "R";
  a.state_name[4] = "F";
  a.state_name[5] = "E";

  a.init_bifurc (0, 1, 3);    // B  -> L R
  a.init_bifurc (1, 2, 2);    // L  -> N N
  a.state_type[2] = Pair_CFG_state_type_enum::Null;  // N
  a.init_bifurc (3, 4, 2);    // R  -> F N
  a.state_type[4] = Pair_CFG_state_type_enum::Null;  // F
  a.init_emit (5, Pair_CFG_state_type_enum::EmitXL, (Score) 0);  // E

  const int s = Grammar_state_enum::Start;
  const int e = Grammar_state_enum::End;

  a.transition(s,0) = Prob2Score(1);    // Start -> B
  a.transition(2,e) = Prob2Score(1);    // N -> End
  a.transition(4,5) = Prob2Score(1);    // F -> E
  a.transition(5,e) = Prob2Score(1);    // E -> end

  // do null state elimination
  vector<Prob> p_empty;
  Concrete_transition_probs tp_orig(1,0.), tp_elim(1,0.);
  vector<int> b2a_state, c2b_state, bifurc_prec_states;

  const Pair_CFG_scores b = a.add_null_bifurc_transitions (p_empty, b2a_state, bifurc_prec_states);
  const Pair_CFG_scores c = b.eliminate_null_states (tp_orig, tp_elim, bifurc_prec_states, c2b_state);

  // log
  CTAG(6,TESTELIM) << "a.state_name: " << a.state_name << '\n';
  a.show (CL);

  CL << "p_empty: " << p_empty << '\n';
  CL << "b2a: " << b2a_state << "\nc2b: " << c2b_state << '\n';
  CL << "bifurc_prec_states: " << bifurc_prec_states << '\n';

  CL << "b.state_name: " << b.state_name << '\n';
  b.show (CL);

  CL << "c.state_name: " << c.state_name << '\n';
  c.show (CL);

  // do DP and traceback
  Named_profile npx;
  Named_profile npy;

  npx.dsq.push_back (0);

  Fold_envelope xenv;
  Fold_envelope yenv;

  xenv.initialise_full (npx.dsq.size());
  yenv.initialise_full (npy.dsq.size());

  xenv.remove_empty_bifurcations();
  yenv.remove_empty_bifurcations();

  Pair_inside_matrix inside (npx, npy, xenv, yenv, c);
  Pair_inside_cell_sorter sorter (inside);
  Pair_CFG_local_path local_path = sorter.traceback_with_coords();

  const vector<int> b_path = Pair_CFG_scores::sample_eliminated_states (tp_orig, tp_elim, c2b_state, b, local_path.path);
  const vector<int> a_path = Pair_CFG_scores::sample_eliminated_bifurcations (p_empty, b2a_state, a, b, b_path);

  // log
  CTAG(6,TESTELIM) << "c_path: " << path2string (local_path.path, c.state_name) << '\n';
  CL << "b_path: " << path2string (b_path, b.state_name) << '\n';
  CL << "a_path: " << path2string (a_path, a.state_name) << '\n';

  // do tests
  test_path (c, local_path.path, "Start E End", ok);
  test_path (b, b_path, "Start B R F E End", ok);
  test_path (a, a_path, "Start B L N End N End R F E End N End", ok);
}

// main
int main (int argc, char** argv)
{
  // test passed?
  bool ok = true;

  // try-catch
  try
    {
      INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
		      "test pairwise SCFG state elimination & restoral");

      bool elim_test, restore_test;
      opts.add ("e -elim", elim_test = true, "do null state/bifurcation elimination test");
      opts.add ("r -restore", restore_test = true, "do null bifurcation restoral test");

      opts.add ("d -debug", "-log TRANSMAT_ELIMINATE -log CFG_ELIM -log TESTELIM", "print lots of highly-relevant debugging stuff");


      opts.parse_or_die();

      // test #1: set up a simple SCFG with null states & bifurcations,
      // and verify that the eliminated SCFG looks as it should.
      if (elim_test)
	do_elim_test (ok);

      // test #2: construct an SCFG that is guaranteed to contain exactly one null subtree.
      // verify (a) state path from traceback to eliminated SCFG and (b) state path with restored null subtree.
      if (restore_test)
	do_restore_test (ok);

    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << (ok ? "ok\n" : "not ok\n");
  return 0;
}
