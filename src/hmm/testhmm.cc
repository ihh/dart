#include "hmm/singlehmm.h"
#include "hmm/pairhmm.h"

void test_pair_HMM (double a, double b, double c, double d, double m, const char* xseq, const char* yseq)
{
  const double match_odds = 4*m;
  const double mismatch_odds = (4./3.) * (1-m);
  cerr << "(testing pair HMM: a=" << a << " b=" << b << " c=" << c << " d=" << d << " m=" << m << " ";
  cerr << "match_odds=" << match_odds << " mismatch_odds=" << mismatch_odds << " ";
  cerr << "xseq=" << xseq << " yseq=" << yseq << ")\n";

  // calculate transition scores
  const double s2i = c;
  const double s2d = (1-c)*(1-a)*b;
  const double s2m = (1-c)*a;
  const double s2e = (1-c)*(1-a)*(1-b);
  const double i2i = d;
  const double i2d = (1-d)*(1-a)*b;
  const double i2m = (1-d)*a;
  const double i2e = (1-d)*(1-a)*(1-b);
  const double d2i = c;
  const double d2d = (1-c)*(1-a)*b;
  const double d2m = (1-c)*a;
  const double d2e = (1-c)*(1-a)*(1-b);
  const double m2i = c;
  const double m2d = (1-c)*(1-a)*b;
  const double m2m = (1-c)*a;
  const double m2e = (1-c)*(1-a)*(1-b);

  // create Pair_HMM_scores
  Pair_HMM_scores pair (3, &DNA_alphabet);
  // set state types
  pair.state_type[0] = Pair_HMM_scores::EmitY;   // ins
  pair.state_type[1] = Pair_HMM_scores::EmitX;   // del
  pair.state_type[2] = Pair_HMM_scores::EmitXY;  // mat
  // set emit scores
  pair.single_emit[0] = pair.single_emit[1] = vector<Score> (4, (Score) 0);
  pair.pair_emit[2] = array2d<Score> (4, 4, Prob2Score (mismatch_odds));
  for (int i = 0; i < 4; ++i) pair.pair_emit[2](i,i) = Prob2Score (match_odds);
  // set transition scores
  pair.start[0] = Prob2Score (s2i);
  pair.start[1] = Prob2Score (s2d);
  pair.start[2] = Prob2Score (s2m);
  pair.start_to_end() = Prob2Score (s2e);
  pair.transition(0,0) = Prob2Score (i2i);
  pair.transition(0,1) = Prob2Score (i2d);
  pair.transition(0,2) = Prob2Score (i2m);
  pair.end[0] = Prob2Score (i2e);
  pair.transition(1,0) = Prob2Score (d2i);
  pair.transition(1,1) = Prob2Score (d2d);
  pair.transition(1,2) = Prob2Score (d2m);
  pair.end[1] = Prob2Score (d2e);
  pair.transition(2,0) = Prob2Score (m2i);
  pair.transition(2,1) = Prob2Score (m2d);
  pair.transition(2,2) = Prob2Score (m2m);
  pair.end[2] = Prob2Score (m2e);

  // create sequences
  Biosequence xbs (xseq);
  Biosequence ybs (yseq);
  Score_profile xprof, yprof;
  Digitized_biosequence xdsq, ydsq;
  DNA_alphabet.seq2score (xbs, xprof);
  DNA_alphabet.seq2score (ybs, yprof);
  DNA_alphabet.seq2dsq (xbs, xdsq);
  DNA_alphabet.seq2dsq (ybs, ydsq);
  const int xlen = xbs.size();
  const int ylen = ybs.size();
  // make probabilistic submat
  array2d<double> submat (4, 4, mismatch_odds);
  for (int i = 0; i < 4; ++i) submat(i,i) = match_odds;

  // create forward matrix
  array2d<double> fmat (xlen+1, ylen+1, 0.);
  array2d<double> fdel (xlen+1, ylen+1, 0.);
  array2d<double> fins (xlen+1, ylen+1, 0.);
  fmat(0,0) = 1.;  // start and match states have same outgoing transitions
  for (int x = 0; x <= xlen; ++x)
    for (int y = 0; y <= ylen; ++y)
      {
	if (x > 0)  // del
	  fdel(x,y) = fmat(x-1,y)*m2d + fdel(x-1,y)*d2d + fins(x-1,y)*i2d;
	if (y > 0)  // ins
	  fins(x,y) = fmat(x,y-1)*m2i + fdel(x,y-1)*d2i + fins(x,y-1)*i2i;
	if (x > 0 && y > 0)  // mat
	  fmat(x,y) = (fmat(x-1,y-1)*m2m + fdel(x-1,y-1)*d2m + fins(x-1,y-1)*i2m) * submat(xdsq[x-1],ydsq[y-1]);
      }
  const double fwd_prob = fmat(xlen,ylen)*m2e + fdel(xlen,ylen)*d2e + fins(xlen,ylen)*i2e;

  // create backward matrix
  array2d<double> bmat (xlen+1, ylen+1, 0.);
  array2d<double> bdel (xlen+1, ylen+1, 0.);
  array2d<double> bins (xlen+1, ylen+1, 0.);
  for (int x = xlen; x >= 0; --x)
    for (int y = ylen; y >= 0; --y)
      {
	if (x == xlen && y == ylen)
	  {
	    bmat(x,y) = m2e;
	    bdel(x,y) = d2e;
	    bins(x,y) = i2e;
	  }
	else
	  {
	    // calculate probabilities in destination cells
	    double pmat = 0, pins = 0, pdel = 0;
	    if (x < xlen)
	      pdel = bdel(x+1,y);
	    if (y < ylen)
	      pins = bins(x,y+1);
	    if (x < xlen && y < ylen)
	      pmat = bmat(x+1,y+1)*submat(xdsq[x],ydsq[y]);
	    // calculate probs. in this cell
	    if (x > 0)
	      bdel(x,y) = d2m*pmat + d2d*pdel + d2i*pins;
	    if (y > 0)
	      bins(x,y) = i2m*pmat + i2d*pdel + i2i*pins;
	    if (x > 0 && y > 0)
	      bmat(x,y) = m2m*pmat + m2d*pdel + m2i*pins;
	  }
      }
  const double back_prob = s2m*bmat(1,1)*submat(xdsq[0],ydsq[0]) + s2i*bins(0,1) + s2d*bdel(1,0);

  // create Pair_forward_backward_DP_matrix
  Pair_forward_backward_DP_matrix fb (pair, xprof, yprof);
  // check scores
  const Score forward_score = fb.forward_score;
  const Score final_score = fb.final_score;
  const Score my_fwd_sc = Prob2Score (fwd_prob);
  const Score my_back_sc = Prob2Score (back_prob);
  if (my_fwd_sc != final_score || my_back_sc != final_score || forward_score != final_score)
    {
      cerr << "fmat:\n" << fmat;
      cerr << "fins:\n" << fins;
      cerr << "fdel:\n" << fdel;
      cerr << "bmat:\n" << bmat;
      cerr << "bins:\n" << bins;
      cerr << "bdel:\n" << bdel;
      cerr << "DART forward matrix:\n";
      fb.dump_forward_cell (cerr);
      cerr << "DART backward matrix:\n";
      fb.show (cerr);
      THROWEXPR ("Forward score: " << my_fwd_sc << " Backward score: " << my_back_sc << " DART forward score: " << forward_score << " DART backward score: " << final_score);
    }
  cerr << "(test-forward, test-backward, DART-forward & DART-backward scores match: " << final_score << ")\n";
}

int main (int argc, char** argv)
{
  try
    {
      test_pair_HMM (.1, .2, .3, .4, .8, "AAC", "AC");
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
