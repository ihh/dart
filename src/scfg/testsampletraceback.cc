#include <math.h>

#include "util/logfile.h"
#include "util/vector_output.h"
#include "scfg/postenv.h"

class QStem : public Pair_CFG_scores
{
public:
  // state indices
  const int loop_states;
  enum { PrePad = 0, PadL = 1, StemBif = 2, PreStem = 3,
	 StemLR = 4, StemL = 5, StemR = 6,
	 Loop = 7 };

  // parameter accessors
  static inline int loop (int i) { return Loop + i; }

  // constructor
  QStem (int min_loop_size = 4);
};

int main (int argc, char** argv)
{
  const int n_iter = 3;
  int min_loop_size = 4;
  int nstruct = 10;
  int max_subseq_len = -1;

  Rnd::set_seed (12345);

  try
    {
      cerr << "(testing " << n_iter << " sequences)\n";
      for (int iter = 0; iter < n_iter; ++iter)
	{
	  // initialise sequence database, alphabet & null model
	  const Alphabet& alphabet = CFG_alphabet;
	  vector<Prob> null_emit_prob (alphabet.size(), 1. / (double) alphabet.size());  // default null model

	  Named_profile nprof;
	  nprof.name = "k10_dmel_utr3_601to800";
	  nprof.seq = "GGCGCGACATAATTGAAACAGAGTGGCCGGATCCTAAGCACTCTGGCATAAATCACACTGTAGGCCTTAGATTACACCACTTGATTGTATTTTTAAATTAATTCTTAAAAACTACAAATTAAGATCACTCTGTGAACGTGTGCTCGATGGTGTGCATCTACGGATTTGTTTTTGTGTTCTTTTTCCCCACCCAACCCCGT";
	  if (iter > 0)
	    for (int n = 0; n < (int) nprof.seq.size(); ++n)
	      nprof.seq[n] = alphabet.int2char (Rnd::rnd_int (alphabet.size()));
	  nprof.seq_update (alphabet);

	  cerr << "(creating local Inside matrix #" << iter+1 << " for sequence '" << nprof.seq.substr(0,10) << "...')\n";
	  
	  // create QStem
	  QStem qs_cfg (min_loop_size);

	  // create x-envelope & empty y-envelope
	  Fold_envelope xenv;
	  const Fold_envelope dummy_yenv;
	  const Named_profile dummy_npy;

	  // initialise fold envelopes, DP matrix, cell sorter
	  xenv.initialise_3prime_banded (nprof.size(), max_subseq_len);
	  const Pair_inside_matrix inside (nprof, dummy_npy, xenv, dummy_yenv, qs_cfg, TRUE);
	  const Pair_inside_cell_sorter cell_sorter (inside);
      
	  // do the tracebacks
	  cerr << " (created matrix; now sampling " << nstruct << " tracebacks...)\n";
	  for (int n = 1; n <= nstruct; ++n)
	    {
	      const Pair_CFG_local_path lp = cell_sorter.traceback_with_coords();
	      if (n % 10 == 0) cerr << " (sampled " << n << " tracebacks)\n";
	    }
	}
    }
  catch (const Dart_exception& e)
    {
      cout << "not ok\n";
      cerr << e.what();
      exit(1);
    }
  cout << "ok\n";
  return 0;
}


QStem::QStem (int min_loop_size)
  : Pair_CFG_scores (Loop + min_loop_size),
    loop_states (min_loop_size)
{
  const Score one = 0;
  const Score half = Prob2Score(0.5);
  const Score quarter = Prob2Score(0.25);

  // set state types
  init_emit (PrePad,  Null,    one);
  init_emit (PadL,    EmitXL,  one);
  init_emit (StemBif, Bifurc,  one);
  init_emit (PreStem, Null,    one);
  init_emit (StemLR,  EmitXLR, one);
  init_emit (StemL,   EmitXL,  one);
  init_emit (StemR,   EmitXR,  one);
  for (int i = 0; i < loop_states; ++i)
    init_emit (loop(i), EmitXL, one);

  // set transitions

  transition (Start,   PrePad)  = one;
  transition (PrePad,  StemBif) = half;
  transition (PrePad,  PreStem) = half;  // this transition allows stems flush to 3' end
  transition (PrePad,  PadL)    = half;
  transition (PrePad,  End)     = one;

  transition (PadL,    PrePad)  = one;

  transition (PreStem, StemLR)  = one;

  transition (StemLR,  StemLR)  = half;
  transition (StemLR,  StemL)   = quarter;
  transition (StemLR,  StemR)   = quarter;
  transition (StemLR,  Loop)    = one;
  
  transition (StemL,   StemLR)  = half;
  transition (StemL,   StemL)   = half;

  transition (StemR,   StemLR)  = half;
  transition (StemR,   StemR)   = half;

  for (int i = 0; i < loop_states; ++i)
    {
      const int next_loop = i+1 < loop_states ? loop(i+1) : End;
      transition (loop(i), loop(i))   = half;
      transition (loop(i), next_loop) = half;
    }

  bifurc[StemBif] = Bifurcation (PreStem, PrePad);

  // set emissions
  const int A = alphabet().size();
  for (int x = 0; x < A; ++x)
    {
      emit[StemL][x] = emit[StemR][x] = one;
      for (int i = 0; i < loop_states; ++i)
	emit[loop(i)][x] = one;
    }
  for (int xl = 0; xl < A; ++xl)
    for (int xr = 0; xr < A; ++xr)
      {
	const int idx = emit_xl_mul(EmitXLR) * xl + emit_xr_mul(EmitXLR) * xr;
	emit[StemLR][idx] = one;
      }
}

