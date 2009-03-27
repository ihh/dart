#include <fstream>
#include <sstream>

#include "evoldoer/tkfstructuretree.h"

// default lambda & mu for TKFST_params
#define INIT_LAMBDA 1.
#define INIT_MU     2.

// name of TKFST pair SCFG in grammar file output
#define TKFST_pair_SCFG_name "TKFST_pair_SCFG"


TKFST_params::TKFST_params()
  : loop_mat (1, CFG_alphabet_size),
    stem_mat (1, CFG_alphabet_size * CFG_alphabet_size),
    loop (loop_mat, INIT_LAMBDA, INIT_MU),
    stem (stem_mat, INIT_LAMBDA, INIT_MU)
{ }

TKFST_single_CFG::TKFST_single_CFG (TKFST_params& params, bool odds_ratio)
  : Pair_CFG_scores (TotalStates),
    params (params)
{
  // set state types
  init_emit (aL, EmitXL);
  init_emit (aS, EmitXLR);

  state_type[L]      = Null;
  state_type[Lprime] = Null;
  state_type[Sprime] = Null;

  state_type[B]      = Bifurc;
  state_type[Bprime] = Bifurc;
  state_type[C]      = Bifurc;

  // make bifurcations
  bifurc[B]      = Bifurcation (Sprime, L);
  bifurc[Bprime] = Bifurcation (Sprime, Lprime);
  bifurc[C]      = Bifurcation (Sprime, Sprime);

  // calculate some score parameters
  const vector<Score> loop_emit_sc = Prob2ScoreVec (params.loop_mat.create_prior());
  const vector<Score> stem_emit_sc = Prob2ScoreVec (params.stem_mat.create_prior());

  const Score kappa1_sc     = Prob2Score (params.loop.kappa());
  const Score not_kappa1_sc = Prob2Score (1. - params.loop.kappa());
  const Score kappa2_sc     = Prob2Score (params.stem.kappa());
  const Score not_kappa2_sc = Prob2Score (1. - params.stem.kappa());

  const Score p1S_sc        = Prob2Score (params.stem_prob);
  const Score not_p1S_sc    = Prob2Score (1. - params.stem_prob);

  // log these scores
  if (CTAGGING(3,TKFST))
    {
      CL << "TKF structure tree singlet scores, odds_ratio=" << odds_ratio << "\n";
      CL << "kappa1_sc=" << kappa1_sc << ",\tnot_kappa1_sc=" << not_kappa1_sc << "\n";
      CL << "p1S_sc=" << p1S_sc << ",\tnot_p1S_sc=" << not_p1S_sc << "\n";
      CL << "kappa2_sc=" << kappa2_sc << ",\tnot_kappa2_sc=" << not_kappa2_sc << "\n";
    }

  // make emissions
  for (int xl = 0; xl < CFG_alphabet_size; ++xl)
    {
      emit[aL][xl] = ScorePMul (not_p1S_sc, odds_ratio ? -kappa1_sc : loop_emit_sc[xl]);
      for (int xr = 0; xr < CFG_alphabet_size; ++xr)
	{
	  const int xlr_idx = emit_xl_mul(EmitXLR) * xl + emit_xr_mul(EmitXLR) * xr;
	  emit[aS][xlr_idx] = ScorePMul (stem_emit_sc[xlr_idx], odds_ratio ? -ScorePMul3(kappa1_sc*2,loop_emit_sc[xl],loop_emit_sc[xr]) : 0);
	}
    }

  // make transitions
  transition (Start, L)    = odds_ratio ? -not_kappa1_sc : 0;

  transition (L, aL)       = kappa1_sc;
  transition (L, B)        = ScorePMul (kappa1_sc, p1S_sc);
  transition (L, End)      = not_kappa1_sc;
  
  transition (aL, L)       = 0;

  transition (aS, aS)      = kappa2_sc;
  transition (aS, Lprime)  = not_kappa2_sc;

  transition (Lprime, aL)     = kappa1_sc;
  transition (Lprime, Bprime) = ScorePMul (kappa1_sc, p1S_sc);
  transition (Lprime, C)      = ScorePMul3 (2*kappa1_sc, 2*p1S_sc, not_kappa1_sc);

  transition (Sprime, aS)  = kappa2_sc;

  // extra transitions that replace prohibited empty-sequence bifurcations
  transition (L, Sprime)   = ScorePMul (transition(L,B), transition(L,End));
}

TKFST_pair_CFG::TKFST_pair_CFG (TKFST_params& params, double time, bool conditional, bool odds_ratio)
  : Pair_CFG_scores (TotalStates),
    params (params),
    time (time)
{
  // set state types
  init_emit (adL1, EmitXLYL);
  init_emit (dL1,  EmitYL);
  init_emit (aL2,  EmitXL);
  init_emit (aL3,  EmitXL);
  init_emit (dL4,  EmitYL);
  init_emit (adS1, EmitXLRYLR);
  init_emit (dS1,  EmitYLR);
  init_emit (aS2,  EmitXLR);
  init_emit (aS3,  EmitXLR);
  init_emit (dS4,  EmitYLR);

  state_type[L1]       = Null;
  state_type[S1]       = Null;
  state_type[L2]       = Null;
  state_type[L3]       = Null;
  state_type[L4]       = Null;
  state_type[L1prime]  = Null;
  state_type[S1prime]  = Null;
  state_type[L3prime]  = Null;
  state_type[S3prime]  = Null;
  state_type[L4prime]  = Null;
  state_type[S4prime]  = Null;

  state_type[B41]      = Bifurc;
  state_type[B32]      = Bifurc;
  state_type[B11]      = Bifurc;
  state_type[B33]      = Bifurc;
  state_type[B44]      = Bifurc;
  state_type[B11prime] = Bifurc;
  state_type[B33prime] = Bifurc;
  state_type[B44prime] = Bifurc;
  state_type[C11]      = Bifurc;
  state_type[C33]      = Bifurc;
  state_type[C44]      = Bifurc;

  // make bifurcations
  bifurc[B41]      = Bifurcation (S4prime, L1);
  bifurc[B32]      = Bifurcation (S3prime, L2);
  bifurc[B11]      = Bifurcation (S1prime, L1);
  bifurc[B33]      = Bifurcation (S3prime, L3);
  bifurc[B44]      = Bifurcation (S4prime, L4);
  bifurc[B11prime] = Bifurcation (S1prime, L1prime);
  bifurc[B33prime] = Bifurcation (S3prime, L3prime);
  bifurc[B44prime] = Bifurcation (S4prime, L4prime);
  bifurc[C11]      = Bifurcation (S1prime, S1prime);
  bifurc[C33]      = Bifurcation (S3prime, S3prime);
  bifurc[C44]      = Bifurcation (S4prime, S4prime);

  // calculate some score parameters
  const vector<Score>  loop_emit_sc = Prob2ScoreVec (params.loop_mat.create_prior());
  const vector<Score>  stem_emit_sc = Prob2ScoreVec (params.stem_mat.create_prior());
  const array2d<Score> loop_pair_sc = Prob2ScoreArray2d (params.loop_mat.create_conditional_substitution_matrix (time));
  const array2d<Score> stem_pair_sc = Prob2ScoreArray2d (params.stem_mat.create_conditional_substitution_matrix (time));
  
  const Score alpha1_sc     = Prob2Score (params.loop.alpha (time));
  const Score not_alpha1_sc = Prob2Score (1. - params.loop.alpha (time));
  const Score alpha2_sc     = Prob2Score (params.stem.alpha (time));
  const Score not_alpha2_sc = Prob2Score (1. - params.stem.alpha (time));

  const Score beta1_sc      = Prob2Score (params.loop.beta (time));
  const Score not_beta1_sc  = Prob2Score (1. - params.loop.beta (time));
  const Score beta2_sc      = Prob2Score (params.stem.beta (time));
  const Score not_beta2_sc  = Prob2Score (1. - params.stem.beta (time));

  const Score gamma1_sc     = Prob2Score (params.loop.gamma (time));
  const Score not_gamma1_sc = Prob2Score (1. - params.loop.gamma (time));
  const Score gamma2_sc     = Prob2Score (params.stem.gamma (time));
  const Score not_gamma2_sc = Prob2Score (1. - params.stem.gamma (time));

  const Score kappa1_sc     = Prob2Score (params.loop.kappa());
  const Score not_kappa1_sc = Prob2Score (1. - params.loop.kappa());
  const Score kappa2_sc     = Prob2Score (params.stem.kappa());
  const Score not_kappa2_sc = Prob2Score (1. - params.stem.kappa());

  const Score p1S_sc        = Prob2Score (params.stem_prob);
  const Score not_p1S_sc    = Prob2Score (1. - params.stem_prob);

  const Score kappa1_anc_sc     = conditional ? 0 : kappa1_sc;
  const Score not_kappa1_anc_sc = conditional ? 0 : not_kappa1_sc;
  const Score kappa2_anc_sc     = conditional ? 0 : kappa2_sc;
  const Score not_kappa2_anc_sc = conditional ? 0 : not_kappa2_sc;

  const Score p1S_anc_sc        = conditional ? 0 : p1S_sc;

  // log these scores
  if (CTAGGING(3,TKFST))
    {
      CL << "TKF structure tree pair scores, time=" << time << ", conditional=" << conditional << ", odds_ratio=" << odds_ratio << "\n";
      CL << "alpha1_sc=" << alpha1_sc << ",\tnot_alpha1_sc=" << not_alpha1_sc << "\n";
      CL << "beta1_sc=" << beta1_sc << ",\tnot_beta1_sc=" << not_beta1_sc << "\n";
      CL << "gamma1_sc=" << gamma1_sc << ",\tnot_gamma1_sc=" << not_gamma1_sc << "\n";
      CL << "kappa1_sc=" << kappa1_sc << ",\tnot_kappa1_sc=" << not_kappa1_sc << "\n";
      CL << "kappa1_anc_sc=" << kappa1_anc_sc << ",\tnot_kappa1_anc_sc=" << not_kappa1_anc_sc << "\n";
      CL << "p1S_sc=" << p1S_sc << ",\tnot_p1S_sc=" << not_p1S_sc << "\n";
      CL << "p1S_anc_sc=" << p1S_anc_sc << "\n";
      CL << "alpha2_sc=" << alpha2_sc << ",\tnot_alpha2_sc=" << not_alpha2_sc << "\n";
      CL << "beta2_sc=" << beta2_sc << ",\tnot_beta2_sc=" << not_beta2_sc << "\n";
      CL << "gamma2_sc=" << gamma2_sc << ",\tnot_gamma2_sc=" << not_gamma2_sc << "\n";
      CL << "kappa2_sc=" << kappa2_sc << ",\tnot_kappa2_sc=" << not_kappa2_sc << "\n";
      CL << "kappa2_anc_sc=" << kappa2_anc_sc << ",\tnot_kappa2_anc_sc=" << not_kappa2_anc_sc << "\n";
    }

  // make emissions
  for (int xl = 0; xl < CFG_alphabet_size; ++xl)
    {
      const Score loopxl_sc = ScorePMul (not_p1S_sc, odds_ratio ? -kappa1_sc : loop_emit_sc[xl]);
      emit[dL1][xl] = loopxl_sc;
      emit[aL2][xl] = conditional ? 0 : loopxl_sc;
      emit[aL3][xl] = conditional ? 0 : loopxl_sc;
      emit[dL4][xl] = loopxl_sc;
      for (int yl = 0; yl < CFG_alphabet_size; ++yl)
	{
	  const int xlyl_idx = emit_xl_mul(EmitXLYL) * xl + emit_yl_mul(EmitXLYL) * yl;
	  emit[adL1][xlyl_idx] = ScorePMul (conditional ? 0 : loopxl_sc,
					    ScorePMul (loop_pair_sc(xl,yl),
						       odds_ratio ? -ScorePMul(kappa1_sc,loop_emit_sc[yl]) : 0));
	}
      for (int xr = 0; xr < CFG_alphabet_size; ++xr)
	{
	  const int xlr_idx = emit_xl_mul(EmitXLR) * xl + emit_xr_mul(EmitXLR) * xr;
	  const Score stemxlr_sc = ScorePMul (stem_emit_sc[xlr_idx], odds_ratio ? -ScorePMul3(kappa1_sc*2,loop_emit_sc[xl],loop_emit_sc[xr]) : 0);
	  emit[dS1][xlr_idx] = stemxlr_sc;
	  emit[aS2][xlr_idx] = conditional ? 0 : stemxlr_sc;
	  emit[aS3][xlr_idx] = conditional ? 0 : stemxlr_sc;
	  emit[dS4][xlr_idx] = stemxlr_sc;
	  for (int yl = 0; yl < CFG_alphabet_size; ++yl)
	    for (int yr = 0; yr < CFG_alphabet_size; ++yr)
	      {
		const int ylr_idx = emit_yl_mul(EmitYLR) * yl + emit_yr_mul(EmitYLR) * yr;
		const int xlrylr_idx = emit_xl_mul(EmitXLRYLR) * xl + emit_xr_mul(EmitXLRYLR) * xr + emit_yl_mul(EmitXLRYLR) * yl + emit_yr_mul(EmitXLRYLR) * yr;
		emit[adS1][xlrylr_idx] = ScorePMul (conditional ? 0 : stemxlr_sc,
						    ScorePMul (stem_pair_sc(xlr_idx,ylr_idx),
							       odds_ratio ? -ScorePMul3(kappa1_sc*2,loop_emit_sc[yl],loop_emit_sc[yr]) : 0));
	      }
	}
    }

  // make transitions
  transition (Start, L1) = odds_ratio ? -2*not_kappa1_sc : 0;

  transition (L1, adL1) = ScorePMul3 (kappa1_anc_sc, not_beta1_sc, alpha1_sc);
  transition (L1, dL1)  = beta1_sc;
  transition (L1, aL2)  = ScorePMul3 (kappa1_anc_sc, not_beta1_sc, not_alpha1_sc);
  transition (L1, B11)  = ScorePMul (ScorePMul (kappa1_anc_sc, p1S_anc_sc), ScorePMul (not_beta1_sc, alpha1_sc));
  transition (L1, B41)  = ScorePMul (beta1_sc, p1S_sc);
  transition (L1, B32)  = ScorePMul (ScorePMul (kappa1_anc_sc, p1S_anc_sc), ScorePMul (not_beta1_sc, not_alpha1_sc));
  transition (L1, End)  = ScorePMul (not_kappa1_anc_sc, not_beta1_sc);

  transition (adL1, L1) = 0;
  transition (dL1,  L1) = 0;

  transition (S1, adS1)    = ScorePMul3 (kappa2_anc_sc, not_beta2_sc, alpha2_sc);
  transition (S1, dS1)     = beta2_sc;
  transition (S1, aS2)     = ScorePMul3 (kappa2_anc_sc, not_beta2_sc, not_alpha2_sc);
  transition (S1, L1prime) = ScorePMul (not_kappa2_anc_sc, not_beta2_sc);

  transition (adS1, S1)    = 0;
  transition (dS1,  S1)    = 0;

  transition (L2, adL1) = ScorePMul3 (kappa1_anc_sc, not_gamma1_sc, alpha1_sc);
  transition (L2, dL1)  = gamma1_sc;
  transition (L2, aL2)  = ScorePMul3 (kappa1_anc_sc, not_gamma1_sc, not_alpha1_sc);
  transition (L2, B11)  = ScorePMul (ScorePMul (kappa1_anc_sc, p1S_anc_sc), ScorePMul (not_gamma1_sc, alpha1_sc));
  transition (L2, B41)  = ScorePMul (gamma1_sc, p1S_sc);
  transition (L2, B32)  = ScorePMul (ScorePMul (kappa1_anc_sc, p1S_anc_sc), ScorePMul (not_gamma1_sc, not_alpha1_sc));
  transition (L2, End)  = ScorePMul (not_kappa1_anc_sc, not_gamma1_sc);

  transition (aL2, L2)  = 0;
  
  transition (aS2, adS1)    = ScorePMul3 (kappa2_anc_sc, not_gamma2_sc, alpha2_sc);
  transition (aS2, dS1)     = gamma2_sc;
  transition (aS2, aS2)     = ScorePMul3 (kappa2_anc_sc, not_gamma2_sc, not_alpha2_sc);
  transition (aS2, L1prime) = ScorePMul (not_kappa2_anc_sc, not_gamma2_sc);

  transition (L3, aL3)      = kappa1_anc_sc;
  transition (L3, B33)      = ScorePMul (kappa1_anc_sc, p1S_anc_sc);
  transition (L3, End)      = not_kappa1_anc_sc;

  transition (aL3, L3)      = 0;

  transition (aS3, aS3)     = kappa2_anc_sc;
  transition (aS3, L3prime) = not_kappa2_anc_sc;

  transition (L4, dL4)      = kappa1_sc;
  transition (L4, B44)      = ScorePMul (kappa1_sc, p1S_sc);
  transition (L4, End)      = not_kappa1_sc;

  transition (dL4, L4)      = 0;

  transition (dS4, dS4)     = kappa2_sc;
  transition (dS4, L4prime) = not_kappa2_sc;

  transition (L1prime, adL1)     = ScorePMul3 (kappa1_anc_sc, not_beta1_sc, alpha1_sc);
  transition (L1prime, dL1)      = beta1_sc;
  transition (L1prime, aL2)      = ScorePMul3 (kappa1_anc_sc, not_beta1_sc, not_alpha1_sc);
  transition (L1prime, B11prime) = ScorePMul (ScorePMul (kappa1_anc_sc, p1S_anc_sc), ScorePMul (not_beta1_sc, alpha1_sc));
  transition (L1prime, B41)      = ScorePMul (beta1_sc, p1S_sc);
  transition (L1prime, B32)      = ScorePMul (ScorePMul (kappa1_anc_sc, p1S_anc_sc), ScorePMul (not_beta1_sc, not_alpha1_sc));
  transition (L1prime, C11)      = ScorePMul (ScorePMul3 (2*kappa1_anc_sc, not_kappa1_anc_sc, 2*p1S_anc_sc),
					      ScorePMul (3*not_beta1_sc, 2*alpha1_sc));

  transition (S1prime, adS1)    = ScorePMul3 (kappa2_anc_sc, not_beta2_sc, alpha2_sc);
  transition (S1prime, dS1)     = beta2_sc;
  transition (S1prime, aS2)     = ScorePMul3 (kappa2_anc_sc, not_beta2_sc, not_alpha2_sc);

  transition (L3prime, aL3)      = kappa1_anc_sc;
  transition (L3prime, B33prime) = ScorePMul (kappa1_anc_sc, p1S_anc_sc);
  transition (L3prime, C33)      = ScorePMul3 (2*kappa1_anc_sc, not_kappa1_anc_sc, 2*p1S_anc_sc);

  transition (S3prime, aS3)      = kappa2_anc_sc;

  transition (L4prime, dL4)      = kappa1_sc;
  transition (L4prime, B44prime) = ScorePMul (kappa1_sc, p1S_sc);
  transition (L4prime, C44)      = ScorePMul3 (2*kappa1_sc, not_kappa1_sc, 2*p1S_sc);

  transition (S4prime, dS4)      = kappa2_sc;

  // extra transitions that replace prohibited empty-sequence bifurcations
  transition (L1, S1prime)      = ScorePMul (transition (L1, B11), transition (L1, End));
  transition (L1, S4prime)      = ScorePMul (transition (L1, B41), transition (L1, End));
  transition (L1, S3prime)      = ScorePMul (transition (L1, B32), transition (L1, End));
  transition (L2, S1prime)      = ScorePMul (transition (L2, B11), transition (L2, End));
  transition (L2, S4prime)      = ScorePMul (transition (L2, B41), transition (L2, End));
  transition (L2, S3prime)      = ScorePMul (transition (L2, B32), transition (L2, End));
  transition (L3, S3prime)      = ScorePMul (transition (L3, B33), transition (L3, End));
  transition (L4, S4prime)      = ScorePMul (transition (L4, B44), transition (L4, End));
  transition (L1prime, S4prime) = ScorePMul (transition (L1prime, B41), transition (L1, End));
  transition (L1prime, S3prime) = ScorePMul (transition (L1prime, B32), transition (L1, End));

  // set state names
  state_name[L1] = "L1";
  state_name[S1] = "S1";
  state_name[L2] = "L2";
  state_name[L3] = "L3";
  state_name[L4] = "L4";
  state_name[L1prime] = "L1prime";
  state_name[S1prime] = "S1prime";
  state_name[L3prime] = "L3prime";
  state_name[S3prime] = "S3prime";
  state_name[L4prime] = "L4prime";
  state_name[S4prime] = "S4prime";
  state_name[B41] = "B41";
  state_name[B32] = "B32";
  state_name[B11] = "B11";
  state_name[B33] = "B33";
  state_name[B44] = "B44";
  state_name[B11prime] = "B11prime";
  state_name[B33prime] = "B33prime";
  state_name[B44prime] = "B44prime";
  state_name[C11] = "C11";
  state_name[C33] = "C33";
  state_name[C44] = "C44";
  state_name[adL1] = "adL1";
  state_name[dL1] = "dL1";
  state_name[aL2] = "aL2";
  state_name[aL3] = "aL3";
  state_name[dL4] = "dL4";
  state_name[adS1] = "adS1";
  state_name[dS1] = "dS1";
  state_name[aS2] = "aS2";
  state_name[aS3] = "aS3";
  state_name[dS4] = "dS4";

  // set grammar name
  name = TKFST_pair_SCFG_name;
}


TKFST_default_params::TKFST_default_params() : TKFST_params()
{
  // Stem rate matrix converted to HSM format from Bjarne Knudsen's files on the PFOLD server website using the following Perl code:
  // perl -e '@i=qw(AA CA GA UA AC CC GC UC AG CG GG UG AU CU GU UU);%i=map(($i[$_]=>$_),0..@i-1);open K,"dbl_mat.txt";$_=<K>;@x=split;@s=map($i{$x[$_]},0..@x-1);while(<K>){($f,@r)=split;$row[$i{$f}]=join("\t",@r[@s])."\n"};close K;open K2,"dbl_p.txt";$_=<K2>;@k=split;while(<K2>){($f,@r)=split;for($i=0;$i<@r;++$i){$eqm{$f.$k[$i]}=$r[$i]}}close K2;print"1 16\n@i\n",join(" ",map($eqm{$_},@i)),"\n",@row,map("0\n",@i)'

  sstring bkstem;
  bkstem << "1 16\n";
  bkstem << "a c g u A C G U b d h v B D H V\n";  // this line changed from the autogenerated perl
  bkstem << "0.001167 0.001806 0.001058 0.177977 0.001806 0.000391 0.266974 0.000763 0.001058 0.266974 0.000406 0.049043 0.177977 0.000763 0.049043 0.002793\n";
  bkstem << "-3.607 0.420 0.589 0.617 0.420 0.000 0.132 0.019 0.589 0.132 -0.000 0.026 0.617 0.019 0.026 0.000\n";
  bkstem << "0.271 -6.070 0.068 2.861 0.024 0.079 0.008 0.010 0.003 2.135 0.000 0.401 0.124 0.008 0.057 0.020\n";
  bkstem << "0.650 0.116 -2.489 0.257 0.006 -0.000 0.734 -0.000 0.097 0.060 -0.000 0.024 0.290 0.019 0.237 -0.000\n";
  bkstem << "0.004 0.029 0.002 -1.163 0.001 0.002 0.115 0.007 0.002 0.501 0.001 0.274 0.185 0.003 0.023 0.016\n";
  bkstem << "0.271 0.024 0.003 0.124 -6.070 0.079 2.135 0.008 0.068 0.008 0.000 0.057 2.861 0.010 0.401 0.020\n";
  bkstem << "0.000 0.365 -0.000 0.799 0.365 -5.666 0.992 0.204 -0.000 0.992 0.000 0.191 0.799 0.204 0.191 0.563\n";
  bkstem << "0.001 0.000 0.003 0.077 0.014 0.001 -0.823 0.004 0.000 0.130 0.004 0.019 0.334 0.001 0.232 0.003\n";
  bkstem << "0.029 0.024 0.000 1.551 0.018 0.105 1.265 -6.126 0.026 0.473 0.000 0.324 0.688 1.105 0.044 0.473\n";
  bkstem << "0.650 0.006 0.097 0.290 0.116 -0.000 0.060 0.019 -2.489 0.734 0.000 0.237 0.257 0.000 0.024 -0.000\n";
  bkstem << "0.001 0.014 0.000 0.334 0.000 0.001 0.130 0.001 0.003 -0.823 0.004 0.232 0.077 0.004 0.019 0.003\n";
  bkstem << "-0.000 -0.000 -0.000 0.252 0.000 -0.000 2.511 -0.000 0.000 2.511 -6.933 0.631 0.252 -0.000 0.631 0.145\n";
  bkstem << "0.001 0.015 0.001 0.996 0.002 0.002 0.101 0.005 0.005 1.262 0.005 -2.554 0.084 0.001 0.037 0.038\n";
  bkstem << "0.004 0.001 0.002 0.185 0.029 0.002 0.501 0.003 0.002 0.115 0.001 0.023 -1.163 0.007 0.274 0.016\n";
  bkstem << "0.029 0.018 0.026 0.688 0.024 0.105 0.473 1.105 -0.000 1.265 0.000 0.044 1.551 -6.126 0.324 0.473\n";
  bkstem << "0.001 0.002 0.005 0.084 0.015 0.002 1.262 0.001 0.001 0.101 0.005 0.037 0.996 0.005 -2.554 0.038\n";
  bkstem << "0.000 0.013 -0.000 0.988 0.013 0.079 0.287 0.129 -0.000 0.287 0.021 0.663 0.988 0.129 0.663 -4.261\n";
  bkstem << "0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n";

  istringstream bkstem_in (bkstem);
  stem_mat.read (bkstem_in);

  // Loop matrix
  // converted manually from Bjarne's
  sstring bkloop;
  bkloop << "1 4\n";
  bkloop << "A C G U\n";
  bkloop << "0.364097 0.151009 0.211881 0.273013\n";
  bkloop << "-0.749  0.164  0.322  0.263\n";
  bkloop << " 0.396 -1.565  0.242  0.927\n";
  bkloop << " 0.553  0.173 -0.964  0.239\n";
  bkloop << " 0.351  0.513  0.185 -1.050\n";
  bkloop << "0\n0\n0\n0\n";

  istringstream bkloop_in (bkloop);
  loop_mat.read (bkloop_in);

  // Stem & loop indel params (almost totally arbitrary)
  stem.lambda = 0.007;
  stem.mu     = 0.01;
  loop.lambda = 0.027;
  loop.mu     = 0.03;
  stem_prob   = 0.01;
}
