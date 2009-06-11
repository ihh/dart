#include "seq/stockholm.h"
#include "tree/tree_alignment.h"
#include "scfg/foldenv.h"
#include "scfg/postenv.h"
#include "scfg/paircfgdp.h"
#include "indiegram/tripletscfgdp.h"

struct TKFST_Params_debug
{
  // branch lengths
  double t, u, v;  // WX, WY, WZ

  // methods & member variables cut-n-pasted from TKFST_Params (inheritance didn't seem to work, sigh)
  double stem_prob;

  double lambda_1;
  double mu_1;

  double lambda_2;
  double mu_2;

  /// singlet transitions probs
 double K1() const { return lambda_1 / mu_1; }
 double K2() const { return lambda_2 / mu_2; }

  /// branch transition probs (alpha, beta, gamma)
 double a1 (double t) const { return exp (-mu_1 * t); }
 double b1 (double t) const { return lambda_1 * (1 - exp((lambda_1-mu_1)*t)) / (mu_1 - lambda_1*exp((lambda_1-mu_1)*t)); }
 double g1 (double t) const { return 1 - mu_1 * (1 - exp((lambda_1-mu_1)*t)) / ((1-exp(-mu_1*t)) * (mu_1 - lambda_1*exp((lambda_1-mu_1)*t))); }
 double a2 (double t) const { return exp (-mu_2 * t); }
 double b2 (double t) const { return lambda_2 * (1 - exp((lambda_2-mu_2)*t)) / (mu_2 - lambda_2*exp((lambda_2-mu_2)*t)); }
 double g2 (double t) const { return 1 - mu_2 * (1 - exp((lambda_2-mu_2)*t)) / ((1-exp(-mu_2*t)) * (mu_2 - lambda_2*exp((lambda_2-mu_2)*t))); }
  
  /// stem construction probability
 double pS() const { return stem_prob; }


  void add_sc (Score& sc, Score delta, const char* term) {
    cout << "lg(" << term << ") = " << Score2Bits (delta) << '\n';
    ScorePMulAcc (sc, delta);
  }


  TKFST_Params_debug()
    {
      // WX, WY, WZ branch lengths; hardcoded
      t = 1;
      u = 1;
      v = .1;

      // TKFST params; hardcoded
      stem_prob = 0.1;
      lambda_1 = 0.025;
      mu_1 = 0.03;
      lambda_2 = 0.007;
      mu_2 = 0.01;
    }

  void run_ref() {
    cout << "\nReference alignment:\n";
    // accumulate scores
    //    grep = alignments/0.1/13.reference.parse.inits | perl -ne 'if(/transition_scores.(.*) = (Prob2Score.*);/){print"add_sc (sc, $2, \"$1\");\n"}'
      Score sc = 0;
      add_sc (sc, Prob2Score ((1)), "start[L_L_L_L]");
      add_sc (sc, Prob2Score ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (L_L_L_L, BiSL_BmSL_BmSL_BmSL)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (S_S_S_S, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score (((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))), "transition (IS_MS_MS_MS, L_L_L_L)");
      add_sc (sc, Prob2Score ((b1(t)) * (1-b1(u)) * (1-b1(v))), "transition (L_L_L_L, L_IL_WL_WL)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v))), "transition (L_IL_WL_WL, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (IL_ML_ML_ML, BiSL_BmSL_BmSL_BmSL)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (S_S_S_S, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score ((K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))), "transition (IS_MS_MS_MS, IS_MS_MS_MS)");
      add_sc (sc, Prob2Score (((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))), "transition (IS_MS_MS_MS, L_L_L_L)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (L_L_L_L, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (IL_ML_ML_ML, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))), "transition (IL_ML_ML_ML, IL_ML_DL_ML)");
      add_sc (sc, Prob2Score ((1-K1()) * ((1-b1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))), "end[IL_ML_DL_ML]");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (L_L_L_L, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (IL_ML_ML_ML, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (IL_ML_ML_ML, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))), "transition (IL_ML_ML_ML, IL_ML_ML_ML)");
      add_sc (sc, Prob2Score ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))), "end[IL_ML_ML_ML]");
      add_sc (sc, Prob2Score ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))), "end[L_L_L_L]");
      cout << "Final score: " << Score2Bits(sc) << " bits\n";
    }

  void run_ig() {
    cout << "\nIndiegram alignment:\n";
    // accumulate scores
    // grep = inferred/0.1/13.indiegram.parse.inits | perl -ne 'if(/transition_scores.(.*) = (Prob2Score.*);/){print"add_sc (sc, $2, \"$1\");\n"}'
      Score sc = 0;
    add_sc (sc, Prob2Score ((1)), "start[L_L_L_L]");
    add_sc (sc, Prob2Score ((b1(u)*pS())), "transition (L_L_L_L, L_L_BiSiL_WL)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_Si_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score (((1-K2())*(1))), "transition (e_e_ISi_e, e_e_Li_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_Li_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((K1()*pS())), "transition (e_e_ILi_e, e_e_BiSiLi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_Si_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_e_ISi_e, e_e_ISi_e)");
    add_sc (sc, Prob2Score (((1-K2())*(1))), "transition (e_e_ISi_e, e_e_Li_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_Li_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_ILi_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((1-K1())), "end[e_e_ILi_e]");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_Li_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_ILi_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_ILi_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_e_ILi_e, e_e_ILi_e)");
    add_sc (sc, Prob2Score ((1-K1())), "end[e_e_ILi_e]");
    add_sc (sc, Prob2Score ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))), "transition (L_L_L_WL, BiSL_BpeL_BpeL_BmSL)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (S_e_e_S, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score (((1-K2())*(1)) * ((1-b2(v))*((1)*(1)))), "transition (IS_e_e_MS, L_e_e_L)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (L_e_e_L, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((K1()*pS()) * ((1-b1(v))*(a1(v)))), "transition (IL_e_e_ML, BiSL_e_e_BmSL)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (S_e_e_S, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score ((K2()) * ((1-b2(v))*(a2(v)))), "transition (IS_e_e_MS, IS_e_e_MS)");
    add_sc (sc, Prob2Score (((1-K2())*(1)) * ((1-b2(v))*((1)*(1)))), "transition (IS_e_e_MS, L_e_e_L)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (L_e_e_L, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (IL_e_e_ML, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (IL_e_e_ML, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((1-K1()) * ((1-b1(v))*(1))), "end[IL_e_e_ML]");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (L_e_e_L, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (IL_e_e_ML, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (IL_e_e_ML, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((K1()) * ((1-b1(v))*(a1(v)))), "transition (IL_e_e_ML, IL_e_e_ML)");
    add_sc (sc, Prob2Score ((1-K1()) * ((1-b1(v))*(1))), "end[IL_e_e_ML]");
    add_sc (sc, Prob2Score ((K2())), "transition (e_Si_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score (((1-K2())*(1))), "transition (e_ISi_e_e, e_Li_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_Li_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_ILi_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1()*pS())), "transition (e_ILi_e_e, e_BiSiLi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_Si_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score ((K2())), "transition (e_ISi_e_e, e_ISi_e_e)");
    add_sc (sc, Prob2Score (((1-K2())*(1))), "transition (e_ISi_e_e, e_Li_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_Li_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_ILi_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_ILi_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((1-K1())), "end[e_ILi_e_e]");
    add_sc (sc, Prob2Score ((K1())), "transition (e_Li_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_ILi_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_ILi_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((K1())), "transition (e_ILi_e_e, e_ILi_e_e)");
    add_sc (sc, Prob2Score ((1-K1())), "end[e_ILi_e_e]");
      cout << "Final score: " << Score2Bits(sc) << " bits\n";
  }

};

int main (int argc, char** argv)
{
  TKFST_Params_debug debug;
  debug.run_ref();
  debug.run_ig();
  return 0;
}

