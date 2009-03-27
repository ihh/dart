#include "ecfg/swexmid.h"
#include "tree/pam.h"
#include "hmm/postpairhmm.h"
#include "handel/tkfhmm.h"
#include "util/maximise.h"
#include "tree/nj.h"
#include "hmm/pairenv.h"

// SWEXMID_params

SWEXMID_params::SWEXMID_params (int iTypes, const Alphabet& alphabet)
  : q (iTypes),
    r (iTypes),
    R (1, alphabet.size())
{
  R.init_alphabet (alphabet);
}

int SWEXMID_params::indelTypes() const
{
  return q.size();
}

int SWEXMID_params::alphabetSize() const
{
  return R.A;
}

double SWEXMID_params::lambda() const
{
  double lambda = 0;
  for (int i = 0; i < indelTypes(); ++i)
    lambda += kappa * mu * q[i] * (1 - r[i]) / (1 - kappa * r[i]);
  return lambda;
}

TKF_params SWEXMID_params::tkf_params() const
{
  return TKF_params ((EM_matrix&) R, lambda(), mu);  // this may not be the best estimate...
}

// Protein_default_SWEXMID_params

Protein_default_SWEXMID_params::Protein_default_SWEXMID_params()
  : SWEXMID_params (2, Protein_alphabet)
{
  // initialise the rate matrix, R, to the PAM matrix
  // eventually use HSM1 for this
  PAM_factory pam;
  R.X[0] = pam.rate_matrix();
  R.update();
  // initialise long indel params
  kappa = .99;  // arbitrary
  q[0] = 0.40;  // alpha from (Miklos, Lunter & Holmes; 2003)
  q[1] = 1. - q[0];
  r[0] = 0.55;  // r_1 from (Miklos, Lunter & Holmes; 2003)
  r[1] = 0.90;  // r_2 from (Miklos, Lunter & Holmes; 2003)
  const double mu_MLH = 0.95;  // mu from (Miklos, Lunter & Holmes; 2003)
  mu = mu_MLH * (q[0]*(1-r[0]) + q[1]*(1-r[1]));  // conversion to SWEXMID mu
}

// SWEXMID_root_params

SWEXMID_root_params::SWEXMID_root_params (PScope& pscope, const Alphabet& alphabet, int indelTypes)
  : q (pscope.new_group (indelTypes, "q")),
    r (indelTypes),
    kappa (pscope.new_boolean_group ("Kappa")),
    pi (pscope.new_alphabet_group (alphabet, "Pi"))
{
  for (int i = 0; i < indelTypes; ++i)
    {
      sstring rName;
      rName << "r" << i+1;
      r[i] = pscope.new_boolean_group (rName.c_str());
    }
}

int SWEXMID_root_params::indelTypes() const
{
  return q.group_size;
}

int SWEXMID_root_params::alphabetSize() const
{
  return pi.group_size;
}

void SWEXMID_root_params::set (PScores& pscore, const SWEXMID_params& params)
{
  if (indelTypes() != params.indelTypes() || alphabetSize() != params.alphabetSize())
    return;
  EM_matrix& R = (EM_matrix&) params.R;  // cast away const
  // calculate stuff
  const vector<Prob> params_pi = R.create_prior();
  // set pscores
  pscore[kappa[0]] = Prob2Score (1. - params.kappa);
  pscore[kappa[1]] = Prob2Score (params.kappa);
  for (int i = 0; i < indelTypes(); ++i)
    {
      pscore[q[i]] = Prob2Score (params.q[i]);
      pscore[r[i][0]] = Prob2Score (1. - params.r[i]);
      pscore[r[i][1]] = Prob2Score (params.r[i]);
    }
  for (int c = 0; c < alphabetSize(); ++c)
    pscore[pi[c]] = Prob2Score (params_pi[c]);
}

// SWEXMID_root_funcs

SWEXMID_root_funcs::SWEXMID_root_funcs (const SWEXMID_root_params& params)
  : notKappaR (params.indelTypes()),
    indelRatioTerm (params.indelTypes()),
    qIns (params.indelTypes())
{
  indelRatio = 0;
  for (int i = 0; i < params.indelTypes(); ++i)
    {
      notKappaR[i] = params.kappa[0] * params.r[i][1] + params.kappa[1] * params.r[i][0] + params.kappa[0] * params.r[i][0];
      indelRatioTerm[i] = params.kappa[1] * params.q[i] * params.r[i][0] / notKappaR[i];
      indelRatio += indelRatioTerm[i];
    }
  for (int i = 0; i < params.indelTypes(); ++i)
    qIns[i] = indelRatioTerm[i] / indelRatio;
}

// SWEXMID_branch_params

SWEXMID_branch_params::SWEXMID_branch_params (PScope& pscope, const Alphabet& alphabet, const char* branchPrefix)
  : M (alphabet.size())
{
  sstring alphaName;
  sstring betaName;
  if (branchPrefix)
    {
      alphaName << branchPrefix;
      betaName << branchPrefix;
    }
  alphaName << "Alpha";
  betaName << "Beta";
  alpha = pscope.new_boolean_group (alphaName.c_str());
  beta = pscope.new_boolean_group (betaName.c_str());

  for (int c = 0; c < alphabet.size(); ++c)
    {
      sstring mName;
      if (branchPrefix)
	mName << branchPrefix;
      mName << "M";
      mName << alphabet.int2char_uc (c);
      M[c] = pscope.new_alphabet_group (alphabet, mName.c_str());
    }
}

int SWEXMID_branch_params::alphabetSize() const
{
  return M.size();
}

void SWEXMID_branch_params::set (PScores& pscore, const SWEXMID_params& params, double time)
{
  if (alphabetSize() != params.alphabetSize())
    return;
  EM_matrix& R = (EM_matrix&) params.R;  // cast away const
  // calculate stuff
  const double muT = params.mu * time;
  const double lambdaT = params.lambda() * time;
  const array2d<Prob> params_M = R.create_conditional_substitution_matrix (time);
  // set pscores
  pscore[alpha[0]] = Prob2Score (1. - Nats2Prob (-muT));
  pscore[alpha[1]] = Nats2Score (-muT);
  pscore[beta[0]] = Nats2Score (-lambdaT);
  pscore[beta[1]] = Prob2Score (1. - Nats2Prob (-lambdaT));
  for (int x = 0; x < alphabetSize(); ++x)
    for (int y = 0; y < alphabetSize(); ++y)
      pscore[M[x][y]] = Prob2Score (params_M(x,y));
}

// SWEXMID_AX_params

SWEXMID_AX_params::SWEXMID_AX_params (PScope& pscope, const Alphabet& alphabet, int indelTypes)
  : root (pscope, alphabet, indelTypes),
    x (pscope, alphabet, "x")
{ }

void SWEXMID_AX_params::set (PScores& pscore, const SWEXMID_params& params, double tx)
{
  root.set (pscore, params);
  x.set (pscore, params, tx);
}

// SWEXMID_AXY_params

SWEXMID_AXY_params::SWEXMID_AXY_params (PScope& pscope, const Alphabet& alphabet, int indelTypes)
  : root (pscope, alphabet, indelTypes),
    x (pscope, alphabet, "x"),
    y (pscope, alphabet, "y")
{ }

void SWEXMID_AXY_params::set (PScores& pscore, const SWEXMID_params& params, double tx, double ty)
{
  root.set (pscore, params);
  x.set (pscore, params, tx);
  y.set (pscore, params, ty);
}

void SWEXMID_AXY_params::find_ancestor (const Score_profile& xprof, const Score_profile& yprof,
					const Pairwise_path& path, const PScores& pscore, 
					Score_profile& a_ret, Alignment_path::Row& apath_ret)
{
  // find out if x is older than y
  const bool x_older = pscore[x.alpha[1]] < pscore[y.alpha[1]];
  // get appropriate alignment path
  apath_ret = x_older ? path[1] : path[0];
  // count columns
  int cols = 0;
  for_const_contents (Alignment_path::Row, apath_ret, c)
    if (*c)
      ++cols;
  // reserve size for a
  a_ret.clear();
  a_ret.reserve (cols);

  // create a
  vector<int> xy_coords = path.create_seq_coords();
  for (int c = 0; c < path.columns(); ++c)
    {
      if (apath_ret[c])
	{
	  a_ret.push_back (Symbol_score_map());
	  Symbol_score_map& ssm = a_ret.back();
	  for (int i = 0; i < root.alphabetSize(); ++i)
	    {
	      Score& asc = ssm[i];
	      if (path[0][c])
		{
		  Score xsc = -InfinityScore;
		  for_const_contents (Symbol_score_map, xprof[xy_coords[0]], xss)
		    ScorePSumAcc (xsc, ScorePMul (xss->second, pscore[x.M[i][xss->first]]));
		  ScorePMulAcc (asc, xsc);
		}
	      if (path[1][c])
		{
		  Score ysc = -InfinityScore;
		  for_const_contents (Symbol_score_map, yprof[xy_coords[1]], yss)
		    ScorePSumAcc (ysc, ScorePMul (yss->second, pscore[y.M[i][yss->first]]));
		  ScorePMulAcc (asc, ysc);
		}
	    }
	}
      path.inc_seq_coords (xy_coords, c);
    }
}


// SWEXMID_branch_funcs

SWEXMID_branch_funcs::SWEXMID_branch_funcs (const SWEXMID_root_params& root, const SWEXMID_branch_params& branch,
					    const SWEXMID_root_funcs& root_funcs)
  : notKappaR (root.indelTypes()),
    WD (root.indelTypes()),
    MI (root.indelTypes()),
    DW (root.indelTypes()),
    DX (root.indelTypes()),
    DWX (root.indelTypes()),
    IW (root.indelTypes()),
    II (root.indelTypes(), root.indelTypes()),
    DI (root.indelTypes(), root.indelTypes()),
    MWD (root.indelTypes()),
    IWM (root.indelTypes()),
    DWM (root.indelTypes()),
    IWD (root.indelTypes(), root.indelTypes()),
    DWXD (root.indelTypes(), root.indelTypes())
{
  // independent of indelType
  gamma = (root_funcs.indelRatio + branch.beta[1]) / (root_funcs.indelRatio + 1);
  notGamma = branch.beta[0] / (root_funcs.indelRatio + 1);
  WM = branch.alpha[1];
  MW = branch.beta[0];
  MWM = MW * WM;

  // dependent on one indelType
  for (int i = 0; i < root.indelTypes(); ++i)
    {
      WD[i] = root.q[i] * branch.alpha[0];
      MI[i] = root_funcs.qIns[i] * branch.beta[1];
      DW[i] = root.r[i][0] * notGamma;
      DX[i] = root.r[i][1];
      DWX[i] = DW[i] + DX[i];
      IW[i] = root_funcs.notKappaR[i] * branch.beta[0];
      MWD[i] = MW * WD[i];
      IWM[i] = IW[i] * WM;
      DWM[i] = DW[i] * WM;
    }
  
  // dependent on two indelTypes
  for (int i = 0; i < root.indelTypes(); ++i)
    for (int j = 0; j < root.indelTypes(); ++j)
      {
	II(i,j) = root_funcs.notKappaR[i] * MI[j];
	DI(i,j) = root.r[i][0] * root_funcs.qIns[j] * gamma;
	IWD(i,j) = IW[i] * WD[j];
	DWXD(i,j) = DW[i] * WD[j];
	if (i == j)
	  {
	    II(i,i) += root.kappa[1] * root.r[i][1];
	    DWXD(i,i) += DX[i];
	  }
      }
}

// SWEXMID_AX_funcs

SWEXMID_AX_funcs::SWEXMID_AX_funcs (const SWEXMID_AX_params& params)
  : root (params.root),
    x (params.root, params.x, root)
{ }

// SWEXMID_AXY_funcs

SWEXMID_AXY_funcs::SWEXMID_AXY_funcs (const SWEXMID_AXY_params& params)
  : root (params.root),
    x (params.root, params.x, root),
    y (params.root, params.y, root)
{ }

// SWEXMID_AX_state_indices

SWEXMID_AX_state_indices::SWEXMID_AX_state_indices (int indelTypes)
  : indelTypes (indelTypes)
{ }

int SWEXMID_AX_state_indices::ss() const { return Start; }
int SWEXMID_AX_state_indices::ww() const { return End; }
int SWEXMID_AX_state_indices::IM() const { return 0; }
int SWEXMID_AX_state_indices::iI (int x) const { return 1 + x; }
int SWEXMID_AX_state_indices::Id (int x) const { return 1 + x + indelTypes; }
int SWEXMID_AX_state_indices::total_EHMM_states() const { return 1 + 2*indelTypes; }

// SWEXMID_AXY_state_indices

SWEXMID_AXY_state_indices::SWEXMID_AXY_state_indices (int indelTypes)
  : indelTypes (indelTypes)
{ }

int SWEXMID_AXY_state_indices::sss() const { return Start; }
int SWEXMID_AXY_state_indices::www() const { return End; }
int SWEXMID_AXY_state_indices::IMM() const { return 0; }
int SWEXMID_AXY_state_indices::iIw (int x) const { return 1 + x; }
int SWEXMID_AXY_state_indices::imI (int y) const { return 1 + y + indelTypes; }
int SWEXMID_AXY_state_indices::IMd (int y) const { return 1 + y + 2*indelTypes; }
int SWEXMID_AXY_state_indices::IdM (int x) const { return 1 + x + 3*indelTypes; }
int SWEXMID_AXY_state_indices::iIx (int x, int y) const { return 1 + x + (y + 4)*indelTypes; }
int SWEXMID_AXY_state_indices::idI (int x, int y) const { return 1 + x + (y + 4 + indelTypes)*indelTypes; }
int SWEXMID_AXY_state_indices::Idd (int x, int y) const { return 1 + x + (y + 4 + 2*indelTypes)*indelTypes; }
int SWEXMID_AXY_state_indices::total_EHMM_states() const { return 1 + (4 + 3*indelTypes)*indelTypes; }

// SWEXMID_AX_PHMM

SWEXMID_AX_PHMM::SWEXMID_AX_PHMM (int indelTypes, const Alphabet& alphabet)
  : SWEXMID_AX_state_indices (indelTypes),
    Pair_PHMM (total_EHMM_states(), alphabet),
    pscore(),
    params (pscore, alphabet, indelTypes),
    funcs (params)
{
  // set state types

  // set transitions
  //   ww=>ee (W->E, W->E) is implicit

	//   ss=>sI (S, S->I)
	//   ss=>ww (S->W, S->W)
	//   ss=>IM (S->I, S->W->M)
	//   ss=>Id (S->I, S->W->D)

	//   sI=>sI (S, I->I)
	//   sI=>ww (S->W, I->W)
	//   sI=>IM (S->I, I->W->M)
	//   sI=>Id (S->I, I->W->D)

	//   IM=>ww (I->W, M->W)
	//   IM=>IM (I->I, M->W->M)
	//   IM=>Id (I->I, M->W->D)
	//   IM=>iI (I, M->I)

	//   Id=>ww (I->W, D->W)
	//   Id=>IM (I->I, D->W->M)
	//   Id=>Id (I->I, D->W->D)
	//   Id=>iI (I, D->I)

	//   iI=>ww (I->W, I->W)
	//   iI=>IM (I->I, I->W->M)
	//   iI=>Id (I->I, I->W->D)
	//   iI=>iI (I, I->I)

  // set emissions
}

Pairwise_path SWEXMID_AX_PHMM::optimal_alignment (const Score_profile& aprof, const Score_profile& xprof,
						  const SWEXMID_params& params, double time)
{
  Pairwise_path path;
  return path;
}

// SWEXMID_AXY_PHMM

SWEXMID_AXY_PHMM::SWEXMID_AXY_PHMM (int indelTypes, const Alphabet& alphabet)
  : SWEXMID_AXY_state_indices (indelTypes),
    Pair_PHMM (total_EHMM_states(), alphabet),
    pscore(),
    params (pscore, alphabet, indelTypes),
    funcs (params)
{
  // print log message
  CTAG(2,SWEXMID) << "Creating AXY Pair HMM in functional form\n";

  // set state types
  const PFunc zero = 0.;
  init_emit (IMM(), EmitXY, zero);
  for (int i = 0; i < indelTypes; ++i)
    {
      init_emit (IMd(i), EmitX, zero);
      init_emit (iIw(i), EmitX, zero);
      init_emit (IdM(i), EmitY, zero);
      init_emit (imI(i), EmitY, zero);
      for (int j = 0; j < indelTypes; ++j)
	{
	  init_emit (iIx(i,j), EmitX, zero);
	  init_emit (idI(i,j), EmitY, zero);
	  init_emit (Idd(i,j), Null, zero);
	}
    }

  // set up some PVar references
  const PVar& IW = params.root.kappa[0];
  const PVar& II = params.root.kappa[1];

  // set transitions
  //   www=>eee (W->E, W->E, W->E) is implicit

  //   sss=>www (S->W, S->W, S->W)
  //   sss=>IMM (S->I, S->W->M, S->W->M)
  //   sss=>sIw (S, S->I, S->W)
  transition (sss(), www()) = IW * funcs.x.MW * funcs.y.MW;
  transition (sss(), IMM()) = II * funcs.x.MWM * funcs.y.MWM;

  //   IMM=>IMM (I->I, M->W->M, M->W->M)
  //   IMM=>www (I->W, M->W, M->W)
  transition (IMM(), IMM()) = II * funcs.x.MWM * funcs.y.MWM;
  transition (IMM(), www()) = IW * funcs.x.MW * funcs.y.MW;
  
  for (int i = 0; i < indelTypes; ++i)
    {
      //   sss=>sIw (S, S->I, S->W)
      transition (sss(), iIw(i)) = 1 * funcs.x.MI[i] * funcs.y.MW;
      
      //   sss=>ssI (S, S, S->I)
      //   sss=>IdM (S->I, S->W->D, S->W->M)
      //   sss=>IMd (S->I, S->W->M, S->W->D)
      transition (sss(), imI(i)) = 1 * 1 * funcs.y.MI[i];
      transition (sss(), IdM(i)) = II * funcs.x.MWD[i] * funcs.y.MWM;
      transition (sss(), IMd(i)) = II * funcs.x.MWD[i] * funcs.y.MWM;

      //   IMM=>iIw (I, M->I, M->W)
      //   IMM=>IdM (I->I, M->W->D, M->W->M)
      //   IMM=>IMd (I->I, M->W->M, M->W->D)
      //   IMM=>imI (I, M, M->I)
      transition (IMM(), iIw(i)) = 1 * funcs.x.MI[i] * funcs.y.MW;
      transition (IMM(), IdM(i)) = II * funcs.x.MWD[i] * funcs.y.MWM;
      transition (IMM(), IMd(i)) = II * funcs.x.MWM * funcs.y.MWD[i];
      transition (IMM(), imI(i)) = 1 * 1 * funcs.y.MI[i];

      //   iIw=>www (I->W, I->W, W)
      //   iIw=>IMM (I->I, I->W->M, W->M)
      transition (iIw(i), www()) = 1 * funcs.x.MI[i] * funcs.y.MW;
      transition (iIw(i), IMM()) = II * funcs.x.IWM[i] * 1;

      //   IdM=>www (I->W, D->W, M->W)
      //   IdM=>IMM (I->I, D->W->M, M->W->M)
      transition (IdM(i), www()) = IW * funcs.x.DWX[i] * funcs.y.MW;
      transition (IdM(i), IMM()) = II * funcs.x.DWM[i] * funcs.y.MWM;

	//   IMd=>www (I->W, M->W, D->W)
	//   IMd=>IMM (I->I, M->W->M, D->W->M)
      transition (IMd(i), www()) = IW * funcs.x.MW * funcs.y.DWX[i];
      transition (IMd(i), IMM()) = II * funcs.x.MWM * funcs.y.DWM[i];

      //   imI=>www (I->W, M->W, I->W)
      //   imI=>IMM (I->I, M->W->M, I->W->M)
      transition (imI(i), www()) = IW * funcs.x.MW * funcs.y.IW[i];
      transition (imI(i), IMM()) = II * funcs.x.MWM * funcs.y.IWM[i];

      for (int j = 0; j < indelTypes; ++j)
	{
	  //   sss=>Idd (S->I, S->W->D, S->W->D)
	  transition (sss(), Idd(i,j)) = II * funcs.x.MWD[i] * funcs.y.MWD[j];

	  //   IMM=>Idd (I->I, M->W->D, M->W->D)
	  transition (IMM(), Idd(i,j)) = II * funcs.x.MWD[i] * funcs.y.MWD[j];

	  //   iIw=>iIw (I, I->I, W)
	  //   iIw=>IdM (I->I, I->W->D, W->M)
	  //   iIw=>IMd (I->I, I->W->M, W->D)
	  transition (iIw(i), iIw(j)) = 1 * funcs.x.II(i,j) * 1;
	  transition (iIw(i), IdM(j)) = II * funcs.x.IWD(i,j) * funcs.y.WM;
	  transition (iIw(i), IMd(j)) = II * funcs.x.IWM[i] * funcs.y.WD[j];

	  //   IdM=>iIw (I, D->I, M->W)
	  //   IdM=>IdM (I->I, D->W->D, M->W->M)
	  //   IdM=>IMd (I->I, D->W->M, M->W->D)
	  //   IdM=>idI (I, D, M->I)
	  transition (IdM(i), iIw(j)) = 1 * funcs.x.DI(i,j) * funcs.y.MW;
	  transition (IdM(i), IdM(j)) = II * funcs.x.DWXD(i,j) * funcs.y.MWM;
	  transition (IdM(i), IMd(j)) = II * funcs.x.DWM[i] * funcs.y.MWD[j];
	  transition (IdM(i), idI(i,j)) = 1 * 1 * funcs.y.MI[j];

	  //   IMd=>iIw (I, M->I, D->W)
	  //   IMd=>IdM (I->I, M->W->D, D->W->M)
	  //   IMd=>IMd (I->I, M->W->M, D->W->D)
	  //   IMd=>imI (I, M, D->I)
	  //   IMd=>iIx (I, M->I, D->X)
	  transition (IMd(i), iIw(j)) = 1 * funcs.x.MI[i] * funcs.y.DW[j];
	  transition (IMd(i), IdM(j)) = II * funcs.x.MWD[i] * funcs.y.DWM[j];
	  transition (IMd(i), IMd(j)) = II * funcs.x.MWM * funcs.y.DWXD(i,j);
	  transition (IMd(i), imI(j)) = 1 * 1 * funcs.y.DI(i,j);
	  transition (IMd(j), iIx(i,j)) = 1 * funcs.x.MI[i] * funcs.y.DX[j];

	  //   Idd=>www (I->W, D->W, D->W)
	  //   Idd=>IMM (I->I, D->W->M, D->W->M)
	  transition (Idd(i,j), www()) = IW * funcs.x.DWX[i] * funcs.y.DWX[j];
	  transition (Idd(i,j), IMM()) = IW * funcs.x.DWM[i] * funcs.y.DWM[j];

	  //   idI=>www (I->W, D->W, I->W)
	  //   idI=>IMM (I->I, D->W->M, I->W->M)
	  transition (idI(i,j), www()) = IW * funcs.x.DWX[i] * funcs.y.IW[j];
	  transition (idI(i,j), IMM()) = II * funcs.x.DWM[i] * funcs.y.IWM[j];

	  //   imI=>iIw (I, M->I, I->W)
	  //   imI=>IdM (I->I, M->W->D, I->W->M)
	  //   imI=>IMd (I->I, M->W->M, I->W->D)
	  //   imI=>imI (I, M, I->I)
	  transition (imI(i), iIw(j)) = 1 * funcs.x.MI[i] * funcs.y.IW[j];
	  transition (imI(i), IdM(j)) = II * funcs.x.MWD[i] * funcs.y.IWM[j];
	  transition (imI(i), IMd(j)) = II * funcs.x.MWM * funcs.y.IWD(i,j);
	  transition (imI(i), imI(j)) = 1 * 1 * funcs.y.II(i,j);

	  //   iIx=>IMd (I->I, I->W->M, W->D)
	  transition (iIx(i,j), IMd(j)) = II * funcs.x.IWM[i] * 1;
	  
	  for (int k = 0; k < indelTypes; ++k)
	    {
	      //   iIw=>Idd (I->I, I->W->D, W->D)
	      transition (iIw(i), Idd(j,k)) = II * funcs.x.IWD(i,j) * funcs.y.WD[k];

	      //   IdM=>Idd (I->I, D->W->D, M->W->D)
	      transition (IdM(i), Idd(j,k)) = II * funcs.x.DWXD(i,j) * funcs.y.MWD[k];
	      
	      //   IMd=>Idd (I->I, M->W->D, D->W->D)
	      transition (IMd(i), Idd(j,k)) = II * funcs.x.MWD[i] * funcs.y.DWXD(j,k);

	      //   Idd=>iIw (I, D->I, D->W)
	      //   Idd=>IdM (I->I, D->W->D, D->W->M)
	      //   Idd=>IMd (I->I, D->W->M, D->W->D)
	      transition (Idd(i,j), iIw(k)) = 1 * funcs.x.DI(i,k) * funcs.y.DW[j];
	      transition (Idd(i,j), IdM(k)) = II * funcs.x.DWXD(i,k) * funcs.y.DWM[j];
	      transition (Idd(i,j), IMd(k)) = II * funcs.x.DWM[i] * funcs.y.DWXD(j,k);

	      //   idI=>iIw (I, D->I, I->W)
	      //   idI=>IdM (I->I, D->W->D, I->W->M)
	      //   idI=>IMd (I->I, D->W->M, I->W->D)
	      //   idI=>idI (I, D, I->I)
	      transition (idI(i,j), iIw(k)) = 1 * funcs.x.DI(i,k) * funcs.y.IW[j];
	      transition (idI(i,j), IdM(k)) = II * funcs.x.DWXD(i,k) * funcs.y.IWM[j];
	      transition (idI(i,j), IMd(k)) = II * funcs.x.DWM[i] * funcs.y.IWD(j,k);
	      transition (idI(i,j), idI(i,k)) = 1 * 1 * funcs.y.II(j,k);
	      
	      //   imI=>Idd (I->I, M->W->D, I->W->D)
	      transition (imI(i), Idd(j,k)) = II * funcs.x.MWD[j] * funcs.y.IWD(i,k);
	      
	      //   iIx=>iIx (I, I->I, X)
	      //   iIx=>Idd (I->I, I->W->D, X->D)
	      transition (iIx(i,j), iIx(k,j)) = 1 * funcs.x.II(i,k) * 1;
	      transition (iIx(i,j), Idd(k,j)) = II * funcs.x.IWD(i,k) * 1;

	      //   Idd=>iIx (I, D->I, D->X)
	      //   Idd=>idI (I, D, D->I)
	      transition (Idd(i,j), iIx(k,j)) = 1 * funcs.x.DI(i,k) * funcs.y.DX[j];
	      transition (Idd(i,j), idI(i,k)) = 1 * 1 * funcs.y.DI(j,k);
	      
	      for (int l = 0; l < indelTypes; ++l)
		{
		  //   Idd=>Idd (I->I, D->W->D, D->W->D)
		  transition (Idd(i,j), Idd(k,l)) = II * funcs.x.DWXD(i,k) * funcs.y.DWXD(j,l);

		  //   idI=>Idd (I->I, D->W->D, I->W->D)
		  transition (idI(i,j), Idd(k,l)) = II * funcs.x.DWXD(i,k) * funcs.y.IWD(j,l);
		}
	    }
	}
    }
  
  // set emissions
  for (int a = 0; a < alphabet.size(); ++a)
    for (int x = 0; x < alphabet.size(); ++x)
      for (int y = 0; y < alphabet.size(); ++y)
	pair_emit[IMM()](x,y) += params.root.pi[a] * params.x.M[a][x] * params.y.M[a][y];
  for (int c = 0; c < alphabet.size(); ++c)
    {
      PFunc px = 0.;
      PFunc py = 0.;
      for (int a = 0; a < alphabet.size(); ++a)
	{
	  px += params.root.pi[a] * params.x.M[a][c];
	  py += params.root.pi[a] * params.y.M[a][c];
	}
      for (int i = 0; i < indelTypes; ++i)
	{
	  single_emit[IMd(i)][c] = px;
	  single_emit[iIw(i)][c] = px;
	  single_emit[imI(i)][c] = py;
	  single_emit[IdM(i)][c] = py;
	  for (int j = 0; j < indelTypes; ++j)
	    {
	      single_emit[iIx(i,j)][c] = px;
	      single_emit[idI(i,j)][c] = py;
	    }
	}
    }
}

SWEXMID_AXY_PHMM::AXY_inference SWEXMID_AXY_PHMM::optimal_alignment (const Score_profile& xprof,
								     const Score_profile& yprof,
								     const SWEXMID_params& swexmid, double tx, double ty)
{
  // it's best to have AX and AY branches the same length for pairwise alignment
  // (actually, SWEXMID_AX_PHMM might be a better class for this, as it has fewer states & so is quicker)
  CTAG(2,SWEXMID) << "Evaluating AXY Pair HMM scores\n";
  const double t = (tx + ty) / 2;
  params.set (pscore, swexmid, t, t);
  Pair_HMM_scores hmm = eval_hmm_scores (pscore);
  CTAG(2,SWEXMID) << "Eliminating null states from Pair HMM\n";
  hmm.eliminate_null_states();
  CTAG(2,SWEXMID) << "Creating forward-backward matrix\n";
  Pair_forward_backward_DP_matrix fb (hmm, xprof, yprof);
  CTAG(2,SWEXMID) << "Getting posterior match probabilities\n";
  Post_pair_HMM post (fb);
  if (CTAGGING(1,SWEXMID DOTPLOT))
    {
      Pair_envelope pair_env;
      pair_env.initialise_from_posterior_matrix (post, 1./8.);  // use min prob of 1/8 for max ANSI color contrast
      Biosequence xseq, yseq;
      alphabet->score2seq (xprof, xseq);
      alphabet->score2seq (yprof, yseq);
      CL << "Posterior matrix:\n";
      pair_env.render_dotplot (CL, xseq, yseq);
    }
  CTAG(2,SWEXMID) << "Creating optimal accuracy matrix\n";
  Optimal_accuracy_DP_matrix optacc (post);
  Pairwise_path xy_path = optacc.traceback();

  AXY_inference axy_inf;
  axy_inf.axy_path = xy_path;
  axy_inf.axy_path.insert_rows (0);

  Score_profile a;
  Alignment_path::Row apath;
  params.set (pscore, swexmid, tx, ty);  // set AX and AY branches back to true length
  params.find_ancestor (xprof, yprof, xy_path, pscore, axy_inf.a, axy_inf.axy_path[0]);
  
  return axy_inf;
}

Pairwise_path SWEXMID_AXY_PHMM::viterbi_alignment (const Score_profile& xprof,
						   const Score_profile& yprof,
						   const SWEXMID_params& swexmid,
						   double txy)
{
  CTAG(2,SWEXMID) << "Evaluating AXY Pair HMM scores\n";
  params.set (pscore, swexmid, txy/2, txy/2);
  Pair_HMM_scores hmm = eval_hmm_scores (pscore);
  CTAG(2,SWEXMID) << "Eliminating null states from Pair HMM\n";
  hmm.eliminate_null_states();
  CTAG(2,SWEXMID) << "Creating Viterbi matrix\n";
  Pair_Viterbi_DP_matrix viterbi (hmm, xprof, yprof);
  CTAG(2,SWEXMID) << "Getting Viterbi traceback alignment\n";
  const vector<int> state_path =  viterbi.optimal_state_path();
  const Pairwise_path xy_path = hmm.convert_state_path_to_alignment (state_path);
  return xy_path;
}

Tree_alignment SWEXMID_AXY_PHMM::progressive_alignment (const Sequence_database_index& seqdb_index, const SWEXMID_params& params)
{
  const TKF_params tkf_params = params.tkf_params();
  const double default_time = .5;
  const double tres = .01;
  const double tmax = 10;
  
  array2d<double> dist (seqdb_index.size(), seqdb_index.size());
  for (int i = 0; i < seqdb_index.size(); ++i)
    for (int j = i + 1; j < seqdb_index.size(); ++j)
      {
	// make copies of some vars
	const Named_profile& npi = *seqdb_index.profile[i];
	const Named_profile& npj = *seqdb_index.profile[j];
	const Score_profile& iprof = npi.prof_sc;
	const Score_profile& jprof = npj.prof_sc;
	// create the distance function object
	const Pairwise_path xy_path = viterbi_alignment (iprof, jprof, params, default_time);
	// print log message
	const Alignment align (xy_path, npi, npj);
	if (CTAGGING(3,SWEXMID))
	  {
	    CL << "Pairwise alignment of '" << align.row_name[0] << "' and '" << align.row_name[1] << "':\n";
	    align.write_MUL (CL, *alphabet);
	  }
	// create distance function object and TKF_functions
	TKF_aligned_counts_function f (iprof, jprof, xy_path, alphabet);
	TKF_functions tkf_funcs (tkf_params, f, FALSE, TRUE, tres);
	// do Brent maximisation
	double t1, t2, t3, tbest, fbest;
	bracket_maximum (tkf_funcs.log_like, t1, t2, t3, 0., tmax);
	brent_deriv (tkf_funcs.log_like, tkf_funcs.log_like_dt, t1, t2, t3, tres, tbest, fbest);
	// set entry
	dist(i,j) = dist(j,i) = tbest;
	// print log message
	CTAG(3,SWEXMID) << "Estimated time " << tbest << " for '" << npi.name << "' and '" << npj.name << "'\n";
      }

  // Build the tree
  // We assume that node == row
  NJ_tree nj;
  nj.build (seqdb_index.name, dist);
  nj.assert_tree_is_binary();

  // print log message
  if (CTAGGING(5,SWEXMID))
    {
      CL << "Estimated tree:\n";
      nj.write (CL);
    }

  // Do progressive alignment
  map<Phylogeny::Node,Score_profile> prof;
  Alignment_path::Decomposition decomp;
  for_rooted_nodes_post (nj, branch_iter)
    {
      const Phylogeny::Node parent = (*branch_iter).first;
      const Phylogeny::Node node = (*branch_iter).second;
      const vector<Phylogeny::Node> child (nj.children_begin(node,parent), nj.children_end(node,parent));
      if (child.size() == 0)
	prof[node] = Score_profile (seqdb_index.profile[node]->prof_sc);  // copy leaf sequence profiles into map
      else if (child.size() == 2)
	{
	  const Phylogeny::Node x = child[0];
	  const Phylogeny::Node y = child[1];
	  const AXY_inference axy = optimal_alignment (prof[x], prof[y], params,
						       nj.branch_length(node,x), nj.branch_length(node,y));
	  prof[node] = axy.a;
	  decomp[Alignment_path::Row_pair (node, x)] = Pairwise_path (axy.axy_path, 0, 1, TRUE);
	  decomp[Alignment_path::Row_pair (node, y)] = Pairwise_path (axy.axy_path, 0, 2, TRUE);

	  // print log message
	  vector<sstring> row_name (3);
	  row_name[0] = nj.node_specifier (node);
	  row_name[1] = nj.node_specifier (x);
	  row_name[2] = nj.node_specifier (y);
	  if (CTAGGING(4,SWEXMID))
	    {
	      CL << "Alignment of node '" << row_name[0] << "' and descendants:\n";
	      axy.axy_path.show (CL, row_name);
	    }
	}
      else
	THROWEXPR ("Not a binary tree");
    }

  // print log message
  CTAG(5,SWEXMID) << "Completed progressive alignment\n";
  // make Tree_alignment
  Tree_alignment tree_align;
  tree_align.set_tree (nj);
  tree_align.make_empty_alignment();
  for (int i = 0; i < seqdb_index.size(); ++i)
    tree_align.align.prof[i] = &seqdb_index.profile[i]->prof_sc;
  tree_align.align.path.compose (decomp, TRUE);

  return tree_align;
}

