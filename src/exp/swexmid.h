#ifndef SWEXMID_INCLUDED
#define SWEXMID_INCLUDED

#include "hmm/pairphmm.h"
#include "hsm/em_matrix.h"
#include "handel/tkfparams.h"

// SWEXMID is a seven-state transducer modeling nonoverlapping single-event gaps for the Long Indel model.
// The D->I transition introduced by Knudsen is used.

// Long indel model parameters
struct SWEXMID_params
{
  // data
  Prob kappa; // equilibrium extend probability
  double mu;  // total deletion rate
  vector<Prob> q;  // deletion type probabilities
  vector<Prob> r;  // deletion extend probability, indexed by type
  EM_matrix R;  // substitution rate matrix

  // constructor
  SWEXMID_params (int indelTypes, const Alphabet& alphabet);

  // accessors
  int indelTypes() const;
  int alphabetSize() const;
  double lambda() const;  // = kappa * mu * \sum_i (q[i] * (1-r[i]) / (1-kappa*r[i]))

  // method to estimate TKF91 parameters
  TKF_params tkf_params() const;
};

// Default long indel model parameters for proteins:
//  1-class amino acid substitution matrix (currently uses PAM; switch to HSM1)
//  2-component indel mixture from (Miklos, Lunter & Holmes, 2004)
struct Protein_default_SWEXMID_params : SWEXMID_params
{
  Protein_default_SWEXMID_params();
};

// PGroups independent of branch length
struct SWEXMID_root_params
{
  // data
  PGroup q;  // indexed by indelType
  vector<Boolean_group> r;  // indexed by indelType
  Boolean_group kappa;  // = lambda/mu
  Alphabet_group pi;
  // constructor
  SWEXMID_root_params (PScope& pscope, const Alphabet& alphabet, int indelTypes);
  // accessors
  int indelTypes() const;
  int alphabetSize() const;
  // set method
  void set (PScores& pscore, const SWEXMID_params& params);
};

// PVars dependent on branch length
struct SWEXMID_branch_params
{
  // data
  Boolean_group alpha, beta;  // alpha=exp(-mu*t), beta=1-exp(-lambda*t)  where lambda = indelRatio * mu
  vector<Alphabet_group> M;  // M[i][j] = exp(R*t)_ij
  // constructor
  SWEXMID_branch_params (PScope& pscope, const Alphabet& alphabet, const char* branchPrefix = 0);
  // accessors
  int alphabetSize() const;
  // set method
  void set (PScores& pscore, const SWEXMID_params& params, double time);
};

// PVars for a two-node AX ETree "((X));"
struct SWEXMID_AX_params
{
  // data
  SWEXMID_root_params root;
  SWEXMID_branch_params x;
  // constructor
  SWEXMID_AX_params (PScope& pscope, const Alphabet& alphabet, int indelTypes);
  // set method
  void set (PScores& pscore, const SWEXMID_params& params, double tx);
};

// PVars for a three-node AXY ETree "((X,Y));"
struct SWEXMID_AXY_params
{
  // data
  SWEXMID_root_params root;
  SWEXMID_branch_params x, y;
  // constructor
  SWEXMID_AXY_params (PScope& pscope, const Alphabet& alphabet, int indelTypes);
  // set method
  void set (PScores& pscore, const SWEXMID_params& params, double tx, double ty);
  // method to create an ancestral Score_profile from a pairwise alignment
  void find_ancestor (const Score_profile& x, const Score_profile& y, const Pairwise_path& path, const PScores& pscore, 
		      Score_profile& a_ret, Alignment_path::Row& apath_ret);
};

// PFuncs independent of branch length
struct SWEXMID_root_funcs
{
  // data
  PFunc indelRatio;  // = \sum_i indelRatioTerm[i]     (==lambda/mu)
  vector<PFunc> notKappaR;  // = (1 - kappa * r[i])
  vector<PFunc> indelRatioTerm;  // = kappa * q[i] * (1-r[i]) / (1-kappa*r[i])
  vector<PFunc> qIns;  // = indelRatioTerm[i] / indelRatio
  // constructor
  SWEXMID_root_funcs (const SWEXMID_root_params& params);
};

// PFuncs dependent on branch length
struct SWEXMID_branch_funcs
{
  // data
  // Intermediate functions
  PFunc gamma;  // = (indelRatio + beta) / (indelRatio + 1)
  PFunc notGamma;  // = (1 - beta) / (indelRatio + 1)
  vector<PFunc> notKappaR;  // = (1 - kappa * r[i])
  // Single-transition paths
  PFunc WM, MW;
  vector<PFunc> WD, MI;  // indexed by [destType]
  vector<PFunc> DW, DX, DWX, IW;  // indexed by [srcType]
  array2d<PFunc> II, DI;  // indexed by [srcType][destType]
  // Two-transition paths
  PFunc MWM;
  vector<PFunc> MWD;  // indexed by [destType]
  vector<PFunc> IWM, DWM;  // indexed by [srcType]
  array2d<PFunc> IWD, DWXD;  // indexed by [srcType][destType]
  // constructor
  SWEXMID_branch_funcs (const SWEXMID_root_params& root, const SWEXMID_branch_params& branch, const SWEXMID_root_funcs& root_funcs);
};

// PVars and PFuncs for a two-node AX ETree "((X));"
struct SWEXMID_AX_funcs
{
  // data
  SWEXMID_root_funcs root;
  SWEXMID_branch_funcs x;
  // constructor
  SWEXMID_AX_funcs (const SWEXMID_AX_params& params);
};

// PFuncs for a three-node AXY ETree "((X,Y));"
struct SWEXMID_AXY_funcs
{
  // data
  SWEXMID_root_funcs root;
  SWEXMID_branch_funcs x, y;
  // constructor
  SWEXMID_AXY_funcs (const SWEXMID_AXY_params& params);
};

// Pair HMM state indices for a two-node AX ETree "((X));"
struct SWEXMID_AX_state_indices : Grammar_state_enum
{
  // data
  int indelTypes;
  // constructor
  SWEXMID_AX_state_indices (int indelTypes);
  // methods
  int ss() const; // Start
  int ww() const; // End
  int IM() const;  // type AX
  int iI (int x) const; // type X
  int Id (int x) const;  // type A
  int total_EHMM_states() const;
};

// Pair HMM state indices for a three-node AXY ETree "((X,Y));"
struct SWEXMID_AXY_state_indices : Grammar_state_enum
{
  // data
  int indelTypes;
  // constructor
  SWEXMID_AXY_state_indices (int indelTypes);
  // methods
  int sss() const; // Start
  int www() const; // End
  int IMM() const;  // type AXY
  int imI (int y) const;  // type Y
  int IMd (int y) const; // type AX
  int IdM (int x) const;  // type AY
  int iIw (int x) const; // type X
  int iIx (int x, int y) const; // type X
  int idI (int x, int y) const;  // type Y
  int Idd (int x, int y) const;  // type A
  int total_EHMM_states() const;
};

// Pair PHMM for a two-node AX ETree "((X));"
struct SWEXMID_AX_PHMM : SWEXMID_AX_state_indices, Pair_PHMM
{
  // data
  PScores pscore;
  SWEXMID_AX_params params;
  SWEXMID_AX_funcs funcs;
  // constructor
  SWEXMID_AX_PHMM (int indelTypes, const Alphabet& alphabet);
  // alignment method
  Pairwise_path optimal_alignment (const Score_profile& aprof, const Score_profile& xprof,
				   const SWEXMID_params& params, double time);
};

// Pair PHMM for a three-node AXY ETree "((X,Y));"
struct SWEXMID_AXY_PHMM : SWEXMID_AXY_state_indices, Pair_PHMM
{
  // Data structure for storing ancestral inferences on a three-node AXY ETree "((X,Y));"
  struct AXY_inference
  {
    Alignment_path axy_path;
    Score_profile a;
  };

  // Data
  PScores pscore;
  SWEXMID_AXY_params params;
  SWEXMID_AXY_funcs funcs;

  // Constructor
  SWEXMID_AXY_PHMM (int indelTypes, const Alphabet& alphabet);

  // Alignment methods.
  // Pairwise Viterbi alignment of X and Y, used to estimate distance.
  Pairwise_path viterbi_alignment (const Score_profile& xprof, const Score_profile& yprof,
				   const SWEXMID_params& params, double txy);
  // Optimal-accuracy pairwise alignment & inference of ancestral profiles.
  AXY_inference optimal_alignment (const Score_profile& xprof, const Score_profile& yprof,
				   const SWEXMID_params& params, double tx, double ty);
  // Progressive multiple alignment using a guide tree estimated by Neighbor-Joining.
  // NJ distance matrix is calculated using substitution model only,
  // for a viterbi_alignment() with txy=1
  // (use TKF_aligned_counts_function, copy code from tkfdistmat.cc)
  Tree_alignment progressive_alignment (const Sequence_database_index& seqdb_index, const SWEXMID_params& params);
};

#endif /* SWEXMID_INCLUDED */
