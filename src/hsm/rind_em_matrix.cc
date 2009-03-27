#include "hsm/rind_em_matrix.h"

RIND_EM_matrix::RIND_EM_matrix (int C, int A, int max_fork,
  const Tree_alignment_database* align_db, double timepoint_res)
  : EM_matrix (C, A, max_fork, align_db, timepoint_res)
{ }

void RIND_EM_matrix::quick_M (const Update_statistics& stats, bool intra, bool inter)	// previously quick_M_rind()
{
  CTAG(5,RATE_EM RATE_EM_RIND) << "Calculating RIND matrix\n";
  // inter-class substitutions
  if (inter && C > 1)
    {
      // count number of times that each class was substituted to from another class
      // also count total wait time for all classes
      vector<double> subst (C, 0.0);
      vector<double> class_pi_total (C, 0.0);
      for (int ci = 0; ci < C; ++ci)
	for (int a = 0; a < A; ++a)
	  {
	    const int i = ca (ci, a);
	    class_pi_total[ci] += pi[i];
	    subst[ci] += stats.s[i];  // each root node counts as one substitution
	    for (int cj = 0; cj < C; ++cj)
	      if (cj != ci)
		{
		  const int j = ca (cj, a);
		  subst[cj] += stats.u(i,j);
		}
	  }
      const double wait = accumulate (stats.w.begin(), stats.w.end(), 0.0);
      const double inter_total = accumulate (subst.begin(), subst.end(), 0.0);
      for (int ci = 0; ci < C; ++ci)
	for (int a = 0; a < A; ++a)
	  {
	    const int i = ca (ci, a);
	    pi[i] = (pi[i] / class_pi_total[ci]) * (subst[ci] / inter_total);
	    for (int cj = 0; cj < C; ++cj)
	      if (ci == cj)
		Y[a](ci,ci) = (subst[ci] - inter_total) / wait;
	      else
		Y[a](ci,cj) = subst[cj] / wait;
	  }
      // output RIND probabilities & rate
      if (CTAGGING(5,RATE_EM RATE_EM_RIND))
	{
	  CL << "Inter-rate= " << inter_total / wait;
	  CL << " inter-counts:";
	  for (int ci = 0; ci < C; ++ci)
	    CL << '\t' << subst[ci];
	  CL << '\n';
	}
    }
  // intra-class substitutions
  if (intra)
    for (int ci = 0; ci < C; ++ci)
      {
	// count number of times that each residue was substituted to from another residue in this class
	// also count total wait time for this class
	vector<double> subst (A, 0.0);
	double wait = 0;
	double class_pi_total = 0;
	for (int ai = 0; ai < A; ++ai)
	  {
	    const int i = ca (ci, ai);
	    class_pi_total += pi[i];
	    wait += stats.w[i];
	    subst[ai] += stats.s[i];  // each root node counts as one substitution
	    for (int aj = 0; aj < A; ++aj)
	      if (aj != ai)
		subst[aj] += stats.u(i,ca(ci,aj));
	  }
	const double intra_total = accumulate (subst.begin(), subst.end(), 0.0);
	// set up RIND-type matrix
	for (int ai = 0; ai < A; ++ai)
	  {
	    pi[ca(ci,ai)] = class_pi_total * subst[ai] / intra_total;  // preserve relative class weights
	    for (int aj = 0; aj < A; ++aj)
	      if (aj == ai)
		X[ci](ai,ai) = (subst[ai] - intra_total) / wait;
	      else
		X[ci](ai,aj) = subst[aj] / wait;
	  }
	// output RIND probabilities & rate
	if (CTAGGING(5,RATE_EM RATE_EM_RIND))
	  {
	    CL << "Class #" << ci+1;
	    CL << " intra-rate= " << intra_total / wait;
	    CL << " intra-counts:";
	    for (int ai = 0; ai < A; ++ai)
	      CL << '\t' << subst[ai];
	    CL << '\n';
	  }
      }
  // let update_S() take care of symmetry hassles
  update();
}

/* void EM_matrix::set_rind_params (const vector<Prob>& rind_probs, double subst_rate)
{
  if ((int) rind_probs.size() != m()) THROWEXPR ("RIND vector size doesn't match matrix");
  pi = rind_probs;
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	X[ci](ai,ai) = Y[ai](ci,ci) = 0;
	for (int cj = 0; cj < C; ++cj)
	  if (cj != ci)
	    Y[ai](ci,ci) -= (Y[ai](ci,cj) = subst_rate * rind_probs[ca(cj,ai)]);
	for (int aj = 0; aj < A; ++aj)
	  if (aj != ai)
	    X[ci](ai,ai) -= (X[ci](ai,aj) = subst_rate * rind_probs[ca(ci,aj)]);
      }
  update();
} */

/* vector<Prob> EM_matrix::Update_statistics::get_rind_counts (const vector<Prob>& rind_probs, double subst_rate) const
{
  vector<Prob> counts (rind_probs.size(), 0.0);
  for (int i = 0; i < states; ++i)
    {
      counts[i] += w[i] * subst_rate * rind_probs[i];
      for (int j = 0; j < states; ++j)
	if (i != j)
	  counts[j] += u(i,j);
    }
  return counts;
} */

