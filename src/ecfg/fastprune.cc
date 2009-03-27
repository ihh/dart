#include "fastprune.h"

void Fast_prune::prepare (const PHYLIP_tree& tree, const vector<const EM_matrix_base*>& submat, const Column_matrix& cm)
{
  root = tree.root;
  nodes = tree.nodes();
  alph_size = nodes ? ((EM_matrix_base*) submat[0])->alphabet().size() : 0;
  branch.clear();
  Q.clear();
  pi.clear();
  for_rooted_nodes_post (tree, b)
    {
      branch.push_back (pair<int,int> ((*b).first, (*b).second));
      pi.push_back (((EM_matrix_base*) submat[(*b).second])->create_prior());  // cast away const
      if ((*b).first >= 0)  // don't do a substitution matrix for branch (-1,0) -- which must be the last one visited
	{
	  const array2d<Loge>& M = submat[(*b).second]->timepoint_data ((*b).length).M;
	  Q.push_back (Nats2ProbArray2d(M));
	}
    }
  F = vector<vector<Prob> > (nodes, vector<Prob> (alph_size));
  colmat = &cm;
  if (CTAGGING(1,FAST_PRUNE))
    {
      CL << "Fast_prune branch data:\n";
      for (int i = 0; i < (int) branch.size(); ++i)
	{
	  CL << "Branch #" << i << " from " << branch[i].first << " to " << branch[i].second;
	  if (i < (int) branch.size() - 1)
	    CL << ":\nPi = (" << pi[i] << ")\n" << Q[i];
	  else
	    CL << " has no matrix\n";
	}
    }
}

void Fast_prune::clear()
{
  for (int n = 0; n < nodes; ++n)
    for (int s = 0; s < alph_size; ++s)
      F[n][s] = 0.;
}

Prob Fast_prune::prune()
{
  const int B = branch.size();

  const vector<int>::const_iterator clique_root = colmat->root.begin();
  const vector<int>::const_iterator gapped = colmat->gapped.begin();
  const vector<vector<int> >::const_iterator allowed_vec = colmat->allowed.begin();

  bool underflow = false;

  Prob final = 1.;
  for (int b = 0; b < B; ++b)
    {
      const pair<int,int>& br = branch[b];
      const int p = br.first;
      const int c = br.second;

      if (!gapped[c])
	{
	  const vector<int>& allowed = allowed_vec[c];
	  const int n_allowed = allowed.size();
	  const vector<int>::const_iterator allowed_begin = allowed.begin();

	  const vector<Prob>::const_iterator Fc = F[c].begin();

	  if (c == clique_root[c])
	    {
	      double prob = 0.;
	      for (int n = 0; n < n_allowed; ++n)
		{
		  const int i = allowed_begin[n];
		  prob += pi[b][i] * Fc[i];
		}
	      final *= prob;
	    }
	  else
	    {
	      const vector<int>& allowed_p = allowed_vec[p];
	      const int n_allowed_p = allowed_p.size();
	      vector<int>::const_iterator allowed_p_begin = allowed_p.begin();

	      const vector<Prob>::iterator Fp = F[p].begin();
	      const array2d<Prob> Qpc = Q[b];

	      for (int m = 0; m < n_allowed_p; ++m)
		{
		  const int i = allowed_p_begin[m];
		  double prob = 0.;
		  for (int n = 0; n < n_allowed; ++n)
		    {
		      const int j = allowed_begin[n];
		      prob += Qpc(i,j) * Fc[j];
		    }
		  Fp[i] *= prob;
		  if (Fp[i] < ZeroProb)
		    underflow = true;
		}
	    }
	}
    }

  if (underflow)
    CLOGERR << "Warning: underflow occurred in Fast_prune\n";

  if (CTAGGING(1,FAST_PRUNE RATE_EM_UP_DOWN RATE_EM_UP_DOWN_MATRIX))
    {
      CL << "F table:\n";
      for (int n = 0; n < nodes; ++n) {
	CL << n << ": log F[" << n << "] =\t(";
	if (gapped[n])
	  CL << "GAPPED";
	else
	  {
	    set<int> allowed_set (allowed_vec[n].begin(), allowed_vec[n].end());
	    for (int s = 0; s < alph_size; ++s)
	      if (allowed_set.find (s) != allowed_set.end())
		CL << Prob2Nats (F[n][s]) << ' ';
	      else
		CL << "DISALLOWED ";
	  }
	CL << ")\n";
      }
    }

  return final;
}
