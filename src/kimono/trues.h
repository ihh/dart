#ifndef TRUES_INCLUDED
#define TRUES_INCLUDED

#include "kimono/secluster.h"

struct True_cluster_set
{
  const Gibbs_dataset& dataset;
  const vector<int>&   model_prior_sc;

  vector<int>     gene2cluster;
  vector<int>     motif_offset;
  vector<int>     motif_reversed;
  vector<sstring> cluster_name;
  vector<int>     cluster_motif_size;
  int             clusters;

  bool            true_enough (const vector<Gibbs_alignment*> gibbs, double truth_level);

  True_cluster_set (const sstring& filename, const Gibbs_alignment_factory& gibbs_factory, const Gaussian_cluster_factory& gauss_factory, const vector<int>& model_prior_sc);
};

#endif
