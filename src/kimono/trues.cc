#include "kimono/trues.h"

True_cluster_set::True_cluster_set (const sstring& filename, const Gibbs_alignment_factory& gibbs_factory, const Gaussian_cluster_factory& gauss_factory, const vector<int>& model_prior_sc)
  : dataset (gibbs_factory.dataset),
    model_prior_sc (model_prior_sc),
    gene2cluster (dataset.profile.size(), -1),
    motif_offset (dataset.profile.size()),
    motif_reversed (dataset.profile.size()),
    clusters (0)
{
  if (filename == sstring("")) return;

  CLOG(6) << "Reading trues file '" << filename << "'\n";
  ifstream trues_stream (filename.c_str());
  if (!trues_stream) { sstring e; e << "Trues file '" << filename << "' not found"; THROW Standard_exception (e); }
  
  sstring s;
  Phonebook cluster_index;
  while (s.getline(trues_stream).size())
    {
      vector<sstring> f = s.split();
      int g = gauss_factory.params.table.probe_idx.lookup (f[0]);
      if (g < 0 ? 1 : !dataset.includes_row(g)) { sstring e; e << "Sequence '" << f[0] << "', listed in trues file '" << filename << "', was not found\n"; THROW Standard_exception (e); }
      int c = cluster_index.lookup (f[1]);
      int m_start = atoi (f[2].c_str());
      int m_end = atoi (f[3].c_str());
      if (c == -1) { c = cluster_index[f[1]] = clusters++; cluster_name.push_back (f[1]); cluster_motif_size.push_back (abs (m_end - m_start) + 1); }
      gene2cluster[g] = c;
      motif_offset[g] = m_start < m_end ? m_start : dataset.profile[g]->size()+1 - m_start;
      motif_reversed[g] = m_start > m_end;
    }
  
  if (CLOGGING(6))
    {
      // calculate trues score for logfile
      //
      
      SE_cluster_model se (gibbs_factory, gauss_factory, clusters);
      se.model_prior_sc = model_prior_sc;

      for (int g = 0; g < (int) gene2cluster.size(); ++g)
	se.assign_gene_to_model (g, gene2cluster[g] + 1);
      se.update_model_multiplicities();
      
      for (int c = 0; c < clusters; ++c)
	{
	  for (int g = 0; g < (int) gene2cluster.size(); ++g)
	    if (gene2cluster[g] == c)
	      {
		se.gibbs[c]->align_row (g, -motif_offset[g] + 1, motif_reversed[g]);
		se.gibbs[c]->motif_start = 0;
		se.gibbs[c]->motif_end = cluster_motif_size[c] - 1;
	      }
	  se.gibbs[c]->recalculate_counts();
	  se.gauss[c]->optimise_parameters();
	  if (CLOGGING(2)) { CL << "True cluster '" << cluster_name[c] << "':\n"; se.gibbs[c]->display(CL); se.gauss[c]->display(CL); }
	  if (CLOGGING(0)) { CL << "Scores for true cluster '" << cluster_name[c] << "':\n"; se.gibbs[c]->display_scores(CL); se.gauss[c]->display_scores(CL); }
	  if (CLOGGING(6))
	    {
	      se.update_scores();
	      double model_prior = Score2Bits (se.model_prior_sc[c+1]);
	      double gibbs_prior = Score2Bits (se.gibbs_param_sc[c]);
	      double gauss_prior = Score2Bits (se.gauss_param_sc[c]);
	      double gibbs_data = 0;
	      double gauss_data = 0;
	      for (int g = 0; g < se.genes(); ++g)
		if (gene2cluster[g] == c)
		  {
		    gibbs_data += Score2Bits (se.gibbs_sc[g][c+1]);
		    gauss_data += Score2Bits (se.gauss_sc[g][c+1]);
		  }
	      double gibbs_total = gibbs_prior + gibbs_data;
	      double gauss_total = gauss_prior + gauss_data;
	      double total = gibbs_total + gauss_total;
	      CL << "Cluster '" << cluster_name[c] << "' score = " << total << " bits";
	      CL << "  = " << model_prior << " (model prior)\n";
	      CL << " + " << gibbs_total << " (alignment = " << gibbs_prior << " (prior) + " << gibbs_data << " (data))";
	      CL << " + " << gauss_total << " (expression = " << gauss_prior << " (prior) + " << gauss_data << " (data))\n";
	    }
	  CL << "True clusters: total score = " << Nats2Bits (se.log_probability()) << " bits\n";
	}
    }
}

bool True_cluster_set::true_enough (const vector<Gibbs_alignment*> gibbs, double truth_level)
{
  if (clusters == 0) return 0;
  
  vector<int> model2cluster (gibbs.size(), -1);
  int clusters_found = 0;
  bool failed = 0;
  for (int m = 0; m < (int) gibbs.size(); ++m)
    {
      vector<int> cluster_invalidated (clusters, 0);
      for (int g = 0; g < (int) dataset.profile.size(); ++g)
	{
	  if (gibbs[m]->in_prob[g] >= truth_level)
	    {
	      int c = gene2cluster[g];
	      if (cluster_invalidated[c]) { failed = 1; break; }                            // fail if cluster is incomplete
	      if (model2cluster[m] != -1 && model2cluster[m] != c) { failed = 1; break; }   // fail if clusters overlap
	      if (model2cluster[m] == -1) { model2cluster[m] = c; ++clusters_found; }
	    }
	  else
	    if (gene2cluster[g] != -1)
	      {
		if (model2cluster[m] == gene2cluster[g]) { failed = 1; break; }             // fail if cluster is incomplete
		cluster_invalidated [gene2cluster[g]] = 1;
	      }
	}
      if (failed) break;   // exit outer loop
    }
  if (clusters_found != clusters) failed = 1;                                               // fail if not all clusters were found
  if (!failed)
    {
      CLOG(5) << "Cluster assignments correct\n";
      for (int m = 0; m < (int) gibbs.size(); ++m)
	if (model2cluster[m] != -1)
	  {
	    const Gibbs_alignment& align = *gibbs[m];
	    int  c = model2cluster[m];
	    typedef map < pair<bool,int>, int > Dir_offset_map;
	    Dir_offset_map dir_offset_count;
	    int total_count = 0;
	    for (int g = 0; g < (int) dataset.profile.size(); ++g)
	      if (gene2cluster[g] == c)
		{
		  bool same_orientation = align.row_reversed[g] == motif_reversed[g];
		  int  relative_offset  = (same_orientation ? align.motif_start-align.start_col(g) : align.end_col(g)-align.motif_end) - motif_offset[g];
		  ++dir_offset_count [pair<bool,int> (same_orientation, relative_offset)];
		  ++total_count;
		}
	    int true_count = 0;
	    for_contents (Dir_offset_map, dir_offset_count, i) true_count = max (true_count, (*i).second);
	    double p = ((double) true_count) / (double) total_count;
	    CLOG(5) << "Alignment #" << m+1 << " (\"" << cluster_name[c] << "\") is " << ((int) (100.0 * p)) << "% correct\n";
	    if (p < truth_level) failed = 1;                                                // fail if too many rows are misaligned (but don't exit the loop)
	  }
    }
  return !failed;
}

