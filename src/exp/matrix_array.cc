#include "exp/matrix_array.h"
#include "util/vector_output.h"

Matrix_array::Matrix_array (int S, int F, int C, int A, int max_fork, const Tree_alignment_database* align_db, double timepoint_res)
  : max_fork (max_fork)
{
  init_array (S, F, C, A, max_fork, align_db, timepoint_res);
}

void Matrix_array::init_array (int new_S, int new_F, int new_C, int new_A, int max_fork, const Tree_alignment_database* new_align_db, double new_timepoint_res)
{
  S = new_S;
  F = new_F;
  C = new_C;
  A = new_A;
  align_db = new_align_db;
  timepoint_res = new_timepoint_res;
  hsm.clear();
  for (int f = 0; f < F; ++f)
    hsm.push_back (EM_matrix (C, A, max_fork, align_db, timepoint_res));
  sprior = vector<Loge> (S);
  fprior = vector<vector<Loge> > (S, vector<Loge> (F));
}

void Matrix_array::init_alphabet (const Alphabet& base_alphabet)
{
  for (int f = 0; f < F; ++f) hsm[f].init_alphabet (base_alphabet);
}

void Matrix_array::update()
{
  for (int f = 0; f < F; ++f) hsm[f].update();
}

Matrix_array::Update_statistics::Update_statistics (int S, int F, int states)
  : S (S),
    F (F),
    states (states)
{
  clear();
}

void Matrix_array::Update_statistics::clear()
{
  log_likelihood = 0;
  scount = vector<double> (S, 0.0);
  fcount = vector<vector<double> > (S, vector<double> (F, 0.0));
  hsm_stats.clear();
  for (int f = 0; f < F; ++f) hsm_stats.push_back (EM_matrix::Update_statistics (states));
}

void Matrix_array::Update_statistics::clear_DJU()
{
  for (int f = 0; f < F; ++f) hsm_stats[f].clear_DJU();
}

void Matrix_array::Update_statistics::transform (const Matrix_array& arr, bool symmetrise)
{
  for (int f = 0; f < F; ++f) hsm_stats[f].transform (arr.hsm[f], symmetrise);
}

ostream& operator<< (ostream& o, const Matrix_array::Update_statistics& stats)
{
  for (int s = 0; s < stats.S; ++s)
    o << "Sequence class #" << s+1 << ": count=" << stats.scount[s] << ", fixed-clique-class counts (" << stats.fcount[s] << ")\n";
  for (int f = 0; f < stats.F; ++f)
    {
      o << "Model #" << f+1 << " of " << stats.F << ":\n";
      o << stats.hsm_stats[f];
    }
  return o;
}

Matrix_array::Alignment_matrix::Alignment_matrix (const Matrix_array& arr, int n_align)
  : arr (arr),
    align (arr.align_db->tree_align[n_align].align),
    S (arr.S),
    F (arr.F),
    total_log_likelihood (-InfinityLoge),
    scount (S, 0.0),
    fcount (S, vector<double> (F, 0.0))
{
  // construct each EM_matrix
  for (int f = 0; f < F; ++f)
    hsm_matrix.push_back (EM_matrix::Alignment_matrix (arr.hsm[f], n_align));
}

int Matrix_array::Alignment_matrix::n_cliques (int column) const
{
  return hsm_matrix.front().colmat[column].clique.size();
}

void Matrix_array::Alignment_matrix::fill_up()
{
  // call fill_up for each EM_matrix
  for (int f = 0; f < arr.F; ++f)
    hsm_matrix[f].fill_up();
}

void Matrix_array::Alignment_matrix::fill_down (Update_statistics& stats)
{
  // add scount & fcount to stats
  for (int s = 0; s < S; ++s)
    {
      stats.scount[s] += scount[s];
      for (int f = 0; f < F; ++f)
	stats.fcount[s][f] += fcount[s][f];
    }
  // call fill_down for each EM_matrix
  for (int f = 0; f < arr.F; ++f)
    hsm_matrix[f].fill_down (stats.hsm_stats[f]);
}

void Matrix_array::Alignment_matrix::calc_post()
{
  // calculate posterior distribution over sequence classes
  vector<Loge> seq_log_post (S, (Loge) 0);  // evidence for each sequence class
  Loge seq_total = -InfinityLoge;  // total evidence
  vector<vector<vector<vector<Loge> > > > clique_log_post (S, vector<vector<vector<Loge> > > (align.columns()));  // posterior distribution over clique class, conditioned on sequence class
  for (int s = 0; s < S; ++s)
    {
      Loge aln_total = 0;
      for (int col = 0; col < align.columns(); ++col)
	{
	  clique_log_post[s][col] = vector<vector<Loge> > (n_cliques(col), vector<Loge> (F));
	  for (int clique = 0; clique < n_cliques(col); ++clique)
	    {
	      vector<Loge>& cur_log_post = clique_log_post[s][col][clique];
	      Loge clique_total = -InfinityLoge;
	      for (int f = 0; f < F; ++f)
		NatsPSumAcc (clique_total, cur_log_post[f] = NatsPMul (hsm_matrix[f].colmat[col].L[hsm_matrix[f].colmat[col].clique[clique]], arr.fprior[s][f]));
	      for (int f = 0; f < F; ++f)
		{
		  NatsPMulAcc (cur_log_post[f], -clique_total);
		  fcount[s][f] += Nats2Prob (cur_log_post[f]);
		}
	      NatsPMulAcc (aln_total, clique_total);
	    }
	}
      NatsPSumAcc (seq_total, seq_log_post[s] = NatsPMul (aln_total, arr.sprior[s]));
    }
  for (int s = 0; s < S; ++s)
    {
      NatsPMulAcc (seq_log_post[s], -seq_total);
      scount[s] += Nats2Prob (seq_log_post[s]);
    }

  // calculate posterior distribution over fixed clique classes
  for (int col = 0; col < align.columns(); ++col)
    for (int clique = 0; clique < n_cliques(col); ++clique)
      for (int f = 0; f < F; ++f)
	{
	  Loge& cur_log_post = hsm_matrix[f].colmat[col].log_post[hsm_matrix[f].colmat[col].clique[clique]];
	  cur_log_post = -InfinityLoge;
	  // sum over sequences
	  for (int s = 0; s < S; ++s)
	    NatsPSumAcc (cur_log_post, NatsPMul (seq_log_post[s], clique_log_post[s][col][clique][f]));
	}

  total_log_likelihood = seq_total;
}

void Matrix_array::up_down (Update_statistics& stats, bool likelihood_only, bool symmetrise) const
{
  if (align_db == 0) THROWEXPR ("No training set!");
  stats.clear_DJU();
  for (int n_align = 0; n_align < align_db->size(); ++n_align)
    {
      Alignment_matrix alnmat (*this, n_align);
      alnmat.fill_up();
      alnmat.calc_post();
      if (!likelihood_only)
	alnmat.fill_down (stats);
      stats.log_likelihood += alnmat.total_log_likelihood;
    }
  if (!likelihood_only)
    stats.transform (*this, symmetrise);
}

Matrix_array::Update_statistics Matrix_array::get_stats_unequilibrated (bool symmetrise) const
{
  Update_statistics stats (S, F, m());
  up_down (stats, 0, symmetrise);
  if (CTAGGING(2,RATE_EM_STATS_UNEQ)) CL << "Unequilibrated update statistics:\n" << stats;
  return stats;
}

Matrix_array::Update_statistics Matrix_array::get_stats() const
{
  Update_statistics stats = get_stats_unequilibrated (1);  // symmetrise
  equilibrate (stats);
  if (CTAGGING(2,RATE_EM_STATS)) CL << "Equilibrated update statistics:\n" << stats;
  return stats;
}

void Matrix_array::equilibrate (Update_statistics& stats) const
{
  for (int f = 0; f < F; ++f)
    hsm[f].equilibrate (stats.hsm_stats[f]);
}

void Matrix_array::randomise (double seq_dev, double fixed_dev, double prior_dev, double intra_min, double intra_dev, double inter_min, double inter_dev)
{
  const double seq_min = 1;
  const double fixed_min = 1;
  double stot = 0;
  vector<double> ps (S);
  vector<double> pf (F);
  for (int s = 0; s < S; ++s)
    {
      stot += ps[s] = seq_min + Rnd::prob() * seq_dev;
      double ftot = 0;
      for (int f = 0; f < F; ++f)
	ftot += pf[f] = fixed_min + Rnd::prob() * fixed_dev;
      for (int f = 0; f < F; ++f)
	fprior[s][f] = Prob2Nats (pf[f] / ftot);
    }
  for (int s = 0; s < S; ++s)
    sprior[s] = Prob2Nats (ps[s] / stot);
  for (int f = 0; f < F; ++f)
    hsm[f].randomise (prior_dev, intra_min, intra_dev, inter_min, inter_dev);
}

Matrix_array::Update_statistics Matrix_array::single_quick_EM (bool rind, bool intra, bool inter)
{
  CTAG(5,RATE_EM) << "Computing update statistics (E-phase)\n";
  const Update_statistics stats = get_stats();
  CTAG(5,RATE_EM) << "Finding nearest reversible rate matrix (M-phase)\n";
  // update sprior & fprior
  const double stot = accumulate (stats.scount.begin(), stats.scount.end(), 0.0);
  for (int s = 0; s < S; ++s)
    {
      sprior[s] = Prob2Nats (stats.scount[s] / stot);
      const double ftot = accumulate (stats.fcount[s].begin(), stats.fcount[s].end(), 0.0);
      for (int f = 0; f < F; ++f)
	fprior[s][f] = Prob2Nats (stats.fcount[s][f] / ftot);
    }
  // update HSM stats
  for (int f = 0; f < F; ++f)
    if (rind)
      hsm[f].quick_M_rind (stats.hsm_stats[f], intra, inter);
    else
      hsm[f].quick_M (stats.hsm_stats[f], intra, inter);
  // return
  return stats;
}

Loge Matrix_array::iterate_quick_EM (bool rind, bool intra, bool inter, int forgive)
{
  Loge best_log_likelihood = 0;
  Matrix_array_params best_params;
  int dec = 0;
  const int em_max_iter = hsm.front().em_max_iter;   // HACKY
  const double em_min_inc = hsm.front().em_min_inc;   // EQUALLY HACKY
  for (int iter = 0; ; ++iter)
    {
      if (iter >= em_max_iter)
	{
	  CTAG(7,RATE_EM) << "EM hit " << em_max_iter << " iterations; stopping\n";
	  break;
	}
      const Matrix_array_params old_params = *this;
      Update_statistics stats = single_quick_EM (rind, intra, inter);
      const Loge prev_best = best_log_likelihood;
      if (iter == 0 || stats.log_likelihood > best_log_likelihood)
	{
	  best_log_likelihood = stats.log_likelihood;
	  best_params = old_params;
	}
      CTAG(6,RATE_EM) << "EM iteration #" << iter+1 << ": log-likelihood = " << Nats2Bits(stats.log_likelihood) << " bits\n";
      if (iter > 0)
	{
	  const double inc = -(stats.log_likelihood - prev_best) / prev_best;
	  if (inc < em_min_inc)
	    {
	      if (stats.log_likelihood < prev_best)
		CTAG(7,RATE_EM) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_best) << " to " << Nats2Bits(stats.log_likelihood) << " bits during EM\n";
	      if (++dec >= forgive)
		{
		  CTAG(7,RATE_EM) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
		  break;
		}
	    }
	  else
	    dec = 0;
	}
    }
  ((Matrix_array_params&) *this) = best_params;
  update();
  return best_log_likelihood;
}

void Matrix_array::write (ostream& out) const
{
  out << S << " " << F << " " << C << " " << A << "\n";
  for (int s = 0; s < S; ++s)
    out << Nats2Prob(sprior[s]) << " ";
  out << "\n";
  for (int s = 0; s < S; ++s)
    {
      for (int f = 0; f < F; ++f)
	out << Nats2Prob(fprior[s][f]) << " ";
      out << "\n";
    }
  out << "\n";
  for (int f = 0; f < F; ++f)
    hsm[f].write (out);
}

void Matrix_array::read (istream& in)
{
  int new_S, new_F, new_C, new_A;
  in >> new_S >> new_F >> new_C >> new_A;
  init_array (new_S, new_F, new_C, new_A, max_fork, align_db, timepoint_res);
  double tmp;
  for (int s = 0; s < S; ++s)
    {
      in >> tmp;
      sprior[s] = Prob2Nats(tmp);
    }
  for (int s = 0; s < S; ++s)
    for (int f = 0; f < F; ++f)
      {
	in >> tmp;
	fprior[s][f] = Prob2Nats (tmp);
      }
  for (int f = 0; f < F; ++f)
    hsm[f].read (in);
}
