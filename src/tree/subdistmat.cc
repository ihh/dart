#include "tree/subdistmat.h"
#include "util/maximise.h"

Substitution_counts::Substitution_counts (const Alignment& align, int row1, int row2)
{
  Alignment_path::Sequence_coords coords = align.path.create_seq_coords();
  for (int col = 0; col < align.columns(); ++col)
    {
      if (!align.is_gap (row1, col) && !align.is_gap (row2, col))
	for_const_contents (Symbol_score_map, (*align.prof[row1]) [coords[row1]], ss1)
	  for_const_contents (Symbol_score_map, (*align.prof[row2]) [coords[row2]], ss2)
	  res_pair_count[Residue_pair (ss1->first, ss2->first)] += Score2Prob (ScorePMul (ss1->second, ss2->second));

      align.path.inc_seq_coords (coords, col);
    }
}

void Substitution_counts::show (ostream& out) const
{
  out << "Substitution counts:\n";
  for_const_contents (Residue_pair_count, res_pair_count, rpc)
    out << "x=" << rpc->first.first << " y=" << rpc->first.second << " count=" << rpc->second << "\n";
}

Subst_log_like::Subst_log_like (Substitution_matrix_factory& submat_factory, const Substitution_counts& subst_counts)
  : submat_factory (submat_factory), subst_counts (subst_counts)
{ }

Loge Subst_log_like::operator() (double t)
{
  const array2d<Prob> pxy = submat_factory.create_joint_substitution_matrix (t);
  Loge L = 0;
  for_const_contents (Substitution_counts::Residue_pair_count, subst_counts.res_pair_count, rpc)
    NatsPMulAcc (L, rpc->second * Prob2Nats (pxy (rpc->first.first, rpc->first.second)));
  return L;
}

Subst_log_like_dt::Subst_log_like_dt (Substitution_matrix_factory& submat_factory, const Substitution_counts& subst_counts)
  : submat_factory (submat_factory), subst_counts (subst_counts)
{ }

double Subst_log_like_dt::operator() (double t)
{
  const array2d<Prob> pxy = submat_factory.create_joint_substitution_matrix (t);
  const array2d<Prob> dpxy = submat_factory.differentiate_joint_substitution_matrix (t);
  double dLdt = 0;
  for_const_contents (Substitution_counts::Residue_pair_count, subst_counts.res_pair_count, rpc)
    dLdt += rpc->second * dpxy (rpc->first.first, rpc->first.second) / pxy (rpc->first.first, rpc->first.second);
  return dLdt;
}

Subst_dist_func_factory::Subst_dist_func* Subst_dist_func_factory::create_dist_func (const Alignment& align)
{
  return new Subst_dist_func (*this, align);
}

Subst_dist_func_factory::Subst_dist_func_factory (Substitution_matrix_factory& submat_factory,
						  double tres, double tmax)
  : submat_factory (submat_factory),
    tres (tres), tmax (tmax)
{ }

double Subst_dist_func_factory::Subst_dist_func::operator() (int i, int j)
{
  // create the function objects
  Substitution_counts subst_counts (align, i, j);
  if (CTAGGING(1,SUBDIST))
    subst_counts.show (CL);
  Subst_log_like f (factory.submat_factory, subst_counts);
  Subst_log_like_dt dfdt (factory.submat_factory, subst_counts);

  // create function caches
  Function_cache<Subst_log_like> f_cached (f, factory.tres);
  Function_cache<Subst_log_like_dt> dfdt_cached (dfdt, factory.tres);

  // do Brent maximisation
  double t1, t2, t3, tbest, fbest;
  bracket_maximum (f_cached, t1, t2, t3, 0., factory.tmax);
  brent_deriv (f_cached, dfdt_cached, t1, t2, t3, factory.tres, tbest, fbest);

  // return the ML distance estimate
  return tbest;
}
