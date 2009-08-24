#include <fstream>
#include <math.h>
#include "newmat/newmatio.h"
#include "tkf/tkfsequence.h"
#include "util/rnd.h"
#include "tree/hasegawa.h"
#include "tree/rate_matrix.h"
#include "util/vector_output.h"
#include "util/logfile.h"

void TKF_sequence::generate (const TKF_params& params)
{
  clear();
  path.reset_rows();
  log_emit_likelihood = log_path_likelihood = 0;

  double e = params.equilibrium_length_prob();
  vector<double> prior = params.submat_factory.create_prior ();

  while (Rnd::decide(e))
    {
      path.append_column (0,1);
      int sym = Rnd::choose(prior);
      push_back (sym);
      log_path_likelihood += log(e);
      log_emit_likelihood += log(prior[sym]);
    }
  log_path_likelihood += log(1-e);
}

void TKF_sequence::generate_predicated_on_length (const TKF_params& params, int len)
{
  clear();
  path.reset_rows();
  log_emit_likelihood = log_path_likelihood = 0;

  double e = params.equilibrium_length_prob();
  vector<double> prior = params.submat_factory.create_prior ();

  for (int i = 0; i < len; ++i)
    {
      path.append_column (0,1);
      int sym = Rnd::choose(prior);
      push_back (sym);
      log_path_likelihood += log(e);
      log_emit_likelihood += log(prior[sym]);
    }
  log_path_likelihood += log(1-e);
}

void TKF_sequence::evolve (const TKF_params& params, const Digitized_biosequence& parent, double time)
{
  clear();
  path.reset_rows();
  log_emit_likelihood = log_path_likelihood = 0;

  double a = params.ancestral_survival_prob (time);
  double m = params.more_descendants_prob (time);
  double o = params.orphan_survival_prob (time);

  array2d<double> submat = params.submat_factory.create_conditional_substitution_matrix (time);
  vector<double>  prior =  params.submat_factory.create_prior ();
  
  while (Rnd::decide(m))
    {
      path.append_column (0,1);
      push_back (Rnd::choose(prior));
      log_path_likelihood += log(m);
    }
  log_path_likelihood += log(1-m);
      
  for_const_contents (Digitized_biosequence, parent, sym)
    {
      double d;
      if (Rnd::decide(a))
	{
	  array2d<double>::Column_view col = submat.column (*sym);
	  vector<double> s (col.begin(), col.end());
	  path.append_column (1,1);
	  int sym = Rnd::choose(s);
	  push_back (sym);
	  log_path_likelihood += log(a);
	  log_emit_likelihood += log(s[sym]);
	  d = m;
	}
      else
	{
	  path.append_column (1,0);
	  log_path_likelihood += log(1-a);
	  d = o;
	}
      while (Rnd::decide(d))
	{
	  path.append_column (0,1);
	  int sym = Rnd::choose(prior);
	  push_back (sym);
	  log_path_likelihood += log(d);
	  log_emit_likelihood += log(prior[sym]);
	  d = m;
	}
      log_path_likelihood += log(1-d);
    }
}


void TKF_emit_align::emit (bool use_expected_len)
{
  make_empty_alignment();    // just to set up the maps

  TKF_sequence tmp_seq;
  Named_profile tmp_np;

  tkf = vector<TKF_sequence> (tree.nodes(), tmp_seq);
  np = vector<Named_profile> (tree.nodes(), tmp_np);
  
  Alignment_path::Decomposition decomp;
  for_iterator (Phylogeny::Branch_iter, b,
		tree.branches_begin (tree.root, -1, 1, 1),
		tree.branches_end ())
    {
      Phylogeny::Node node = (*b).second, parent = (*b).first;
      if (parent == -1)
	{
	  if (use_expected_len)
	    tkf[node].generate_expected_length (params);
	  else
	    tkf[node].generate (params);
	}
      else
	{
	  tkf[node].evolve (params, tkf[parent], (*b).length);
	  decomp[Alignment_path::Row_pair (node2row[parent], node2row[node])] = tkf[node].path;
	}
    }

  if (CLOGGING(8))
    CL << "Alignment generated with log-likelihood " << log_likelihood() / log(2.) << " bits\n";

  if (CLOGGING(5))
    CL << "Log path likelihood = " << log_path_likelihood() / log(2.) << " bits, log emission likelihood = " << log_emit_likelihood() / log(2.) << " bits\n";

  align.path.compose_and_log (decomp, 1);

  for (int i = 0; i < (int) tkf.size(); i++)
    { np[i].prof_sc = Score_profile (tkf[i]); align.prof[i] = &np[i].prof_sc; }

  update_alignment_row_names_from_maps();
  for (int i = 0; i < (int) np.size(); i++) np[i].name = align.row_name[i];
}

double TKF_emit_align::log_path_likelihood() const
{
  double loglike = 0;
  for_const_contents (vector<TKF_sequence>, tkf, t) loglike += (*t).log_path_likelihood;
  return loglike;
}

double TKF_emit_align::log_emit_likelihood() const
{
  double loglike = 0;
  for_const_contents (vector<TKF_sequence>, tkf, t) loglike += (*t).log_emit_likelihood;
  return loglike;
}

int TKF_emit_align::score() const
{
  return Nats2Score (log_likelihood());
}

