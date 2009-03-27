#include "hsm/phylo_em.h"
#include "hsm/em_matrix_base.h"
#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
#include "util/newmat_adaptors.h"
#include "util/vector_output.h"

#define Default_alphabet_name "HSM_alphabet"

#define HSM_MINIMUM_WAIT_TIME    .0001 /* negative wait times are set to this value */
#define HSM_MINIMUM_NONZERO_RATE TINY  /* rates smaller than this value are treated as zero */

typedef EM_matrix_base::Update_statistics Update_statistics;

EM_matrix_base::Timepoint_data::Timepoint_data (int m)
  : M (m, m),
    J (m, m)
{ }

EM_matrix_base::EM_matrix_base (int C, int A, int max_fork,
	const Tree_alignment_database* align_db, double timepoint_res)
  : Piper (max_fork),
    hidden_alphabet ("", 0)
{
  init_matrix (C, A, align_db, timepoint_res);
  em_min_inc = .001;
  em_max_iter = m() * 20;
}

void EM_matrix_base::init_matrix (int new_C, int new_A, const Tree_alignment_database* new_align_db, double new_timepoint_res)
{
  const int states = new_C * new_A;

  init_matrix_params (new_C, new_A);
  init_update_flags (true);
  init_matrix_factory (states);
  
  R = Matrix (states, states);
  log_pi = vector<double> (states);

  // set up the cache
  init_cache (new_align_db, new_timepoint_res);
}

void EM_matrix_base::init_update_flags (bool default_flag)
{
  X_update_flag = vector<array2d<short> > (C, array2d<short> (A, A, default_flag));
  Y_update_flag = vector<array2d<short> > (A, array2d<short> (C, C, default_flag));
}

void EM_matrix_base::init_cache (const Tree_alignment_database* new_align_db, const double new_timepoint_res)
{
  // update the cache parameters
  align_db = new_align_db;
  timepoint_res = new_timepoint_res;
  // update the cache
  timepoint_cache.clear();
  if (align_db != 0)
    for (int n_align = 0; n_align < align_db->size(); ++n_align)
      {
	const Tree_alignment& tree_align = *align_db->tree_align[n_align];
	tree_align.assert_leaves_equal_rows();
	add_branch_lengths_to_cache (tree_align.tree);
      }
}

void EM_matrix_base::add_branch_lengths_to_cache (const PHYLIP_tree& tree)
{
  for_rooted_branches_post (tree, b)
    {
      const int tp = discrete_time ((*b).length);
      Timepoint_cache::iterator map_iter = timepoint_cache.find (tp);
      if (map_iter == timepoint_cache.end())
	timepoint_cache[tp] = Timepoint_data (m());
    }
}

void EM_matrix_base::init_alphabet (const Alphabet& base_alphabet)
{
  if ((int) base_alphabet.size() != A)
    THROWEXPR ("Size of alphabet '" << base_alphabet.name << "' != " << A);
  hidden_alphabet.init_hidden (base_alphabet, C);
}

const Alphabet& EM_matrix_base::alphabet()
{
  if (hidden_alphabet.size() != m()) THROWEXPR ("Alphabet uninitialised");
  return hidden_alphabet;
}

void EM_matrix_base::assign (Substitution_matrix_factory& submat_factory)
{
  const array2d<double> R = submat_factory.rate_matrix();
  init_matrix (1, R.xsize());
  init_alphabet (submat_factory.alphabet());
  X[0] = R;
  pi = submat_factory.create_prior();
  update();
}

const EM_matrix_base::Timepoint_data& EM_matrix_base::timepoint_data (double t) const
{
  const int tp = discrete_time (t);
  Timepoint_cache::const_iterator map_iter = timepoint_cache.find (tp);
  if (map_iter == timepoint_cache.end())
    {
      CTAG(4,RATE_EM) << "Adding timepoint " << t << " to cache\n";
      ((Timepoint_cache&) timepoint_cache)[tp] = Timepoint_data (m());  // cast away const
      map_iter = timepoint_cache.find (tp);
	  // current implementation recalculates all timepoints - could easily be optimized
	  ((EM_matrix_base&) (*this)).update_timepoint_cache(); // cast away const again
    }
  return map_iter->second;
}

void EM_matrix_base::update_R()
{
  for (int i = 0; i < m(); ++i)
    for (int j = 0; j < m(); ++j)
      R(i+1,j+1) = 0;

  for (int c = 0; c < C; ++c)
    for (int ai = 0; ai < A; ++ai)
      for (int aj = 0; aj < A; ++aj)
	R (ca(c,ai) + 1, ca(c,aj) + 1) = X[c] (ai, aj);

  for (int a = 0; a < A; ++a)
    for (int ci = 0; ci < C; ++ci)
      for (int cj = 0; cj < C; ++cj)
	if (ci == cj)
	  R (ca(ci,a) + 1, ca(cj,a) + 1) += Y[a] (ci, cj);
	else
	  R (ca(ci,a) + 1, ca(cj,a) + 1) = Y[a] (ci, cj);
}


void EM_matrix_base::update_log_pi()
{
  for (int i = 0; i < m(); ++i)
    log_pi[i] = Prob2Nats (pi[i]);

  irrev_prior = pi;

  // print out debugging message
  // this is only for logging - not using the QRZ decomposition here
  if (CTAGGING(-1,RATE_EM_CHECK_EQM))
  {
    // decompose R.t() by QR
    Matrix Q = R.t();
    UpperTriangularMatrix upper;
    QRZ (Q, upper);
    vector_to_RowVector_adaptor pi_row (pi);
    CL << "EM_matrix_base::update_sqrt_pi_transform(): QR transform not used here\n";
    CL << "QR upper diagonal of R.t():\n" << upper;
    CL << "pi:\n" << pi_row;
    CL << "pi * R:\n" << pi_row * R;
  }
}

void EM_matrix_base::update_timepoint_cache()
{
  CTAG(3,RATE_EM RATE_EM_CACHE) << "Updating cached M() & J()\n";
  // loop over timepoints
  for_contents (Timepoint_cache, timepoint_cache, tp_data)
    {
      const int tp = tp_data->first;
      Timepoint_data& dat = tp_data->second;
      
      const double T = tp * timepoint_res;

      if (CTAGGING(2,RATE_EM_CACHE))
	CL << "Calculating M() & J() for time T=" << T << "\n";

      // calculate exp_mu_T vector
      vector<Complex> exp_mu_T (m());
      for (int k = 0; k < m(); ++k)
	exp_mu_T[k] = exp (mu[k] * T);
      
      // fill M()
      for (int i = 0; i < m(); ++i)
	for (int j = 0; j < m(); ++j)
	  {
	    double sum = 0;
	    for (int k = 0; k < m(); ++k)
	      sum += std::real (U(i,k) * exp_mu_T[k] * Uinv(k,j));
	    sum = max (sum, tp==0 ? 0.0 : min_prob);
	    dat.M(i,j) = Prob2Nats (sum);
	  }

      // fill J()
      for (int k = 0; k < m(); ++k)
	for (int l = 0; l < m(); ++l)
	  {
	    const Complex mu_diff = mu[k] - mu[l];
	    dat.J(k,l) =
	      (k == l || abs(mu_diff) <= EM_matrix_eigenvalue_tolerance)   // EQUIVALENT EIGENVALUE HACK: IH, 4/20/2005
	      ? (T * exp_mu_T[k])
	      : (exp_mu_T[k] - exp_mu_T[l]) / mu_diff;
	  }

      // log
      if (CTAGGING(1,RATE_EM_CACHE))
	{
	  CL << "M(i,j) at T=" << T << ":\n" << dat.M;
	  CL << "J(k,l) at T=" << T << ":\n" << dat.J;
	}
    }
}

void EM_matrix_base::update()
{
  update_R();
  update_log_pi();
  diagonalize();
  sanitize_eigenvalues();
  update_timepoint_cache();
}

void EM_matrix_base::sanitize_eigenvalues()
{
  // find max negative eigenvalue
  double max_neg_eval = 1.;
  for (int k = 0; k < m(); ++k)
    if (std::real(mu[k]) < 0 && (max_neg_eval > 0 || std::real(mu[k]) > max_neg_eval))
      max_neg_eval = std::real (mu[k]);

  // find max eigenvalue
  double max_eval = 0;
  int max_eval_idx = -1;
  for (int k = 0; k < m(); ++k)
    if (k == 0 ? 1 : std::real(mu[k]) > max_eval)
      max_eval = std::real (mu[max_eval_idx = k]);

  if (CTAGGING(-1,RATE_EM_CHECK_EQM))
    {
      vector<Complex> pi_minus_eqm_evec;
      for (int i = 0; i < m(); ++i)
	pi_minus_eqm_evec.push_back (pi[i] - U(i,max_eval_idx) * Uinv(max_eval_idx,i));
      CL << "pi - U(0)^2: (" << pi_minus_eqm_evec << ")\n";
    }
  
  // round down positive eigenvalues to 0 (changed from max_neg_eval/100 in response to Richard Goldstein's bug report -- IH, 10/14/2008)
  int negval = 0;
  for (int k = 0; k < m(); ++k)
    if (std::real(mu[k]) > 0) { ++negval; mu[k] = 0.; }
  //  mu[max_eval_idx] = 0.;  // set equilibrium eigenvalue to zero (commented out because we now set ALL +ve evals to 0 -- IH, 10/14/2008)
  if (negval)
    CTAG(4,RATE_EM RATE_EM_NEGVAL) << "Warning: zeroed " << negval << " positive eigenvalues; largest was " << max_eval << "\n";
  else
    CTAG(4,RATE_EM) << "Largest eigenvalue was " << max_eval << "; set to 0\n";
  
  // EQUIVALENT EIGENVALUE HACK: IH, 4/20/2005
  // if any two eigenvalues are closer than EM_matrix_eigenvalue_tolerance, then make them equal
  for (int k = 0; k < m(); ++k)
  {
    if (abs(mu[k]) < TINY)
      continue;
    for (int l = k+1; l < m(); ++l)
      if (abs (mu[k] - mu[l]) <= EM_matrix_eigenvalue_tolerance)
	{
	  CTAG(3,RATE_EM) << "Setting eigenvalue #" << l << " (" << mu[l] << ") equal to eigenvalue #" << k << " (" << mu[k] << ")\n";
	  mu[l] = mu[k];
	}
  }

  // recalculate R after messing with the eigenvalues, so everything's consistent
  for (int i = 0; i < m(); ++i)
    for (int j = 0; j < m(); ++j)
      {
	double R_sum = 0;
	for (int k = 0; k < m(); ++k)
	  R_sum += std::real (U(i,k) * mu[k] * Uinv(k,j));
	R(i+1,j+1) = R_sum;
      }
  
  // log
  if (CTAGGING(1,RATE_EM RATE_EM_EIGEN))
    CL << "Eigenvalues: " << mu << '\n';
  if (CTAGGING(1,RATE_EM_EIGEN))
    CL << "Right eigenvectors:\n" << U << "Left eigenvectors:\n" << Uinv;
  if (CTAGGING(1,RATE_EM_EQ))
    CL << "R post-sanitization:\n" << R;
}

EM_matrix_base::Update_statistics::Update_statistics (int states) : states (states)
{
  clear();
}

void EM_matrix_base::Update_statistics::clear (double pseud_init, double pseud_mutate, double pseud_wait)
{
  log_likelihood = 0.0;
  s = vector<double> (states, pseud_init);
  w = vector<double> (states, pseud_wait);
  u = array2d<double> (states, states, pseud_mutate);
  for (int i = 0; i < states; ++i)
    u(i,i) = 0.;  // zero diagonal of u
  clear_DJU();
}

void EM_matrix_base::Update_statistics::clear_DJU()
{
  DJU = array2d<Complex> (states, states, Complex(0.,0.));
  DJU_rev = array2d<Complex> (states, states, Complex(0.,0.));
}

void EM_matrix_base::Update_statistics::check_waits_transitions (const EM_matrix_base& hsm)
{
  // get alphabet
  const Alphabet& alphabet = ((EM_matrix_base&) hsm).alphabet();  // cast away const
  const int C = hsm.C;
  const int A = hsm.A;

  // zero counts that have somehow become negative
  // also zero counts where the corresponding rate is effectively zero (but have become nonzero due to precision errors)
  // (this latter problem reported by Carolin Kosiol, 10/5/2005)
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
    {
      const int i = hsm.ca (ci, ai);
      if (s[i] < 0)
      {
        CTAG(4,RATE_EM_WARNING) << "Warning: zeroing negative start count (s[" << alphabet.int2char(ai) << ci+1 << "] = " << s[i] << ")\n";
        s[i] = 0;
      }
      if (w[i] <= 0)
      {
        CTAG(4,RATE_EM_WARNING) << "Warning: setting negative wait time (w[" << alphabet.int2char(ai) << ci+1 << "] = " << w[i] << ") to "
          << HSM_MINIMUM_WAIT_TIME << "\n";
        w[i] = HSM_MINIMUM_WAIT_TIME;
	// ih, 7/8/05: removed the "short_wait" hack which was over-complicated & causing problems
	// (in particular, no reason to expect that first cycle of EM is exempt from negative wait time problems)
      }
      for (int aj = 0; aj < A; ++aj)
        if (ai != aj)
        {
          const int j = hsm.ca (ci, aj);
          if (u(i,j) < 0)
          {
            CTAG(4,RATE_EM_WARNING) << "Warning: zeroing negative transition usage (u[" << alphabet.int2char(ai) << ci+1 << "->" << alphabet.int2char(aj) << ci+1 << "] = " << u(i,j) << ")\n";
            u(i,j) = 0;
          }
        }
      for (int cj = 0; cj < C; ++cj)
        if (ci != cj)
        {
          const int j = hsm.ca (cj, ai);
          if (u(i,j) < 0)
          {
            CTAG(4,RATE_EM_WARNING) << "Warning: zeroing negative transition usage (u[" << alphabet.int2char(ai) << ci+1 << "->" << alphabet.int2char(ai) << cj+1 << "] = " << u(i,j) << ")\n";
            u(i,j) = 0;
          }
        }
    }
}

void EM_matrix_base::Update_statistics::transform (const EM_matrix_base& hsm, bool symmetrise)
{
  // derived class specific update of wait counts, transition usages
  hsm.transform_symmetrised_waits_transitions (*this, symmetrise);

  // check for negative or zero wait times, transition usages
  // also zero any transition usages if the corresponding rates are zero
  check_waits_transitions (hsm);
}

void EM_matrix_base::transform_waits_transitions (Update_statistics& stats, bool symmetrise) const
{
  CTAG(4,RATE_EM RATE_EM_EIGCOUNT) << "Transforming update estimates from eigenvector basis\n";

  const int states = m();
  if (CTAGGING(1,RATE_EM_EIGCOUNT))
  {
    CL << "Eigencounts:\n";
    for (int j = 0; j < states; ++j)
      {
	for (int i = 0; i < states; ++i)
	  CL << (symmetrise ? (.5 * (stats.DJU(i,j) + stats.DJU_rev(i,j))) : stats.DJU(i,j)) << ' ';
	CL << '\n';
      }
  }

  for (int k = 0; k < states; ++k)
    for (int l = 0; l < states; ++l)
    {
      const Complex DJU_kl = symmetrise ? ((stats.DJU(k,l) + stats.DJU_rev(k,l)) / 2.) : stats.DJU(k,l);
      for (int ci = 0; ci < C; ++ci)
        for (int ai = 0; ai < A; ++ai)
        {
          const int i = ca (ci, ai);
          const Complex Uinv_ki = Uinv(k,i);
	  // update w_i
	  const Complex wi = Uinv_ki * U(i,l) * DJU_kl;
          stats.w[i] += std::real (wi);

	  // loop over j
	  // ci=cj, ai!=aj
          for (int aj = 0; aj < A; ++aj)
            if (ai != aj)
            {
              const int j = ca (ci, aj);
	      // update u_ij
	      const Complex uij = R(i+1,j+1) * Uinv_ki * U(j,l) * DJU_kl;
              stats.u(i,j) += std::real (uij);
            }

	  // cj!=cj, ai=aj
          for (int cj = 0; cj < C; ++cj)
            if (ci != cj)
            {
              const int j = ca (cj, ai);
	      // update u_ij
	      const Complex uij = R(i+1,j+1) * Uinv_ki * U(j,l) * DJU_kl;
              stats.u(i,j) += std::real (uij);
            }
        }
    }
}

void EM_matrix_base::display_classes (ostream& out) const
{
  for (int n_align = 0; n_align < align_db->size(); ++n_align)
    {
      if (CTAGGING(4,RATE_EM RATE_EM_INIT RATE_EM_PROGRESS))
	CL << "Initialising up-down matrix for alignment #" << n_align+1 << ": '" << align_db->name[n_align] << "'\n";
      Alignment_matrix alnmat (*this, *align_db->tree_align[n_align], TRUE);
      alnmat.display_class_labels (out);  // calls fill_up_down
      out << "//\n";
    }
}

EM_matrix_base::Update_statistics EM_matrix_base::get_stats_unequilibrated (bool symmetrise, bool infer_class_labels) const
{
  Update_statistics stats (m());
  up_down (stats, 0, symmetrise, infer_class_labels);
  if (CTAGGING(2,RATE_EM_STATS_UNEQ)) CL << "Unequilibrated update statistics:\n" << stats;
  return stats;
}

Loge EM_matrix_base::log_likelihood() const
{
  Update_statistics stats (m());
  up_down (stats, 1, 0);
  return stats.log_likelihood;
}

void EM_matrix_base::Column_matrix::alloc (int _nodes, int _states, bool alloc_class_labels)
{
  nodes = _nodes;
  states = _states;

  // U, D, L tables
  U = vector<vector<Loge> > (nodes, vector<Loge> (states));
  D = vector<vector<Loge> > (nodes, vector<Loge> (states));
  L = vector<Loge> (nodes);

  // allowed states for each node; clique roots
  gapped = vector<int> (nodes, (int) 0);
  allowed = vector<vector<int> > (nodes, vector<int>());
  root = vector<int> (nodes, (int) -1);
  leaf = vector<int> (nodes, (int) 1);

  // allocate space for class labels, if mandated
  if (alloc_class_labels) class_label = sstring (nodes, '-');
}

void EM_matrix_base::Column_matrix::initialise (const Tree_alignment& tree_align, int column, const vector<int>& seq_coords, const Symbol_score_map* wildcard)
{
  CTAG(2,RATE_EM_INIT) << "Initialising position #" << column << "\n";

  const PHYLIP_tree& tree = tree_align.tree;
  const Alignment& align = tree_align.align;
  vector<const Symbol_score_map*> node_scores (nodes, (const Symbol_score_map*) 0);

  for_rooted_nodes_post (tree, b)
    {
      const int p = (*b).first;
      const int n = (*b).second;
      const int r = tree_align.node2row[n];

      if (r < 0)  // row has no data; allocate '*'s generously
	{
	  for_children (tree, n, p, c)
	    if (node_scores[*c] != 0)
	      {
		node_scores[n] = wildcard;
		break;
	      }
	}
      else if (align.path(r,column))  // allocate whatever's at this row, if it's not a gap
	if (align.prof[r])
	  node_scores[n] = &(*align.prof[r])[seq_coords[r]];
    }

  // prune long lineages of '*'s
  for_rooted_nodes_pre (tree, b)
    {
      const int p = (*b).first;
      const int n = (*b).second;
      if (p < 0 ? TRUE : node_scores[p] == 0)  // only proceed if parent is gapped
	if (node_scores[n] == wildcard && wildcard != 0)  // only proceed if node is a wildcard
	  {
	    int n_children = 0;
	    int n_ungapped_children = 0;
	    for_children (tree, n, p, c)
	      {
		++n_children;
		if (node_scores[*c] != 0)
		  ++n_ungapped_children;
	      }
	    if (n_ungapped_children < 2 && n_children > 1)
	      node_scores[n] = 0;
	  }
    }
  
  initialise (tree_align.tree, node_scores);
}

void EM_matrix_base::Column_matrix::initialise (const PHYLIP_tree& tree, const Alphabet& alphabet, const char* node_chars)
{
  if (!node_chars) THROWEXPR ("Null pointer to column string");
  if (CTAGGING(2,RATE_EM_INIT))
    {
      CL << "Initialising column '";
      for (int n = 0; n < nodes; ++n) CL << node_chars[n];
      CL << "'\n";
    }

  vector<const Symbol_score_map*> node_scores (nodes);

  for (int n = 0; n < nodes; ++n)
    node_scores[n] = Alignment::is_gapspacenull_char (node_chars[n])
      ? (Symbol_score_map*) 0
      : new Symbol_score_map (alphabet.char2score (node_chars[n]));  // uck, since Alphabet returns Symbol_score_map on the stack...

  initialise (tree, node_scores);

  for (int n = 0; n < nodes; ++n)
    if (node_scores[n])
      delete node_scores[n];
}

void EM_matrix_base::Column_matrix::init_signposts (const PHYLIP_tree& tree)
{
  for (int n = 0; n < nodes; ++n)
    {
      gapped[n] = 0;
      root[n] = -1;
      leaf[n] = 1;
    }

  clique.clear();
  for_rooted_nodes_pre (tree, b)
    {
      const Phylogeny::Node p = (*b).first;
      const Phylogeny::Node c = (*b).second;

      // get allowed symbols
      allowed[c].clear();
      for (int s = 0; s < states; ++s)
	if (U[c][s] > -InfinityLoge)
	  allowed[c].push_back (s);

      // if no symbols are allowed, this is a gap
      gapped[c] = (allowed[c].size() == 0);
      
      // if child row has a gap, it's not a leaf
      if (gapped[c])
	leaf[c] = 0;
      
      // execute the following code only if child row is ungapped
      if (!gapped[c])
	{
	  if (p >= 0) leaf[p] = 0;  // if child row is ungapped, parent is not a leaf
      
	  // if parent row is ungapped, child is not a clique root
	  if (p >= 0 ? !gapped[p] : 0)
	    root[c] = root[p];
	  else
	    {
	      root[c] = c;  // start a new clique with this as the root
	      clique.push_back (c);
	    }
	}
    }
}

void EM_matrix_base::Column_matrix::initialise (const PHYLIP_tree& tree, const vector<const Symbol_score_map*>& node_scores)
{
  // clear U, D, L, tll
  clear();

  // set the relevant U's
  for (int n = 0; n < nodes; ++n)
    if (node_scores[n])
      for_const_contents (Symbol_score_map, *node_scores[n], ss)
	U[n][ss->first] = Score2Nats (ss->second);

  // initialise signposts
  init_signposts (tree);
}

void EM_matrix_base::Column_matrix::fill_up (const EM_matrix_base& hsm, const PHYLIP_tree& tree, int column)
{
  sstring position_descriptor;
  if (column >= 0)
    position_descriptor << "column #" << column;
  fill_up (hsm, tree, position_descriptor);
}

void EM_matrix_base::Column_matrix::fill_up (const EM_matrix_base& hsm, const PHYLIP_tree& tree, const sstring& position_descriptor)
{
  const vector<const EM_matrix_base*> hsm_vec (tree.nodes(), &hsm);
  fill_up (hsm_vec, tree, position_descriptor);
}

void EM_matrix_base::Column_matrix::fill_up (const vector<const EM_matrix_base*>& hsm, const PHYLIP_tree& tree, const sstring& position_descriptor)
{
  // fill U table
  for_rooted_nodes_post (tree, b)
    {
      const Phylogeny::Node p = (*b).first;
      const Phylogeny::Node c = (*b).second;
      
      if (!gapped[c])
	{
	  // If the child node, c, is the root of its clique, then
	  // don't bother filling out node p (because there's nothing there).
	  // Instead, just update L[c] = log [sum(j) (pi[c][j] * exp(U[c][j]))].
	  // If, otoh, c is not a root, then fill out U[p][].
	  if (c == root[c])
	    {
	      for_const_contents (vector<int>, allowed[c], j)
		NatsPSumAcc (L[c], NatsPMul (hsm[c]->log_pi[*j], U[c][*j]));
	    }
	  else
	    {
	      const Timepoint_data& pc = hsm[c]->timepoint_data ((*b).length);
	      for_const_contents (vector<int>, allowed[p], i)
		{
		  Loge csum = -InfinityLoge;
		  for_const_contents (vector<int>, allowed[c], j)
		    NatsPSumAcc (csum, NatsPMul (pc.M(*i,*j), U[c][*j]));
		  NatsPMulAcc (U[p][*i], csum);
		}
	    }
	}
    }
  
  calculate_total_log_likelihood();
  
  if (CTAGGING(1,RATE_EM_UP_DOWN RATE_EM_UP_DOWN_MATRIX))
    {
      CL << "U table";
      if (position_descriptor.size())
	CL << " for " << position_descriptor;
      CL << ":\n";
      for (int n = 0; n < tree.nodes(); ++n) {
	CL << n << ": log U[" << tree.node_specifier(n) << "] =\t(";
	for (int s = 0; s < hsm[n]->m(); ++s) CL << U[n][s] << ' ';
	CL << ")\n";
      }
    }
  
  if (CTAGGING(2,RATE_EM_UP_DOWN))
    {
      CL << "Filled U table";
      if (position_descriptor.size())
	CL << " for " << position_descriptor;
      for_const_contents (vector<int>, clique, c)
	CL << " (clique " << *c << ": " << Nats2Bits(L[*c]) << " bits)";
      CL << '\n';
    }
}

void EM_matrix_base::Column_matrix::fill_down (const EM_matrix_base& hsm,
					       const PHYLIP_tree& tree,
					       Update_statistics& stats,
					       int column,
					       double weight)
{
  // prepare position descriptor
  sstring position_descriptor;
  if (column >= 0)
    position_descriptor << "column #" << column;
  fill_down (hsm, tree, stats, position_descriptor, weight);
}

void EM_matrix_base::Column_matrix::fill_down (const EM_matrix_base& hsm,
					       const PHYLIP_tree& tree,
					       Update_statistics& stats,
					       const sstring& position_descriptor,
					       double weight)
{
  const vector<const EM_matrix_base*> hsm_vec (tree.nodes(), &hsm);
  const vector<Update_statistics*> stats_vec (tree.nodes(), &stats);
  fill_down (hsm_vec, tree, stats_vec, position_descriptor, weight);
}

void EM_matrix_base::Column_matrix::fill_down (const vector<const EM_matrix_base*>& hsm_vec,
					       const PHYLIP_tree& tree,
					       const vector<Update_statistics*>& stats_vec,
					       const sstring& position_descriptor,
					       double weight)
{
  // fill the D table

  // precalculate some logging flags for speed
#define CUMULATIVE_BRANCH_LOG_TAGS RATE_EM_CUMULATIVE_BRANCH_EIGCOUNT PHYLO_EM
  const bool branch_log = CTAGGING(-2,CUMULATIVE_BRANCH_LOG_TAGS);
  if (branch_log)
    {
      CL << "Calculating per-branch DJU";
      if (position_descriptor.size())
	CL << " for " << position_descriptor;
      CL << "\n";
    }

  // do we want stats? setting stats.states=0 is a hacky way of signaling bypass of EM statistic collection during down-fill
  bool want_stats = true;
  for_const_contents (vector<Update_statistics*>, stats_vec, s)
    if ((*s)->states == 0)
      {
	want_stats = false;
	break;
      }

  // calculate update counts for start states
  if (want_stats)
    for_const_contents (vector<int>, clique, c)
      for_const_contents (vector<int>, allowed[*c], j)
      stats_vec[*c]->s[*j] += weight * Nats2Prob (NatsPMul3 (hsm_vec[*c]->log_pi[*j], U[*c][*j], -L[*c]));

  // for clique roots, fill D table with zeroes for allowed states
  for_const_contents (vector<int>, clique, r)
    for_const_contents (vector<int>, allowed[*r], i)
    D[*r][*i] = 0.;

  // fill remainder of D table, calculating relevant branch posterior probabilities & updating counts
  for_rooted_branches_pre (tree, b)
    {
      const Phylogeny::Node p  = (*b).first;
      const Phylogeny::Node n  = (*b).second;
      const Phylogeny::Node g  = tree.parent[p];

      const EM_matrix_base& hsm_n = *hsm_vec[n];

      if (!gapped[n])
	if (n != root[n])
	  {
	    const double root_likelihood = Nats2Prob (L[root[n]]);
	    if (root_likelihood == 0.)
	      {
		if (branch_log)
		  CTAG(-2,CUMULATIVE_BRANCH_LOG_TAGS) << "(clique root likelihood is too small; skipping)\n";
	      }
	    else
	      {

		// do g-p branch
		if (p == root[p])
		  for_const_contents (vector<int>, allowed[p], i)
		    D[n][*i] = hsm_vec[p]->log_pi[*i];
		else
		  {
		    const EM_matrix_base::Timepoint_data& gp = hsm_vec[p]->timepoint_data (tree.branch_length (g, p));
		    for_const_contents (vector<int>, allowed[g], i)
		      for_const_contents (vector<int>, allowed[p], j)
		      NatsPSumAcc (D[n][*j], NatsPMul (D[p][*i], gp.M(*i,*j)));
		  }

		// do p-c branches, excluding p-n branch
		for_rooted_children (tree, p, c)
		  if (!gapped[*c] && *c != n)
		    {
		      const EM_matrix_base::Timepoint_data& pc = hsm_vec[*c]->timepoint_data (tree.branch_length (p, *c));
		      for_const_contents (vector<int>, allowed[p], i)
			{
			  Loge csum = -InfinityLoge;
			  for_const_contents (vector<int>, allowed[*c], j)
			    NatsPSumAcc (csum, NatsPMul (pc.M(*i,*j), U[*c][*j]));
			  NatsPMulAcc (D[n][*i], csum);
			}
		    }

		// project D[n] & U[n] onto eigenvector basis for p-n branch
		vector<Complex> D_basis (states, Complex (0., 0.));
		for_const_contents (vector<int>, allowed[p], a)
		  {
		    const double D_na = Nats2Prob (D[n][*a]);
		    for (int k = 0; k < states; ++k)
		      D_basis[k] += D_na * hsm_n.U(*a,k);
		  }

		vector<Complex> U_basis (states, Complex (0., 0.));
		for_const_contents (vector<int>, allowed[n], b)
		  {
		    const double U_nb = Nats2Prob (U[n][*b]);
		    for (int l = 0; l < states; ++l)
		      U_basis[l] += U_nb * hsm_n.Uinv(l,*b);
		  }

		// get DJU stats for p-n branch
		const double pn_branch_length = tree.branch_length (p, n);
		const EM_matrix_base::Timepoint_data& pn = hsm_n.timepoint_data (pn_branch_length);
		if (want_stats)
		  {
		    EM_matrix_base::Update_statistics& stats = *stats_vec[n];
		    for (int k = 0; k < states; ++k)
		      for (int l = 0; l < states; ++l)
			{
			  const Complex weighted_Jkl_over_r = weight * pn.J(k,l) / root_likelihood;
			  stats.DJU(k,l) += D_basis[k] * weighted_Jkl_over_r * U_basis[l];
			  stats.DJU_rev(k,l) += U_basis[k] * weighted_Jkl_over_r * D_basis[l];
			}
		  }

		// log, if being hyper-verbose
		if (branch_log)
		  {
		    // the following perl snippet extracts the weighted counts from the logfile
		    // cat LOG | perl -ne 'if(/^ P.*= (\S+)/){$w=$1;$f=0;$p=1}elsif($p && /^\s+([0-9\-\+\.]+)/){@g=split;for$i(0..@g-1){$tot[$f]->[$i]+=$g[$i]*$w}++$f}elsif($p && $f){$p=0}elsif(/end phylo/){for$row(@tot){print"@$row\n"}@tot=();print"\n"}'
		    EM_matrix_base::Update_statistics tmp_stats = *stats_vec[n];
		    tmp_stats.transform (hsm_n, true);

		    const array2d<double> rmx = hsm_n.rate_matrix();

		    CTAG(-2,CUMULATIVE_BRANCH_LOG_TAGS) << "(begin phylo-EM block)\n";
		    CL << "Branch length " << pn_branch_length << " from node " << p << " '" << tree.node_specifier(p) << "' (parent) to node " << n << " '" << tree.node_specifier(n) << "' (child)\n";
		    CL << "U_basis=(" << U_basis << "), D_basis=(" << D_basis << "), root_likelihood=" << root_likelihood << ", cumulative eigencounts matrix:\n" << stats_vec[n]->DJU;
		    CL << "Cumulative wait times: (" << tmp_stats.w << ")\n";
		    CL << "Cumulative transition counts:\n" << tmp_stats.u;
		    CL << "Rate matrix (reconstructed from eigenvalues/eigenvectors):\n";
		    for (int ri = 0; ri < rmx.xsize(); ++ri)
		      {
			for (int rj = 0; rj < rmx.ysize(); ++rj)
			  CL << rmx(ri,rj) << ' ';
			CL << '\n';
		      }
		    CL << "Posterior branch state probabilities:\n";
		    for_const_contents (vector<int>, allowed[p], i)
		      for_const_contents (vector<int>, allowed[n], j)
		      {
			CL << " P(" << *i << "->" << *j << ") = " << Nats2Prob (NatsPMul (NatsPMul (D[n][*i], U[n][*j]), NatsPMul (pn.M(*i,*j), -L[root[n]]))) << ", counts:\n";

			const array2d<double> counts = hsm_n.clean_phylo_EM (*i, *j, pn_branch_length);
			for (int ci = 0; ci < counts.xsize(); ++ci)
			  {
			    CTAG(-2,CUMULATIVE_BRANCH_LOG_TAGS) << ' ';
			    for (int cj = 0; cj < counts.ysize(); ++cj)
			      CL << ' ' << counts(ci,cj);
			    CL << '\n';
			  }
		      }
		    CL << "(end phylo-EM block)\n\n";  // leave this line here to mark the end of the large PHYLO_EM block
		  }
	      }

	    // get class labels
	    if (!gapped[n] && estimate_class_labels())
	      {
		vector<Loge> log_class_post (hsm_n.C, -InfinityLoge);
		const Loge root_loglike = L[root[n]];
		if (n == root[n])  // clique root
		  for_const_contents (vector<int>, allowed[n], j)
		    NatsPSumAcc (log_class_post[hsm_n.state_class(*j)], NatsPMul3 (U[n][*j], hsm_n.log_pi[*j], -root_loglike));
		else  // not clique root
		  {
		    const EM_matrix_base::Timepoint_data& pc = hsm_n.timepoint_data (tree.branch_length (p, n));
		    for_const_contents (vector<int>, allowed[n], j)
		      {
			Loge nsum = -InfinityLoge;
			for_const_contents (vector<int>, allowed[p], i)
			  NatsPSumAcc (nsum, NatsPMul (pc.M(*i,*j), D[n][*i]));
			NatsPMulAcc (nsum, NatsPMul (U[n][*j], -root_loglike));
			NatsPSumAcc (log_class_post[hsm_n.state_class(*j)], nsum);
		      }
		  }
		if (CTAGGING(1,RATE_EM_CLASS_POST)) CL << "Class log-posteriors (nats): (" << log_class_post << ")\n";
		const int best_class = max_element (log_class_post.begin(), log_class_post.end()) - log_class_post.begin();
		char& best_class_char = class_label[n];
		if (log_class_post[best_class] >= -Loge2)
		  best_class_char = 'A' + best_class;
		else
		  best_class_char = 'a' + best_class;
	      }
	  }
    }

  if (CTAGGING(1,RATE_EM_UP_DOWN RATE_EM_UP_DOWN_MATRIX))
    {
      CL << "D table";
      if (position_descriptor.size())
	CL << " for " << position_descriptor;
      CL << ":\n";
      for (int n = 0; n < tree.nodes(); ++n) {
	CL << n << ": log D[" << tree.node_specifier(n) << "] =\t(";
	for (int s = 0; s < hsm_vec[n]->m(); ++s) CL << D[n][s] << ' ';
	CL << ")\n";
      }
    }

  if (CTAGGING(2,RATE_EM_UP_DOWN))
    {
      CL << "Filled D-table";
      if (position_descriptor.size())
	CL << " for " << position_descriptor;
      CL << '\n';
    }

  if (estimate_class_labels() && CTAGGING(5,RATE_EM_CLASS_LABELS))
    CL << "Class labels: " << class_label << "\n";
}

void EM_matrix_base::Column_matrix::calculate_total_log_likelihood()
{
  Loge loglike = 0;
  for_const_contents (vector<int>, clique, c)
    NatsPMulAcc (loglike, L[*c]);
  // check for infinity
  if (loglike <= -InfinityLoge)
    loglike = -InfinityLoge;
  // store
  tll = loglike;
}

EM_matrix_base::Alignment_matrix::Alignment_matrix (const EM_matrix_base& hsm, const Tree_alignment& ta, bool alloc_class_labels)
  : hsm (hsm),
    tree_align (ta),
    tree (tree_align.tree),
    align (tree_align.align),
    states (hsm.m())
{
  // create wildcard map
  for (int i = 0; i < hsm.m(); ++i) wildcard_ssm[i] = 0;
  colmat.alloc (tree.nodes(), states, alloc_class_labels);
}

void EM_matrix_base::Alignment_matrix::fill_initialised_colmat (int col)
{
  CTAG(3,RATE_EM_COLUMN) << "Doing up-down algorithm for column #" << col << "\n";
  colmat.fill_up (hsm, tree, col);
  const Loge col_loglike = colmat.total_log_likelihood();
  if (col_loglike <= -InfinityLoge)
    CTAG(5,RATE_EM) << "Warning: likelihood for column #" << col << " is zero\n";
  NatsPMulAcc (total_loglike, col_loglike);
}

void EM_matrix_base::Alignment_matrix::fill_up()
{
  // loop over columns
  total_loglike = 0.;
  vector<int> seq_coords = align.path.seq_coords_begin();
  for (int col = 0; col < align.columns(); ++col)
    {
      colmat.initialise (tree_align, col, seq_coords, &wildcard_ssm);
      fill_initialised_colmat (col);
      align.path.inc_seq_coords (seq_coords, col);
    }
  check_finite();
}

void EM_matrix_base::Alignment_matrix::check_finite() const
{
  // check for infinity
  if (total_loglike <= -InfinityLoge)
    CTAG(5,RATE_EM) << "Warning: alignment likelihood is zero (perhaps due to rounding error?) (alignment has "
		    << align.rows() << " rows and " << align.columns() << " columns)\n";
  // print log message
  if (CTAGGING(3,RATE_EM RATE_EM_PROGRESS RATE_EM_ALIGN_LOG_LIKE))
    CL << "Alignment log-likelihood = " << Nats2Bits(total_log_likelihood()) << " bits\n";
}

void EM_matrix_base::Alignment_matrix::fill_up_down (Update_statistics& stats)
{
  // loop over columns
  total_loglike = 0.;
  vector<int> seq_coords = align.path.seq_coords_begin();
  for (int col = 0; col < align.columns(); ++col)
    {
      colmat.initialise (tree_align, col, seq_coords, &wildcard_ssm);
      fill_initialised_colmat (col);
      colmat.fill_down (hsm, tree, stats, col);
      align.path.inc_seq_coords (seq_coords, col);
    }
}

EM_matrix_base::Alignment_matrix& EM_matrix_base::Alignment_matrix::operator= (const EM_matrix_base::Alignment_matrix& aln)
{
  colmat = aln.colmat;
  return *this;
}

vector<sstring> EM_matrix_base::Alignment_matrix::get_class_labels (Update_statistics& stats)
{
  total_loglike = 0.;
  vector<int> seq_coords = align.path.seq_coords_begin();
  vector<sstring> class_label (tree.nodes());
  for (int col = 0; col < align.columns(); ++col)
    {
      colmat.initialise (tree_align, col, seq_coords, &wildcard_ssm);
      fill_initialised_colmat (col);
      colmat.fill_down (hsm, tree, stats, col);
      for (int n = 0; n < tree.nodes(); ++n)
	class_label[n] << colmat.class_label[n];
      align.path.inc_seq_coords (seq_coords, col);
    }
  return class_label;
}

void EM_matrix_base::Alignment_matrix::display_class_labels (ostream& out)
{
  Update_statistics stats (hsm.m());  // dummy
  const vector<sstring> class_label = get_class_labels (stats);

  sstring::size_type max_name_len = 0;
  for (int n = 0; n < tree.nodes(); ++n)
    {
      const sstring node_name = tree.node_specifier (n);
      max_name_len = max (max_name_len, node_name.size());
    }

  for (int n = 0; n < tree.nodes(); ++n)
    {
      const sstring node_name = tree.node_specifier (n);
      out << node_name;
      for (int pad = node_name.size(); pad < (int) max_name_len + 1; ++pad) out << ' ';
      out << class_label[n];
      out << "\n";
    }
}

void EM_matrix_base::up_down (Update_statistics& stats, bool likelihood_only, bool symmetrise, bool infer_class_labels) const
{
  if (align_db == 0) THROWEXPR ("No training set!");
  const int max_n = align_db->size();
  stats.clear_DJU();
  // decide whether or not to fork
  if (max_fork < 2)
    {
      for (int n_align = 0; n_align < max_n; ++n_align)
	{
	  const sstring& align_name = align_db->name[n_align];
	  if (CTAGGING(4,RATE_EM RATE_EM_INIT RATE_EM_PROGRESS))
	    CL << "Initialising up-down matrix for alignment #" << n_align+1 << ": '" << align_name << "'\n";
	  Alignment_matrix alnmat (*this, *align_db->tree_align[n_align], infer_class_labels);
	  if (likelihood_only)
	    alnmat.fill_up();
	  else
	    alnmat.fill_up_down (stats);
	  stats.log_likelihood += alnmat.total_log_likelihood();
	}
    }
  else
    {
      // open some pipes
      EM_matrix_base& noconst = *((EM_matrix_base*)this);  // cast away const
      noconst.open_pipes();
      // fill the profile, forkily
      for (int n_start = 0; n_start < max_n; n_start += max_fork)
	{
	  const int max_offset = min (max_fork, max_n - n_start);
	  for (int offset = 0; offset < max_offset; ++offset)
	    if (noconst.fork_child(offset))  // child
	      {
		ostream_cfile pipe_out (write_fd (offset));
		const sstring& align_name = align_db->name[n_start + offset];
		if (CTAGGING(4,RATE_EM RATE_EM_INIT RATE_EM_PROGRESS))
		  CL << "Initialising up-down matrix for alignment #" << n_start+offset+1 << ": '" << align_name << "'\n";
		Alignment_matrix alnmat (*this, *align_db->tree_align[n_start + offset], infer_class_labels);
		alnmat.fill_up();
		if (likelihood_only)
		  alnmat.fill_up();
		else
		  alnmat.fill_up_down (stats);
		stats.log_likelihood = alnmat.total_log_likelihood();
		stats.send (pipe_out.out);
		pipe_out.close();
		exit(0);
	      }
	  
	  // parent
	  for (int offset = 0; offset < max_offset; ++offset)
	    {
	      istream_cfile pipe_in (read_fd (offset));
	      stats.receive (pipe_in.in);
	      wait_for_child(offset);  // wait for each child to exit
	    }
	}
    }
  if (!likelihood_only)
    stats.transform (*this, symmetrise);
  // check for infinity
  if (stats.log_likelihood <= -InfinityLoge)
    CLOGERR << "Warning: likelihood for alignment database is zero\n";
}

ostream& operator<< (ostream& o, const EM_matrix_base::Update_statistics& stats)
{
  o << "Start count = (" << stats.s << ")\n";
  o << "Wait time = (" << stats.w << ")\n";
  o << "Transition usage:\n" << stats.u.transpose();
  return o;
}

void EM_matrix_base::Update_statistics::send (ostream& out) const
{
  out << log_likelihood << "\n";
  out << s << "\n";
  out << DJU << "\n";
  out << DJU_rev << "\n";
}

void EM_matrix_base::Update_statistics::receive (istream& in)
{
  double real_inval;
  Complex cplx_inval;
  in >> real_inval;
  log_likelihood += real_inval;
  const int n = s.size();
  for (int i = 0; i < n; ++i) { in >> real_inval; s[i] += real_inval; }
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      { in >> cplx_inval; DJU(i,j) += cplx_inval; }
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      { in >> cplx_inval; DJU_rev(i,j) += cplx_inval; }
}

void EM_matrix_base::randomise (double prior_dev, double intra_min, double intra_dev, double inter_min, double inter_dev)
{
  const double prior_min = 1;
  double beta = 0;
  for (int ci = 0; ci < C; ++ci)
    for (int ai = 0; ai < A; ++ai)
      {
	const int i = ca(ci,ai);
	pi[i] = prior_min + Rnd::prob() * prior_dev;
	beta += pi[i];

	double& X_ii = X[ci](ai,ai);
	for (int aj = 0; aj < A; ++aj)
	  if (aj != ai)
	    {
	      X[ci](ai,aj) = aj < ai ? X[ci](aj,ai) : (intra_min + Rnd::prob() * intra_dev);
	      X_ii -= X[ci](ai,aj);
	    }

	double& Y_ii = Y[ai](ci,ci);
	for (int cj = 0; cj < C; ++cj)
	  if (cj != ci)
	    {
	      Y[ai](ci,cj) = cj < ci ? Y[ai](cj,ci) : (inter_min + Rnd::prob() * inter_dev);
	      Y_ii -= Y[ai](ci,cj);
	    }
      }
  for (int i = 0; i < m(); ++i)
    pi[i] /= beta;
  update();
}

void EM_matrix_base::write (ostream& out) const
{
  // write C and A
  out << C << " " << A << "\n";
  // write alphabet
  for (int i = 0; i < A; ++i)
    out << hidden_alphabet.int2char(i) << ' ';
  out << '\n';
  // write initial residue distribution
  for (int i = 0; i < m(); ++i)
    out << " " << pi[i];
  out << "\n";
  // write intra-class substitution rates
  for (int c = 0; c < C; ++c)
    for (int ai = 0; ai < A; ++ai)
      {
	for (int aj = 0; aj < A; ++aj)
	  out << X[c](ai,aj) << " ";
	out << "\n";
      }
  // write inter-class substitution rates
  for (int a = 0; a < A; ++a)
    for (int ci = 0; ci < C; ++ci)
      {
	for (int cj = 0; cj < C; ++cj)
	  out << Y[a](ci,cj) << " ";
	out << "\n";
      }
}

void EM_matrix_base::read (istream& in)
{
  // read C and A
  int new_C, new_A;
  in >> new_C >> new_A;
  init_matrix (new_C, new_A, align_db, timepoint_res);
  // read alphabet
  sstring tok, tokstr;
  for (int i = 0; i < new_A; ++i)
    {
      in >> tok;
      if (tok.size() != 1)
	THROWEXPR ("Can't currently handle multi-character tokens in alphabets");
      tokstr.push_back (tok[0]);
    }
  Alphabet base_alph (Default_alphabet_name, A);
  base_alph.init_chars (tokstr.c_str());
  init_alphabet (base_alph);
  // read initial residue distribution
  for (int i = 0; i < m(); ++i)
    in >> pi[i];
  // read intra-class substitution rates & set/clear update flags
  for (int c = 0; c < C; ++c)
    for (int ai = 0; ai < A; ++ai)
      for (int aj = 0; aj < A; ++aj)
	{
	  in >> X[c](ai,aj);
	  X_update_flag[c](ai,aj) = X[c](ai,aj) > HSM_MINIMUM_NONZERO_RATE;
	}
  // read inter-class substitution rates & set/clear update flags
  for (int a = 0; a < A; ++a)
    for (int ci = 0; ci < C; ++ci)
      for (int cj = 0; cj < C; ++cj)
	{
	  in >> Y[a](ci,cj);
	  Y_update_flag[a](ci,cj) = Y[a](ci,cj) > HSM_MINIMUM_NONZERO_RATE;
	}
  // update eigenvectors etc
  update();
}

void EM_matrix_base::guess_alphabet()
{
  switch (A)
    {
    case 4: init_alphabet (DNA_alphabet); break;
    case 20: init_alphabet (Protein_alphabet); break;
    default: THROWEXPR ("Don't know any alphabets with " << A << " symbols\n");
    }
}

Loge EM_matrix_base::iterate_quick_EM (bool intra, bool inter, int forgive, bool infer_class_labels)
{
  Loge best_log_likelihood = 0;
  EM_matrix_params best_params;
  int dec = 0;
  for (int iter = 0; ; ++iter)
  {
    if (iter >= em_max_iter)
	{
	  CTAG(7,RATE_EM RATE_EM_PROGRESS) << "EM hit " << em_max_iter << " iterations; stopping\n";
	  break;
	}
    const EM_matrix_params old_params = *this;
    Update_statistics stats = single_quick_EM (intra, inter, infer_class_labels);
    const Loge prev_best = best_log_likelihood;
    if (iter == 0 || stats.log_likelihood > best_log_likelihood)
	{
	  best_log_likelihood = stats.log_likelihood;
	  best_params = old_params;
	}
    CTAG(6,RATE_EM RATE_EM_PROGRESS) << "EM iteration #" << iter+1 << ": log-likelihood = " << Nats2Bits(stats.log_likelihood) << " bits\n";
    if (iter > 0)
	{
	  const double inc = (stats.log_likelihood - prev_best) / (abs(prev_best) < TINY ? 1. : -prev_best);  // IH, 4/20/2005
	  CTAG(3,RATE_EM RATE_EM_PROGRESS) << "(previous best = " << Nats2Bits(prev_best) << " bits; fractional improvement = " << inc << ")\n";
	  if (inc < em_min_inc)
	  {
	    if (stats.log_likelihood < prev_best)
	      CTAG(7,RATE_EM RATE_EM_PROGRESS) << "Warning: log-likelihood dropped from " << Nats2Bits(prev_best) << " to " << Nats2Bits(stats.log_likelihood) << " bits during EM\n";
	    if (++dec > forgive)
		{
		  CTAG(7,RATE_EM RATE_EM_PROGRESS) << "Failed EM improvement threshold for the " << dec << "th time; stopping\n";
		  break;
		}
	  }
	  else
	    dec = 0;
	}
  }
  // since we only update the best log-likelihood the round *after* the new params are calculated, we need to check for an increase one last time
  CTAG(6,RATE_EM RATE_EM_PROGRESS) << "Checking post-iteration log-likelihood\n";
  const Loge final_log_likelihood = log_likelihood();
  CTAG(6,RATE_EM RATE_EM_PROGRESS) << "Post-iteration log-likelihood = " << Nats2Bits(final_log_likelihood) << " bits\n";
  if (em_max_iter != 0 && final_log_likelihood < best_log_likelihood)
  {
    CTAG(6,RATE_EM RATE_EM_PROGRESS) << "Restoring previous best paramaters\n";
    ((EM_matrix_params&) *this) = best_params;
    update();
  }
  else
    best_log_likelihood = final_log_likelihood;
  // and return
  return best_log_likelihood;
}


bool EM_matrix_base::hsm_help (Opts_list* ol)
{
  cout << "\n";
  cout << "File format for hidden substitution matrix\n";
  cout << "==========================================\n";
  cout << "Line 0 has two integer parameters \"C A\"\n";
  cout << " where C is the number of hidden classes\n";
  cout << "       A is the number of alphabet symbols (e.g. 20 for proteins)\n";
  cout << "Line 1 has A whitespace-separated alphabet tokens\n";
  cout << "Line 2 has C*A floating-point parameters representing probabilities at equilibrium\n";
  cout << "Lines    3..A+2     are the intra-class substitution matrix for class #1\n";
  cout << "                     (A*A reversible rate matrix, R_ij = rate from i to j)\n";
  cout << "Lines  A+3..2A+2    are the intra-class matrix for class #2\n";
  cout << "  ...etc...\n";
  cout << "Lines    X..X+C-1  (where X=CA+3) are the inter-class matrix for residue #1\n";
  cout << "Lines  X+C..X+2C-1   are the inter-class matrix for residue #2\n";
  cout << "  ...etc...\n";
  cout << "\n";
  cout << "All numeric parameters are delimited by whitespace.\n";
  cout << "\n";
  exit(0);
  return 0;
}

void EM_matrix_base::normalise_substitution_rate()
{
  const double r = expected_substitution_rate();

  if (CTAGGING(2,RATE_EM_NORM))
    {
      CL << "Matrix before normalisation:\n";
      write (CL);
    }
  CTAG(5,RATE_EM_NORM RATE_EM_NORM_FACTOR) << "Expected substitution rate/site = " << r << "\n";

  for (int c = 0; c < C; ++c)
    for (int ai = 0; ai < A; ++ai)
      for (int aj = 0; aj < A; ++aj)
	X[c] (ai, aj) /= r;

  for (int a = 0; a < A; ++a)
    for (int ci = 0; ci < C; ++ci)
      for (int cj = 0; cj < C; ++cj)
	Y[a] (ci, cj) /= r;

  update();

  if (CTAGGING(2,RATE_EM_NORM))
    {
      CL << "Matrix after normalisation:\n";
      write (CL);
    }
}

array2d<double> EM_matrix_base::clean_phylo_EM (int src, int dest, double T) const
{
  Phylo_EM<Complex> phylo_EM (mu, U, Uinv);
  array2d<Complex> complex_counts = phylo_EM.get_counts (src, dest, T);
 
  array2d<double> real_counts (m(), m());
  double largest_imag = 0.;
  int n_nonzero_imag = 0;
  for (int i = 0; i < m(); ++i)
    for (int j = 0; j < m(); ++j)
      {
 	real_counts(i,j) = std::real (complex_counts(i,j));
	const double im = abs(std::imag (complex_counts(i,j)));
 	if (im > TINY)
	  {
	    ++n_nonzero_imag;
	    if (im > largest_imag)
	      largest_imag = im;
	  }
      }
  if (n_nonzero_imag)
    CLOGERR << "Warning: complex_counts had " << n_nonzero_imag << " entries with nonzero imaginary component (largest was " << largest_imag << ")\n";
 
  return real_counts;
}

