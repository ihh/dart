#include "hsm/em_matrix.h"
#include "util/vector_output.h"

#define CMP_TOL .001
#define grad_xdelta .000000001
#define grad_tol .1

void test_gradient (Grad_function& gf, const vector<double>& x)
{
  double f;
  vector<double> grad_f (gf.dim());
  gf.eval (x, f, grad_f);
  vector<double> xnew (x);
  vector<double> dfdx (gf.dim());
  for (int i = 0; i < gf.dim(); ++i)
    {
      xnew[i] += grad_xdelta;
      double fnew;
      vector<double> tmp (gf.dim());
      gf.eval (xnew, fnew, tmp);
      dfdx[i] = (fnew - f) / grad_xdelta;
      xnew[i] = x[i];
    }
  for (int i = 0; i < gf.dim(); ++i)
    {
      const double abserr = abs (dfdx[i] - grad_f[i]);
      if (abserr >= grad_tol) THROWEXPR ("f: " << f << "\nx: (" << x << ")\nReturned gradient: (" << grad_f << ")\nActual gradient:  (" << dfdx << ")\n");
    }
}


bool cmp_stats (double a, double b)
{ return abs(a-b) < CMP_TOL; }

void test_str_eq (const char* a, const char* b)
{
  const int alen = strlen (a);
  const int blen = strlen (b);
  for (int i = 0; i < alen; ++i)
    if (i >= blen || b[i] != a[i])
      {
	sstring asub;
	sstring bsub;
	const int lookback = 10;
	if (i < lookback) asub = bsub = "...";
	for (int j = max (i - 10, -1); j <= i; ++j)
	  {
	    asub << a[j];
	    bsub << (j < blen ? b[j] : '$');
	  }
	THROWEXPR ("String mismatch at position " << i << ": \"" << asub << "\" != \"" << bsub << "\"\nFull strings:\na = \"" << a << "\"\nb = \"" << b << "\"\n");
      }
}

int main (int argc, char** argv)
{
  try
    {
      Opts_list opts (argc, argv);
      Log_stream::add_opts (opts);
      int max_fork;
      opts.add ("-fork fork f", max_fork = 1, "\tnumber of processes to fork [BUGGY under gcc 3.3.1; BEWARE]");

      if (!opts.parse()) THROWEXPR ("bad command-line options");

      cerr << "(calculating update stats for ((A:B):(C:D)) tree, one-column alignment)\n";

      Tree_alignment tree_align;
      Sequence_database seq_db;

      ifstream tree_stream ("testemmatrix.phylip");
      ifstream align_stream ("testemmatrix.mul");

      if (!tree_stream || !align_stream) THROWEXPR ("Missing test files");

      tree_align.read_PHYLIP (tree_stream);
      tree_align.read_MUL (align_stream, seq_db);
      
      tree_align.build_maps_from_names();
      seq_db.seqs2scores (DNA_alphabet);
      
      Tree_alignment_database tree_align_db (seq_db);
      tree_align_db.tree_align_list.push_back (tree_align);
      tree_align_db.tree_align.push_back (&(*tree_align_db.tree_align_list.rbegin()));
      tree_align_db.name.push_back (sstring ("test"));
      
      Simple_alphabet binary_alphabet ("Binary", "AB");

      EM_matrix em_mat (1, 2, max_fork, &tree_align_db, .01);
      em_mat.init_alphabet (binary_alphabet);

      em_mat.X[0](0,0) = -3;
      em_mat.X[0](0,1) = 3;
      em_mat.X[0](1,0) = 4;
      em_mat.X[0](1,1) = -4;

      em_mat.pi[0] = 4./7.;
      em_mat.pi[1] = 3./7.;
      
      em_mat.update();

      EM_matrix::Update_statistics stats = em_mat.get_stats_unequilibrated(0);

      EM_matrix::Update_statistics targ (2);  // target stats
      // the following numbers courtesy of Maple
      targ.log_likelihood = Prob2Nats (.04444543807);
      targ.s[0] = .826866;
      targ.s[1] = .173134;
      targ.w[0] = .162978;
      targ.w[1] = .047022;
      targ.u(0,1) = .893225;
      targ.u(1,0) = .582443;

      /* for reference, U & D tables should look like this (command line '-log RATE_EM_UP_DOWN'):
# U table:
# log L[root] =	(-2.74443 -4.02002)
# log U[A::B] =	(-0.0870199 -5.84858)
# log U[A] =	(0 -9.87654e+99)
# log U[B] =	(0 -9.87654e+99)
# log U[C::D] =	(-2.62167 -2.08387)
# log U[C] =	(-9.87654e+99 0)
# log U[D] =	(0 -9.87654e+99)
# D table:
# log D[root] =	(-9.87654e+99 -9.87654e+99)
# log D[A::B] =	(-3.08171 -3.01621)
# log D[A] =	(-3.0875 -5.66286)
# log D[B] =	(-3.05928 -6.32162)
# log D[C::D] =	(-0.781954 -2.69841)
# log D[C] =	(-1.0189 -4.07886)
# log D[D] =	(-3.41946 -2.2241)
       */

      // compare target & actual update stats
      bool fail = 0;
      for (int i = 0; i < 2; ++i)
	{
	  if (!cmp_stats (stats.s[i], targ.s[i])) fail = 1;
	  if (!cmp_stats (stats.w[i], targ.w[i])) fail = 1;
	  for (int j = 0; j < 2; ++j)
	    if (!cmp_stats (stats.u(i,j), targ.u(i,j))) fail = 1;
	}
      if (!cmp_stats (stats.log_likelihood, targ.log_likelihood)) fail = 1;
      if (fail) THROWEXPR ("Bad update stats.\nTarget:\n" << targ << "Actual:\n" << stats);

      cerr << "(testing grad(L) indices)\n";
      int C = 3;
      int A = 5;
      EM_matrix im (C, A, max_fork, &tree_align_db, .01);
      int idx = 0;
      for (int c = 0; c < C; ++c)
	for (int a = 0; a < A; ++a)
	  if (im.pi_idx(c,a) != idx++) fail = 1;
      for (int c = 0; c < C; ++c)
	for (int ai = 0; ai < A; ++ai)
	  for (int aj = ai + 1; aj < A; ++aj)
	    if (im.X_idx(c,ai,aj) != idx++) fail = 1;
      for (int a = 0; a < A; ++a)
	for (int ci = 0; ci < C; ++ci)
	  for (int cj = ci + 1; cj < C; ++cj)
	    if (im.Y_idx(ci,cj,a) != idx++) fail = 1;
      if (im.beta_idx() != idx++) fail = 1;
      for (int c = 0; c < C; ++c)
	for (int a = 0; a < A; ++a)
	  if (im.kappa_idx(c,a) != idx++) fail = 1;
      if (im.n_params() != idx) fail = 1;
      if (fail) THROWEXPR ("Indexing error");

      cerr << "(checking Lagrangian)\n";

      C = 2;
      A = 2;
      EM_matrix lm (C, A, max_fork, &tree_align_db, .01);
      lm.init_alphabet (binary_alphabet);

      EM_matrix::Update_statistics lstat (lm.m());
      for (int i = 0; i < 4; ++i)
	{
	  lstat.s[i] = .1 * (i + 1);
	  lstat.w[i] = 11 * (i + 1);
	  for (int j = 0; j < 4; ++j)
	    lstat.u(i,j) = .1 * (j + 1) + .01 * (i + 1);
	}
      
      Diff lag = lm.lagrangian (lstat);

      sstring lagstr = lag.to_string();
      test_str_eq (lagstr.c_str(), "0.1*log(PiA0)+11*(-1*XAB0-1*YA01)+0.21*log(XAB0)+0.31*log(YA01)+0.2*log(PiB0)+22*(-1*XAB0*PiA0*PiB0^(-1)-1*YB01)+0.12*log(XAB0*PiA0*PiB0^(-1))+0.42*log(YB01)+0.3*log(PiA1)+33*(-1*XAB1-1*YA01*PiA0*PiA1^(-1))+0.43*log(XAB1)+0.13*log(YA01*PiA0*PiA1^(-1))+0.4*log(PiB1)+44*(-1*XAB1*PiA1*PiB1^(-1)-1*YB01*PiB0*PiB1^(-1))+0.34*log(XAB1*PiA1*PiB1^(-1))+0.24*log(YB01*PiB0*PiB1^(-1))-1*log(PiA0+PiB0+PiA1+PiB1)+kappaA0*(PiA0*(-1*XAB0-1*YA01)+PiB0*XAB0*PiA0*PiB0^(-1)+PiA1*YA01*PiA0*PiA1^(-1))+kappaB0*(PiA0*XAB0+PiB0*(-1*XAB0*PiA0*PiB0^(-1)-1*YB01)+PiB1*YB01*PiB0*PiB1^(-1))+kappaA1*(PiA0*YA01+PiA1*(-1*XAB1-1*YA01*PiA0*PiA1^(-1))+PiB1*XAB1*PiA1*PiB1^(-1))+kappaB1*(PiB0*YB01+PiA1*XAB1+PiB1*(-1*XAB1*PiA1*PiB1^(-1)-1*YB01*PiB0*PiB1^(-1)))+beta*(PiA0+PiB0+PiA1+PiB1-1)");

      cerr << "(checking Laplacian gradient at initial seed position)\n";
      
      const vector<double> lseed = lm.unconstrained_optimum (lstat);
      // const double lag_val = lag.eval (lseed);
      Laplacian lap (lag);
      test_gradient (lap, lseed);

      cerr << "(HSM training code good. Let's kick some substitution model ass)\n";
    }
  catch (const Dart_exception& e)
    {
      cerr << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << "ok\n";
  return 0;
}
