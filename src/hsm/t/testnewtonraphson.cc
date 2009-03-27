#include "hsm/newtonraphson.h"
#include "util/vector_output.h"
#include "util/math_fn.h"
#include "util/logfile.h"

#define CMP_TOL .00001

bool cmp_val (double a, double b)
{ return abs(a-b) < CMP_TOL; }

int main (int argc, char** argv)
{
  try
    {
      Opts_list opts (argc, argv);
      Log_stream::add_opts (opts);
      if (!opts.parse()) THROWEXPR ("bad command-line options");
      
      /* Maple:
	 
f = -(a-2*b)^2+exp(-(c-d)^2)-(e-4)^4+log(a)/(a+1)^2+(c+.1*b-2)^2;

L = diff(f,a)^2+diff(f,b)^2+diff(f,c)^2+diff(f,d)^2+diff(f,e)^2;

(a,b,c,d,e) = 2.093495237, 1.046747619, 1.895325238, 1.895325238, 4;

f = 1.077205560;

L = .233600e-16;
			    
*/
      Diff_param a (0);
      Diff_param b (1);
      Diff_param c (2);
      Diff_param d (3);
      Diff_param e (4);

      Diff_const one (1);
      Diff_const two (2);
      Diff_const minus_one (-1);
      Diff_const minus_two (-2);
      Diff_const minus_four (-4);
      Diff_const point_one (.1);

      Diff two_b = two * b;
      Diff a_minus_two_b = a - two_b;
      Diff a_minus_two_b_squared = a_minus_two_b.pow(2);
      Diff minus_a_minus_two_b_squared = minus_one * a_minus_two_b_squared;
      
      Diff c_minus_d = c - d;
      Diff c_minus_d_squared = c_minus_d.pow(2);
      Diff minus_c_minus_d_squared = minus_one * c_minus_d_squared;
      Diff exp_minus_c_minus_d_squared = minus_c_minus_d_squared.exp();
      
      Diff e_minus_four = e + minus_four;
      Diff e_minus_four_pow4 = e_minus_four.pow(4);
      Diff minus_e_minus_four_pow4 = minus_one * e_minus_four_pow4;
      
      Diff a_plus_one = a + one;
      Diff a_plus_one_squared = a_plus_one.pow(2);
      Diff log_a = a.log();
      Diff log_a_over_a_plus_one_squared = log_a / a_plus_one_squared;

      Diff point_one_b = point_one * b;
      Diff point_one_b_minus_two = point_one_b + minus_two;
      Diff c_plus_point_one_b_minus_two = c + point_one_b_minus_two;
      Diff c_plus_point_one_b_minus_two_squared = c_plus_point_one_b_minus_two.pow(2);
      
      vector<Diff> D_series;
      D_series.push_back (minus_a_minus_two_b_squared);
      D_series.push_back (exp_minus_c_minus_d_squared);
      D_series.push_back (minus_e_minus_four_pow4);
      D_series.push_back (log_a_over_a_plus_one_squared);
      D_series.push_back (c_plus_point_one_b_minus_two_squared);
      
      vector<double> abcde_seed (5, 1.0);
      //      cerr << "Series: ";
      //      for (int i = 0; i < D_series.size(); ++i)
      //	cerr << " term[" << i << "]=" << D_series[i]->eval(abcde_seed);
      //      cerr << "\n";

      /* Maple: */
      /* D(1,1,1,1,1)=-80.19, dD/d{a,b,c,d,e} = 9/4, -4.18, -1.8, 0, 108 */

      vector<sstring> abcde_text (5);
      for (int i = 0; i < 5; ++i) abcde_text[i] << (char) (i + 'a');

      Diff D = Diff::sum (D_series);
      D.param_name = &abcde_text;

      cerr << "(created function \"" << D.to_string() << "\")\n";

      cerr << "(trying to maximize D=-(a-2*b)^2+exp(-(c-d)^2)-(e-4)^4+log(a)/(a+1)^2+(c+.1*b-2)^2)\n";
      cerr << "(by finding zeroes of L=grad(D)^2)\n";
      cerr << "(testing seed a=b=c=d=e=1)\n";

      const double maple_D = -80.19;
      const double maple_D_deriv[5] = { 2.25, -4.18, -1.8, 0, 108 };

      bool fail = 0;
      const double D_val = D.eval(abcde_seed);
      cerr << "(D=" << D_val;
      if (!cmp_val (D_val, maple_D)) fail = 1;
      for (int i = 0; i < 5; ++i)
	{
	  Diff D_deriv = D.derivative (i);
	  const double D_deriv_val = D_deriv.eval(abcde_seed);
	  cerr << " dD/d" << (char) ('a' + i) << "=" << D_deriv_val;
	  if (!cmp_val (D_deriv_val, maple_D_deriv[i])) fail = 1;
	}
      cerr << ")\n";

      if (fail) THROWEXPR ("Bad D");

      /* Maple: */
      /* L=11689.77490, dL/d{a,b,c,d,e} = {-45.81500000, 83.9928, -1.672, -7.2, -23328 */

      Laplacian laplace (D);
      const double maple_L = 11689.77490;
      const double maple_L_deriv[5] = { -45.81500000, 83.9928, -1.672, -7.2, -23328 };

      double L;
      vector<double> L_deriv (5);
      laplace.eval (abcde_seed, L, L_deriv);
      cerr << "(L=" << L;
      if (!cmp_val (L, maple_L)) fail = 1;
      for (int i = 0; i < 5; ++i)
	{
	  cerr << " dL/d" << (char) ('a' + i) << "=" << L_deriv[i];
	  if (!cmp_val (L_deriv[i], maple_L_deriv[i])) fail = 1;
	}
      cerr << ")\n";

      if (fail) THROWEXPR ("Bad L");

      cerr << "(seed values OK; trying Newton-Raphson. Expect about 2625 iterations...)\n";

      const double nr_tol = .0000001;
      const int nr_max_iter = 10000;
      map<int,double> no_positives;
      vector<double> abcde = Newton_Raphson::iterate (laplace, abcde_seed, no_positives, nr_tol, nr_max_iter);
      laplace.eval (abcde, L, L_deriv);
      cerr << "(L=" << L;
      if (L > nr_tol) fail = 1;
      const double maple_abcde[5] = { 2.093495237, 1.046747619, 1.895325238, 1.895325238, 4 };
      for (int i = 0; i < 5; ++i)
	{
	  cerr << " " << (char) ('a' + i) << "=" << abcde[i];
	  if (abs (abcde[i] - maple_abcde[i]) > .05) fail = 1;
	}
      cerr << ")\n";
      if (fail) THROWEXPR ("Newton-Raphson failed");

      cerr << "(you are SOOO optimal)\n";
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << "ok\n";
}
