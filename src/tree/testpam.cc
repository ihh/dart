#include "tree/pam.h"

int main (int argc, char** argv)
{
  // doesn't really do any numerical checks, but prints OK if it makes it through without dying
  try
    {
      //      if (argc != 3) THROWEXPR ("Usage: " << argv[0] << " <time> <score to bits ratio>");

      //      double t = atof (argv[1]);
      //      double b = atof (argv[2]);
      
      double t = 1.20;  // PAM120
      double b = 3;  // score of 3 == 1 bit

      PAM_factory pam;
      array2d<double> prob_mx = pam.create_conditional_substitution_matrix (t);
      cerr << "Probability matrix:\n";
      for (int i = 0; i < 20; ++i)
	for (int j = 0; j < 20; ++j)
	  cerr << prob_mx(i,j) << (j == 19 ? "\n" : " ");
      cerr << "\n";

      array2d<double> odds_mx = pam.create_substitution_odds_ratio_matrix (t);
      cerr << "Score matrix:\n";
      for (int i = 0; i < 20; ++i)
	for (int j = 0; j < 20; ++j)
	  cerr << (int) (Prob2Bits(odds_mx(i,j)) * b) << (j == 19 ? "\n" : " ");
      cerr << "\n";
    }
  catch (const Dart_exception& e)
    {
      cout << "not ok\n";
      CLOGERR << e.what();
      exit(1);
    }
  cout << "ok\n";
}

