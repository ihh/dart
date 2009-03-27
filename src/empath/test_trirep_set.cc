#include "empath/trirep.h"

int main (int argc, char** argv)
{
  try
    {
      Trirep_set trirep_set (3, 3);
      if (trirep_set.reps.size() != 30) THROWEXPR ("Should be 30 reps");
      cerr << trirep_set;
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cout << "ok\n";
}

