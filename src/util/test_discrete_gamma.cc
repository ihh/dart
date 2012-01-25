#include "util/discrete_gamma.h"
#include "util/vector_output.h"
#include "util/math_fn.h"

void check_4dp (double a, double b) {
  if (abs(a-b) >= .00005) {
    cout << "not ok\n";
    exit(1);
  }
}

void check_int (int a, int b) {
  if (a != b) {
    cout << "not ok\n";
    exit(1);
  }
}

void check_4 (vector<double>& x, double y0, double y1, double y2, double y3) {
  check_int (x.size(), 4);
  cerr << "Figure 1 of Yang Z, JME, 1994: (" << y0 << ' ' << y1 << ' ' << y2 << ' ' << y3 << ")\n";
  cerr << "Actual values: (" << x << ")\n";

  check_4dp (x[0], y0);
  check_4dp (x[1], y1);
  check_4dp (x[2], y2);
  check_4dp (x[3], y3);
}

int main (int argc, char** argv) {
  cerr << "Medians:\n";
  Discrete_gamma dg_median (.5, .5, 4, true);
  check_4 (dg_median, 0.0291, 0.2807, 0.9248, 2.7654);

  cerr << "Means:\n";
  Discrete_gamma dg_mean (.5, .5, 4, false);
  check_4 (dg_mean, 0.0334, 0.2519, 0.8203, 2.8944);

  cout << "ok\n";

  return 1;
}
