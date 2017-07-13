#include "calc.hpp"
#include <random>
#include <ctime>

using namespace std;
int main() {
  mt19937_64 foo(time(NULL));
  normal_distribution<> bar;
  uniform_real_distribution<> baz(-1, 1);
  arma::vec xs(3);
  xs[0] = 1;
  arma::rowvec xp = {1, 2, 3};
  setup(3);
  for(int i = 0; i < 1000; i++) {
    xs[1] = bar(foo);
    xs[2] = bar(foo);
    updateTheta(1, xs, as_scalar(xp * xs) + baz(foo));
  }
  cout << theta;
}
