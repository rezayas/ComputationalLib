#include "linreg.hpp"
#include <random>
#include <ctime>
#include <iostream>

using namespace std;
int main() {
  Eigen::MatrixX3d xs(1000, 3); // Why do I still need to pass the 3?
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(-1,1);
  std::mt19937_64 rands(std::time(NULL));
  Eigen::Vector3d xp(1, 2, 3);
  Eigen::VectorXd ys(1000);
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = 1;
    xs(i,1) = norm(rands);
    xs(i,2) = norm(rands);
    ys(i) = uni(rands);
  }
  ys += xs * xp;
  linreg::LinearRegression lr(3);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i));
  cout << lr.getCoefficients() << endl << endl;
  cout << linreg::LinearRegression::runRegression(xs, ys) << endl << endl;

  linreg::LinearRegression lrb(3, .5);
  for(int i = 0; i < 2; i++)
    lrb.updateCoefficients(xs.row(i).transpose(), ys(i));
  cout << lrb.getCoefficients() << endl << endl;
  cout << linreg::LinearRegression::runRegression(1, .5, xs.topRows(2), ys.head(2)) << endl;
  // Test 1: regular LS 
  // Test 2: updating algorithm with omega = 0, lambda = 1 (we should get the same coefficient as Test 1)

  // Test 3: LS + omega > 0 on a X matrix that is not full rank
  // Test 4: updating algorithm with omega > 0 on a X matrix that is not full rank (same results as Test 3)
}
