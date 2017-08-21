#include "linreg.hpp"
#include <random>
#include <ctime>
#include <iostream>
// Issue 13
using namespace std;
int main() {
  
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(-1,1);
  std::mt19937_64 rands(std::time(NULL));
  Eigen::VectorXd xp(4); // fitting y = 1 + 2x1 - 3x3
  xp << 1, 2, 0, -3;
  Eigen::MatrixXd xs(1000, 4); // X matrix
  Eigen::VectorXd ys(1000);	// y vector
  // populate X matrix and y vector
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = 1;
    xs(i,1) = norm(rands);
    xs(i,2) = norm(rands);
    xs(i,3) = norm(rands);
    ys(i) = uni(rands);
  }
  ys += xs * xp;
  auto co = complib::LinearRegression::runRegression(xs, ys);
  // Test 1: regular LS 
  cout << complib::LinearRegression::runRegression(xs, ys) << endl << endl;
  // Test 5: updating algorithm with omega = 0, lambda = 1 (we should get the same coefficient as Test 1)
  complib::LinearRegression lr(4);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i));
  cout << lr.getCoefficients() << endl << endl;

  // Test 2: w = .001, l = 1
  cout << complib::LinearRegression::runRegression(1, .001, xs, ys)
       << endl << endl;
  // Test 6: same, online
  lr.reset(4, .001);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i));
  cout << lr.getCoefficients() << endl << endl;

  // Test 3: w = .001, l = .95
  cout << complib::LinearRegression::runRegression(.95, .001, xs, ys)
       << endl << endl;
  // Test 7: same, online
  lr.reset(4, .001);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i), .95);
  cout << lr.getCoefficients() << endl << endl;

  // Test 4: w = 0, l = .95
  cout << complib::LinearRegression::runRegression(.95, xs, ys)
       << endl << endl;
  // Test 8: same, online
  lr.reset(4);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i), .95);
  cout << lr.getCoefficients() << endl;
}
