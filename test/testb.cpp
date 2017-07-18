#include "linreg.hpp"
#include <random>
#include <ctime>
#include <iostream>

using namespace std;
int main() {
  Eigen::MatrixX2d xs(1000, 2); // Why do I still need to pass the 2?
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(0,2);
  std::mt19937_64 rands(std::time(NULL));
  Eigen::Vector2d xp(2, 3);
  Eigen::VectorXd ys(1000);
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = norm(rands);
    xs(i,1) = norm(rands);
    ys(i) = uni(rands);
  }
  ys += xs * xp;
  linreg::LinearRegression lr(2);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i));
  cout << lr.getCoefficients() << endl;
  cout << lr.getConstant() << endl << endl;
  cout << linreg::LinearRegression::runRegression(xs, ys) << endl;
}
