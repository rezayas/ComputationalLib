#include "linreg.hpp"
#include <random>
#include <ctime>
#include <iostream>

using namespace std;
int main() {
  
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(-1,1);
  std::mt19937_64 rands(std::time(NULL));

  Eigen::Vector3d xp(1, 2, 3); // fitting y = 1 + 2x1 + 3x3 
  Eigen::MatrixX3d xs(1000, 3); // X matrix
  Eigen::VectorXd ys(1000);	// y vector
  Eigen::VectorXd zs(1000);
  // populate X matrix and y vector
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = 1;
    xs(i,1) = norm(rands);
    xs(i,2) = norm(rands);
    ys(i) = zs(i) = uni(rands);
  }
  ys += xs * xp;

  // Test 1: regular LS 
  cout << linreg::LinearRegression::runRegression(xs, ys) << endl << endl;

  // Test 2: updating algorithm with omega = 0, lambda = 1 (we should get the same coefficient as Test 1)
  linreg::LinearRegression lr(3);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).transpose(), ys(i));
  cout << lr.getCoefficients() << endl << endl;


  // Test 3: LS + omega > 0 on a X matrix that is not full rank
  // fiting y = 2 + 3x1 + 4x3 where x2 is always 2*x1
  xs.col(2) = xs.col(1) * 2;
  zs += xs * Eigen::Vector3d(2, 3, 4)  
  
  linreg::LinearRegression lrb(3, .5);
  cout << linreg::LinearRegression::runRegression(1, .5, xs, zs) << endl << endl;
  
  // Test 4: updating algorithm with omega > 0 on a X matrix that is not full rank (same results as Test 3)
  for(int i = 0; i < 1000; i++)
    lrb.updateCoefficients(xs.row(i).transpose(), zs(i));
  cout << lrb.getCoefficients() << endl ;
  

}
