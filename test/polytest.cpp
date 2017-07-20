#include "polyreg.hpp"
#include <random>
#include <ctime>
#include <iostream>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;
int main() {
  Eigen::MatrixX2d xs(1000, 2); // Why do I still need to pass the 2?
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(-1,1);
  std::mt19937_64 rands(std::time(NULL));
  Eigen::Vector3d fp(1, 2, 3);
  Eigen::Matrix<double, 6, 1> gp;
  gp << 1, 2, 5, 3, 4, 6;
  Eigen::VectorXd ys(1000);
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = norm(rands);
    xs(i,1) = norm(rands);
    ys(i) = uni(rands);
  }
  Eigen::VectorXd zs = ys;
  Eigen::MatrixX3d fs(1000, 3);
  Eigen::MatrixXd gs(1000, 6);
  
  linreg::PolynomialRegression prf(2, 1);
  linreg::PolynomialRegression prg(2,2);
  for(int i = 0; i < 1000; i++) {
    fs.row(i) = prf.polynomial(xs.row(i).head<1>()).transpose();
    gs.row(i) = prg.polynomial(xs.row(i).transpose()).transpose();
  }
  ys += fs * fp;
  zs += gs * gp;

  for(int i = 0; i < 1000; i++) {
    prf.updateCoefficients(xs.row(i).head<1>(), ys(i));
    prg.updateCoefficients(xs.row(i).transpose(), zs(i));
  }

  cout << prf.coefficients << endl << endl;
  cout << linreg::LinearRegression::runRegression(fs, ys) << endl << endl;
  cout << prg.coefficients << endl << endl;
  cout << linreg::LinearRegression::runRegression(gs, zs) << endl;
}
