#include "polyreg.hpp"
#include <random>
#include <ctime>
#include <iostream>
#include <boost/math/special_functions/binomial.hpp>
#include <sstream>

typedef Eigen::Matrix<double, 1, 1> Matrix1d;

using namespace std;
int main() {
  Eigen::MatrixX2d xs(1000, 2); // Why do I still need to pass the 2?
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(-1,1);
  std::mt19937_64 rands(std::time(NULL));
  Eigen::VectorXd ys(1000), zs(1000);
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = norm(rands);
    xs(i,1) = norm(rands);
    ys(i) = uni(rands);
    zs(i) = uni(rands);
  }

  complib::Polynomial f(Eigen::VectorXi::LinSpaced(3, 0, 2),
			(Eigen::VectorXd(3) << 1, 2, -3).finished());

  complib::Polynomial g((Eigen::MatrixX2i(6,2) << 0,0, 1,0, 0,1,
			1,1, 2,0, 0,2).finished(),
			(Eigen::VectorXd(6) << 1, 2, -3, 4, -5, 6).finished());
  
  cout << f.evaluate(Matrix1d(0)) << ' ' << f.evaluate(Matrix1d(-1)) << endl;
  cout << g.evaluate(Eigen::Vector2d(0, 0)) << ' '
       << g.evaluate(Eigen::Vector2d(-1, 1)) << endl;
  cout << f.derivative(Matrix1d(0)) << ' '
       << f.derivative(Matrix1d(-1)) << endl;
  cout << g.derivative(Eigen::Vector2d(0,0)).transpose() << endl
       << g.derivative(Eigen::Vector2d(-1,1)).transpose() << endl << endl;
    
  for(int i = 0; i < 1000; i++) {
    ys(i) += f.evaluate(xs.row(i).head<1>());
    zs(i) += g.evaluate(xs.row(i).transpose());
  }
  
  for(double omega = 0; omega < .002; omega += .001)
    for(double lambda = 1; lambda > .94; lambda -= .05) {
      complib::PolynomialRegression prf(f, omega);
      complib::PolynomialRegression prg(g, omega);

      // Estimate with the online algorithm.
      // First line is estimating the parameters of f(x) = 1 + 2x - 3x^2.
      // Second is g(x) = 1 + 2x - 3y + 4xy - 5x^2 + 6y^2.
      for(int i = 0; i < 1000; i++) {
	prf.updateCoefficients(xs.row(i).head<1>(), ys(i), lambda);
	prg.updateCoefficients(xs.row(i).transpose(), zs(i), lambda);
      }

      // Estimate with the offline algorithm.
      // First line is estimating the parameters of f(x) = 1 + 2x - 3x^2.
      // Second is g(x) = 1 + 2x - 3y + 4xy - 5x^2 + 6y^2.
      f.fitToData(xs.col(0), ys, lambda, omega);
      g.fitToData(xs, zs, lambda, omega);

      /* Print out (first two should be about 1,2,-3;
	 last two should be 1,2,-3,4,-5,6 or so.) */
      cout << prf.poly().coefficients << endl << endl;
      cout << f.coefficients << endl << endl;
      cout << prg.poly().coefficients << endl << endl;
      cout << g.coefficients << endl << endl;
    }
}
