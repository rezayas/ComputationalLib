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
  Eigen::VectorXd ys(1000), zs(1000);
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = norm(rands);
    xs(i,1) = norm(rands);
    ys(i) = uni(rands);
    zs(i) = uni(rands);
  }

  complib::Polynomial f(Eigen::VectorXi::LinSpaced(3, 0, 2),
		       Eigen::VectorXd::LinSpaced(3, 1, 3));

  complib::Polynomial g((Eigen::MatrixX2i(6,2) << 0,0, 1,0, 0,1,
			1,1, 2,0, 0,2).finished(),
		       Eigen::VectorXd::LinSpaced(6, 1, 6));

  complib::PolynomialRegression prf(f);
  complib::PolynomialRegression prg(g);

  for(int i = 0; i < 1000; i++) {
    ys(i) += f.evaluate(xs.row(i).head<1>());
    zs(i) += g.evaluate(xs.row(i).transpose());
  }

  // Estimate with the online algorithm.
  for(int i = 0; i < 1000; i++) {
    prf.updateCoefficients(xs.row(i).head<1>(), ys(i));
    prg.updateCoefficients(xs.row(i).transpose(), zs(i));
  }

  // Estimate with the offline algorithm
  f.fitToData(xs.col(0), ys);
  g.fitToData(xs, zs);

  /* Print out (first two should be about 1,2,3;
     last two should be 1,2,3,4,5,6 or so.) */
  cout << prf.poly().coefficients << endl << endl;
  cout << f.coefficients << endl << endl;
  cout << prg.poly().coefficients << endl << endl;
  cout << g.coefficients << endl;
}
