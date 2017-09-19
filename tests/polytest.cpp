#include "../include/ComputationalLib/polyreg.hpp"
#include <random>
#include <ctime>
#include <iostream>
#include <boost/math/special_functions/binomial.hpp>
#include <sstream>

using namespace std;
int main() {
  // xs is a matrix of 1000 columns, 2 rows of doubles
  Eigen::MatrixX2d xs(1000, 2); // Why do I still need to pass the 2?
  std::normal_distribution<> norm;
  std::uniform_real_distribution<> uni(-1,1);
  std::mt19937_64 rands(std::time(NULL));
  // ys and zs are vectors of length 1000
  Eigen::VectorXd ys(1000), zs(1000);
  for(int i = 0; i < 1000; i++) {
    xs(i,0) = norm(rands);
    xs(i,1) = norm(rands);
    ys(i) = uni(rands);
    zs(i) = uni(rands);
  }

// this is a polynomial equivalent to f(x) = 1 + 2x + 3x^2
// this is a int (i) dynamic-sized (X) vector (1 column matrix) of size 3: 0, 1, 2
  ComputationalLib::Polynomial f(Eigen::VectorXi::LinSpaced(3, 0, 2),
		       Eigen::VectorXd::LinSpaced(3, 1, 3));
// ^ this is a double (d) dynamic-sized (X) vector (1 column matrix) of size 3: 1, 2, 3

// this is a polynomial equivalent to g(x,y) = 1 + 2x + 3y + 4xy + 5x^2 + 6y^2
// this is a matrix of with dynamic (X) number of rows and 2 columns:
// exs col1 (x exponents): [0,1,0,1,2,0]ᵀ, exs col2 (y exponents): [0,0,1,1,0,2]ᵀ
  ComputationalLib::Polynomial g((Eigen::MatrixX2i(6,2) << 0,0, 1,0, 0,1,
			1,1, 2,0, 0,2).finished(),
		       Eigen::VectorXd::LinSpaced(6, 1, 6));
  // this is a vector D of: 1, 2, 3, 4, 5, 6

// prf is polynomial regression of f.
  ComputationalLib::PolynomialRegression prf(f, 0.001);
  ComputationalLib::PolynomialRegression prg(g);

// ys start with uni(rands)
// then you evaluate f(x) at the appropriate value of x to get each y.
// zs also starts with uni(rands)
// but you evaluate g(x,y) instead
  for(int i = 0; i < 1000; i++) {
    // returns a 1x1 block containing just the i-th row, 1-st element
    ys(i) += f.evaluate(xs.row(i).head<1>());
    // returns a column vector equivalent to the transpose of the i-th row
    zs(i) += g.evaluate(xs.row(i).transpose());
  }

  // Estimate with the online algorithm.
  // First line is estimating the parameters of f(x) = 1 + 2x + 3x^2.
  // Second is g(x) = 1 + 2x + 3y + 4xy + 5x^2 + 6y^2.
  for(int i = 0; i < 1000; i++) {
    prf.updateCoefficients(xs.row(i).head<1>(), ys(i), 0.99);
    prg.updateCoefficients(xs.row(i).transpose(), zs(i));
  }

  // Estimate with the offline algorithm.
  // First line is estimating the parameters of f(x) = 1 + 2x + 3x^2.
  // Second is g(x) = 1 + 2x + 3y + 4xy + 5x^2 + 6y^2.
  f.fitToData(xs.col(0), ys);
  g.fitToData(xs, zs);

  /* Print out (first two should be about 1,2,3;
     last two should be 1,2,3,4,5,6 or so.) */
  cout << prf.poly().coefficients << endl << endl;
  cout << f.coefficients << endl << endl;
  cout << prg.poly().coefficients << endl << endl;
  cout << g.coefficients << endl << endl;


  // Literally no idea what is going on after this point:
  ostringstream crtstr;
  crtstr << f << g;
  istringstream readr(crtstr.str());
  ComputationalLib::Polynomial a = ComputationalLib::Polynomial::readIn(readr);
  ComputationalLib::Polynomial b = ComputationalLib::Polynomial::readIn(readr);

  cout << (a.exponents() - f.exponents()).array().abs().maxCoeff() << endl;
  cout << (a.coefficients - f.coefficients).array().abs().maxCoeff() << endl;
  cout << (b.exponents() - g.exponents()).array().abs().maxCoeff() << endl;
  cout << (b.coefficients - g.coefficients).array().abs().maxCoeff() << endl;

}
