#include "linreg.hpp"

using namespace std;
int main() {
  arma::arma_rng::set_seed_random();
  arma::mat xs(1000, 2, arma::fill::randn);
  arma::vec xp = {2, 3};
  arma::vec ys(1000, arma::fill::randu);
  ys *= 2;
  ys += xs * xp;
  linreg::LinearRegression lr(2);
  for(int i = 0; i < 1000; i++)
    lr.updateCoefficients(xs.row(i).t(), ys(i));
  cout << lr.getCoefficients();
  cout << "   " << lr.getConstant() << endl;
  cout << linreg::LinearRegression::runRegression(xs, ys);
}
