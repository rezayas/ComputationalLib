#include "linreg.hpp"

namespace linreg {
  class PolynomialRegression {
  public:
    const arma::umat &exponents;
    PolynomialRegression(int order, int nvars);
    const arma::vec &coefficients;
    bool updateCoefficients(const arma::vec &x, double y);
    arma::vec polynomial(const arma::vec &);
  private:
    LinearRegression lin;
    arma::umat pexp;
    arma::vec pcoef;
  };
}
