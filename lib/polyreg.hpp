#include "linreg.hpp"

namespace linreg {
  class PolynomialRegression {
  public:
    const Eigen::MatrixXi &exponents;
    PolynomialRegression(int order, int nvars, double omega = 0);
    const vec &coefficients;
    bool updateCoefficients(const vec &x, double y, double lambda = 1);
    vec polynomial(const vec &);
  private:
    LinearRegression lin;
    Eigen::MatrixXi pexp;
    vec pcoef;
  };
}
