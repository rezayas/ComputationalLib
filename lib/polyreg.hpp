#include "linreg.hpp"

namespace linreg {
  class PolynomialRegression {
  public:
    const Eigen::MatrixXi &exponents;
    PolynomialRegression(int order, int nvars);
    const vec &coefficients;
    bool updateCoefficients(const vec &x, double y);
    vec polynomial(const vec &);
  private:
    LinearRegression lin;
    Eigen::MatrixXi pexp;
    vec pcoef;
  };
}
