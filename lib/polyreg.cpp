#include "polyreg.hpp"
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>

namespace complib {
  using namespace Eigen;
  bool PolynomialRegression::updateCoefficients(const vec &x, double y, double lambda) {
    bool b = plin.updateCoefficients(ppoly.expand(x), y, lambda);
    if(b)
      ppoly.coefficients = plin.getCoefficients();
    return(b);
  }
}
