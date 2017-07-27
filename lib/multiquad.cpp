#include "multiquad.hpp"
#include <stdexcept>

namespace linreg {
  void Multiquad::addDataPoint(const vec &x, double y, double lambda) {
    const_cast<PolynomialRegression &>(getQR(x))
      .updateCoefficients(x, y, lambda);
  }
  // Get by indices
  const PolynomialRegression
  &Multiquad::getQRix(const Eigen::VectorXi &ix) const {
    int posn = 0;
    for(int i = brackets.size() - 1, step = 1; i >= 0; i--) {
      if(ix(i) < 0) {
	std::ostringstream msg;
	msg << "Index " << i << " (value: " << ix(i) << ") is negative";
	throw std::out_of_range(msg.str());
      } else if(ix(i) >= brackets(i)) {
	std::ostringstream msg;
	msg << "Index " << i << " (value: " << ix(i)
	    << ") is too big (dimension: " << brackets(i) << ')';
        throw std::out_of_range(msg.str());
      }
      posn += ix(i) * step;
      step *= brackets(i);
    }
    return(lrs[posn]);
  }
  // Get by coordinates
  const PolynomialRegression
  &Multiquad::getQR(const Eigen::VectorXd &x) const {
    try {
      return(getQRix(((x - lows).array() / steps.array()).cast<int>()));
    } catch(const std::out_of_range &e) {
      int i;
      for(i = 0; x(i) > lows(i)
	    && x(i) < steps(i) * brackets(i) + lows(i); i++);
      std::ostringstream msg;
      msg << "Coordinate " << i << "(value: " << x(i) << ") is too ";
      if(x(i) < lows(i))
	msg << "low (minimum: " << lows(i);
      else
	msg << "high (maximum: " << steps(i) * brackets(i) + lows(i);
      msg << ')';
      throw std::out_of_range(msg.str());
    }
  }
}
