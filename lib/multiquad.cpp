#include "multiquad.hpp"
#include <stdexcept>

namespace complib {
  void Multiquad::addDataPoint(const vec &x, double y, double lambda) {
    // Basically, for each i,
    // test if brackets(i) is 0 (in which case it's fine),
    // or if x(i) < mins(i) or x(i) > mins(i) + steps(i) * brackets(i)
    // (in which case it's not fine).
    // If all are fine, then continue.
    if((!brackets.cast<bool>().array()
	|| ((x.array() >= mins.array())
	    && ((x - mins).array() <=
		(stps.cwiseProduct(brackets.cast<double>())).array())))
       .all())
      // This is basically a hack.
      // For technical reasons, it's really hard to avoid code duplication
      // without this cast.
      ((PolynomialRegression &)getQR(x)).updateCoefficients(x, y, lambda);
  }
  // Get by indices
  const PolynomialRegression
  &Multiquad::getQRix(const Eigen::VectorXi &ix) const {
    int posn = 0;
    // Basically, work backwards.
    // The formula is posn = ∑_i=0^n-1 ix(i)∏_j=i+1^n-1 brackets(j).
    for(int i = brackets.size() - 1, step = 1; i >= 0; i--) {
      if(brackets(i)) {
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
    }
    return(regs[posn]);
  }
  // Get by coordinates
  const PolynomialRegression
  &Multiquad::getQR(const Eigen::VectorXd &x) const {
    try {
      // The method here is to take (x-lows)/stps
      // (interpreted coefficientwise), and then take the floors.
      return(getQRix(((x - mins).array() / stps.array()).cast<int>()));
    } catch(const std::out_of_range &e) {
      // If that throws an exception,
      // we translate it to what the actual problem was.
      int i;
      for(i = 0; !brackets(i) ||
	    (x(i) > mins(i) && x(i) < stps(i) * brackets(i) + mins(i)); i++);
      std::ostringstream msg;
      msg << "Coordinate " << i << "(value: " << x(i) << ") is too ";
      if(x(i) < mins(i))
	msg << "low (minimum: " << mins(i);
      else
	msg << "high (maximum: " << stps(i) * brackets(i) + mins(i);
      msg << ')';
      throw std::out_of_range(msg.str());
    }
  }
}
