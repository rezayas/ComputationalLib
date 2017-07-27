#pragma once
#include "polyreg.hpp"
#include <vector>

namespace linreg {
  class Multiquad {
  public:
    inline Multiquad(const vec &lows, const vec &highs,
		     const Eigen::VectorXi &brackets, double omega = 0)
      : lows(lows), brackets(brackets),
	stps((highs - lows).array() / brackets.cast<double>().array()),
	lrs(brackets.prod(), PolynomialRegression(2, lows.size(), omega)) {
      for(int i = 0; i < lows.size(); i++)
	if(stps(i) != stps(i))
	  stps(i) = INFINITY;
    }

    const vec lows;
    inline const vec &steps() const {
      return(stps);
    }
    const Eigen::VectorXi brackets;
    void addDataPoint(const vec &x, double y, double lambda = 1);
    // Get by indices
    const PolynomialRegression &getQRix(const Eigen::VectorXi &) const;
    // Get by coordinates
    const PolynomialRegression &getQR(const Eigen::VectorXd &) const;
  private:
    std::vector<PolynomialRegression> lrs;
    vec stps;
  };
}
