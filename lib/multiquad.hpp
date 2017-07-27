#pragma once
#include "polyreg.hpp"
#include <vector>

namespace linreg {
  class Multiquad {
  public:
    inline Multiquad(vec lows, vec highs,
		     Eigen::VectorXi brackets, double omega = 0)
      : lows(lows), brackets(brackets),
	steps((highs - lows).array() / brackets.cast<double>().array()),
	lrs(brackets.prod(), PolynomialRegression(2, lows.size(), omega)) {}

    const vec lows, steps;
    const Eigen::VectorXi brackets;
    void addDataPoint(const vec &x, double y, double lambda = 1);
    // Get by indices
    const PolynomialRegression &getQRix(const Eigen::VectorXi &) const;
    // Get by coordinates
    const PolynomialRegression &getQR(const Eigen::VectorXd &) const;
  private:
    std::vector<PolynomialRegression> lrs;
  };
}
