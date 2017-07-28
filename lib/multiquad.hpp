#pragma once
#include "polyreg.hpp"
#include <vector>

namespace complib {
  class Multiquad {
  public:
    /*************************************************************************\
     * Pass in a vector of low points, a vector of high points,		     *
     * and a vector of integers that gives the number of brackets.	     *
     * So for example, passing in the vectors (0, 0, 0),		     *
     * (1, 3, 5), and (1, 2, 4) means that we will cover the space	     *
     * with 0 < x₀ < 1, 0 < x₁ < 3, 0 < x₂ < 5,				     *
     * and that the boxes break at x₀ = 1.5, x₂ ϵ {1.25, 2.5, 3.75}.	     *
     * Passing in a zero for the number of brackets means that this	     *
     * coordinate should not be considered for purposes of determining which *
     * box the point lies in.						     *
     * 									     *
     * Explanation of the stps line:					     *
     * (highs - lows) is a vector subtraction.				     *
     * brackets.cast<double> simply turns the integers into floating-point   *
     * numbers.								     *
     * .array() turns it into arrays,					     *
     * for which mathematical operations are coordinate-wise.		     *
     * 									     *
     * Explanation of the regs line;					     *
     * We compute the number of polynomial regressions we want by	     *
     * turning brackets into an array and computing the element-wise minimum *
     * of its elements with 1, then multiplying all of these.		     *
     * Then we create a quadratic regression of the proper dimension,	     *
     * and then make regs a list of that many quadratic regressions.	     *
    \*************************************************************************/
    inline Multiquad(const vec &lows, const vec &highs,
		     const Eigen::VectorXi &brackets, double omega = 0)
      : mins(lows.array().min(highs.array())), brackets(brackets),
	stps((highs - lows).array().abs() / brackets.cast<double>().array()),
	regs(brackets.array().min(1).prod(),
	    PolynomialRegression(2, lows.size(), omega)) {
      // Basically, if there are no brackets for this coordinate
      // so we should ignore it,
      // set the step size to infinity.
      // This prevents float-to-int errors.
      for(int i = 0; i < mins.size(); i++)
	if(!brackets(i))
	  stps(i) = INFINITY;
    }

    const vec mins;
    inline const vec &steps() const {
      return(stps);
    }
    const Eigen::VectorXi brackets;
    // Add a data point to the pertinent regression.
    // Does nothing if the data point is outside the limits.
    void addDataPoint(const vec &x, double y, double lambda = 1);
    // Get by indices.
    // Throws std::out_of_range if any ix(i) is negative or at least
    // brackets(i).
    // But if brackets(i) is 0, then ix(i) isn't considered.
    const PolynomialRegression &getQRix(const Eigen::VectorXi &ix) const;
    // Get by coordinates.
    // Throws std::out_of_range if any coordinate is outside the limits.
    // If brackets(i) is 0, then the limits for xᵢ do not exist.
    const PolynomialRegression &getQR(const Eigen::VectorXd &) const;
  private:
    std::vector<PolynomialRegression> regs;
    vec stps;
  };
}
