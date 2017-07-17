#pragma once
#include <armadillo>

namespace linreg {

  class LinearRegression {
  public:
    static arma::vec runRegression(double, const arma::mat &, const arma::vec &);
    inline static arma::vec runRegression(const arma::mat &xs,
					  const arma::vec &ys) {
      return(runRegression(1, xs, ys));
    }
    inline LinearRegression(int dim, double lambda = 1) : theta(dim + 1),
							  b(dim + 1, dim + 1),
							  yvec(dim + 1),
							  xmat(dim+1, dim+1),
							  points(dim + 1),
							  dim(dim),
							  lambda(lambda) {
      xmat.col(dim).ones();
    }
    bool updateCoefficients(const arma::vec &, double);
    inline const arma::vec getCoefficients() const {
      return(theta.head(dim));
    }
    inline double getConstant() const {
      return(theta(dim));
    }
    const int dim;
    const double lambda;
  private:
    arma::vec theta, yvec;
    arma::mat b, xmat;
    int points;
  };
}
