#pragma once
#include <Eigen/Core>

namespace linreg {
  typedef Eigen::MatrixXd mat;
  typedef Eigen::VectorXd vec;
  class LinearRegression {
  public:
    static vec runRegression(double, const mat &, const vec &);
    inline static vec runRegression(const mat &xs,
				    const vec &ys) {
      return(runRegression(1, xs, ys));
    }
    inline LinearRegression(int dim, double lambda = 1) : theta(dim + 1),
							  b(dim + 1, dim + 1),
							  yvec(dim + 1),
							  xmat(dim+1, dim+1),
							  points(dim + 1),
							  dim(dim),
							  lambda(lambda) {
      xmat.col(dim).setOnes();
    }
    bool updateCoefficients(const vec &, double);
    inline const vec getCoefficients() const {
      return(theta.head(dim));
    }
    inline double getConstant() const {
      return(theta(dim));
    }
    const int dim;
    const double lambda;
  private:
    vec theta, yvec;
    mat b, xmat;
    int points;
  };
}
