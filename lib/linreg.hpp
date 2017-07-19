#pragma once
#include <Eigen/Core>

namespace linreg {
  typedef Eigen::MatrixXd mat;
  typedef Eigen::VectorXd vec;

  class LinearRegression {
  public:
	
	// theta = (X'WX)^-1 X'WY, where the (i,i)th element of W is lambda^(n-i-1); n is the number of observations
    static vec runRegression(double, const mat &, const vec &);

	// theta = (X'X)^-1 X'Y
	inline static vec runRegression(const mat &xs, const vec &ys) {
		return(runRegression(1, xs, ys));
	}

	// theta = (X'WX + wI)^-1 X'WY, where w > 0 is the L2 regularization factor
	static vec runRegression(double, double, const mat &, const vec &);

	// constructor to use the training algorithm to estimate coefficient
	// omega > 0 is the L2 regularization factor
    inline LinearRegression(int dim, double omega = 0) : 
							  theta(dim + 1),
							  b(dim + 1, dim + 1),
							  yvec(dim + 1),
							  xmat(dim+1, dim+1),
							  points(dim + 1),
							  dim(dim),
							  omega(omega) {
		xmat.col(dim).setOnes();
    }

	// update theta based on new x' and y using forgetting factor lambda
    bool updateCoefficients(const vec &, double, double = 1);

	// get theta
    inline const vec getCoefficients() const {
      return(theta.head(dim));
    }

	// TODO: REMOVE
    inline double getConstant() const {
      return(theta(dim));
    }

    const int dim; // mumber of columns of X
	const double omega; // L2 regularization factor

  private:
    vec theta, yvec;
    mat b, xmat;
    int points;
  };
}
