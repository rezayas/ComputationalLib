#include "linreg.hpp"
#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace linreg {
  using namespace Eigen;
  vec LinearRegression::runRegression(double lambda,
				      const mat &xs, const vec &ys) {
    mat xmat(xs.rows(), xs.cols() + 1);
    xmat.col(xs.cols()).fill(1);
    xmat.leftCols(xs.cols()) = xs;
    // If x is square, then it turns out we might as well
    // do x⁻¹y
    if(xmat.rows() == xmat.cols())
      return(xmat.fullPivLu().solve(ys));
    mat w = MatrixXd::Zero(ys.size(), ys.size());
    for(int i = ys.size() - 1, p = 1; i >= 0; i--, p /= lambda)
      w(i,i) = p;
    // (xᵀwx)⁻¹(xᵀwy)
    return((xmat.transpose() * w * xmat).ldlt().solve(xmat.transpose() * w * ys));
  }
  
  bool LinearRegression::updateCoefficients(const vec &x, double y) {
    if(points) {
      xmat.row(--points).head(dim) = x.transpose();
      yvec(points) = y;
      if(!points) {
	// Translation: Use LDLT to construct the inverse.
	b = (xmat.transpose() * xmat).ldlt().solve(mat::Identity(dim + 1, dim + 1));
	theta = runRegression(lambda, xmat.leftCols(dim), yvec);
      }
    } else {
      vec xv(dim + 1);
      xv.head(dim) = x;
      xv(dim) = 1;
      // The .eval() function does nothing, but it prevents issues
      // due to modifying what we are using.
      // γ = λ² + xᵀbx
      // θ' = θ + bx(y-θᵀx)/γ
      // b' = (b - bxxᵀb/γ)/λ²
      double gamma = lambda * lambda + (xv.transpose() * b * xv)(0,0);
      theta += (b * xv * (y - theta.transpose() * xv) / gamma).eval();
      b -= (b * xv * xv.transpose() * b / gamma).eval();
      b /= lambda * lambda;
    }
    return(!points);
  }
}
