#include "linreg.hpp"
#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace linreg {
  using namespace Eigen;
  vec LinearRegression::runRegression(double lambda, double omega,
				      const mat &xs, const vec &ys) {
    // If x is square, then it turns out we might as well
    // do x⁻¹y
    if(xs.rows() == xs.cols() && omega != 0)
      return(xs.fullPivLu().solve(ys));
    mat w = mat::Zero(ys.size(), ys.size());
    for(int i = ys.size() - 1, p = 1; i >= 0; i--, p /= lambda)
      w(i,i) = p;
    // (xᵀwx)⁻¹(xᵀwy)
    return((xs.transpose() * w * xs
	    + mat::Identity(xs.cols(), xs.cols()) * omega * w(0,0))
	   .ldlt().solve(xs.transpose() * w * ys));
  }
  
  bool LinearRegression::updateCoefficients(const vec &x, double y,
					    double lambda) {
    if(points) {
      int r = dim - points--;
      xmat.row(r) = x.transpose();
      yvec(r) = y;
      b *= lambda;
      b(r,r) = 1;
      if(!points) {
	// b = (xᵀwx + ωIλⁿ)⁻¹, where n is the number of xs so far.
	// w is the matrix that contains weighting.
	// We never need both b and w, so we do it like this.
	b = (xmat.transpose() * b * xmat
	     + mat::Identity(xmat.rows(), xmat.rows()) * omega * b(0,0))
	  .ldlt().solve(mat::Identity(dim, dim)).eval();
	theta = b * yvec;
      }
    } else {
      // The .eval() function does nothing, but it prevents issues
      // due to modifying what we are using.
      // γ = λ + xᵀbx
      // θ' = θ + bx(y-θᵀx)/γ
      // b' = (b - bxxᵀb/γ)/λ
      double gamma = lambda + (x.transpose() * b * x)(0,0);
      theta += (b * x * (y - theta.transpose() * x) / gamma).eval();
      b -= (b * x * x.transpose() * b / gamma).eval();
      b /= lambda;
    }
    return(!points);
  }
}
