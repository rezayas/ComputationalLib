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
	    + mat::Identity(xs.cols(), xs.cols()) * omega * w(0,0) * lambda)
	   .ldlt().solve(xs.transpose() * w * ys));
  }
  
  void LinearRegression::reset(int idim, double iomega) {
    points = 0;
    pdim = idim;
    omega = iomega;
    theta.setZero(dim);
    if(omega > 0) {
      ready = true;
      b = Eigen::MatrixXd::Identity(dim, dim) / omega;
    } else {
      xmat.setZero(dim, dim);
      yvec.setZero(dim);
      wvec.setZero(dim);
    }
  }

  
  bool LinearRegression::updateCoefficients(const vec &x, double y,
					    double lambda) {
    if(!ready) {
      if(points == xmat.rows()) {
	xmat.conservativeResize(points * 2, NoChange);
	yvec.conservativeResize(points * 2);
	wvec.conservativeResize(points * 2);
	wvec.tail(points).setZero();
	yvec.tail(points).setZero();
	xmat.bottomRows(points).setZero();
      }
      xmat.row(points) = x.transpose();
      yvec(points) = y;
      wvec.head(points) *= lambda;
      wvec(points) = 1;
      omega *= lambda;
      if(points >= dim) {
		// b = (xᵀwx + ωIλⁿ)⁻¹, where n is the number of xs so far.
		// w is the matrix that contains weighting.
	auto c = (xmat.transpose() * wvec.asDiagonal() * xmat
		  + mat::Identity(dim, dim) * omega).ldlt();
	if(c.rcond() > .0001) {
	  b = c.solve(mat::Identity(dim, dim));
	  theta = b * xmat.transpose() * wvec.asDiagonal() * yvec;
	  xmat.resize(0,0);
	  yvec.resize(0);
	  ready = true;
	}
      }
      points++;
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
    return(ready);
  }
}
