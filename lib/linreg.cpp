#include "linreg.hpp"
#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace complib {
  using namespace Eigen;
  vec LinearRegression::runRegression(double lambda, double omega,
				      const mat &xs, const vec &ys) {
    // If x is square and ω = 0, then it turns out we might as well
    // do x⁻¹y (because (xᵀwx)⁻¹(xᵀwy) = x⁻¹w⁻¹xᵀ⁻¹xᵀwy = x⁻¹w⁻¹wy = x⁻¹y)
    if(xs.rows() == xs.cols() && omega == 0)
      return(xs.fullPivLu().solve(ys));
    mat w = mat::Zero(ys.size(), ys.size());
    for(int i = ys.size() - 1, p = 1; i >= 0; i--, p /= lambda)
      w(i,i) = p;
    return((xs.transpose() * w * xs
	    + mat::Identity(xs.cols(), xs.cols()) * omega * w(0,0) * lambda)
	   .ldlt().solve(xs.transpose() * w * ys));
  }
  
  void LinearRegression::reset(int idim, double iomega) {
    points = 0;
    pdim = idim;
    omega = iomega;
    theta.setZero(idim);
    if((ready = omega > 0))
      b = Eigen::MatrixXd::Identity(idim, idim) / omega;
    else {
      xmat.setZero(idim, idim);
      yvec.setZero(idim);
      wvec.setZero(idim);
    }
  }

  
  bool LinearRegression::updateCoefficients(const vec &x, double y,
					    double lambda) {
    if(!ready) {
      if(points == xmat.rows()) {
	// If we're out of room, double the space
	// (this is the right thing to do)
	xmat.conservativeResize(points * 2, NoChange);
	yvec.conservativeResize(points * 2);
	wvec.conservativeResize(points * 2);
	wvec.tail(points).setZero();
	yvec.tail(points).setZero();
	xmat.bottomRows(points).setZero();
      }
      // Insert the new x and y in our list
      xmat.row(points) = x.transpose();
      yvec(points) = y;
      // Decrease the early weights by a factor of lambda
      wvec.head(points) *= lambda;
      // And set this point's weight to 1.
      wvec(points) = 1;
      // (ω is like point 0)
      omega *= lambda;
      // If we have too few points, why bother checking if we can invert?
      if(points >= pdim) {
		// b = (xᵀwx + ωIλⁿ)⁻¹, where n is the number of xs so far.
		// w is the matrix that contains weighting.
	auto c = (xmat.transpose() * wvec.asDiagonal() * xmat
		  + mat::Identity(pdim, pdim) * omega).ldlt();
	if(c.rcond() > .0001) {
	  // If our matrix isn't almost singular, we now set b and θ.
	  b = c.solve(mat::Identity(pdim, pdim));
	  // (θ = bxᵀwy)
	  theta = c.solve(xmat.transpose() * wvec.cwiseProduct(yvec));
	  // We won't need these anymore.
	  xmat.resize(0,0);
	  yvec.resize(0);
	  wvec.resize(0);
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
      double gamma = lambda + x.dot(b * x);
      theta += (b * x * (y - theta.dot(x)) / gamma).eval();
      b -= (b * x * x.transpose() * b / gamma).eval();
      b /= lambda;
    }
    return(ready);
  }
}
