#include "linreg.hpp"

namespace linreg {
  using namespace arma;
  vec LinearRegression::runRegression(double lambda,
				      const mat &xs, const vec &ys) {
    mat xmat(xs.n_rows, xs.n_cols + 1);
    xmat.unsafe_col(xs.n_cols).fill(1);
    xmat.head_cols(xs.n_cols) = xs;
    if(xmat.n_rows == xmat.n_cols)
      return(solve(xmat, ys));
    mat w(ys.n_elem, ys.n_elem, fill::zeros);
    for(int i = ys.n_elem - 1, p = 1; i >= 0; i--, p /= lambda)
      w(i,i) = p;
    return(solve(xmat.t() * w * xmat, xmat.t() * w * ys));
  }
  
  bool LinearRegression::updateCoefficients(const vec &x, double y) {
    if(points) {
      xmat.row(--points).head(dim) = x.t();
      yvec(points) = y;
      if(!points) {
	b = (xmat.t() * xmat).i();
	theta = runRegression(lambda, xmat.head_cols(dim), yvec);
      }
    } else {
      vec xv(dim + 1);
      xv.head(dim) = x;
      xv(dim) = 1;
      double gamma = lambda + as_scalar(xv.t() * b * xv);
      theta += b * xv * (y - theta.t() * xv) / gamma;
      b -= b * xv * xv.t() * b / gamma;
      b /= lambda;
    }
    return(!points);
  }
}
