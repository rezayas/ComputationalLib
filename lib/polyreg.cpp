#include "polyreg.hpp"
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>

namespace linreg {
  using namespace Eigen;
  PolynomialRegression::PolynomialRegression(int order, int nvars, double omega)
    : exponents(pexp), coefficients(pcoef),
      lin((int)boost::math::binomial_coefficient<double>(order + nvars,
							 nvars), omega) {
    pexp.resize(lin.dim, nvars);
    pcoef.resize(pexp.rows());
    // The algorithm here is as follows:
    // We keep track of the sum of the powers, S.
    // Then, if S < O (the order of our polynomial),
    // the next row is the current row + (0,0,…,0,1),
    // and we set S := S + 1.
    // (We store the current row in es.)
    // Otherwise, we set i such that es(i) ≠ 0,
    // but if j > i, es(i) = 0.
    // Then we set S := S - es(i) + 1,
    // then es(i) := 0,
    // and es(i - 1) := es(i - 1) + 1.
    RowVectorXi es = RowVectorXi::Zero(nvars);
    for(int crow = 0, total = 0; crow < pexp.rows(); crow++) {
      pexp.row(crow) = es;
      if(total < order) {
	es(nvars - 1)++;
	total++;
      } else {
	int i = nvars;
	while(!es(--i));
	total -= es(i) - 1;
	es(i) = 0;
	if(i)
	  es(i - 1)++;
      }
    }
  }

  bool PolynomialRegression::updateCoefficients(const vec &x, double y, double lambda) {
    bool b = lin.updateCoefficients(polynomial(x), y, lambda);
    if(b)
      pcoef = lin.getCoefficients();
    return(b);
  }
  
  vec PolynomialRegression::polynomial(const vec &x) {
    vec ans(pexp.rows());
    for(int i = 0; i < ans.size(); i++) {
      double ai = 1;
      for(int j = 0; j < pexp.cols(); j++) {
	unsigned v = pexp(i, j);
	double d = x(j);
	while(v) {
	  if(v & 1)
	    ai *= d;
	  v >>= 1;
	  d *= d;
	}
      }
      ans(i) = ai;
    }
    return(ans);
  }
}
