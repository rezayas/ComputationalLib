#include "polyreg.hpp"
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>

namespace linreg {
  using namespace arma;
  PolynomialRegression::PolynomialRegression(int order, int nvars)
    : exponents(pexp), coefficients(pcoef),
      lin((uword)boost::math::binomial_coefficient<double>(order + nvars,
							   nvars) - 1) {
    pexp.set_size(lin.dim + 1, nvars);
    pcoef.set_size(pexp.n_rows);
    urowvec es(nvars, fill::zeros);
    for(int crow = 0, total = 0; crow < pexp.n_rows; crow++) {
      pexp.row(crow) = es;
      if(total < order) {
	es(nvars)++;
	total++;
      } else {
	int i = nvars;
	while(!es(--i));
	total -= es(i);
	es(i) = 0;
	if(i)
	  es(i - 1)++;
      }
    }
  }

  bool PolynomialRegression::updateCoefficients(const vec &x, double y) {
    bool b = lin.updateCoefficients(polynomial(x).tail(lin.dim), y);
    if(b) {
      pcoef.tail(lin.dim) = lin.getCoefficients();
      pcoef(0) = lin.getConstant();
    }
    return(b);
  }
  
  vec PolynomialRegression::polynomial(const vec &x) {
    vec ans(pexp.n_rows);
    for(int i = 0; i < ans.n_elem; i++) {
      double ai = 1;
      for(int j = 0; j < pexp.n_cols; j++) {
	uword v = pexp(i, j);
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
