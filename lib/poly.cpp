#include "poly.hpp"
#include <boost/math/special_functions/binomial.hpp>

using namespace Eigen;
namespace linreg {
  Polynomial::Polynomial(int order, int nvars) :
    monomials((int)boost::math::binomial_coefficient<double>(order + nvars,
							     nvars)),
    variables(nvars), exponents(exps) {
    coefficients = VectorXd::Zero(monomials);
    exps.resize(monomials, nvars);
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
    for(int crow = 0, total = 0; crow < monomials; crow++) {
      exps.row(crow) = es;
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
  VectorXd Polynomial::expand(const VectorXd &x) const {
    VectorXd ans(monomials);
    for(int i = 0; i < monomials; i++) {
      double ai = 1;
      for(int j = 0; j < variables; j++) {
	unsigned v = exps(i, j);
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
  
  double Polynomial::evaluate(const VectorXd &x) const {
    double ans;
    for(int i = 0; i < monomials; i++) {
      double ai = 1;
      for(int j = 0; j < variables; j++) {
	unsigned v = exps(i, j);
	double d = x(j);
	while(v) {
	  if(v & 1)
	    ai *= d;
	  v >>= 1;
	  d *= d;
	}
      }
      ans += ai;
    }
    return(ans);
  } 
}