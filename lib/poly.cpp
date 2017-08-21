#include "poly.hpp"
#include <boost/math/special_functions/binomial.hpp>
#include "linreg.hpp"

using namespace Eigen;
namespace complib {
  Polynomial::Polynomial(int order, int nvars) :
    monomials((int)boost::math::binomial_coefficient<double>(order + nvars,
							     nvars)),
    variables(nvars) {
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
    // In other words, we find the last nonzero element of es.
    // Then we set S := S - es(i) + 1,
    // then es(i) := 0,
    // and es(i - 1) := es(i - 1) + 1.
    // So if order is 2 and nvars is 2,
    // we do as follows (recording only at the time of writing to the matrix)
    // | es(0) | es(1) | total | notes |
    // |-------+-------+-------+-------|
    // |   0   |   0   |   0   |  t<o  |
    // |   0   |   1   |   1   |  t<o  |
    // |   0   |   2   |   2   |  t=o  |
    // |   1   |   0   |   1   |  t<o  |
    // |   1   |   1   |   2   |  t=o  |
    // |   2   |   0   |   2   | done  |

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
      // For each monomial:
      // Set ai to 1.
      double ai = 1;
      // For each coordinate:
      for(int j = 0; j < variables; j++) {
	// Set v to the exponent of the coordinate in this monomial
	unsigned v = exps(i, j);
	// Set d to the value of the coordinate
	double d = x(j);
	// From now on, let w be the number of times x(j) has been multiplied
	// into ai, and let n be such that d = x(j)^n.
	// We mantain the invariant that nv + w = exps(i,j).
	// We start with n=1, v=exps(i,j), w = 0.
	// If v is odd, we set v' = (v-1)/2, w' = w + n, n' = 2n,
	// in which case
	// n'v'+w' = 2n(v-1)/2 + w + n = n(v-1) + w + n
	// = nv - n + w + n = nv + w;
	// if v is even, we set v' = v/2, w' = w, n' = 2n,
	// in which case n'v'+w' = 2nv/2 + w = nv+w.
	while(v) {
	  if(v & 1)
	    ai *= d;
	  v >>= 1;
	  d *= d;
	}
      }
      // Now store ai into the answer.
      ans(i) = ai;
    }
    return(ans);
  }

  // The algorithm here is the same as in expand,
  // except that we evaluate the polynomial completely.
  double Polynomial::evaluate(const VectorXd &x) const {
    double ans=0;
    // For each monomial
    for(int i = 0; i < monomials; i++) {
      // Take the coefficient
      double ai = coefficients(i);
      // Multiply it by the proper powers of the variables
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
      // And add it to the answer.
      ans += ai;
    }
    return(ans);
  }

  void Polynomial::fitToData(const MatrixXd &xs, const VectorXd &ys,
			     double lambda, double omega) {
    // Create a matrix to fit the monomials
    MatrixXd longs(xs.rows(), monomials);
    // Now for each row in the xs,
    for(int i = xs.rows() - 1; i >= 0; i--)
      // Set the corresponding row of the monomial matrix properly
      longs.row(i) = expand(xs.row(i).transpose()).transpose();
    // We're done with that.
    // Do linear regression.
    coefficients = LinearRegression::runRegression(lambda, omega, longs, ys);
  }

  std::tuple<MatrixXd, VectorXd, double> Polynomial::quaddec() const {
    MatrixXd g = MatrixXd::Zero(variables, variables);
    VectorXd v = VectorXd::Zero(variables);
    double d = 0;
    // For each monomial:
    for(int i = 0; i < monomials; i++) {
      int v1 = -1, v2 = -1, order = 0;
      // Go through the variables. Set order to the total order of this
      // monomial, unless that's higher than 2, in which case we don't care.
      // If the order is no more than 2,
      // set v1 to the index of the first variable that is included,
      // and v2 to that of the second.
      // If a variable is included twice, set both v1 and v2 to it.
      // (If there aren't enough variables included, we don't care.)
      for(int j = 0; j < variables && order <= 2; j++) {
	char c = exps(i, j);
	order += c;
	if(c == 1) {
	  if(v1 < 0)
	    v1 = j;
	  else
	    v2 = j;
	} else if(c == 2)
	  v1 = v2 = j;
      }
      // If no variables are included, this is a constant term,
      // so add its coefficient to d.
      if(order == 0)
	d += coefficients(i);
      // If the term is first-order, including variable i,
      // add its coefficient to v(i).
      else if(order == 1)
	v(v1) += coefficients(i);
      // If it's second order, and is of the form cx_ix_j,
      // add c to g(i, j) and g(j, i).
      else if(order == 2) {
	g(v1, v2) += coefficients(i);
	g(v2, v1) += coefficients(i);
      }
    }
    // And return.
    return(std::make_tuple(g, v, d));
  }
  Polynomial Polynomial::readIn(std::istream &is) {
    int monomials, vars;
    is >> monomials >> vars;
    MatrixXi exs(monomials, vars);
    VectorXd coes(monomials);
    for(int i = 0; i < monomials; i++) {
      for(int j = 0; j < vars; j++)
	is >> exs(i, j);
      is >> coes(i);
    }
    return(Polynomial(exs, coes));
  }
  std::ostream &operator<<(std::ostream &os, const Polynomial &p) {
    os << p.monomials << ' ' << p.variables << std::endl;
    for(int i = 0; i < p.monomials; i++)
      os << p.exponents().row(i) << ' ' << p.coefficients(i) << std::endl;
    return(os);
  }

  Eigen::VectorXd Polynomial::derivative(Eigen::VectorXd xs) {
    Eigen::VectorXd ret;
    ret.setZero(vars);
    for(int i = 0; i < monomials; i++) {
      double ai = coes(i);
      for(int j = 0; j < vars; j++) {
	unsigned v = exs(i, j);
	double d = xs(j);
	while(v) {
	  if(v & 1)
	    ai *= d;
	  v >>= 1;
	  d *= d;
	}
      }
      for(int j = 0; j < vars; j++)
	if(exs(i, j))
	  ret(j) += v * exs(i, j) / xs(j);
      
    }
  }
}
