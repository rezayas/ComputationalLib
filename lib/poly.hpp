#pragma once
#include <Eigen/Core>

namespace linreg {
  class Polynomial {
  public:
    const Eigen::MatrixXi &exponents;
    Eigen::VectorXd coefficients;
    Polynomial(int order, int nvars);
    inline Polynomial(const Eigen::MatrixXi &exs, const Eigen::VectorXd &coes)
      : exps(exs), exponents(exps), coefficients(coes), variables(exs.cols()),
	monomials(exs.rows()) {}
    double evaluate(const Eigen::VectorXd &) const;
    // Returns the products in the proper order.
    Eigen::VectorXd expand(const Eigen::VectorXd &) const;
    void fitToData(const Eigen::MatrixXd &, const Eigen::VectorXd &,
		   double = 1, double = 0);
    const int variables, monomials;
  private:
    Eigen::MatrixXi exps;
  };
}
