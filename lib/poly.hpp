#pragma once
#include <Eigen/Core>
#include <tuple>

namespace linreg {
  class Polynomial {
  public:
    inline const Eigen::MatrixXi &exponents() const {
      return(exps);
    }
    Eigen::VectorXd coefficients;
    Polynomial(int order, int nvars);
    inline Polynomial(const Eigen::MatrixXi &exs, const Eigen::VectorXd &coes)
      : exps(exs), coefficients(coes), variables(exs.cols()),
	monomials(exs.rows()) {}
    double evaluate(const Eigen::VectorXd &) const;
    // Returns the products in the proper order.
    Eigen::VectorXd expand(const Eigen::VectorXd &) const;
    void fitToData(const Eigen::MatrixXd &, const Eigen::VectorXd &,
		   double = 1, double = 0);
    // Returns (g, v, c) such that f(x) = xᵀgx/2 + xᵀv + c,
    // plus possibly terms that are third-order or higher in x.
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd, double> quaddec() const;
    const int variables, monomials;
  private:
    Eigen::MatrixXi exps;
  };
}
