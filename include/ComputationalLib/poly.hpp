#pragma once
#include <Eigen/Core>
#include <tuple>
#include <ostream>
#include <istream>
#include <vector>

namespace ComputationalLib {
  class Polynomial {
  public:

    // There were no line breaks in this file, so I (Nathan) tried my 
    // best to associate the comments with the correspdoning functions



    inline const Eigen::MatrixXi &exponents() const {
      return(exps);
    }


    Eigen::VectorXd coefficients;

    
    // Returns the zero polynomial with space to be any polynomial of
    // the proper order on that many variables.
    // You can later edit the coefficients.
    static Polynomial readIn(std::istream &is);
    Polynomial(int order, int nvars);




    // Let m be the number of rows in exs and coes.
    // Let n be the number of columns in exs.
    // Then the polynomial represented here is:
    // ∑_i=1^m coes(i)*∏_j=1^n x_j^exs(i,j).
    // For example, if we want a + bx + cx^2,
    // then setting exs to [0, 1, 2]ᵀ and coes to [a, b, c]ᵀ
    // works quite well.
    // Permuting the rows of both does not change the polynomial.
    // If we want a + bx + cy + dx^2 + exy + fy^2,
    // we set exs to [[0,0],[1,0],[0,1],[2,0],[1,1],[0,2]]
    // and coes to [a,b,c,d,e,f]ᵀ.
    inline Polynomial(const Eigen::MatrixXi &exs, const Eigen::VectorXd &coes)
      : exps(exs), coefficients(coes), variables(exs.cols()),
	monomials(exs.rows()) {}



    double evaluate(const Eigen::VectorXd &) const;
    // Returns the products in the proper order.
    // Ignores the coefficients.
    Eigen::VectorXd expand(const Eigen::VectorXd &) const;
    Polynomial &fitToData(const Eigen::MatrixXd &, const Eigen::VectorXd &,
			  double = 1, double = 0);

    // Returns (g, v, c) such that f(x) = xᵀgx/2 + xᵀv + c,
    // plus possibly terms that are third-order or higher in x.
    // In other words,
    // we express f(x) as ∑_i=0^n f_i(x) where f_i is homogenous of degree i
    // (g is homogenous of degree n if g(λx) = λⁿg(x) for any x),
    // and return a matrix G such that f_2(x) = xᵀGx/2,
    // a vector v such that f_1(x) = vᵀx,
    // and a scalar c such that f_0(x) = c.
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd, double> quaddec() const;
    const int variables, monomials;
    Eigen::VectorXd derivative(const Eigen::VectorXd &) const;
    std::vector<Polynomial> formalDerivatives() const;
    Polynomial formalDerivative(int) const;

  private:
    Eigen::MatrixXi exps;
  };
  std::ostream &operator<<(std::ostream &os, const Polynomial &p);
}
