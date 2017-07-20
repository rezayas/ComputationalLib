#include "linreg.hpp"
#include "poly.hpp"

namespace linreg {
  class PolynomialRegression {
  public:
    inline PolynomialRegression(int order, int nvars, double omega = 0) :
      PolynomialRegression(Polynomial(order, nvars), omega) {}
    inline PolynomialRegression(const Polynomial &p, double omega = 0) :
      lin(p.monomials), ppoly(p), poly(ppoly) {}
    const Polynomial &poly;
    bool updateCoefficients(const vec &x, double y, double lambda = 1);
  private:
    LinearRegression lin;
    Polynomial ppoly;
  };
}
