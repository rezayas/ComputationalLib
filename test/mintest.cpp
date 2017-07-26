#include <iostream>
#include "poly.hpp"
#include <Eigen/Core>
#include <Eigen/Cholesky>

using namespace std;
int main() {
  // 3x^2 + 2x + 1
  linreg::Polynomial f(Eigen::VectorXi::LinSpaced(3, 0, 2),
		       Eigen::VectorXd::LinSpaced(3, 1, 3));

  // 1 + 2x + 3y + 4xy + 5x^2 + 6y^2
  linreg::Polynomial g((Eigen::MatrixX2i(6,2) << 0,0, 1,0, 0,1,
			1,1, 2,0, 0,2).finished(),
		       Eigen::VectorXd::LinSpaced(6, 1, 6));
  
  Eigen::MatrixXd fmat, gmat;
  Eigen::VectorXd fvec, gvec, ftheta, gtheta;
  // How to get out the matrix and vector.
  // To minimize, you do not need the scalar.
  tie(fmat, fvec, ignore) = f.quaddec();
  tie(gmat, gvec, ignore) = g.quaddec();
  
  // How to find the absolute minimum.
  // If we wanted the minimum in some convex area,
  // we'd want to use linreg::solve_quadprog.
  // This also requires that we have the absolute minimum and the inverse
  // of the G-matrix at hand.
  ftheta = fmat.ldlt().solve(fvec);
  gtheta = gmat.ldlt().solve(gvec);

  cout << ftheta << endl << f.evaluate(ftheta) << endl << endl;
  cout << gtheta << endl << g.evaluate(gtheta) << endl;
}
