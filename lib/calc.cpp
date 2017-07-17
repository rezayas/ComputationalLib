#include <armadillo>
#include "no_lapack_workaround.hpp"

using namespace arma;
vec thetaOLS(double lambda, mat x, vec y) {
  mat w(x.n_rows, x.n_rows, fill::zeros);
  double p;
  for(int i = x.n_rows - 1, p = 1; i >= 0; i--, p /= lambda)
    w(i, i) = p;
  return(SOLVE(x.t() * w * x, x.t() * w * y));
}

vec theta;
mat b;
static mat xmat;
static vec yvec;
static int points;
static int dim;
bool ready;

void setup(int dimension) {
  xmat.set_size(dimension, dimension);
  theta.set_size(dimension);
  yvec.set_size(dimension);
  dim = dimension;
  points = 0;
  ready = false;
}

bool updateTheta(double lambda, vec x, double y) {
  if(!ready) {
    int i;
    for(i = 0; i < dim; i++)
      xmat(points, i) = x(i);
    yvec(points) = y;
    points++;
    if(points == dim) {
      ready = true;
      b = INVERT(xmat.t() * xmat);
      theta = thetaOLS(lambda, xmat, yvec);
      xmat.reset();
      yvec.reset();
    }
  } else {
    double gamma = lambda + as_scalar(x.t() * b * x);
    theta += b * x * (y - theta.t() * x) / gamma;
    b -= b * x * x.t() * b / gamma;
    b /= lambda;
  }
  return(ready);
}
