#pragma once
#include <Eigen/Core>

namespace complib {
  extern Eigen::VectorXd regressFile(const char *fname, int dim,
				     double lambda = 1, double omega = 0);
}
