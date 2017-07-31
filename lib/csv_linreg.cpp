#include "csv_linreg.hpp"
#include "linreg.hpp"
#include <cstdio>

namespace complib {
  Eigen::VectorXd regressFile(const char *fname, int dim,
			      double lambda, double omega) {
    LinearRegression lr(dim, omega);
    FILE *f = fopen(fname, "r");
    Eigen::VectorXd x(dim);
    double y;
    int c;
    while(!feof(f)) {
      for(int i = 0; i < dim; i++)
	fscanf(f, "%lf,", &x(i));
      fscanf(f, "%lf", &y);
      while((c = getc(f)) != EOF && c != '\n');
      lr.updateCoefficients(x, y, lambda);
    }
    return(lr.getCoefficients());
  }
}
