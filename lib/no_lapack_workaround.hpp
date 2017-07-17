#pragma once

// NOTE: If LAPACK is not installed, then this breaks down if the matrix is not
// symmetric.
#ifdef ARMA_USE_LAPACK
#define SOLVE solve
#define INVERT inv_sympd
#else
#define SOLVE workaround_solve
#define INVERT workaround_invert
#endif
#include <armadillo>
extern arma::vec &&workaround_solve(arma::mat, arma::vec);
extern arma::mat &&workaround_invert(arma::mat);
