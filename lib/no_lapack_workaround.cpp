#include <armadillo>

using namespace arma;

// Stores L in the lower left, and D in the diagonal.
static void LDL_decompose(mat &a) {
  for(int i = 0; i < a.n_rows; i++)
    for(int j = 0; j <= a.n_rows; j++) {
      for(int k = 0; k < j; k++)
	a(i, j) -= a(i, k) * a(j, k) * a(k, k);
      if(i != j)
	a(i, j) /= a(j, j);
    }
}

// Lower triangular.
static void forsub(const mat &a, vec &b) {
  for(int i = 1; i < a.n_rows; i++)
    for(int j = 0; j < i; j++)
      b(i) -= b(j) * a(i, j);
}

// Upper triangular, but works with the lower triangle.
static void revsub(const mat &a, vec &b) {
  for(int i = a.n_rows - 2; i >= 0; i--)
    for(int j = i + 1; j < a.n_rows; j++)
      b(i) -= a(j, i) * b(j);
}

// Inverts the input lower-triangular matrix, with ones on the diagonal.
static mat ltinv(const mat &a) {
  mat b(a.n_rows, a.n_rows, fill::eye);
  b.each_col([&](vec &bc) {forsub(a, bc); });
  return(b);
}

vec workaround_solve(mat A, vec B) {
  LDL_decompose(A);
  forsub(A, B);
  B /= A.diag();
  revsub(A, B);
  return(B);
}
mat workaround_invert(mat A) {
  LDL_decompose(A);
  mat lowinv = ltinv(A);
  return(lowinv.t() * diagmat(vec(A.n_rows, fill::ones) / A.diag()) * lowinv);
}

