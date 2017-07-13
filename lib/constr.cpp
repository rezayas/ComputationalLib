#include <armadillo>
#include <list>
#include <utility>
#include <algorithm>
#include <cassert>

using namespace arma;
using namespace std;  

vec theta;
static mat G, Gi;
static int dim, equal;
static vec c;
static list<pair<vec, double> > equalities, inequalities, w;

void setupAb(mat &A, vec &b);

void setup(int dimension, list<pair<vec, double> > e,
	   list<pair<vec, double> > i) {
  c.set_size(dimension);
  theta.set_size(dimension);
  G.set_size(dimension, dimension);
  Gi.set_size(dimension, dimension);
  dim = dimension;
  equalities = e;
  inequalities = i;
  w = e;
}

void addConstr(vec a, mat &h, mat &hi, int size) { // See paper for details
  a.resize(h.n_rows); // Lengthen it (so as to avoid throwing errors)
  vec u = a - h * hi * a; // u = (1-H H†)a
  if(any(abs(u) > 1e-4)) { // u != 0
    vec k = hi * a; // k = -H†a
    rowvec ui = u.t() / as_scalar(u.t() * u); // u† = uT/uuT
    double zi = as_scalar(a.t() * hi * a); // z = -1/(aT H† a)
    hi -= k * ui + ui.t() * k.t(); // H†' = [[A v][vT 0]]; A = H† - k u - u† kT
    hi += ui.t() * ui * zi; // - (u uT)†/z
    hi.row(size) = ui;
    hi.col(size) = ui.t();
  } else {
    double z = -1 / as_scalar(a.t() * hi * a); // z = -1/(aT H† a)
    vec w = -z * hi * a; // w = -H† a z
    hi -= hi * a * w.t(); // H†' = [[A w][wT z]]; A = H†(1 - a wT)
    hi.row(size) = w.t();
    hi.col(size) = w;
    hi(size, size) = z;
  }
  h.row(size) = a.t();
  h.col(size) = a;
}

void addConstr(const mat &A, mat &h, mat &hi) {
  int size = h.n_rows + equalities.size() + inequalities.size();
  int hsize = h.n_rows;
  h.resize(size, size); // Increase size of H to handle it
  hi.resize(size, size); // And H†
  int i;
  for(i = 0; i < A.n_rows; i++, hsize++)
    addConstr(A.row(i).t(), h, hi, hsize); // Recursively add constraints.
}

void setupAb(mat &A, vec &b) {
  A.set_size(w.size(), dim);
  b.set_size(w.size());
  int i = 0;
  for_each(begin(w), end(w), [&](pair<vec, double> p) {
      A.row(i) = p.first.t();
      b(i) = p.second;
    }); // Sets up A, b to contain what they should from w.
}

void fulladdConstr(mat &A, vec &b, mat &h, mat &hi, pair<vec, double> a) {
  // Adds a constraint to A, b, w.
  int size = A.n_rows + dim;
  A.resize(A.n_rows + 1, dim); // Add a row to A;
  b.resize(A.n_rows); // And to b.
  A.row(A.n_rows - 1) = a.first.t(); // Set that row correctly
  b(b.n_elem - 1) = a.second;
  addConstr(a.first, h, hi, size); // Add the constraint
  w.push_back(a); // And add it to w.
}

void removeConstr(mat &A, vec &b, mat &h, mat &hi, pair<vec, double> a) {
  bool ret = false;
  for_each(begin(equalities), end(equalities), [a, &ret](pair<vec, double> p) {
      ret |= p.second == a.second && all(p.first == a.first);
    });
  if(ret) return;
  w.remove_if([&a](pair<vec, double> p) {
      return(p.second == a.second && all(p.first == a.first));
      });  // Throw out the constraint from w
  h = G; // Reset H to G
  hi = Gi;
  setupAb(A, b); // Set up the constraints again
  addConstr(A, h, hi); // And put them back in H
}

void solveUp() { // Standard quadratic programming.
  mat A, H = G, Hi = Gi;
  vec b;
  setupAb(A, b);
  addConstr(A, H, Hi);
  while(true) {
    vec g = G * theta + c; // g = Gθ + c
    g.resize(dim + equalities.size() + inequalities.size());
    // lengthen that
    vec pmu = Hi * g; // [[G AT][A 0]][[-p][μ]] = [[g][0]]
    vec p = -pmu;
    p.resize(dim); // p now just holds the p vector
    if(any(abs(p) > 1e-4)) { // It isn't 0
      pair<vec, double> block;
      bool isblock = false;
      double alpha = 1;
      for_each(begin(inequalities), end(inequalities),
	       [&](pair<vec, double> pr) {
		 bool fnd = false;
		 for_each(begin(w), end(w), [pr, &fnd](pair<vec, double> p) {
		     fnd |= p.second == pr.second && all(p.first == pr.first);
		   });
		 if(!fnd) {
		   double poss = pr.second - as_scalar(pr.first.t() * theta);
		   double div = as_scalar(pr.first.t() * p);
		   poss /= div;// See eqn 29.
		   if(poss < alpha && div < 0 && poss > 0) {
		     alpha = poss;
		     isblock = true;
		     block = pr;
		   }
		 }
	       });
      theta += alpha * p;
      if(isblock) {
	fulladdConstr(A, b, H, Hi, block); // Add a blocking constraint
      }
      //if there is one
    } else {
      int es = equalities.size();
      vec mu(A.n_rows - es);
      int i;
      for(i = 0; i < A.n_rows - es; i++)
	mu(i) = pmu(i + dim + es); // Set mu to hold the μ vector
      if(any(mu < 0) && mu.n_rows) {
      	uword j;
	mu.min(j);
	removeConstr(A, b, H, Hi, make_pair(A.row(j + es).t(), b(j + es)));
	//Remove said constraint
      } else return;
    }
  }
}

inline mat newG(vec x, double lambda, mat g) {
  return(lambda * g + x * x.t()); // G' = λG + x xT
}

inline vec newC(vec x, double y, double lambda, vec c) {
  return(lambda * c - y * x); // c' = λc - xy
}

mat newGi(vec &x, const mat &g, const mat &gi, double lambda) {
  // See eqn. 34
  vec u = x - g * gi * x; // u = (1 - GG†)x
  double gamma = lambda + as_scalar(x.t() * gi * x); // γ = λ + xT G† x
  mat ngi;
  if(!any(abs(u) > 1e-4))
    ngi = gi - gi * x * x.t() * gi / gamma; // G†' = G† - G† x xT G† / γ
  else {
    rowvec ui = u.t() / as_scalar(u.t() * u); // u† = uT/(uT u)
    ngi = gi * x * ui; // G†' = -G† x u†
    ngi += ngi.t(); // - u†T xT G†
    ngi = gi - ngi + gamma * ui.t() * ui; // + G† + (u uT)† γ
  }
  return(ngi / lambda); // but we must divide by λ.
}

void addObs(vec x, double y, double lambda) {
  static bool foo = false;
  if(foo) {
    Gi = newGi(x, G, Gi, lambda); // Correct first G†,
    G = newG(x, lambda, G); // then G and c.
    c = newC(x, y, lambda, c);
    solveUp();
  } else {
    foo = true;
    G = x * x.t();
    rowvec xi;
    pinv(xi, x);
    Gi = xi.t() * xi;
    c = x * y;
    mat H = G;
    mat Hi = Gi;
    vec b;
    mat A;
    setupAb(A, b);
    addConstr(A, H, Hi);
    b.insert_rows(0, -c);
    b.resize(Hi.n_rows);
    theta = Hi * b;
    theta.resize(dim);
  }
}
