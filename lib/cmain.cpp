#include <armadillo>
#include <list>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>

using namespace arma;
using namespace std;


extern void setup(int, list<pair<vec, double> >, list<pair<vec, double> >);
extern void addObs(vec, double, double);
extern vec theta;

int main(int argc, char **argv) {
  int dim = atoi(argv[1]) + 1;
  float lambda = stof(argv[4]);
  ifstream equalities(argv[2]);
  ifstream inequalities(argv[3]);
  list<pair<vec, double> > e, i;
  vec x(dim);
  double y;
  int j;
  while(!equalities.eof()) {
    for(j = 0; j < dim; j++)
      equalities >> x(j);
    equalities >> y;
    if(!equalities.eof())
      e.push_front(make_pair(x, y));
  }
  while(!inequalities.eof()) {
    for(j = 0; j < dim; j++)
      inequalities >> x(j);
    inequalities >> y;
    if(!inequalities.eof())
      i.push_front(make_pair(x, y));
  }
  setup(dim, e, i);
  x(0) = 1;
  while(!cin.eof()) {
    for(j = 1; j < dim; j++)
      cin >> x(j);
    cin >> y;
    addObs(x, y, lambda);
#ifdef ALL
    cout << theta.t();
#endif
  }
#ifndef ALL
  cout << theta.t();
#endif
}
