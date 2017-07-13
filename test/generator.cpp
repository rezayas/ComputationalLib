#include <armadillo>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <random>

using namespace arma;
using namespace std;

void read(vec &means, mat &sigmas, istream &desc, int dim) {
  int i, j;
  for(i = 0; i < dim; i++)
    desc >> means(i);
  for(i = 0; i < dim; i++)
    for(j = 0; j <= i; j++) {
      desc >> sigmas(i, j);
      sigmas(j, i) = sigmas(i, j);
    }
}

void readb(rowvec &theta, double &cterm, double &sy, istream &desc, int dim) {
  int i;
  desc >> cterm;
  for(i = 0; i < dim; i++)
    desc >> theta(i);
  desc >> sy;
}

int main(int argc, char **argv) {
  int i = 0;
  bool constxdist = false;
  while(!strcmp(argv[i + 1], "-n")) {
    i++;
    constxdist = true;
  }
  int n = atoi(argv[i++]);
  while(!strcmp(argv[i + 1], "-n")) {
    i++;
    constxdist = true;
  }
  ifstream desc(argv[i + 1]);
  int dim;
  desc >> dim;
  vec means(dim), means_errors(dim), means_lambdas(dim);
  rowvec theta(dim), theta_errors(dim), theta_lambdas(dim);
  double cterm, cterm_error, cterm_lambda;
  double sy, sy_error, sy_lambda;
  mat sigmas(dim, dim), sigma_errors(dim, dim), sigma_lambdas(dim, dim);
  mat a;
  normal_distribution norm;
  read(means, sigmas, desc, dim);
  readb(theta, cterm, sy, desc, dim);
  if(!constxdist) read(means_errors, sigma_errors, desc, dim);
  readb(theta_errors, cterm_error, sy_error, desc, dim);
  if(!constxdist) read(means_lambdas, sigma_lambdas, desc, dim);
  readb(theta_lambdas, cterm_lambda, sy_lambda, desc, dim);
  if(constxdist)
    a = chol(sigmas).t();
  for(i = 0; i < n; i++) {
    if(!constxdist)
      a = chol(sigmas + sigma_errors).t();
    vec rand(dim);
    vec xs;
    double y;
    rand.randn();
    cout << (xs = (a * rand + means + means_errors));
    y = as_scalar((theta + theta_errors) * xs) + cterm + cterm_error;
    cout << y + norm() * (sy + sy_error) << endl;
    if(!constxdist) {
      sigma_errors %= sigma_lambdas;
      means_errors %= means_lambdas;
    }
    sy_error /= sy_lambda;
    cterm_error /= cterm_lambda;
    theta_errors %= theta_lambda;
  }
}
    
