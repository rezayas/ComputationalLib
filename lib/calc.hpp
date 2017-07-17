#pragma once
#include <armadillo>

extern void setup(int);
extern bool updateTheta(double lambda, arma::vec x, double y);
extern arma::vec theta;
extern arma::mat b;
extern bool ready;

