#include "../include/ComputationalLib/polyreg.hpp"

#include <cmath>

#include <Bound.h>
#include <TimeSeries.h>
#include <Normal.h>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;
using namespace StatisticalDistributions;

using TWO_PI_ratio = std::ratio_divide<PI_ratio, std::ratio<1,2>>;

// for now BoundDouble is hard coded with bounds, 0, 2PI
using BoundDouble = Bound<double, std::ratio<0,1>, TWO_PI_ratio>;

double standardDeviation = 1;
StatisticalDistributions::Normal nDist(0, standardDeviation);

// A Mersenne Twister pseudo-random generator of 64-bit 
// numbers with a state size of 19937 bits.
std::mt19937_64 rands(std::time(NULL));

using F = function<double(BoundDouble)>;
F f = [] (BoundDouble x) {
	return std::sin(x()) + nDist(rands); 
}; 

BoundDouble x{RatioToNumeric<double>(PI_ratio{})};

BoundDouble CalibrateSin(BoundDouble _x, F _f) {
	BoundDouble x_i {0};
	BoundDouble x_prev {0};
	F f;
	double lower, upper, omega, alpha, epsilon, y_i, learningWeight_i;
	int maxIterations, b, idx;

	x_i = _x;
	f = _f;

	// omega is the penalty; it is one of the variables used to
	// construct the poly reg
	omega = 0.001;

	// alpha is the step size
	alpha = 0.05;

	// epsilon is the distance from prev to current used to 
	// check if the step was small enough to end the regression
	epsilon = 0.001;

	// b is the constant used to define the learning weight
	// idx is the current iteration index, starts at 1
	// the learning weight is: idx/(b + idx)
	b = 50;
	idx = 1;

	maxIterations = 1000;

	// use the bounds from the BoundDouble x to set the bounds for the regression 
	lower = x.Lower;
	upper = x.Upper;


	// First create a quadratic polynomial with coefficients of all 0
	// Eigen::VectorXi::LinSpaced(3, 0, 0) = [0,0,0]
	// Eigen::VectorXd::LinSpaced(3, 1, 3) = [1.0, 2.0, 3.0]
	Polynomial g(Eigen::VectorXi::LinSpaced(3, 0, 0),
		         Eigen::VectorXd::LinSpaced(3, 1, 3));

	// Initialize the poly regression with the polynomial and omega
	PolynomialRegression PR_g(g, omega);

	// do 
	do {
		// get the value at x_i
		y_i = f(x_i());

		// calculate the learning weight
		learningWeight_i = idx / (b + idx);

		// x value need to be converted to matrix to use
		// updateCoefficients
		Eigen::VectorXd x_matrix(1);
		x_matrix(0) = x_i();

		// update the coefficients
		PR_g.updateCoefficients(x_matrix, y_i,learningWeight_i);

		// calculate the derivative, for now manually
		Eigen::VectorXd coefficients = PR_g.poly().coefficients;


		idx++;
	} while (idx < maxIterations || abs(x_i() - x_prev()) < epsilon);



	// for each value of x do the following :

	// - fit points to get values from function

	// - give points, values to regression fit

	// - get constants from fit

	// - take derivative

	// - evaluate slope

	// - choose new value for x

	// break if new x is less than epsilon away from old x






	

}



