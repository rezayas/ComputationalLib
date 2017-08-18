#include "../include/ComputationalLib/CalibrateSin.hpp"
#include <iostream>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

// BoundDouble x{RatioToNumeric<double>(PI_ratio{})};

BoundDouble CalibrateSin(BoundDouble _x, F _f) {
	BoundDouble x_i {0};
	BoundDouble x_prev {0};
	F f;
	double lower, upper, omega, alpha, epsilon;
	double y_i, lambda_i, slope_i;
	int maxIterations, b, idx;

	vector<double> xs;
	vector<double> ys;

	x_i = _x;
	f = _f;

	// omega is the penalty; it is one of the variables used to
	// construct the polynomial regression
	omega = 0.001;

	// alpha is the step size
	alpha = 0.01;

	// epsilon is the distance from prev to current used to 
	// check if the step was small enough to end the regression
	epsilon = 0.00001;

	// b is the constant used to define the learning weight, lambda
	// idx is the current iteration index, starts at 1
	// the learning weight is: idx/(b + idx)
	b = 50;
	idx = 1;

	maxIterations = 1000;

	// use the bounds from the BoundDouble x to set the bounds for the regression 
	lower = x_i.Lower;
	upper = x_i.Upper;


	// First create a quadratic polynomial with coefficients of all 0
	// exs: Eigen::VectorXi::LinSpaced(3, 0, 2) = [0,1,2]
	// coefs: Eigen::VectorXd::LinSpaced(3, 0, 0) = [0.0, 0.0, 0.0]
	Polynomial g(Eigen::VectorXi::LinSpaced(3, 0, 2),
		         Eigen::VectorXd::LinSpaced(3, 1, 1));

	// Initialize the poly regression with the polynomial and omega
	PolynomialRegression PR_g(g, omega);

	// Attempted to add preliminary values; However, causes the regression to be
	// skewed because of equal weight
	double r = 0.1;
	for (int i = 0; i < 5; ++i) {
		// randomDouble is a random double between -1.0 and 1.0
		int randomInt = rand() % 10000;
		double randomDouble = randomInt/5000 - 1.0;

		double x_r, y_r;
		x_r = x_i() + randomDouble * r;

		if (x_r < lower || x_r > upper) 
			--i;
		else {
			y_r = f(x_r);
			Eigen::VectorXd x_matrix_r(1);
			x_matrix_r(0) = x_r;
			PR_g.updateCoefficients(x_matrix_r, y_r);
		}
	}

	// do 
	do {
		// get the value at x_i
		y_i = f(x_i());

		// calculate the learning weight
		lambda_i = idx / (b + idx);

		// x value need to be converted to matrix to use
		// updateCoefficients
		Eigen::VectorXd x_matrix(1);
		x_matrix(0) = x_i();

		// update the coefficients, ideally you would use lambda, but
		// the code for the lambda isn't working
		PR_g.updateCoefficients(x_matrix, y_i);
		// PR_g.updateCoefficients(x_matrix, y_i, lambda_i);

		// calculate the derivative, for now manually
		Eigen::VectorXd coefficients = PR_g.poly().coefficients;

		// coefficients will be a vector : [a b c]
		// correspdonding to: a + bx + cx^2 

		// derivative of a + bx + cx^2 is : b + 2cx
		// find slope by evaluating derivative at x = x_i
		slope_i = coefficients(1) + 2 * coefficients(2) * x_i();

		// set x_prev to previous value of x_i, update x_i
		x_prev = x_i();
		x_i = max(lower,
				  min(upper, x_i() + alpha * slope_i));


		// push values of x and y onto vectors for data output
		// although because this isn't a real class, 
		// there isn't any way to actually access these values
		// when the function resolves
		xs.push_back(x_prev());
		ys.push_back(y_i);

		idx++;

		// Debugging print statements
		printf("n: %d x_i: %f x_prev: %f slope: %f\n", idx, x_i(), x_prev(), slope_i);

		cout << "Coefficents:\n" << coefficients << endl;
	} while (idx < maxIterations && !(abs(x_i() - x_prev()) < epsilon));

	return x_i;
}



