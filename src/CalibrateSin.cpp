#include "../include/ComputationalLib/CalibrateSin.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

UnitDouble CalibrateSin(UnitDouble _x, F _f) {
	// This model takes in a UnitDouble _x and a function _f
	// There are several parameters to this model:
	//  - omega:   this is the penalty for steeper slopes; it is
	//			   used in the construction of the Polynomial
	//             Regression
	//  - alpha:   alpha is the step size; it is multiplied
	//             by the slope to calculate the step to the 
	//             next value of x_i
	//  - epsilon: epsilon is baseline minimum distance for each 
	//             step; if a step is smaller than epsilon, then the
	//             regression stops
	//  - b:       b is a constant that is used to calculate the 
	//             learning weight, lambda. lambda approaches 1, as 
	//             b increases. Thus, the values of b and maxIterations
	//             should be chosen such that b is very close to 1
	//             as b appraoches maxIterations


	UnitDouble x_i    {_x};
	UnitDouble x_prev {0};

	F f = _f;

	double lower, upper;

    double omega, alpha, epsilon;

	double y_i, lambda_i, slope_i;

	int maxIterations, b, idx;

	// Will write out data into csv file
	ofstream fout("CalibrateSin.csv");
	fout << "Iteration,x,f(x)" << endl;

	f = _f;

	// These are the paramters for the model:
	omega   = 0.001;
	alpha   = 0.05;
	epsilon = 0.00001;

	maxIterations = 1000;
	b             = 50;

	// This is the lambda index used with b to calculate learning weight
	idx = 1;

	// Use the bounds from the UnitDouble x to set the bounds
	// for this regression
	lower = x_i.Lower;
	upper = x_i.Upper;

	// First create a quadratic polynomial with coefficients of all 0
	// exs:   Eigen::VectorXi::LinSpaced(3, 0, 2) = [0,1,2]
	// coefs: Eigen::VectorXd::LinSpaced(3, 0, 0) = [0.0, 0.0, 0.0]
	Polynomial g( Eigen::VectorXi::LinSpaced(3, 0, 2),
				  Eigen::VectorXd::LinSpaced(3, 0, 0) );

	// Initialize the poly regression with the polynomial and omega
	PolynomialRegression PR_g(g, omega);

	// Add preliminary values to get preliminary coefficients
	// r is the radius for these preliminary random values
	double r = 0.1;
	for (int i = 0; i < 5; ++i) {
		// there is probably a better way to randomize
		// randomDouble is a random double between -1.0 and 1.0
		int randomInt = rand() % 10000;
		double randomDouble = (double)randomInt/5000.0 - 1.0;

		double x_r, y_r;
		// x_r is some random value near x
		x_r = x_i() + randomDouble * r;

		if (x_r < lower || x_r > upper)
			--i;
		else {
			// Update the coefficients, updateCoefficients
			// function needs to be given a 1x1 matrix
			y_r = f(x_r);
			Eigen::VectorXd x_matrix_r(1);
			x_matrix_r(0) = x_r;
			PR_g.updateCoefficients(x_matrix_r, y_r, 0.018);

			// cout << x_r << ": " << y_r << endl;
		}
	}

	// Do the polynomial regression
	do {
		// Get the value at x_i
		y_i = f(x_i);

		// Calculate the learning weight
		lambda_i = (double) idx / (double) (b + idx);

		// x value need to be converted to matrix to use
		// updateCoefficients function so we create a 1x1 matrix
		Eigen::VectorXd x_matrix(1);
		x_matrix(0) = x_i();

		// Update the coefficients
		// PR_g.updateCoefficients(x_matrix, y_i);
		PR_g.updateCoefficients(x_matrix, y_i, lambda_i);

		// Calculate the derivative, for now manually
		Eigen::VectorXd coefficients = PR_g.poly().coefficients;

		// Coefficients will be a vector : [a b c]
		// correspdonding to: (a + bx + cx^2)
		// Derivative of (a + bx + cx^2) is : (b + 2cx)
		// Find slope by evaluating derivative at x = x_i
		slope_i = coefficients(1) + 2 * coefficients(2) * x_i();

		// Set x_prev to previous value of x_i, update x_i
		x_prev = x_i();
		x_i = max(lower,
				  min(upper, x_i() + alpha * slope_i));

		// Write out data to csv:
		fout << idx << "," << x_prev() << "," << y_i << endl;

		// // Debugging print statements
		// printf("n: %d x_i: %f x_prev: %f slope: %f\n", idx, x_i(), x_prev(), slope_i);
		// cout << "Coefficents:\n" << coefficients << endl;

		idx++;
	} while (idx < maxIterations && !(abs(x_i() - x_prev()) < epsilon));

	// Flush and close output buffer
	fout.flush();
	fout.close();

	return x_i;
}



