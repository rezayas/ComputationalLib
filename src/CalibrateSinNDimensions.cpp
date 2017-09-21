#include "../include/ComputationalLib/CalibrateSinNDimensions.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace ComputationalLib;
using namespace SimulationLib;

using std::get, std::min, std::max;

// Comment the following line to disable logging
#define DEBUG

// Constants for algorithm
#define OMEGA    0.001
#define ALPHA    0.05
#define EPSILON  0.00001

#define B_       50
#define MAXITERS 1000

// Are all arguments truthy?
template <typename... InTs>
bool all(InTs... args) { return (... && args); }

// Are any arguments truthy?
template <typename... InTs>
bool any(InTs... args) { return (... || args); }

// Returns a function 'f: Ts... => I' whose output
// is always 'value'
template <typename I,
          typename... Ts>
std::function<I(Ts...)>
IdentityFnGen(I value)
{
    return [=value] (Ts&&... params) -> I {
        return value;
    };
}

template <typename... Xts,
          size_t...     I,
          typename     Xs = std::tuple<Xts...>,
          typename     Fn = std::function<double(Xts...)>,
Xs CalibrateSin(Xs x_i,
                const Fn &f,
                const size_t &d,
                std::index_sequence<I...>)
{
    // Create a function 'Zero: size_t => double' which always
    // returns zero
    auto Zero = IdentityFnGen<double, size_t>(0);

    // Create a tuple of 'previous' values, which in this case
    // we assume to be zero
    Xs x_prev {Zero(I)...};

#ifdef DEBUG
    // Open a CSV file to record the progress of the calibration
    // algorithm
    ofstream fout("CalibrateSin.csv");
    fout << "Iteration, [Xs], f(Xs...)" << endl;
#endif

    // Omega: The penalty for regressing towards steeper slopes,
    //   used in the construction of the polynomial used in regression
    // Alpha: The step size. It is multiplied by the slope to calculate
    //   the distance to the next vector x_i
    // Epsilon: The minimum distance for each step. When a step becomes
    //   smaller than epsilon, the algorithm terminates
    // 'b': A constant used to calculate Lambda, the learning weight.
    //   As 'b' increases, lambda approaches 1. Therefore, the values of
    //   'b' and 'maxIters' should be chosen such that "b is very
    //   close to 1 as b approaches 'maxIters'??
    double omega   {OMEGA};
    double alpha   {ALPHA};
    double epsilon {EPSILON};

    int b          {B_};
    int maxIters   {MAXITERS};
    int idx        {1};

    // Form two tuples of lower and upper bounds on each x in X
    Xs lowers {get<I>(x_i).Lower...};
    Xs uppers {get<I>(x_i).Upper...};

    // Set up a regression using a 2nd-degree polynomial for a
    // function f(X) where the set X has 'd' members
    PolynomialRegression Regression(2, d, omega);

    // The algorithm must stop iterating when we hit 'maxIters', or
    // when one member of x_i fails to change value by more than
    // 'epsilon'
    auto terminationPred = [&] (const auto &idx,
                                const auto &x_i,
                                const auto &x_prev,
                                const auto &epsilon) -> bool {

        bool hitMaxIters {idx >= maxIters};

        bool insignificantDifference {
          any(
          // all( Do we halt on the first parameter satisfying the
          //      condition or do we wait until all of them satisfy?
            (abs(get<I>(x_i)() - get<I>(x_prev)()) < epsilon)...
             )
        };

        // Terminate if either of these conditions is satisfied
        return hitMaxIters || insignificantDifference;
    };

    do {
        // Calculate learning weight for regression
        double lambda_i { (double)idx / (double)(b + idx) };
        double y_i      { f(x_i) };

        // Represent x_i as an Eigen vector for use with Regression obj.
        // Comma-initializer syntax described at:
        //   https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
        Eigen::VectorXd x_i_E(d);
        x_i_E << get<I>(x_i)...;

        // Add a new point to the regression
        Regression.updateCoefficients(x_i_E, y_i, lambda_i);

        // Numerically calculate the derivative of each value x in X at the
        // point X
        Eigen::VectorXd slope_i_E = Regression.poly().derivative(x_i_E);

        // Prepare for modification of 'x_i'
        x_prev = x_i;

        // For each x in X, travel along the slope a little bit in
        // the positive direction
        x_i = { max(get<I>(lowers),
                    min(get<I>(uppers),
                        get<I>(x_i)() + (alpha*slope_i_E[I])) )... };

#ifdef DEBUG
        fout << idx << "," << "[" << (string(get<I>(x_prev)) + ",")... << "]," << y_i << endl;
#endif

        idx += 1;
    } while (!terminationPred(idx, x_i, x_prev, epsilon));

#ifdef DEBUG
    fout.flush();
    fout.close();
#endif

    return x_i;
}

template <typename... Xts,
          typename     Xs = std::tuple<Xts...>,
          typename     Fn = std::function<double(Xts...)>,
          size_t        d = sizeof...(Xts),
          typename     IS = std::make_index_sequence<Xts...> >
Xs CalibrateSin(const Xs &x_i, const Fn &f)
{
    return CalibrateSin(std::forward<decltype(x_i)>(x_i),
                        std::forward<decltype(f)>(f),
                        d
                        IS{});
}